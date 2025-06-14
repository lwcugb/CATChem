!> \file catchem_nuopc_cap.F90
!! \brief NUOPC cap for CATChem atmospheric chemistry model
!!
!! \defgroup catchem_nuopc_group CATChem NUOPC Interface
!! \brief NUOPC interface drivers and utilities for CATChem
!! \ingroup catchem
!!
!! This group contains all NUOPC-compliant interface modules and utilities
!! for integrating CATChem with the NUOPC framework, including the main
!! cap module, data transformation utilities, and I/O capabilities.
!!
!! \details
!! This module provides the NUOPC (National Unified Operational Prediction
!! Capability) cap for the CATChem (Configurable ATmospheric Chemistry) model.
!! The cap enables CATChem to run within the NUOPC/ESMF framework as a
!! component in coupled Earth system models.
!!
!! The cap implements the standard NUOPC phases:
!! - Initialize Phase 1: Advertise import and export fields
!! - Initialize Phase 2: Realize fields and initialize the CATChem model
!! - Run: Execute chemistry calculations and data exchange
!! - Finalize: Clean up resources and finalize CATChem processes
!!
!! Key features:
!! - Standard NUOPC/ESMF compliance for easy integration
!! - Flexible field mapping and data exchange
!! - Support for various grid configurations
!! - Configurable chemistry processes and diagnostics
!! - Parallel execution capabilities
!! - Error handling and logging
!!
!! \note This cap follows NUOPC conventions and requires ESMF/NUOPC libraries
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module catchem_nuopc_cap

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS        => SetServices, &
    model_label_Advance => label_Advance, &
    model_label_CheckImport => label_CheckImport, &
    model_label_SetRunClock => label_SetRunClock

  use CATChem
  use catchem_types, only: catchem_container_type
  use catchem_nuopc_interface
  use catchem_nuopc_utils

  implicit none

  private

  public :: SetServices

  !> \brief CATChem component data containers
  !!
  !! These module-level variables maintain the state of the CATChem component
  !! throughout the simulation, including all chemistry, meteorological,
  !! emissions, and diagnostic state information.
  !! \{
  type(catchem_container_type), save :: catchem_states  !< Main CATChem state container
  type(ConfigType), save :: config                      !< CATChem configuration object
  type(DustStateType), save :: dustState               !< Dust emission and transport state
  type(SeaSaltStateType), save :: seaSaltState         !< Sea salt emission and transport state
  type(DryDepStateType), save :: dryDepState           !< Dry deposition process state
  !! \}

  !> \brief Component configuration parameters
  !! \{
  character(len=256), save :: config_file = 'CATChem_config.yml' !< Configuration file path
  logical, save :: do_chemistry = .true.                         !< Enable chemistry calculations
  !! \}

contains

  !> Set services for the CATChem NUOPC cap
  !!
  !! This is the main entry point called by NUOPC to set up the CATChem component.
  !! It registers the initialize, run, and finalize phase entry points and sets
  !! component metadata attributes.
  !!
  !! @param model NUOPC model component to configure
  !! @param rc ESMF return code (ESMF_SUCCESS on success)
  !!
  !! This routine performs the following setup operations:
  !! - Derives the component from the NUOPC model template
  !! - Sets component metadata (name, version)
  !! - Registers entry points for all required NUOPC phases:
  !!   - IPDv00p1: Initialize Phase 1 (advertise fields)
  !!   - IPDv00p2: Initialize Phase 2 (realize fields and initialize model)
  !!   - RunPhase1: Model advance (execute chemistry)
  !!   - FinalizePhase1: Model cleanup
  !!
  !! @note This routine must be called before any other component operations
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Set the model services
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Set component metadata
    call ESMF_AttributeSet(model, name="model_name", value="CATChem", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(model, name="model_version", value="1.0", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Set entry points for standard NUOPC phases
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_RUN, &
      phaseLabelList=(/"RunPhase1"/), userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_FINALIZE, &
      phaseLabelList=(/"FinalizePhase1"/), userRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  !> \brief Initialize Phase 1 - Advertise import and export fields
  !!
  !! In this phase, the component advertises what fields it can import
  !! (receive from other components) and export (provide to other components).
  !! This establishes the interface contract without creating actual field objects.
  !!
  !! \param[inout] model NUOPC model component
  !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
  !!
  !! \details
  !! This routine performs the following operations:
  !! - Retrieves the component's import and export states
  !! - Calls the field advertisement routine to declare all required fields
  !! - Sets up the field interface for meteorological inputs and chemistry outputs
  !! - Enables field matching and connection with other components
  !!
  !! Fields advertised typically include:
  !! - Import: Temperature, humidity, pressure, winds, surface properties
  !! - Export: Chemical species concentrations, deposition fluxes, emissions
  !!
  !! \note This is a standard NUOPC initialization phase that must complete
  !!       successfully before Phase 2 can proceed
  !!
  !! \ingroup catchem_nuopc_group
  subroutine InitializeP1(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    character(len=*), parameter :: routine = 'InitializeP1'

    rc = ESMF_SUCCESS

    ! Get import and export states
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Advertise fields
    call catchem_advertise_fields(importState, exportState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Log successful completion
    call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeP1

  !> \brief Initialize Phase 2 - Realize fields and initialize CATChem model
  !!
  !! In this phase, the component creates the actual ESMF fields for the
  !! advertised imports and exports, and performs complete initialization
  !! of the CATChem model including memory allocation and configuration.
  !!
  !! \param[inout] model NUOPC model component
  !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
  !!
  !! \details
  !! This routine performs comprehensive model initialization:
  !! - Retrieves component information (local PET, PET count)
  !! - Gets the computational grid from the driver
  !! - Determines grid dimensions for memory allocation
  !! - Creates actual ESMF field objects on the grid
  !! - Reads CATChem configuration from file
  !! - Initializes all CATChem state containers and processes
  !! - Sets up chemistry, meteorology, emissions, and diagnostic systems
  !! - Prepares the model for time stepping
  !!
  !! The grid is typically provided by the parent driver or mediator
  !! component and defines the spatial discretization for all field
  !! operations and data exchange.
  !!
  !! \note This phase must complete successfully before any model
  !!       advance operations can be performed
  !!
  !! \ingroup catchem_nuopc_group
  subroutine InitializeP2(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    type(ESMF_Grid) :: grid
    type(ESMF_Config) :: config_esmf
    character(len=*), parameter :: routine = 'InitializeP2'
    character(len=512) :: errmsg
    integer :: localPet, petCount
    integer :: im, jm  ! Grid dimensions

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get import and export states
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get or create grid (this should be provided by the driver)
    call ESMF_StateGet(importState, "grid", grid, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      ! Create a simple grid if not provided
      grid = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), maxIndex=(/144,91/), &
        regDecomp=(/petCount,1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    ! Get grid dimensions
    call ESMF_GridGet(grid, tile=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
      computationalCount=(/im, jm/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Realize fields
    call catchem_realize_fields(importState, exportState, grid, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Read configuration
    config_esmf = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigLoadFile(config_esmf, filename=config_file, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite("CATChem: Could not load config file, using defaults", &
        ESMF_LOGMSG_WARNING, rc=rc)
      rc = ESMF_SUCCESS
    end if

    ! Initialize CATChem using the interface
    call catchem_nuopc_init(config, catchem_states, dustState, seaSaltState, &
                           dryDepState, im*jm, config_file, grid, rc, errmsg)
    if (rc /= CC_SUCCESS) then
      call ESMF_LogWrite("CATChem: Failed to initialize - " // trim(errmsg), &
        ESMF_LOGMSG_ERROR, rc=rc)
      rc = ESMF_FAILURE
      return
    end if

    ! Clean up
    call ESMF_ConfigDestroy(config_esmf, rc=rc)

    ! Log successful completion
    call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeP2

  !> \brief Model advance routine - Execute one time step of chemistry calculations
  !!
  !! This routine is called at each time step to advance the chemistry model.
  !! It handles the complete workflow of importing meteorological data,
  !! running chemistry calculations, and exporting the computed results.
  !!
  !! \param[inout] model NUOPC model component
  !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
  !!
  !! \details
  !! This routine performs the following operations for each time step:
  !! - Retrieves the component clock and current simulation time
  !! - Calculates the time step duration for chemistry integration
  !! - Imports meteorological fields from other model components
  !! - Transforms NUOPC field data into CATChem state format
  !! - Executes all enabled chemistry processes (if do_chemistry=true):
  !!   - Dust emission and transport
  !!   - Sea salt emission and transport
  !!   - Gas-phase and aerosol chemistry
  !!   - Dry and wet deposition processes
  !! - Transforms computed results back to NUOPC field format
  !! - Exports chemistry fields to other model components
  !! - Updates diagnostic outputs and logging
  !!
  !! The routine handles both sequential and parallel execution,
  !! with appropriate logging and error checking throughout.
  !!
  !! \note This routine is called repeatedly during model integration
  !!       and must maintain consistent state between calls
  !!
  !! \ingroup catchem_nuopc_group
  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    type(ESMF_Clock) :: clock
    type(ESMF_Time) :: currTime
    type(ESMF_TimeInterval) :: timeStep
    character(len=*), parameter :: routine = 'ModelAdvance'
    character(len=512) :: errmsg
    integer :: localPet, im, jm, kme
    real(ESMF_KIND_R8) :: dt_seconds

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get states and clock
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get current time and time step
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_TimeIntervalGet(timeStep, s_r8=dt_seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Running chemistry for dt = " // &
        trim(adjustl(real_to_string(dt_seconds))) // " seconds", &
        ESMF_LOGMSG_INFO, rc=rc)
    end if

    ! Get grid dimensions (placeholder - should come from grid)
    im = catchem_states%im
    jm = 1  ! Assume 1D for now
    kme = 1  ! Will be updated based on actual grid

    ! Import meteorological data from other components
    call transform_nuopc_to_catchem(importState, catchem_states, im, kme, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (do_chemistry) then
      ! Run chemistry processes using new interface with current time
      call catchem_nuopc_run(config, catchem_states, dustState, seaSaltState, &
                            dryDepState, dt_seconds, currTime, rc, errmsg)
      if (rc /= CC_SUCCESS) then
        call ESMF_LogWrite("CATChem: Failed to run chemistry - " // trim(errmsg), &
          ESMF_LOGMSG_ERROR, rc=rc)
        rc = ESMF_FAILURE
        return
      end if
    end if

    ! Export results to other components
    call transform_catchem_to_nuopc(exportState, catchem_states, im, kme, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Log successful completion
    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
    end if

  end subroutine ModelAdvance

  !> \brief Execute all chemistry processes for one time step
  !!
  !! This internal routine runs all enabled atmospheric chemistry processes
  !! for each horizontal grid point, including emissions, transport, and
  !! deposition calculations.
  !!
  !! \param[in] dt Time step duration in seconds
  !! \param[out] rc Return code (CC_SUCCESS on success)
  !!
  !! \details
  !! This routine executes the following chemistry processes sequentially:
  !! - Dust emission and transport calculations
  !! - Sea salt emission and transport calculations
  !! - Dry deposition processes for all species
  !! - Additional chemistry processes as configured
  !!
  !! Each process is called for all horizontal grid points in sequence,
  !! with error checking after each major process group.
  !!
  !! \note This routine is called internally by ModelAdvance
  !!
  !! \ingroup catchem_nuopc_group
  subroutine run_chemistry_processes(dt, rc)
    real(ESMF_KIND_R8), intent(in) :: dt
    integer, intent(out) :: rc

    integer :: i

    rc = CC_SUCCESS

    do i = 1, catchem_states%im
      ! Run dust process
      call cc_dust_run(dustState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run seasalt process
      call cc_seasalt_run(seaSaltState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run dry deposition process
      call cc_drydep_run(dryDepState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run other chemistry processes as needed
      call cc_run_process(config, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return
    end do

  end subroutine run_chemistry_processes

  !> \brief Finalize the CATChem model component
  !!
  !! This routine performs cleanup operations and finalizes all CATChem
  !! processes, deallocating memory and closing any open resources.
  !!
  !! \param[inout] model NUOPC model component to finalize
  !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
  !!
  !! \details
  !! This routine performs the following cleanup operations:
  !! - Finalizes all CATChem chemistry and emission processes
  !! - Deallocates memory for all state containers
  !! - Closes any open files or external resources
  !! - Performs diagnostic output if enabled
  !! - Logs completion status and any warnings
  !!
  !! This is a standard NUOPC finalization phase that should be called
  !! when the component is no longer needed or at the end of simulation.
  !!
  !! \note Failure to call this routine may result in memory leaks
  !!       or incomplete output files
  !!
  !! \ingroup catchem_nuopc_group
  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    character(len=*), parameter :: routine = 'ModelFinalize'
    character(len=512) :: errmsg
    integer :: localPet

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Finalize CATChem using the interface
    call catchem_nuopc_finalize(config, catchem_states, dustState, seaSaltState, &
                               dryDepState, rc, errmsg)
    if (rc /= CC_SUCCESS) then
      call ESMF_LogWrite("CATChem: Warning - " // trim(errmsg), &
        ESMF_LOGMSG_WARNING, rc=rc)
    end if

    ! Log successful completion
    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
    end if

    rc = ESMF_SUCCESS

  end subroutine ModelFinalize

  !> \brief Convert real number to string for logging purposes
  !!
  !! This utility function converts a real number to a formatted string
  !! representation suitable for logging and diagnostic messages.
  !!
  !! \param[in] val Real value to convert
  !! \return str String representation of the value
  !!
  !! \details
  !! The function formats the real number with appropriate precision
  !! and removes leading/trailing whitespace for clean output in log
  !! messages and diagnostic information.
  !!
  !! \ingroup catchem_nuopc_group
  function real_to_string(val) result(str)
    real(ESMF_KIND_R8), intent(in) :: val
    character(len=32) :: str

    write(str, '(f0.2)') val
    str = adjustl(str)
  end function real_to_string

end module catchem_nuopc_cap
