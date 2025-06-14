!> \file catchem_nuopc_interface.F90
!! \brief NUOPC interface for CATChem atmospheric chemistry model
!!
!! \details
!! This module provides the core interface routines for the CATChem NUOPC cap,
!! handling data transformations between ESMF/NUOPC fields and CATChem states.
!! It follows similar patterns as the CCPP interface but is adapted for NUOPC
!! framework requirements including ESMF grids, fields, and parallel operations.
!!
!! Key functionalities include:
!! - CATChem model initialization and finalization within NUOPC framework
!! - Data transformation between ESMF fields and CATChem state objects
!! - Field mapping configuration and management
!! - Integration with CF-compliant input and NetCDF output systems
!! - Support for flexible grid configurations and parallel decomposition
!! - Error handling and logging for NUOPC applications
!!
!! The interface supports both sequential and parallel execution modes
!! and provides standardized data exchange capabilities for coupling
!! with other Earth system model components.
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module catchem_nuopc_interface

  use ESMF
  use NUOPC
  use CATChem
  use catchem_types, only: catchem_container_type
  use catchem_nuopc_cf_input
  use catchem_nuopc_netcdf_out
  use machine, only: kind_phys

  implicit none

  private

  public :: catchem_nuopc_init
  public :: catchem_nuopc_run
  public :: catchem_nuopc_finalize
  public :: transform_nuopc_to_catchem
  public :: transform_catchem_to_nuopc
  public :: load_field_config

  !> \brief Field mapping configuration structure
  !!
  !! Defines the mapping between NUOPC standard names and CATChem variables,
  !! including metadata for proper data transformation and validation.
  !! \{
  type :: field_mapping_type
    character(len=128) :: standard_name !< NUOPC/CF standard field name
    character(len=128) :: catchem_var   !< Corresponding CATChem variable path
    integer :: dimensions               !< Number of spatial dimensions (2D/3D)
    character(len=64) :: units          !< Physical units for conversion
    logical :: optional = .false.       !< Whether field is required or optional
  end type field_mapping_type
  !! \}

  !> \brief Module-level field configuration storage
  !!
  !! These variables maintain the field mapping configuration loaded
  !! from external YAML files and used throughout the interface.
  !! \{
  type(QFYAML_t) :: field_config_yaml                    !< YAML configuration parser
  type(field_mapping_type), allocatable :: import_fields(:) !< Import field mapping array
  type(field_mapping_type), allocatable :: export_fields(:) !< Export field mapping array
  integer :: n_import_fields = 0                         !< Number of import fields
  integer :: n_export_fields = 0                         !< Number of export fields
  !! \}

contains

  !> Initialize CATChem model for NUOPC interface
  !!
  !! This routine performs comprehensive initialization of the CATChem model
  !! within the NUOPC framework, including configuration loading, memory
  !! allocation, and setup of I/O systems for external data and diagnostics.
  !!
  !! @param config CATChem configuration object to initialize
  !! @param catchem_states Main container for all CATChem state variables
  !! @param dustState Dust emission and transport state object
  !! @param seaSaltState Sea salt emission and transport state object
  !! @param dryDepState Dry deposition process state object
  !! @param im Number of horizontal grid points
  !! @param config_file Path to CATChem configuration file
  !! @param grid ESMF grid for regridding and I/O operations
  !! @param errflg Error flag (CC_SUCCESS on success)
  !! @param errmsg Error message string if errflg indicates failure
  !!
  !! This routine performs the following initialization steps:
  !! - Reads and validates the CATChem configuration file
  !! - Initializes all CATChem state objects and process modules
  !! - Allocates memory for MetState, EmisState, ChemState, and DiagState arrays
  !! - Sets up field mapping configuration for NUOPC data exchange
  !! - Initializes CF-compliant input system for external data
  !! - Initializes NetCDF output system for diagnostic data
  !! - Validates grid compatibility and spatial dimensions
  !!
  !! The routine is similar to catchem_init in the CCPP interface but includes
  !! additional NUOPC-specific functionality such as ESMF grid integration
  !! and YAML-based field configuration management.
  !!
  !! @note This routine must be called before any CATChem calculations
  !!       and requires valid ESMF grid and configuration files
  !!
  !! @warning Proper error checking should be performed on errflg after calling
  subroutine catchem_nuopc_init(config, catchem_states, dustState, seaSaltState, &
                               dryDepState, im, config_file, grid, errflg, errmsg)

    type(ConfigType), intent(inout) :: config
    type(catchem_container_type), intent(inout) :: catchem_states
    type(DustStateType), intent(inout) :: dustState
    type(SeaSaltStateType), intent(inout) :: seaSaltState
    type(DryDepStateType), intent(inout) :: dryDepState
    integer, intent(in) :: im
    character(len=*), intent(in) :: config_file
    type(ESMF_Grid), intent(in) :: grid
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    ! Local variables
    type(EmisStateType) :: EmisState
    type(ChemStateType) :: ChemState
    type(MetStateType) :: MetState
    type(DiagStateType) :: DiagState
    type(GridStateType) :: GridState
    integer :: i

    ! Initialize
    errflg = CC_SUCCESS
    errmsg = ''

    ! Read configuration
    call cc_read_config(config, GridState, EmisState, ChemState, errflg, config_file)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error reading CATChem configuration file: '//trim(config_file)
      return
    end if

    ! Initialize states and processes
    call cc_init_met(GridState, MetState, errflg)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error initializing MetState'
      return
    end if

    call cc_init_emis(GridState, EmisState, errflg)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error initializing EmisState'
      return
    end if

    call cc_init_chem(GridState, ChemState, errflg)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error initializing ChemState'
      return
    end if

    call cc_init_diag(config, DiagState, ChemState, errflg)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error initializing DiagState'
      return
    end if

    call cc_init_process(config, ChemState, EmisState, dustState, seaSaltState, dryDepState, errflg)
    if (errflg /= CC_SUCCESS) then
      errmsg = 'Error initializing processes'
      return
    end if

    ! Store the grid state (not an array)
    catchem_states%GridState = GridState
    catchem_states%im = im

    ! Allocate and initialize state arrays
    if (.not. allocated(catchem_states%MetState)) then
      allocate(catchem_states%MetState(im), stat=errflg)
      if (errflg /= 0) then
        errmsg = 'Error allocating MetState array'
        errflg = CC_FAILURE
        return
      end if
      do i = 1, im
        catchem_states%MetState(i) = MetState
      end do
    end if

    if (.not. allocated(catchem_states%EmisState)) then
      allocate(catchem_states%EmisState(im), stat=errflg)
      if (errflg /= 0) then
        errmsg = 'Error allocating EmisState array'
        errflg = CC_FAILURE
        return
      end if
      do i = 1, im
        catchem_states%EmisState(i) = EmisState
      end do
    end if

    if (.not. allocated(catchem_states%ChemState)) then
      allocate(catchem_states%ChemState(im), stat=errflg)
      if (errflg /= 0) then
        errmsg = 'Error allocating ChemState array'
        errflg = CC_FAILURE
        return
      end if
      do i = 1, im
        catchem_states%ChemState(i) = ChemState
      end do
    end if

    if (.not. allocated(catchem_states%DiagState)) then
      allocate(catchem_states%DiagState(im), stat=errflg)
      if (errflg /= 0) then
        errmsg = 'Error allocating DiagState array'
        errflg = CC_FAILURE
        return
      end if
      do i = 1, im
        catchem_states%DiagState(i) = DiagState
      end do
    end if

    ! Load field mapping configuration
    call load_field_mappings('catchem_field_config.yml', errflg, errmsg)
    if (errflg /= CC_SUCCESS) return

    ! Initialize CF input system
    call cf_input_init('catchem_input_config.yml', grid, errflg)
    if (errflg /= ESMF_SUCCESS) then
      errmsg = 'Error initializing CF input system'
      errflg = CC_FAILURE
      return
    end if

    ! Initialize NetCDF output system
    call output_diagnostics_init('catchem_output_config.yml', grid, errflg)
    if (errflg /= ESMF_SUCCESS) then
      errmsg = 'Error initializing NetCDF output system'
      errflg = CC_FAILURE
      return
    end if

  end subroutine catchem_nuopc_init

  ! Run CATChem processes for NUOPC interface
  !!
  !! \param config         CATChem configuration
  !! \param catchem_states CATChem container
  !! \param dustState      Dust process state
  !! \param seaSaltState   Sea salt process state
  !! \param dryDepState    Dry deposition process state
  !! \param    dt             Time step (seconds)
  !! \param    current_time   Current model time
  !! \param   errflg         Error flag
  !! \param   errmsg         Error message
  !!
  subroutine catchem_nuopc_run(config, catchem_states, dustState, seaSaltState, &
                              dryDepState, dt, current_time, errflg, errmsg)

    type(ConfigType), intent(inout) :: config
    type(catchem_container_type), intent(inout) :: catchem_states
    type(DustStateType), intent(inout) :: dustState
    type(SeaSaltStateType), intent(inout) :: seaSaltState
    type(DryDepStateType), intent(inout) :: dryDepState
    real(kind_phys), intent(in) :: dt
    type(ESMF_Time), intent(in) :: current_time
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    integer :: i, rc

    errflg = CC_SUCCESS
    errmsg = ''

    ! Update CF input data if needed
    call cf_input_update(current_time, catchem_states, rc)
    if (rc /= ESMF_SUCCESS) then
      errmsg = 'Error updating CF input data'
      errflg = CC_FAILURE
      return
    end if
    errmsg = ''

    ! Run chemistry processes for each grid point
    do i = 1, catchem_states%im

      ! Run dust process
      call cc_dust_run(dustState, catchem_states%GridState, &
                      catchem_states%MetState(i), catchem_states%ChemState(i), &
                      catchem_states%EmisState(i), catchem_states%DiagState(i), &
                      dt, errflg)
      if (errflg /= CC_SUCCESS) then
        write(errmsg, '(A,I0)') 'Error in dust process at grid point ', i
        return
      end if

      ! Run sea salt process
      call cc_seasalt_run(seaSaltState, catchem_states%GridState, &
                         catchem_states%MetState(i), catchem_states%ChemState(i), &
                         catchem_states%EmisState(i), catchem_states%DiagState(i), &
                         dt, errflg)
      if (errflg /= CC_SUCCESS) then
        write(errmsg, '(A,I0)') 'Error in sea salt process at grid point ', i
        return
      end if

      ! Run dry deposition process
      call cc_drydep_run(dryDepState, catchem_states%GridState, &
                        catchem_states%MetState(i), catchem_states%ChemState(i), &
                        catchem_states%EmisState(i), catchem_states%DiagState(i), &
                        dt, errflg)
      if (errflg /= CC_SUCCESS) then
        write(errmsg, '(A,I0)') 'Error in dry deposition process at grid point ', i
        return
      end if

      ! Run main chemistry solver
      call cc_run_process(config, catchem_states%GridState, &
                         catchem_states%MetState(i), catchem_states%ChemState(i), &
                         catchem_states%EmisState(i), catchem_states%DiagState(i), &
                         dt, errflg)
      if (errflg /= CC_SUCCESS) then
        write(errmsg, '(A,I0)') 'Error in chemistry process at grid point ', i
        return
      end if

    end do

    ! Write NetCDF output diagnostics if needed
    call output_diagnostics_write(current_time, catchem_states, rc)
    if (rc /= ESMF_SUCCESS) then
      errmsg = 'Error writing NetCDF output diagnostics'
      errflg = CC_FAILURE
      return
    end if

  end subroutine catchem_nuopc_run

  ! Finalize CATChem for NUOPC interface
  !!
  !! \param catchem_states CATChem container
  !! \param dustState      Dust process state
  !! \param seaSaltState   Sea salt process state
  !! \param dryDepState    Dry deposition process state
  !! \param   errflg         Error flag
  !! \param   errmsg         Error message
  !!
  subroutine catchem_nuopc_finalize(catchem_states, dustState, seaSaltState, &
                                   dryDepState, errflg, errmsg)

    type(catchem_container_type), intent(inout) :: catchem_states
    type(DustStateType), intent(inout) :: dustState
    type(SeaSaltStateType), intent(inout) :: seaSaltState
    type(DryDepStateType), intent(inout) :: dryDepState
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    errflg = CC_SUCCESS
    errmsg = ''

    ! Finalize CF input system
    call cf_input_finalize()

    ! Finalize NetCDF output system
    call output_diagnostics_finalize()

    ! Finalize processes
    call cc_dust_finalize(dustState, errflg)
    call cc_seasalt_finalize(seaSaltState, errflg)
    call cc_drydep_finalize(dryDepState, errflg)

    ! Deallocate state arrays
    if (allocated(catchem_states%MetState)) deallocate(catchem_states%MetState)
    if (allocated(catchem_states%ChemState)) deallocate(catchem_states%ChemState)
    if (allocated(catchem_states%EmisState)) deallocate(catchem_states%EmisState)
    if (allocated(catchem_states%DiagState)) deallocate(catchem_states%DiagState)

    ! Deallocate field mappings
    if (allocated(import_fields)) deallocate(import_fields)
    if (allocated(export_fields)) deallocate(export_fields)

    ! Clean up YAML parser
    call cc_yaml_cleanup(field_config_yaml)

  end subroutine catchem_nuopc_finalize

  ! Transform NUOPC import fields to CATChem states
  !!
  !! \param    importState    NUOPC import state
  !! \param catchem_states CATChem container
  !! \param    im             Horizontal dimension
  !! \param    kme            Vertical dimension
  !! \param   rc             ESMF return code
  !!
  subroutine transform_nuopc_to_catchem(importState, catchem_states, im, kme, rc)

    type(ESMF_State), intent(in) :: importState
    type(catchem_container_type), intent(inout) :: catchem_states
    integer, intent(in) :: im, kme
    integer, intent(out) :: rc

    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:), fptr2d(:,:)
    integer :: i, j, k, n
    logical :: isPresent

    rc = ESMF_SUCCESS

    ! Loop through all import fields and transform to CATChem states
    do n = 1, n_import_fields

      ! Check if field is present in import state
      call ESMF_StateGet(importState, trim(import_fields(n)%standard_name), &
                        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      if (.not. isPresent) then
        if (.not. import_fields(n)%optional) then
          call ESMF_LogWrite("Required field not found: "// &
            trim(import_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
          rc = ESMF_FAILURE
          return
        else
          cycle  ! Skip optional fields that are not present
        end if
      end if

      ! Get the field
      call ESMF_StateGet(importState, trim(import_fields(n)%standard_name), field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! Transform based on field type and dimensions
      call transform_field_to_catchem(field, import_fields(n), catchem_states, im, kme, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

    end do

  end subroutine transform_nuopc_to_catchem

  ! Transform CATChem states to NUOPC export fields
  !!
  !! \param exportState    NUOPC export state
  !! \param    catchem_states CATChem container
  !! \param    im             Horizontal dimension
  !! \param    kme            Vertical dimension
  !! \param   rc             ESMF return code
  !!
  subroutine transform_catchem_to_nuopc(exportState, catchem_states, im, kme, rc)

    type(ESMF_State), intent(inout) :: exportState
    type(catchem_container_type), intent(in) :: catchem_states
    integer, intent(in) :: im, kme
    integer, intent(out) :: rc

    type(ESMF_Field) :: field
    integer :: n
    logical :: isPresent

    rc = ESMF_SUCCESS

    ! Loop through all export fields and transform from CATChem states
    do n = 1, n_export_fields

      ! Check if field is present in export state
      call ESMF_StateGet(exportState, trim(export_fields(n)%standard_name), &
                        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      if (.not. isPresent) then
        if (.not. export_fields(n)%optional) then
          call ESMF_LogWrite("Required export field not found: "// &
            trim(export_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
          rc = ESMF_FAILURE
          return
        else
          cycle  ! Skip optional fields that are not present
        end if
      end if

      ! Get the field
      call ESMF_StateGet(exportState, trim(export_fields(n)%standard_name), field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! Transform from CATChem to field
      call transform_catchem_to_field(catchem_states, export_fields(n), field, im, kme, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

    end do

  end subroutine transform_catchem_to_nuopc

  ! Load field mappings from configuration file using the utilities module
  !!
  !! \param  config_file Configuration file path
  !! \param errflg      Error flag
  !! \param errmsg      Error message
  !!
  subroutine load_field_mappings(config_file, errflg, errmsg)

    character(len=*), intent(in) :: config_file
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    ! Simply call the local load_field_config routine
    call load_field_config(config_file, errflg, errmsg)

  end subroutine load_field_mappings

  ! Transform individual field to CATChem state
  !!
  !! \param    field          ESMF field
  !! \param    field_map      Field mapping information
  !! \param catchem_states CATChem container
  !! \param    im             Horizontal dimension
  !! \param    kme            Vertical dimension
  !! \param   rc             ESMF return code
  !!
  subroutine transform_field_to_catchem(field, field_map, catchem_states, im, kme, rc)

    type(ESMF_Field), intent(in) :: field
    type(field_mapping_type), intent(in) :: field_map
    type(catchem_container_type), intent(inout) :: catchem_states
    integer, intent(in) :: im, kme
    integer, intent(out) :: rc

    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:), fptr2d(:,:)
    integer :: i, j, k

    rc = ESMF_SUCCESS

    ! Transform based on field mapping
    select case (trim(field_map%catchem_var))

      ! 3D meteorological fields
      case ("MetState%Temperature")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%Temperature(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do
      case ("MetState%Pressure")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%Pressure(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do

      case ("MetState%U")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%U(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do

      case ("MetState%V")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%V(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do

      case ("MetState%W")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%W(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do

      ! 2D meteorological fields
      case ("MetState%SurfacePressure")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%SurfacePressure(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%SkinTemperature")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%SkinTemperature(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%SolarRadiation")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%SolarRadiation(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%PBL_Height")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%PBL_Height(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%FrictionVelocity")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%FrictionVelocity(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%SnowDepth")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%SnowDepth(j) = fptr2d(i,j)
          end do
        end do

      ! Grid state fields
      case ("GridState%LandMask")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%GridState%LandMask(i,j) = fptr2d(i,j)
          end do
        end do

      case ("GridState%VegetationFraction")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%GridState%VegetationFraction(i,j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%SpecificHumidity")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              catchem_states%MetState(i)%SpecificHumidity(j,k) = fptr3d(i,j,k)
            end do
          end do
        end do

      case ("MetState%SensibleHeatFlux")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%SensibleHeatFlux(j) = fptr2d(i,j)
          end do
        end do

      case ("MetState%LatentHeatFlux")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            catchem_states%MetState(i)%LatentHeatFlux(j) = fptr2d(i,j)
          end do
        end do

      case default
        call ESMF_LogWrite("Unknown field mapping: " // trim(field_map%catchem_var), &
                          ESMF_LOGMSG_WARNING, rc=rc)

    end select

  end subroutine transform_field_to_catchem

  ! Transform CATChem state to ESMF field
  !!
  !! \param    catchem_states CATChem container
  !! \param    field_map      Field mapping information
  !! \param field          ESMF field
  !! \param    im             Horizontal dimension
  !! \param    kme            Vertical dimension
  !! \param   rc             ESMF return code
  !!
  subroutine transform_catchem_to_field(catchem_states, field_map, field, im, kme, rc)

    type(catchem_container_type), intent(in) :: catchem_states
    type(field_mapping_type), intent(in) :: field_map
    type(ESMF_Field), intent(inout) :: field
    integer, intent(in) :: im, kme
    integer, intent(out) :: rc

    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:), fptr2d(:,:)
    integer :: i, j, k, species_id

    rc = ESMF_SUCCESS

    ! Transform based on field mapping
    select case (trim(field_map%catchem_var))

      ! 3D chemical species fields
      case ("ChemState%Species%O3")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "O3", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%NO2")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "NO2", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%CO")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "CO", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%NO")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "NO", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%SO2")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "SO2", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%DUST")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "DUST", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      case ("ChemState%Species%SEAS")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call cc_find_species_by_name(catchem_states%ChemState(1), "SEAS", species_id)
        if (species_id > 0) then
          do i = 1, im
            do k = 1, kme
              do j = 1, size(fptr3d, 2)
                fptr3d(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
              end do
            end do
          end do
        else
          fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
        end if

      ! 2D diagnostic and emission fields
      case ("DiagState%DustDryDep")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = catchem_states%DiagState(i)%DustDryDep(j)
          end do
        end do

      case ("DiagState%SeasDryDep")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = catchem_states%DiagState(i)%SeasDryDep(j)
          end do
        end do

      case ("DiagState%O3DryDep")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = catchem_states%DiagState(i)%O3DryDep(j)
          end do
        end do

      case ("EmisState%DustEmission")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = catchem_states%EmisState(i)%DustEmission(j)
          end do
        end do

      case ("EmisState%SeasEmission")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = catchem_states%EmisState(i)%SeasEmission(j)
          end do
        end do

      case default
        call ESMF_LogWrite("Unknown export field mapping: " // trim(field_map%catchem_var), &
                          ESMF_LOGMSG_WARNING, rc=rc)

    end select
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              fptr3d(i,j,k) = 0.0_ESMF_KIND_R8  ! Placeholder - need species lookup
            end do
          end do
        end do

      case ("ChemState%Species%CO")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              fptr3d(i,j,k) = 0.0_ESMF_KIND_R8  ! Placeholder - need species lookup
            end do
          end do
        end do

      case ("ChemState%Species%DUST")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              fptr3d(i,j,k) = 0.0_ESMF_KIND_R8  ! Placeholder - need species lookup
            end do
          end do
        end do

      case ("ChemState%Species%SEAS")
        call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do k = 1, kme
            do j = 1, size(fptr3d, 2)
              fptr3d(i,j,k) = 0.0_ESMF_KIND_R8  ! Placeholder - need species lookup
            end do
          end do
        end do

      ! 2D diagnostic and emission fields
      case ("DiagState%DustDryDep")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = 0.0_ESMF_KIND_R8  ! Placeholder - need actual diagnostic data
          end do
        end do

      case ("DiagState%SeasDryDep")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = 0.0_ESMF_KIND_R8  ! Placeholder - need actual diagnostic data
          end do
        end do

      case ("EmisState%DustEmission")
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do i = 1, im
          do j = 1, size(fptr2d, 2)
            fptr2d(i,j) = 0.0_ESMF_KIND_R8  ! Placeholder - need actual emission data
          end do
        end do

      case default
        call ESMF_LogWrite("Unknown export field mapping: "//trim(field_map%catchem_var), &
                          ESMF_LOGMSG_WARNING, rc=rc)

    end select

  end subroutine transform_catchem_to_field

  ! Load field configuration from YAML file
  !!
  !! \param  config_file Configuration file name
  !! \param errflg      Error flag
  !! \param errmsg      Error message
  !!
  subroutine load_field_config(config_file, errflg, errmsg)

    character(len=*), intent(in) :: config_file
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    integer :: rc, i, n_fields
    character(len=128) :: field_key, temp_str
    character(len=512) :: err_msg

    errflg = CC_SUCCESS
    errmsg = ''

    ! Initialize YAML parser
    call cc_yaml_init(field_config_yaml, trim(config_file), rc)
    if (rc /= QFYAML_Success) then
      errflg = CC_FAILURE
      write(errmsg, '(A,A)') 'Failed to initialize YAML parser for: ', trim(config_file)
      return
    end if

    ! Get number of import fields
    call cc_yaml_get(field_config_yaml, 'import_fields%n_fields', n_import_fields, rc)
    if (rc /= QFYAML_Success) then
      errflg = CC_FAILURE
      errmsg = 'Failed to read import_fields%n_fields from config'
      return
    end if

    ! Allocate import fields array
    if (allocated(import_fields)) deallocate(import_fields)
    allocate(import_fields(n_import_fields))

    ! Read import field configurations
    do i = 1, n_import_fields
      write(field_key, '(A,I0,A)') 'import_fields%field_', i, '%standard_name'
      call cc_yaml_get(field_config_yaml, trim(field_key), import_fields(i)%standard_name, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read import field standard_name for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'import_fields%field_', i, '%catchem_var'
      call cc_yaml_get(field_config_yaml, trim(field_key), import_fields(i)%catchem_var, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read import field catchem_var for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'import_fields%field_', i, '%dimensions'
      call cc_yaml_get(field_config_yaml, trim(field_key), import_fields(i)%dimensions, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read import field dimensions for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'import_fields%field_', i, '%units'
      call cc_yaml_get(field_config_yaml, trim(field_key), import_fields(i)%units, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read import field units for field ', i
        errflg = CC_FAILURE
        return
      end if

      ! Optional field (default to false)
      write(field_key, '(A,I0,A)') 'import_fields%field_', i, '%optional'
      call cc_yaml_get(field_config_yaml, trim(field_key), temp_str, rc)
      if (rc == QFYAML_Success) then
        import_fields(i)%optional = (trim(temp_str) == 'true' .or. trim(temp_str) == 'True' .or. trim(temp_str) == 'TRUE')
      else
        import_fields(i)%optional = .false.
      end if
    end do

    ! Get number of export fields
    call cc_yaml_get(field_config_yaml, 'export_fields%n_fields', n_export_fields, rc)
    if (rc /= QFYAML_Success) then
      errflg = CC_FAILURE
      errmsg = 'Failed to read export_fields%n_fields from config'
      return
    end if

    ! Allocate export fields array
    if (allocated(export_fields)) deallocate(export_fields)
    allocate(export_fields(n_export_fields))

    ! Read export field configurations
    do i = 1, n_export_fields
      write(field_key, '(A,I0,A)') 'export_fields%field_', i, '%standard_name'
      call cc_yaml_get(field_config_yaml, trim(field_key), export_fields(i)%standard_name, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read export field standard_name for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'export_fields%field_', i, '%catchem_var'
      call cc_yaml_get(field_config_yaml, trim(field_key), export_fields(i)%catchem_var, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read export field catchem_var for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'export_fields%field_', i, '%dimensions'
      call cc_yaml_get(field_config_yaml, trim(field_key), export_fields(i)%dimensions, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read export field dimensions for field ', i
        errflg = CC_FAILURE
        return
      end if

      write(field_key, '(A,I0,A)') 'export_fields%field_', i, '%units'
      call cc_yaml_get(field_config_yaml, trim(field_key), export_fields(i)%units, rc)
      if (rc /= QFYAML_Success) then
        write(errmsg, '(A,I0)') 'Failed to read export field units for field ', i
        errflg = CC_FAILURE
        return
      end if

      ! Optional field (default to false)
      write(field_key, '(A,I0,A)') 'export_fields%field_', i, '%optional'
      call cc_yaml_get(field_config_yaml, trim(field_key), temp_str, rc)
      if (rc == QFYAML_Success) then
        export_fields(i)%optional = (trim(temp_str) == 'true' .or. trim(temp_str) == 'True' .or. trim(temp_str) == 'TRUE')
      else
        export_fields(i)%optional = .false.
      end if
    end do

  end subroutine load_field_config

end module catchem_nuopc_interface
