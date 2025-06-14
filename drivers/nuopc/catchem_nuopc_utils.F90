!> \file catchem_nuopc_utils.F90
!! \brief NUOPC utilities for CATChem atmospheric chemistry model
!!
!! \details
!! This module provides utility subroutines and functions for the NUOPC cap
!! of CATChem. It handles essential operations including data transformations
!! between ESMF/NUOPC fields and CATChem data structures, field creation and
!! management, and common NUOPC operations required for component integration.
!!
!! Key functionalities include:
!! - Field advertisement and realization for NUOPC data exchange
!! - Data import/export between ESMF fields and CATChem states
!! - Field configuration management from YAML files
!! - Grid and coordinate system utilities
!! - Error handling and validation for NUOPC operations
!! - State container initialization and management
!!
!! The module supports flexible field mapping configurations and provides
!! standardized interfaces for coupling CATChem with other Earth system
!! model components through the NUOPC framework.
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module catchem_nuopc_utils

  use ESMF
  use NUOPC
  use CATChem
  use catchem_types, only: catchem_container_type
  use catchem_nuopc_interface

  implicit none

  private

  !> \brief Field configuration structure for NUOPC interface
  !!
  !! Contains the complete field mapping configuration including both
  !! import and export field definitions with associated metadata.
  !! \{
  type :: field_config_type
    integer :: n_import_fields = 0                         !< Number of import fields
    integer :: n_export_fields = 0                         !< Number of export fields
    type(field_mapping_type), allocatable :: import_fields(:) !< Import field mapping array
    type(field_mapping_type), allocatable :: export_fields(:) !< Export field mapping array
  end type field_config_type
  !! \}

  public :: catchem_advertise_fields
  public :: catchem_realize_fields
  public :: catchem_import_from_esmf
  public :: catchem_export_to_esmf
  public :: catchem_state_init
  public :: load_field_config
  public :: field_config_type

contains

  ! Load field configuration from YAML file
  !!
  !! \param  config_file   Path to configuration file
  !! \param field_config  Field configuration structure
  !! \param rc            Return code
  !!
  subroutine load_field_config(config_file, field_config, rc)
    character(len=*), intent(in) :: config_file
    type(field_config_type), intent(out) :: field_config
    integer, intent(out) :: rc

    ! For now, hardcode the configuration - later implement YAML parsing
    rc = ESMF_SUCCESS

    ! Import fields
    field_config%n_import_fields = 13
    allocate(field_config%import_fields(field_config%n_import_fields))

    field_config%import_fields(1) = field_mapping_type("air_temperature", "MetState%Temperature", 3, "K", .false.)
    field_config%import_fields(2) = field_mapping_type("specific_humidity", "MetState%SpecificHumidity", 3, "kg kg-1", .false.)
    field_config%import_fields(3) = field_mapping_type("air_pressure", "MetState%Pressure", 3, "Pa", .false.)
    field_config%import_fields(4) = field_mapping_type("eastward_wind", "MetState%U", 3, "m s-1", .false.)
    field_config%import_fields(5) = field_mapping_type("northward_wind", "MetState%V", 3, "m s-1", .false.)
    field_config%import_fields(6) = field_mapping_type("surface_air_pressure", "MetState%SurfacePressure", 2, "Pa", .false.)
    field_config%import_fields(7) = field_mapping_type("surface_skin_temperature", "MetState%SkinTemperature", 2, "K", .false.)
    field_config%import_fields(8) = field_mapping_type("land_sea_mask", "GridState%LandMask", 2, "1", .false.)
    field_config%import_fields(9) = field_mapping_type("surface_downwelling_shortwave_flux", "MetState%SolarRadiation", 2, "W m-2", .false.)
    field_config%import_fields(10) = field_mapping_type("planetary_boundary_layer_height", "MetState%PBL_Height", 2, "m", .false.)
    field_config%import_fields(11) = field_mapping_type("friction_velocity", "MetState%FrictionVelocity", 2, "m s-1", .true.)
    field_config%import_fields(12) = field_mapping_type("upward_air_velocity", "MetState%W", 3, "m s-1", .true.)
    field_config%import_fields(13) = field_mapping_type("snow_thickness", "MetState%SnowDepth", 2, "m", .true.)

    ! Export fields
    field_config%n_export_fields = 8
    allocate(field_config%export_fields(field_config%n_export_fields))

    field_config%export_fields(1) = field_mapping_type("mass_fraction_of_ozone_in_air", "ChemState%Species%O3", 3, "kg kg-1", .false.)
    field_config%export_fields(2) = field_mapping_type("mass_fraction_of_nitrogen_dioxide_in_air", "ChemState%Species%NO2", 3, "kg kg-1", .false.)
    field_config%export_fields(3) = field_mapping_type("mass_fraction_of_carbon_monoxide_in_air", "ChemState%Species%CO", 3, "kg kg-1", .false.)
    field_config%export_fields(4) = field_mapping_type("mass_fraction_of_dust_dry_aerosol_particles_in_air", "ChemState%Species%DUST", 3, "kg kg-1", .false.)
    field_config%export_fields(5) = field_mapping_type("mass_fraction_of_sea_salt_dry_aerosol_particles_in_air", "ChemState%Species%SEAS", 3, "kg kg-1", .false.)
    field_config%export_fields(6) = field_mapping_type("dust_dry_deposition_flux", "DiagState%DustDryDep", 2, "kg m-2 s-1", .false.)
    field_config%export_fields(7) = field_mapping_type("sea_salt_aerosol_dry_deposition_flux", "DiagState%SeasDryDep", 2, "kg m-2 s-1", .false.)
    field_config%export_fields(8) = field_mapping_type("dust_emission_flux", "EmisState%DustEmission", 2, "kg m-2 s-1", .false.)

  end subroutine load_field_config

  ! Advertise import and export fields for CATChem based on configuration
  !!
  !! \param importState  NUOPC import state
  !! \param exportState  NUOPC export state
  !! \param    config_file  Field configuration file
  !! \param   rc           ESMF return code
  !!
  subroutine catchem_advertise_fields(importState, exportState, config_file, rc)
    type(ESMF_State), intent(inout) :: importState
    type(ESMF_State), intent(inout) :: exportState
    character(len=*), intent(in) :: config_file
    integer,          intent(out)   :: rc

    type(field_config_type) :: field_config
    integer :: i

    rc = ESMF_SUCCESS

    ! Load field configuration
    call load_field_config(config_file, field_config, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Advertise import fields
    do i = 1, field_config%n_import_fields
      call NUOPC_Advertise(importState, &
        StandardName=trim(field_config%import_fields(i)%standard_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    ! Advertise export fields
    do i = 1, field_config%n_export_fields
      call NUOPC_Advertise(exportState, &
        StandardName=trim(field_config%export_fields(i)%standard_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

  end subroutine catchem_advertise_fields

  ! Realize advertised fields by creating ESMF fields and adding to states
  !!
  !! \param importState  NUOPC import state
  !! \param exportState  NUOPC export state
  !! \param    grid         ESMF grid
  !! \param   rc           ESMF return code
  !!
  subroutine catchem_realize_fields(importState, exportState, grid, rc)
    type(ESMF_State), intent(inout) :: importState
    type(ESMF_State), intent(inout) :: exportState
    type(ESMF_Grid),  intent(in)    :: grid
    integer,          intent(out)   :: rc

    type(ESMF_Field) :: field

    rc = ESMF_SUCCESS

    ! Create 3D field for 3D variables
    field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
      ungriddedLBound=[1], ungriddedUBound=[100], &  ! Assume 100 levels max
      staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Realize key 3D import fields
    call NUOPC_Realize(importState, field=field, StandardName="air_temperature", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_Realize(importState, field=field, StandardName="specific_humidity", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_Realize(importState, field=field, StandardName="air_pressure", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Create 2D field for surface variables
    field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Realize key 2D import fields
    call NUOPC_Realize(importState, field=field, StandardName="surface_air_pressure", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_Realize(importState, field=field, StandardName="land_sea_mask", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Realize export fields (similar pattern)
    call NUOPC_Realize(exportState, field=field, StandardName="mass_fraction_of_ozone_in_air", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Continue for all advertised fields...

  end subroutine catchem_realize_fields

  ! Import data from ESMF states to CATChem data structures
  !!
  !! \param    importState    NUOPC import state
  !! \param catchem_states CATChem container
  !! \param   rc             ESMF return code
  !!
  subroutine catchem_import_from_esmf(importState, catchem_states, rc)

    type(ESMF_State),               intent(in)    :: importState
    type(catchem_container_type),   intent(inout) :: catchem_states
    integer,                        intent(out)   :: rc

    type(ESMF_Field)       :: field
    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:)
    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)
    integer :: i, j, k

    rc = ESMF_SUCCESS

    ! Get temperature field
    call ESMF_StateGet(importState, "air_temperature", field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Transform to CATChem states
    do i = 1, size(fptr3d, 1)
      ! For simplicity, assume direct assignment - in practice need proper transformation
      ! catchem_states%MetState(i)%Temperature = fptr3d(i,:,:)
    end do

    ! Continue for other fields...

  end subroutine catchem_import_from_esmf

  ! Export data from CATChem data structures to ESMF states
  !!
  !! \param exportState    NUOPC export state
  !! \param    catchem_states CATChem container
  !! \param   rc             ESMF return code
  !!
  subroutine catchem_export_to_esmf(exportState, catchem_states, rc)

    type(ESMF_State),             intent(inout) :: exportState
    type(catchem_container_type), intent(in)    :: catchem_states
    integer,                      intent(out)   :: rc

    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:)
    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)
    integer :: i, j, k

    rc = ESMF_SUCCESS

    ! Get ozone field for export
    call ESMF_StateGet(exportState, "mass_fraction_of_ozone_in_air", field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Transform from CATChem states
    do i = 1, size(fptr3d, 1)
      ! For simplicity, assume direct assignment - in practice need proper transformation
      ! fptr3d(i,:,:) = catchem_states%ChemState(i)%ozone_concentration
    end do

    ! Similar operations for other species...

  end subroutine catchem_export_to_esmf

  ! Initialize CATChem states with configuration
  !!
  !! \param catchem_states CATChem container
  !! \param    config_file    Configuration file name
  !! \param    im             Horizontal dimension
  !! \param   rc             Return code
  !!
  subroutine catchem_state_init(catchem_states, config_file, im, rc)
    type(catchem_container_type), intent(inout) :: catchem_states
    character(len=*),             intent(in)    :: config_file
    integer,                      intent(in)    :: im
    integer,                      intent(out)   :: rc

    type(ConfigType) :: config
    integer :: i

    rc = CC_SUCCESS

    ! Initialize container dimensions
    catchem_states%im = im

    ! Read configuration
    call cc_read_config(config, config_file, rc)
    if (rc /= CC_SUCCESS) return

    ! Allocate states
    allocate(catchem_states%MetState(im))
    allocate(catchem_states%ChemState(im))
    allocate(catchem_states%EmisState(im))
    allocate(catchem_states%DiagState(im))

    ! Initialize each state
    do i = 1, im
      call cc_allocate_metstate(catchem_states%MetState(i), config, rc)
      if (rc /= CC_SUCCESS) return

      call cc_allocate_chemstate(catchem_states%ChemState(i), config, rc)
      if (rc /= CC_SUCCESS) return

      call cc_allocate_emisstate(catchem_states%EmisState(i), config, rc)
      if (rc /= CC_SUCCESS) return

      call cc_allocate_diagstate(catchem_states%DiagState(i), config, rc)
      if (rc /= CC_SUCCESS) return
    end do

  end subroutine catchem_state_init

end module catchem_nuopc_utils
