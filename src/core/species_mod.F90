!> \file species_mod.F90
!! \brief Species definition and management for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2023
!!
!! This module contains the SpeciesType derived type and related routines
!! for managing chemical species in the CATChem atmospheric chemistry model.
!! It handles species properties, concentrations, and classification.
!!

module species_mod

   use precision_mod
   implicit none

   !> Derived type for chemical species
   !!
   !! This type contains all properties and data for a chemical species
   !! in the CATChem atmospheric chemistry model, including physical
   !! properties, classification flags, and concentration data.
   !!
   !! @param long_name Long descriptive name for species (NetCDF attribute)
   !! @param short_name Short identifier name for species
   !! @param description Detailed description of the species
   !! @param is_gas Logical flag: true if species is gaseous
   !! @param is_aerosol Logical flag: true if species is an aerosol
   !! @param is_tracer Logical flag: true if species is a passive tracer
   !! @param is_advected Logical flag: true if species undergoes advection
   !! @param is_drydep Logical flag: true if species undergoes dry deposition
   !! @param is_photolysis Logical flag: true if species undergoes photolysis
   !! @param is_gocart_aero Logical flag: true if species is a GOCART aerosol
   !! @param is_dust Logical flag: true if species is dust
   !! @param is_seasalt Logical flag: true if species is sea salt
   !! @param mw_g Gaseous molecular weight [g/mol]
   !! @param density Particle density [kg/m³]
   !! @param radius Mean molecular diameter [m]
   !! @param lower_radius Lower radius bound [m]
   !! @param upper_radius Upper radius bound [m]
   !! @param viscosity Kinematic viscosity [m²/s]
   !! @param BackgroundVV Background concentration [v/v]
   !! @param species_index Index in species array
   !! @param drydep_index Index in dry deposition array
   !! @param photolysis_index Index in photolysis array
   !! @param gocart_aero_index Index in GOCART aerosol array
   !! @param conc Species concentration [v/v] or [kg/kg]
   type, public :: SpeciesType

      ! Names
      character(len=30) :: long_name   !< Long name for species used for NetCDF attribute "long_name"
      character(len=30) :: short_name  !< Short name for species
      character(len=50) :: description !< Description of species

      ! Logical switches
      logical :: is_gas               !< If true, species is a gas and not an aerosol
      logical :: is_aerosol           !< If true, species is aerosol and not a gas
      logical :: is_tracer            !< If true, species is a tracer and not an aerosol or gas that undergoes chemistry or photolysis
      logical :: is_advected          !< If true, species is advected
      logical :: is_drydep            !< If true, species undergoes dry deposition
      logical :: is_photolysis        !< If true, species undergoes photolysis
      logical :: is_gocart_aero       !< If true, species is a GOCART aerosol species
      logical :: is_dust              !< If true, species is dust
      logical :: is_seasalt           !< If true, species is sea salt

      ! Numerical properties
      real(kind=fp) :: mw_g                 !< Gaseous molecular weight [g/mol]
      real(kind=fp) :: density              !< Particle density [kg/m³]
      real(kind=fp) :: radius               !< Mean molecular diameter [m]
      real(kind=fp) :: lower_radius         !< Lower radius [m]
      real(kind=fp) :: upper_radius         !< Upper radius [m]
      real(kind=fp) :: viscosity            !< Kinematic viscosity [m²/s]

      ! Default background concentration
      real(kind=fp) :: BackgroundVV        !< Background concentration [v/v]

      ! Indices
      integer :: species_index        !< Species index in species array
      integer :: drydep_index         !< Dry deposition index in drydep array
      integer :: photolysis_index     !< Photolysis index in photolysis array
      integer :: gocart_aero_index    !< GOCART aerosol index in gocart_aero array

      ! Concentration
      real(kind=fp), ALLOCATABLE :: conc(:)             !< Species concentration [v/v] or [kg/kg]

   end type SpeciesType

   !
   ! !DEFINED PARAMETERS:
   !
   !=========================================================================
   ! Missing species concentration value if not in restart file and special
   ! background value not defined
   !=========================================================================
   REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp !< Missing species concentration value

contains

   !> Initialize a species object
   !!
   !! This subroutine initializes a SpeciesType object with basic properties.
   !!
   !! @param Species_State The species object to initialize
   !! @param species_name Short name identifier for the species
   !! @param atomic_num Atomic number or molecular weight
   subroutine init(Species_State, species_name, atomic_num)
      type(SpeciesType), intent(inout) :: Species_State
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: atomic_num

      Species_State%short_name = species_name
      Species_State%mw_g = atomic_num
   end subroutine init

   ! function get_name(this) result(species_name)
   !    character(len=30) :: species_name

   !    species_name = this%short_name
   ! end function get_name

   ! function get_atomic_number(this) result(atomic_num)
   !    class(Species), intent(in) :: this
   !    integer :: atomic_num

   !    atomic_num = this%atomic_number
   ! end function get_atomic_number

end module species_mod
