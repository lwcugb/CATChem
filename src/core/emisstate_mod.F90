!>
!! \file emisstate_mod.F90
!! \brief Emission Species State
!!
!! \details The Emission Species States contain information about the emitted species and
!!          the emitted fluxes, scale factors, and speciation to concentration mappings.
!!
!! \ingroup core_modules
!!!>
module EmisState_Mod

   ! Uses

   USE Error_Mod
   USE Precision_Mod

   IMPLICIT NONE
   PRIVATE

   !> \brief Emission Species State
   !!
   !! \details The Emission Species States contain information about the emitted species and
   !!          the emitted fluxes, scale factors, and speciation to concentration mappings.
   !!
   !! \ingroup core_modules
   !!
   !! \param name Name of the species
   !! \param long_name Long name of the species
   !! \param units Units of the species
   !! \param nEmisMap Number of Emission mappings per emitted species
   !! \param EmisMapIndex Emission mapping to concentration index
   !! \param EmisScale Scale factor
   !! \param Flux Emission flux
   !!
   !!!>
   TYPE :: EmisSpeciesType
      ! Character
      character(len=10)               :: name             !< Name of the species
      character(len=50)               :: long_name        !< Long name of the species
      character(len=10)               :: units            !< Units of the species
      character(len=100), allocatable :: EmisMapName(:)   !< Emission mapping to concentration species

      ! Integers
      integer              :: nEmisMap         !< Number of Emission mappings per emitted species
      integer, ALLOCATABLE :: EmisMapIndex(:)  !< Emission mapping to concentration index
      integer              :: plumerise        !< Plumerise option (default = 0 for no plumerise) sofiev = 1 briggs = 2 simple linear weighted thickness = 3
      integer              :: EmisLayer        !< Emission layer

      ! Real
      real(fp), ALLOCATABLE :: Scale(:)        !< Scale factor
      real(fp), ALLOCATABLE :: Flux(:)         !< Emission flux
      real(fp)              :: EmisHeight      !< Emission Height [m] - Simple emission height or -1 for PBLH
      real(fp), ALLOCATABLE :: PlmSrcFlx(:)    !< Plumerise source emission flux [kg/m2/s]
      real(fp), ALLOCATABLE :: FRP(:)          !< Fire Radiative Power (W/m^2)
      real(fp), ALLOCATABLE :: STKDM(:)        !< Briggs stack diameter [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKHT(:)        !< Briggs stack thickness [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKTK(:)        !< Briggs stack thickness [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKVE(:)        !< Briggs stack velocity [m/s] (array of all point sources in grid cell)

   END TYPE EmisSpeciesType

   !> \brief Emission Category State
   !!
   !! \details The Emission Category State contains all the Species types for a given emission
   !!          category. It includes the name of the category, the number of emitted species, and
   !!          the emitted species containers.
   !!
   !! \param name Name of the category
   !! \param nSpecies Number of emitted species
   !! \param Species Emitted Species container
   !!
   !! \ingroup core_modules
   !!!>
   TYPE :: EmisCategoryType

      ! Character
      character(len=20)    :: name                         !< Name of the species

      ! Integers
      integer              :: nSpecies                 !< Number of emitted species per process
      integer              :: nPlumerise               !< Number of Species undergoing a Plumerise

      ! Types
      type(EmisSpeciesType), ALLOCATABLE :: Species(:) !< Emitted species container

   END TYPE EmisCategoryType

   !>
   !! \brief Emission State
   !!
   !! \details The Emission State contains all the emission category containers for the simulation and
   !!          the number of emission categories.
   !!
   !! \ingroup core_modules
   !!
   !! \param State Name of the state
   !! \param nCats Number of emission categories
   !! \param Cats Emission category containers
   !!
   !!!>
   TYPE, PUBLIC :: EmisStateType

      ! Character
      CHARACTER(LEN=4) :: State = 'Emis' !< Name of this state
      CHARACTER(LEN=20), ALLOCATABLE :: TotEmisNames(:)         !< Name of the state

      ! Integers
      integer :: nCats         !< Number of emission categories
      integer :: nEmisTotal
      integer :: nEmisCats

      ! Types
      type(EmisSpeciesType), ALLOCATABLE :: TotSpecies(:) !< Emitted species container
      type(EmisCategoryType), ALLOCATABLE :: Cats(:) !< Emission categories container

   END TYPE EmisStateType

CONTAINS

   !> \brief Allocate the emission state
   !!
   !! \details Allocate the emission state.
   !!
   !! \ingroup core_modules
   !!
   !! \param GridState Grid State
   !! \param EmisState Emission State
   !! \param RC Return code
   !!
   !!!>
   subroutine Emis_Allocate(GridState, EmisState, RC)

      ! Uses
      USE GridState_Mod, ONLY : GridStateType

      IMPLICIT NONE

      !Input
      TYPE(GridStateType), INTENT(IN) :: GridState

      ! Input/Output
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState

      ! Output
      INTEGER, INTENT(OUT) :: RC

      ! Locals
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species


      ! Initialize return code
      RC = CC_SUCCESS

      ! Initialize local variables
      ErrMsg = ''
      ThisLoc = ' -> at Emis_Allocate (in core/emisstate_mod.F90)'

      if (EmisState%nCats > 0) then
         do c = 1, EmisState%nCats
            do s = 1, EmisState%Cats(c)%nSpecies

               ALLOCATE(EmisState%Cats(c)%Species(s)%Flux(GridState%number_of_levels), STAT=RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = 'Error allocating "EmisState%Cats%Species%Flux"!'
                  call Handle_Error(ErrMsg, RC, ThisLoc)
                  return
               endif

               ALLOCATE(EmisState%Cats(c)%Species(s)%EmisMapIndex(EmisState%Cats(c)%Species(s)%nEmisMap), STAT=RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = 'Error allocating "EmisState%Cats%Species%EmisMapIndex"!'
                  call Handle_Error(ErrMsg, RC, ThisLoc)
                  return
               endif

            end do
         end do
      end if

   end subroutine Emis_Allocate

   !> \brief Allocate the total emitted species
   !!
   !! \details Allocate the total emitted species
   !!
   !! \ingroup core_modules
   !!
   !! \param GridState Grid State
   !! \param EmisState Emission State
   !! \param RC Return code
   !!
   !!!>
   subroutine TotEmisSpecies_Allocate(GridState, EmisState, RC)

      USE GridState_Mod, ONLY : GridStateType

      IMPLICIT NONE

      TYPE(GridStateType), INTENT(IN) :: GridState

      TYPE(EmisStateType), INTENT(INOUT) :: EmisState

      INTEGER, INTENT(OUT) :: RC

      integer :: s ! Loop counter for emitted species
      integer :: c ! Loop counter for emission Cats
      integer :: t ! loop counter for total emitted species

      integer :: nT ! total number of emitted species counter
      integer :: nS ! number of the same emitted species in different categories

      logical :: IsIn

      character(len=255) :: ErrMsg
      character(len=255) :: ThisLoc

      character(len=10) ::currName

      ! Initialize return code
      RC = CC_SUCCESS

      ! Initialize local variables
      ErrMsg = ''
      ThisLoc = ' -> at TotEmisSpecies_Allocate (in core/emisstate_mod.F90)'

      ALLOCATE(EmisState%TotEmisNames(0), STAT=RC)
      if (EmisState%nCats > 0) then
         EmisState%nEmisTotal = 0
         do c = 1, EmisState%nCats
            do s = 1, EmisState%Cats(c)%nSpecies
               currName = EmisState%Cats(c)%Species(s)%name
               if ( ANY( EmisState%TotEmisNames == TRIM(currName) ) .eqv. .False. ) THEN
                  EmisState%nEmisTotal = EmisState%nEmisTotal + 1
                  EmisState%TotEmisNames = [EmisState%TotEmisNames, TRIM(currName)]
               endif
            end do
         end do

         ALLOCATE(EmisState%TotSpecies(EmisState%nEmisTotal), STAT=RC)
         IF (RC /= CC_SUCCESS) THEN
            ErrMsg = 'Error allocating "EmisState%TotSpecies"!'
            call Handle_Error(ErrMsg, RC, ThisLoc)
            return
         ENDIF

         do t = 1, EmisState%nEmisTotal
            ! Fill Total Species fluxes and properties
            EmisState%TotSpecies(t)%name = EmisState%TotEmisNames(s)
            EmisState%TotSpecies(t)%Flux = 0. ! set all fluxes to zero
            currName = EmisState%TotEmisNames(s)
            do c = 1, EmisState%nCats
               do s = 1, EmisState%Cats(c)%nSpecies
                  if (currName == EmisState%Cats(c)%Species(s)%name) then
                     EmisState%TotSpecies(t)%Flux(:) = EmisState%TotSpecies(t)%Flux(:) +  &
                        EmisState%Cats(c)%Species(s)%Flux(:)
                     EmisState%TotSpecies(t)%units = EmisState%Cats(c)%Species(s)%units
                     EmisState%TotSpecies(t)%long_name = EmisState%Cats(c)%Species(s)%long_name
                  endif
               end do
            end do
         end do
      endif

   end subroutine TotEmisSpecies_Allocate

   !> \brief Find the index of the mapped species
   !!
   !! \details Loops through the emitted species and finds the index to all of the mapped emitted species to the chemical species
   !!
   !! \ingroup core_modules
   !!
   !! \param EmisState Emission State
   !! \param ChemState Chemical State
   !! \param RC Return code
   !!
   !!!>
   subroutine Emis_Find_Chem_Map_Index(EmisState, ChemState, RC)

      USE ChemState_Mod, ONLY : ChemStateType, FindSpecByName

      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      TYPE(ChemStateType), INTENT(INOUT) :: ChemState
      INTEGER, INTENT(OUT) :: RC

      ! Local
      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species
      integer :: n ! Loop counter for mapped species

      integer :: index ! index of the chemical species species

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      ThisLoc = ' -> at Emis_find_chem_map_indexs (in core/emisstate_mod.F90)'

      do c = 1, EmisState%nCats
         do s = 1, EmisState%Cats(c)%nSpecies
            do n = 1, EmisState%Cats(c)%Species(s)%nEmisMap
               call FindSpecByName(ChemState, EmisState%Cats(c)%Species(s)%EmisMapName(n), index, RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = 'Error in find_species_by_name'
                  call Handle_Error(ErrMsg, RC, ThisLoc)
                  return
               endif
               EmisState%Cats(c)%Species(s)%EmisMapIndex(n) = index
            enddo
         enddo
      enddo

   end subroutine Emis_find_chem_map_index

   !> \brief Apply emissions to the chemical state
   !!
   !! \details TODO: Apply emissions to the chemical state
   !! \ingroup core_modules
   !!
   !! \param EmisState Emission State
   !! \param ChemState Chemical State
   !! \param RC Return code
   !!
   !!!>
   subroutine Apply_Emis_to_Chem(Species, MetState, ChemState, RC)
      USE ChemState_Mod, ONLY : ChemStateType
      USE MetState_Mod, ONLY : MetStateType

      TYPE(EmisSpeciesType), INTENT(IN)    :: Species
      TYPE(MetStateType),    INTENT(IN)    :: MetState
      TYPE(ChemStateType),   INTENT(INOUT) :: ChemState

      INTEGER, INTENT(OUT) :: RC

      ! Local
      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      ThisLoc = ' -> at Apply_Emis_to_Chem (in core/emisstate_mod.F90)'

      ! TODO: implement this

   end subroutine Apply_Emis_to_Chem

   subroutine EmisState_CleanUp(EmisState, RC)

      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      INTEGER, INTENT(OUT) :: RC

      character(len=255) :: ErrMsg
      character(len=255) :: ThisLoc

      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species

      ! Initialize return code
      RC = CC_SUCCESS

      ! Initialize local variables
      ErrMsg = ''
      ThisLoc = ' -> at EmisState_CleanUp (in core/emisstate_mod.F90)'

      ! Deallocate total variables
      if (allocated(EmisState%TotEmisNames)) deallocate(EmisState%TotEmisNames)
      do c = 1, EmisState%nEmisTotal
         if (allocated(EmisState%TotSpecies(c)%Flux)) deallocate(EmisState%TotSpecies(c)%Flux)
      end do

      ! Deallocate emission variables in each category
      do c = 1, EmisState%nCats
         do s = 1, EmisState%Cats(c)%nSpecies
            if (allocated(EmisState%Cats(c)%Species(s)%Flux)) &
               deallocate(EmisState%Cats(c)%Species(s)%Flux)
            if (allocated(EmisState%Cats(c)%Species(s)%EmisMapIndex)) &
               deallocate(EmisState%Cats(c)%Species(s)%EmisMapIndex)
         end do
      end do
      if (allocated(EmisState%Cats)) deallocate(EmisState%Cats)

   end subroutine EmisState_CleanUp

END MODULE EmisState_Mod
