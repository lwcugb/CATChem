! \file CCPr_Plumerise_mod.F90
!! \brief Plume rise process module for CATChem
!!
!! This module implements plume rise calculations for point source emissions
!! in the CATChem atmospheric chemistry modeling system. It provides various
!! schemes for calculating the vertical distribution of emissions from
!! elevated point sources such as smokestacks and fires.
!!
!! The module includes multiple plume rise schemes:
!! - Briggs scheme for buoyant plumes
!! - Sofiev scheme for fire emissions
!!
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!! \ingroup catchem_plumerise_process
!!
!! \defgroup catchem_plumerise_process CATChem Plume Rise Process
!! \brief Point source emission vertical distribution
!! \ingroup process_modules
!!
!! This group contains modules and functions for calculating the vertical
!! distribution of point source emissions through plume rise calculations.
!! Supports multiple schemes for different types of emission sources.
!!
!! \par Available Plume Rise Schemes:
!! - **Briggs Scheme**: Standard buoyant plume calculations
!! - **Sofiev Scheme**: Specialized for fire and wildfire emissions
!!!>
MODULE CCPr_Plumerise_mod

   ! USES:
   USE State_Mod
   USE Error_MOD

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CCPR_Plumerise_Run
   PUBLIC :: CCPR_Plumerise_Init
   PUBLIC :: CCPR_Plumerise_Finalize

   TYPE, PUBLIC :: PlumeRiseStateType
      Logical :: activate      ! Activate Plume Rise
   END TYPE PlumeRiseStateType

CONTAINS

   !>
   !! \brief Initialize the CATCHem Plumerise Process
   !!
   !! \param PlumeriseState CATCHem Plumerise state
   !! \param RC Error return code
   !!
   !! \ingroup catchem_plumerise_process
   !!!>
   SUBROUTINE CCPR_Plumerise_Init(PlumeriseState, RC)
      ! USES

      IMPLICIT NONE

      ! INPUT PARAMETERS
      !-----------------
      TYPE(PlumeriseStateType), intent(inout) :: PlumeriseState  ! Nullify Plumerise State During INIT

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      INTEGER,          INTENT(INOUT) :: RC

      ! LOCAL VARIABLES
      !----------------
      ! INTEGER :: k ! Loop Counter

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! Initialize Error handling
      !--------------------------
      RC = CC_SUCCESS
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_Plumerise_INIT (in process/Plumerise/ccpr_Plumerise_mod.F90)'

      ! Initialize
      !-----------
      if (Config%Plumerise_activate) then

         ! Activate Plumerise Process
         !----------------------
         Plumerisestate%Activate = .true.

      else

         PlumeriseState%Activate = .false.

      endif


   END SUBROUTINE CCPR_Plumerise_INIT

   !>
   !! \brief Run the Plumerise scheme
   !!
   !! \param PlumeriseState The PlumeriseState object
   !! \param RC Return code
   !! \param verbose Optional verbose output flag
   !!
   !! \ingroup catchem_plumerise_process
   !!!>
   SUBROUTINE CCPr_Plumerise_Run(PlumeriseState, RC, verbose )

      ! USE
      USE ccpr_scheme_sofiev_mod, ONLY : CCPr_Sofiev_Plmrise
      use CCPR_Scheme_Briggs_Mod, only : CCPr_Briggs_Plumerise

      IMPLICIT NONE

      ! INPUT PARAMETERS
      !-----------------
      ! TYPE(MetStateType),  INTENT(IN) :: MetState       ! MetState Instance
      ! TYPE(GridStateType), INTENT(IN) :: GridState      ! GridState Instance
      ! TYPE(EmisStateType), INTENT(IN) :: EmisState      ! EmisState Instance

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(PlumeriseStateType), INTENT(INOUT) :: PlumeriseState   ! PlumeriseState Instance
      ! TYPE(ChemStateType), INTENT(INOUT) :: ChemState  ! ChemState Instance

      ! OUTPUT PARAMETERS
      !------------------
      INTEGER, INTENT(OUT) :: RC                         ! Return Code

      LOGICAL, OPTIONAL, INTENT(IN) :: verbose

      ! LOCAL VARIABLES
      !----------------
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      Integer :: c ! Emission Categories
      Integer :: s ! Emitted Species
      Integer :: p ! Plumerise source length counter
      Integer :: z ! Vertical counter
      ! integer :: nPlumes ! temporary variable for number of plumes
      Logical :: verbose_

      REAL(fp) :: plmHGT                              ! Plumerise Height [m]
      REAL(fp) :: EFRAC(GridState%number_of_levels)   ! Emission fraction
      real(fp) :: ColEmis(GridState%number_of_levels) ! Column emission rate [kg/m2/s]

      ! Initialize
      !-----------
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_Plumerise_Run (in process/Plumerise/ccpr_Plumerise_mod.F90)'

      if (present(verbose) .eqv. .true.) then
         if (verbose .eqv. .true.) then
            verbose_=.TRUE.
         else
            verbose_=.FALSE.
         endif
      else
         verbose_=.FALSE.
      endif

      if (PlumeriseState%Activate) then
         if (EmisState%nEmisTotalPlumerise == 0) RETURN  ! no plumerise species listed in CATCHem_emission.yml
         cats: do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%nPlumerise /= 0) then  ! in EmisState%Cats(c) plumerise options activated
               species: do s = 1, EmisState%Cats(c)%nSpecies ! loop over emitted species
                  print*, 'Running Plumerise for ', EmisState%Cats(c)%Species(s)%name
                  plume: do p = 1, EmisState%Cats(c)%Species(s)%nPlmSrc ! loop over plume sources
                     select case (EmisState%Cats(c)%Species(s)%plumerise)
                      case (1) ! Sofiev Plumerise

                        call CCPr_Sofiev_Plmrise(MetState%Z,     &
                           MetState%T,                           &
                           MetState%PMID,                        &
                           MetState%PBLH,                        &
                           MetState%PS,                          &
                           EmisState%Cats(c)%Species(s)%frp(p),  &
                           plmHGT,                               &
                           EFRAC,                                &
                           RC)
                        if (RC /= CC_SUCCESS) then
                           errMsg = 'Error in CCPr_Sofiev_Plmrise'
                           CALL CC_Error( errMsg, RC, thisLoc )
                        endif

                        ! Add emission to ColEmis for total Emission in grid cell due to plumerise
                        ! Will Speciate out afterwards to concentration at end
                        do z = 1, GridState%number_of_levels
                           ColEmis(z) = ColEmis(z) + EmisState%Cats(c)%Species(s)%PlmSrcFlx(p) * EFRAC(z)
                        end do

                        ! Add Plume Rise height to species
                        EmisState%Cats(c)%Species(s)%PlmRiseHgt(p) = plmHGT

                      case (2) ! Brigg's Plumerise
                        call CCPr_Briggs_Plumerise(MetState%Z,          &
                           MetState%ZMID,                               &
                           MetState%T,                                  &
                           MetState%QV,                                 &
                           MetState%U,                                  &
                           MetState%V,                                  &
                           MetState%PMID,                               &
                           MetState%HFLUX,                              &
                           MetState%PBLH,                               &
                           MetState%USTAR,                              &
                           MetState%T2M,                                &
                           MetState%PS,                                 &
                           SIZE(MetState%T),                            &
                           EmisState%Cats(c)%Species(s)%STKDM(p),   &
                           EmisState%Cats(c)%Species(s)%STKHT(p),   &
                           EmisState%Cats(c)%Species(s)%STKTK(p),   &
                           EmisState%Cats(c)%Species(s)%STKVE(p),   &
                           plmHGT,                                      &
                           EFRAC)
                        if (RC /= CC_SUCCESS) then
                           errMsg = 'Error in CCPr_Briggs_Plumerise'
                           CALL CC_Error( errMsg, RC, thisLoc )
                        endif

                        do z = 1, GridState%number_of_levels
                           ColEmis(z) = ColEmis(z) + EmisState%Cats(c)%Species(s)%PlmSrcFlx(p) * EFRAC(z)
                        end do

                        ! Add Plume Rise height to species
                        EmisState%Cats(c)%Species(s)%PlmRiseHgt(p) = plmHGT

                      case default

                     end select

                     if ((verbose_ .eqv. .true.) .and. (EmisState%Cats(c)%Species(s)%plumerise > 0)) then
                        write(*,*) '------------ PLUMERISE ---------------'
                        write(*,*) 'Category: ', EmisState%Cats(c)%Name
                        write(*,*) 'Species: ', EmisState%Cats(c)%Species(s)%Name
                        write(*,*) 'PlmHGT: ', plmHGT
                        write(*,*) 'PlmSrcFlx: ', EmisState%Cats(c)%Species(s)%PlmSrcFlx(p)
                        if (EmisState%Cats(c)%Species(s)%plumerise == 2) then
                           write(*,*) 'PlumriseOpt: BRIGGS'
                           write(*,*) '  STKDM: ', EmisState%Cats(c)%Species(s)%STKDM(p)
                           write(*,*) '  STKHT: ', EmisState%Cats(c)%Species(s)%STKHT(p)
                           write(*,*) '  STKTK: ', EmisState%Cats(c)%Species(s)%STKTK(p)
                           write(*,*) '  STKVE: ', EmisState%Cats(c)%Species(s)%STKVE(p)
                        else if (EmisState%Cats(c)%Species(s)%plumerise == 1) then
                           write(*,*) 'PlumriseOpt: SOFIEV'
                           write(*,*) '  FRP: ', EmisState%Cats(c)%Species(s)%frp(p)
                        endif
                     endif
                  end do plume

                  ! Add emission to ColEmis to Species total flux in grid cell
                  EmisState%Cats(c)%Species(s)%Flux = ColEmis

               end do species
            end if ! End Plumerise Condition within EmisState%Cats(c)
         end do cats

      endif ! Activate

   END SUBROUTINE CCPr_Plumerise_Run

   !>
   !! \brief Finalize the Plumerise scheme
   !!
   !! \param RC Return code
   !!
   !! \ingroup catchem_plumerise_process
   !!!>
   SUBROUTINE CCPr_Plumerise_Finalize( RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      ! TYPE(PlumeriseStateType), INTENT(INOUT) :: PlumeriseState ! PlumeriseState Instance
      INTEGER, INTENT(OUT) :: RC                       ! Return Code

      ! LOCAL VARIABLES
      !----------------
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      !-----------
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_Plumerise_Finalize (in process/plumerise/CCPr_Plumerise.F90)'

      ! Nothing to do yet

      RETURN

   END SUBROUTINE CCPr_Plumerise_Finalize

END MODULE CCPr_Plumerise_Mod
