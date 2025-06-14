!> \file CCPr_SeaSalt_mod.F90
!! \brief Driver for the CATChem Process: SeaSalt
!! \ingroup catchem_seasalt_process
!!
!! \author CATChem Development Team
!! \date 2024
!!
!! This module provides the main driver for sea salt emission processes in CATChem,
!! including initialization, execution, and finalization of sea salt calculations.
!! The module supports multiple sea salt emission schemes and handles the
!! integration with the atmospheric chemistry system.
!!
MODULE CCPr_SeaSalt_mod

   ! USES:
   USE Precision_Mod
   USE Error_MOD
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod, Only : MetStateType
   USE Config_Opt_Mod, Only : ConfigType
   USE ChemState_Mod, Only : ChemStateType
   USE EmisState_Mod, Only : EmisStateType
   USE CCPr_SeaSalt_Common_Mod, Only : SeaSaltStateType

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: CCPR_SeaSalt_Run
   PUBLIC :: CCPR_SeaSalt_Init
   PUBLIC :: CCPR_SeaSalt_Finalize

CONTAINS

   !>
   !! \brief Initialize the CATChem SeaSalt Process
   !!
   !! Initializes the sea salt aerosol emission process with default or configured parameters.
   !! Sets up sea salt bin properties, size distributions, and links to chemical species.
   !!
   !! \param Config CATChem configuration options
   !! \param SeaSaltState CATChem sea salt state to be initialized
   !! \param ChemState CATChem chemical state
   !! \param EmisState CATChem emission state
   !! \param RC Error return code
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   SUBROUTINE CCPR_SeaSalt_Init( Config, SeaSaltState, ChemState, EmisState, RC)
      ! USES

      IMPLICIT NONE

      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType),    intent(in)    :: Config     ! Config options
      TYPE(SeaSaltStateType), intent(inout) :: SeaSaltState  ! Nullify SeaSalt State During INIT
      TYPE(ChemStateType), intent(in)    :: ChemState  ! Chemical State
      TYPE(EmisStateType), intent(in)    :: EmisState  ! Emission State

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      INTEGER,          INTENT(INOUT) :: RC

      ! LOCAL VARIABLES
      !----------------
      Integer, parameter :: nSeaSaltBinsDefault = 5
      REAL(fp), DIMENSION(nSeaSaltBinsDefault), Parameter :: DefaultSeaSaltDensity  = (/ 2200., &
         2200., &
         2200., &
         2200., &
         2200. /)
      REAL(fp), DIMENSION(nSeaSaltBinsDefault), Parameter :: DefaultEffectiveRadius = (/ 0.079, &
         0.316, &
         1.119, &
         2.818, &
         7.772 /)
      REAL(fp), DIMENSION(nSeaSaltBinsDefault), Parameter :: DefaultLowerBinRadius  = (/ 0.03, &
         0.1,  &
         0.5,  &
         1.5,  &
         5.0  /)
      REAL(fp), DIMENSION(nSeaSaltBinsDefault), Parameter :: DefaultUpperBinRadius  = (/ 0.1, &
         0.5, &
         1.5, &
         5.0, &
         10. /)

      INTEGER :: c, k ! Loop Counter
      INTEGER, DIMENSION(1) :: min_ind
      REAL(fp)  :: radius_temp(10)   ! radius of sea salt bin holder
      LOGICAL   :: mask(10) = .FALSE. ! flag for sorting bins by radius

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! Initialize Error handling
      !--------------------------
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_SeaSalt_INIT (in process/SeaSalt/ccpr_SeaSalt_mod.F90)'

      ! Initialize
      !-----------
      if (Config%seasalt_activate) then

         ! Activate SeaSalt Process
         !-------------------------
         SeaSaltState%Activate = .true.

         ! Set number of seasalt species
         !------------------------------
         SeaSaltState%nSeaSaltSpecies = ChemState%nSpeciesSeaSalt

         ! Set Scheme Options | Default GEOS 2012 scheme
         !----------------------------------------------
         if (Config%seasalt_scheme < 0) then
            SeaSaltState%SeaSaltScaleFactor = 3
         else
            SeaSaltState%SeaSaltScaleFactor = Config%seasalt_scheme
         endif
         SeaSaltState%SchemeOpt = Config%seasalt_scheme

         ! Set Tuning Scale Factor | Default = 1 if not set
         !-------------------------------------------------
         if (Config%seasalt_scalefactor < 0) then
            SeaSaltState%SeaSaltScaleFactor = 1
         else
            SeaSaltState%SeaSaltScaleFactor = Config%seasalt_scalefactor
         endif

         ! Set Weibull Distribution flag following Fan and Toon 2011 | Default = .true.
         SeaSaltState%WeibullFlag = Config%seasalt_weibull

         ! Set Hoppel Correction flag following Fan and Toon 2011 | Default = .true.
         SeaSaltState%HoppelFlag = Config%seasalt_hoppel

         !Find emission caterory index in EmisState for future use
         !--------------------------------------------
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'seasalt') then
               SeaSaltState%CatIndex = c
               exit
            endif
         end do

         if (SeaSaltState%nSeaSaltSpecies == 0) then

            ! Set default bin properties for schemes that need them
            !------------------------------------------------------
            ALLOCATE(SeaSaltState%LowerBinRadius(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%LowerBinRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%LowerBinRadius(k) = DefaultLowerBinRadius(k)
            end do

            ALLOCATE(SeaSaltState%UpperBinRadius(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%UpperBinRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%UpperBinRadius(k) = DefaultUpperBinRadius(k)
            enddo

            ALLOCATE(SeaSaltState%EffectiveRadius(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%EffectiveRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%EffectiveRadius(k) = DefaultEffectiveRadius(k)
            end do

            ALLOCATE(SeaSaltState%SeaSaltDensity(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%SeaSaltDensity', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%SeaSaltDensity(k) = DefaultSeaSaltDensity(k)
            end do

            ALLOCATE(SeaSaltState%EmissionPerSpecies(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('EmissionPerSpecies', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%EmissionPerSpecies(k) = ZERO
            end do

            ALLOCATE(SeaSaltState%NumberEmissionBin(nSeaSaltBinsDefault), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%NumberEmissionBin', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN
            do k = 1, nSeaSaltBinsDefault
               SeaSaltState%NumberEmissionBin(k) = ZERO
            end do

         else

            ! seasalt Aerosols are present in ChmState
            !--------------------------------------
            ALLOCATE(SeaSaltState%LowerBinRadius(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%LowerBinRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ALLOCATE(SeaSaltState%UpperBinRadius(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%UpperBinRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ALLOCATE(SeaSaltState%EffectiveRadius(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%EffectiveRadius', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ALLOCATE(SeaSaltState%SeaSaltDensity(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%SeaSaltDensity', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ALLOCATE(SeaSaltState%EmissionPerSpecies(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%EmissionPerSpecies', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ALLOCATE(SeaSaltState%NumberEmissionBin(SeaSaltState%nSeaSaltSpecies), STAT=RC)
            CALL CC_CheckVar('SeaSaltState%NumberEmissionBin', 0, RC)
            IF (RC /= CC_SUCCESS) RETURN

            ! Set the default values for the SeaSalt species from lower to upper bins
            ! TODO: here we assume every sea salt emission species is mapped to the only one species in concentration
            radius_temp(1:SeaSaltState%nSeaSaltSpecies) = ChemState%ChemSpecies(ChemState%SeaSaltIndex(:))%radius
            mask(1:SeaSaltState%nSeaSaltSpecies) = .TRUE.
            do k = 1, SeaSaltState%nSeaSaltSpecies
               min_ind = MINLOC(radius_temp, mask)  ! Find the index of the minimum radius in the mask
               SeaSaltState%LowerBinRadius(k) = ChemState%ChemSpecies(ChemState%SeaSaltIndex(min_ind(1)))%lower_radius !TODO: keep um, not m
               SeaSaltState%UpperBinRadius(k) = ChemState%ChemSpecies(ChemState%SeaSaltIndex(min_ind(1)))%upper_radius
               SeaSaltState%EffectiveRadius(k) = ChemState%ChemSpecies(ChemState%SeaSaltIndex(min_ind(1)))%radius
               SeaSaltState%SeaSaltDensity(k) = ChemState%ChemSpecies(ChemState%SeaSaltIndex(min_ind(1)))%density
               mask(min_ind) = .FALSE. ! Set the minimum to false so it won't be selected again
               SeaSaltState%EmissionPerSpecies(k) = 0.0_fp ! Initialize to zero
               SeaSaltState%NumberEmissionBin(k) = 0.0_fp ! Initialize to zero
            end do
            SeaSaltState%TotalEmission = 0.0_fp ! Initialize to zero

         endif

      else
         SeaSaltState%TotalEmission = ZERO
         SeaSaltState%TotalNumberEmission = ZERO
         SeaSaltState%SeaSaltScaleFactor = ZERO
         SeaSaltState%SchemeOpt = 3
         SeaSaltState%Activate = .false.

      endif


   END SUBROUTINE CCPR_SeaSalt_INIT

   !>
   !! \brief Run the sea salt emission scheme
   !!
   !! Executes sea salt aerosol emission calculations using the selected scheme
   !! (Gong 1997, Gong 2003, or GEOS-12). Computes sea salt fluxes based on
   !! wind speed and sea surface conditions.
   !!
   !! \param MetState The meteorological state containing atmospheric conditions
   !! \param SeaSaltState The sea salt state containing process parameters
   !! \param EmisState The emission state for storing sea salt emission fluxes
   !! \param RC Return code indicating success or failure
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   SUBROUTINE CCPr_SeaSalt_Run( MetState, SeaSaltState, EmisState, RC )

      ! USE
      USE CCPr_Scheme_Gong03_mod,  ONLY: CCPr_Scheme_Gong03   ! Gong2003 SeaSalt Scheme
      USE CCPr_Scheme_Gong97_mod,  ONLY: CCPr_Scheme_Gong97   ! Gong1997 SeaSalt Scheme
      USE CCPr_Scheme_GEOS12_mod,  ONLY: CCPr_Scheme_GEOS12   ! Gong1997 SeaSalt Scheme

      IMPLICIT NONE

      ! INPUT PARAMETERS
      !-----------------
      TYPE(MetStateType),  INTENT(IN) :: MetState       ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      ! TYPE(DiagStateType), INTENT(INOUT)    :: DiagState   ! DiagState Instance
      TYPE(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState   ! SeaSaltState Instance
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState   ! Emission State
      ! TYPE(ChemStateType), INTENT(INOUT)    :: ChemState  ! ChemState Instance

      ! OUTPUT PARAMETERS
      !------------------
      INTEGER, INTENT(OUT) :: RC                         ! Return Code


      ! LOCAL VARIABLES
      !----------------
      INTEGER :: i ! Loop Counter
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      !-----------
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_SeaSalt_Run (in process/seasalt/ccpr_seasalt_mod.F90)'

      if (SeaSaltState%Activate) then

         ! Run the SeaSalt Scheme
         !--------------------
         if (SeaSaltState%SchemeOpt == 1) then ! Gong2003
            call CCPr_Scheme_Gong03(MetState%FROCEAN,                 &
               MetState%FRSEAICE,                &
               MetState%U10M,                    &
               MetState%V10M,                    &
               MetState%SST,                     &
               SeaSaltState%WeibullFlag,         &
               SeaSaltState%SeaSaltScaleFactor,  &
               SeaSaltState%UpperBinRadius,      &
               SeaSaltState%LowerBinRadius,      &
               SeaSaltState%EffectiveRadius,     &
               SeaSaltState%SeaSaltDensity,      &
               SeaSaltState%EmissionPerSpecies,  &
               SeaSaltState%NumberEmissionBin,   &
               SeaSaltState%TotalEmission,       &
               SeaSaltState%TotalNumberEmission, &
               RC)
            if (RC /= CC_SUCCESS) then
               errMsg = 'Error in CCPr_Scheme_Gong03'
               CALL CC_Error( errMsg, RC, thisLoc )
            endif
         else if (SeaSaltState%SchemeOpt == 2) then ! Gong1997
            ! call CCPr_Scheme_Gong97( MetState, DiagState, SeaSaltState, RC )
            call CCPr_Scheme_Gong97(MetState%FROCEAN,                 &
               MetState%FRSEAICE,                &
               MetState%U10M,                    &
               MetState%V10M,                    &
               MetState%SST,                     &
               SeaSaltState%WeibullFlag,         &
               SeaSaltState%SeaSaltScaleFactor,  &
               SeaSaltState%UpperBinRadius,      &
               SeaSaltState%LowerBinRadius,      &
               SeaSaltState%EffectiveRadius,     &
               SeaSaltState%SeaSaltDensity,      &
               SeaSaltState%EmissionPerSpecies,  &
               SeaSaltState%NumberEmissionBin,   &
               SeaSaltState%TotalEmission,       &
               SeaSaltState%TotalNumberEmission, &
               RC)
            if (RC /= CC_SUCCESS) then
               errMsg = 'Error in CCPr_Scheme_Gong97'
               CALL CC_Error( errMsg, RC, thisLoc )
            endif
         else if (SeaSaltState%SchemeOpt == 3) then ! GEOS2012
            call CCPr_Scheme_GEOS12(MetState%FROCEAN,                 &
               MetState%FRSEAICE,                &
               MetState%USTAR,                   &
               MetState%SST,                     &
               SeaSaltState%SeaSaltScaleFactor,  &
               SeaSaltState%UpperBinRadius,      &
               SeaSaltState%LowerBinRadius,      &
               SeaSaltState%EffectiveRadius,     &
               SeaSaltState%SeaSaltDensity,      &
               SeaSaltState%EmissionPerSpecies,  &
               SeaSaltState%NumberEmissionBin,   &
               SeaSaltState%TotalEmission,       &
               SeaSaltState%TotalNumberEmission, &
               RC)
            if (RC /= CC_SUCCESS) then
               errMsg = 'Error in CCPr_Scheme_GEOS12'
               CALL CC_Error( errMsg, RC, thisLoc )
            endif
         else
            errMsg =  'ERROR: Unknown seasalt scheme option'
            RC = CC_FAILURE
            CALL CC_Error( errMsg, RC, thisLoc )
            return
         endif

         !Fill Emission State.
         do i = 1, EmisState%Cats(SeaSaltState%CatIndex)%nSpecies
            EmisState%Cats(SeaSaltState%CatIndex)%Species(i)%Flux(1) = SeaSaltState%EmissionPerSpecies(i)
         end do

      endif



   END SUBROUTINE CCPr_SeaSalt_Run

   !>
   !! \brief Finalize the sea salt emission process
   !!
   !! Cleans up and deallocates memory used by the sea salt emission process.
   !! Frees arrays and resets state variables.
   !!
   !! \param SeaSaltState The sea salt state to be finalized
   !! \param RC Return code indicating success or failure
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   SUBROUTINE CCPr_SeaSalt_Finalize( SeaSaltState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState ! SeaSaltState Instance
      INTEGER, INTENT(OUT) :: RC                       ! Return Code

      ! LOCAL VARIABLES
      !----------------
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      !-----------
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_SeaSalt_Finalize (in process/seasalt/ccpr_SeaSalt.F90)'

      DEALLOCATE( SeaSaltState%LowerBinRadius, STAT=RC )
      CALL CC_CheckVar('SeaSaltState%LowerBinRadius', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN

      DEALLOCATE( SeaSaltState%UpperBinRadius, STAT=RC )
      CALL CC_CheckVar('SeaSaltState%UpperBinRadius', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN

      DEALLOCATE( SeaSaltState%EffectiveRadius, STAT=RC )
      CALL CC_CheckVar('SeaSaltState%EffectiveRadius', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN

      DEALLOCATE( SeaSaltState%SeaSaltDensity, STAT=RC )
      CALL CC_CheckVar('SeaSaltState%SeaSaltDensity', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN

      DEALLOCATE( SeaSaltState%EmissionPerSpecies, STAT=RC )
      CALL CC_CheckVar('SeaSaltState%EmissionPerSpecies', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN

   END SUBROUTINE CCPr_SeaSalt_Finalize

END MODULE CCPr_SeaSalt_Mod
