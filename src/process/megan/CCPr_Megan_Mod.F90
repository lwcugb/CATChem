!> \brief Driver for the CATCHem Process: Megan
!!
!!
!! \defgroup catchem_megan_process
!!
!! \author Wei Li
!! \date 07/2024
!!!>
MODULE CCPR_Megan_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Mod,    Only : ConfigType
   USE CCPr_Megan_Common_Mod, Only : MeganStateType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_Megan_Init
   PUBLIC :: CCPR_Megan_Run
   PUBLIC :: CCPR_Megan_Final

CONTAINS

   !>
   !! \brief Initialize the CATChem Megan module
   !!
   !! \param Config_Opt       CATCHem configuration options
   !! \param MeganState       CATCHem Megan state
   !! \param ChmState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !!!>
   SUBROUTINE CCPR_Megan_Init( Config, ChemState, MeganState, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigOptType), POINTER    :: Config    ! Module options
      TYPE(ChemStateType), POINTER    :: ChemState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(MeganStateType), POINTER   :: MeganState ! <PROCESS> state
      INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------

      ! Put any local variables here

      !=================================================================
      ! CCPR_Megan_Init begins here!
      !=================================================================
      RC = CC_SUCCESS
      ThisLoc = ' -> at CCPR_Megan_INIT (in process/megan/ccpr_megan_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%megan_activate) then

         ! Activate Process
         !------------------
         MeganState%Activate = .true.

         ! Set number of species
         !----------------------
         MeganState%nMeganSpecies = 21

         ! CO2 inhibition option
         !TODO: what if it is not given in the configuration file properly
         !------------------
         MeganState%CO2Inhib = Config%CO2_Inhib_Opt

         ! Set CO2 concentration (ppm)
         !----------------------------
         if (Config%CO2_conc < 0) then !TODO: DO we give it a negative value if it is missing in config
            MeganState%CO2conc = 390.0_fp
         else
            MeganState%CO2conc = Config%CO2_conc_ppm
         endif

         ! Check GLOBCO2 if CO2 inhibition is turned on (LISOPCO2 = .TRUE.)
         ! GLOBCO2 should be between 150-1250 ppmv. Isoprene response to
         ! CO2 outside this range has no empirical basis.
         if ( MeganState%CO2Inhib ) then
            if ( MeganState%CO2conc <  150.0_fp .or. &
               MeganState%CO2conc > 1250.0_fp     ) then
               RC = CC_FAILURE
               ErrMsg = 'Global CO2 outside valid range of 150-1250 ppmv!'
               call CC_Error( errMsg, RC, thisLoc )
               return
            endif
         endif

         ! Allocate emission species index
         ALLOCATE( MeganState%MeganSpeciesIndex(MeganState%nMeganSpecies) )
         CALL CC_CheckVar('MeganState%MeganSpeciesIndex', 0, RC)  
         IF (RC /= CC_SUCCESS) RETURN
         
         ! Allocate emission speceis names
         ALLOCATE( MeganState%MeganSpeciesName(MeganState%nMeganSpecies) )
         CALL CC_CheckVar('MeganState%MeganSpeciesName', 0, RC)  
         IF (RC /= CC_SUCCESS) RETURN

         ! Allocate emission flux
         ALLOCATE( MeganState%EmissionPerSpecies(MeganState%nMeganSpecies) )
         CALL CC_CheckVar('MeganState%EmissionPerSpecies', 0, RC)  
         IF (RC /= CC_SUCCESS) RETURN

         ! Allocate normalized factor
         ! There should be a different normalization factor for each compound, but
         ! we calculate only 1 normalization factor for all compounds
         ALLOCATE( MeganState%EmisNormFactor(1) )
         CALL CC_CheckVar('MeganState%EmisNormFactor', 0, RC)  
         IF (RC /= CC_SUCCESS) RETURN

         !TODO: emission factor from 7 speceis are read from files and put in MetSate by now
         !      Others are calculated using 'PFT_16', which is also added in MetState
         !      Some met values (last 15 day T average) may need another function and be saved to restart file.

         !TODO: emission species name and ID should read from a namelist. Give them values for now
         MeganState%MeganSpeciesIndex = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/)
         MeganState%MeganSpeciesName(/'ISOP','APIN','BPIN','LIMO','SABI','MYRC','CARE',
                                      'OCIM','OMON','ALD2','MOH', 'EOH', 'MBOX','FAXX',
                                      'AAXX','ACET','PRPE','C2H4','FARN','BCAR','OSQT' /)

      else

         MeganState%Activate = .false.

      endif

   end subroutine CCPR_Megan_Init

   !>
   !! \brief Run the Megan process
   !!
   !! \param [IN] MetState The MetState object
   !! \param [INOUT] DiagState The DiagState object
   !! \param [INOUT] MeganState The MeganState object
   !! \param [INOUT] ChemState The ChemState object
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_Megan_Run( MetState, EmisState, DiagState, MeganState, ChemState, RC )

      ! USE
      USE CCPr_Scheme_Megan_Mod, ONLY: CCPr_Scheme_Megan  ! Megan scheme
      USE EmisState_Mod

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetState_type),  INTENT(IN) :: MetState       ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(EmisStateType), INTENT(INOUT)  :: EmisState   ! Emission Instance
      TYPE(DiagState_type), INTENT(INOUT) :: DiagState   ! DiagState Instance
      TYPE(MeganStateType), INTENT(INOUT) :: MeganState  ! Megan State Instance
      TYPE(ChemState_type), INTENT(INOUT) :: ChemState   ! ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                         ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      integer            :: c, s, cat_megan_idx
      logical, save      :: FIRST = .TRUE.

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_Megan_Run (in process/megan/ccpr_megan_mod.F90)'

      ! Run the Megan Scheme if activated
      !----------------------------------
      if (MeganState%Activate) then
         
         if (FIRST) then
            ! Calculate normalization factor
            ! Not really used now based on Sam's method in which 0.21 is used 
            CALL CALC_NORM_FAC( MetState%D2RAD, MeganState%EmisNormFactor(1), RC )
            if (RC /= CC_SUCCESS) then
               MSG = 'call on CALC_NORM_FAC failed!'
               call CC_Error( MSG, RC , thisLoc)
            endif
            FIRST = .FALSE.
         endif
         
         ! Run the megan Scheme
         ! TODO: 1, if keep current structure (CCpr_Schem_Megan once), need to create a bool array based on EmisState
         !          to decide which species to calculated and then map back to EmisState%Cat%species%flux
         !       2, Or we can put 'CCPr_Scheme_Megan' in a loop based on EmisSate%cat%nSpecies and calculate
         !          the flux we need directly and do not need to map the flux back to EmisState (use this for now)
         !-------------------------
         !find megan caterory index
         !TODO: maybe this can be added to MeganState in case of future usage
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'megan') then
               cat_megan_idx = c
               exit
            endif          
         end do

         do s = 1, EmisState$Cats(cat_megan_idx)%nSpecies

            call CCPr_Scheme_Megan(             &
                  !MeganState%nMeganSpecies,     &
                  !MeganState%MeganSpeciesName,  &   
                  !MeganState%EmissionPerSpecies,&
                  !MeganState%TotalEmission,     &
                  EmisState%Cats(cat_megan_idx)%Species(s)%name, &
                  !TODO: give Flux(1) directly since megan is only at thhe surface
                  !      TotalEmission is not used since it can be calculated through EmisState later
                  EmisState%Cats(cat_megan_idx)%Species(s)%Flux(1), &   
                  MetState%LAI,                 &
                  MetState%PFT_16,              &
                  MetState%PMISOLAI,            &   
                  MetState%Q_DIR_2,             &
                  MetState%Q_DIFF_2,            &   
                  MetState%PARDR_LASTXDAYS,     &
                  MetState%PARDF_LASTXDAYS,     &   
                  MetState%TS,                  &
                  MetState%T_LASTXDAYS,         &
                  MetState%T_LAST24H,           &   
                  MetState%GWETROOT,            &
                  MetState%CO2Inhib,            &
                  MetState%CO2conc,             &   
                  MetState%SUNCOS,              &
                  MetState%LAT,                 &
                  MetState%DOY,                 &
                  MetState%LocalHour,           &
                  MetState%D_BTW_M,             &   
                  MetState%AEF_ISOP,            &
                  MetState%AEF_MBOX,            &
                  MetState%AEF_BPIN,            &
                  MetState%AEF_CARE,            &
                  MetState%AEF_LIMO,            &
                  MetState%AEF_OCIM,            &
                  MetState%AEF_SABI             &   
                  MetState%D2RAD,               &
                  MetState%RAD2D,               &   
                  RC)
            if (RC /= CC_SUCCESS) then
               errMsg = 'Error in CCPr_Scheme_Megan'
               CALL CC_Error( errMsg, RC, thisLoc )
            endif
         end do ! for each species requested in EmisState
      endif

   end subroutine CCPr_Megan_Run

   !>
   !! \brief Finalize Megan
   !!
   !! \param [INOUT] MeganState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_Megan_Final( MeganState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(MeganStateType), INTENT(INOUT) :: MeganState  ! MeganState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_Megan_Final (in process/megan/ccpr_megan_mod.F90)'

      ! Deallocate any arrays here
      IF ( ASSOCIATED( MeganState%SpcIDs ) ) THEN
         DEALLOCATE( MeganState%SpcIDs, STAT=RC )
         CALL CC_CheckVar('MeganState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( MeganState%MeganSpeciesIndex ) ) THEN
         DEALLOCATE( MeganState%MeganSpeciesIndex, STAT=RC )
         CALL CC_CheckVar('MeganState%MeganSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( MeganState%MeganSpeciesName ) ) THEN
         DEALLOCATE( MeganState%MeganSpeciesName, STAT=RC )
         CALL CC_CheckVar('MeganState%MeganSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( MeganState%EmissionPerSpecies ) ) THEN
         DEALLOCATE( MeganState%EmissionPerSpecies, STAT=RC )
         CALL CC_CheckVar('MeganState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( MeganState%EmisNormFactor ) ) THEN
         DEALLOCATE( MeganState%EmisNormFactor, STAT=RC )
         CALL CC_CheckVar('MeganState%EmisNormFactor', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

   end subroutine CCPr_Megan_Final

END MODULE CCPR_Megan_Mod
