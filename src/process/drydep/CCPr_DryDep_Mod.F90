!> \brief CCPR drydep state types
!!
!! \defgroup catchem_drydep_process
!!
!! \author Lacey Holland
!! \date 07/2024
!!!>
MODULE CCPR_DryDep_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod, Only : ConfigType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_DryDep_Init
   PUBLIC :: CCPR_DryDep_Run
   PUBLIC :: CCPR_DryDep_Finalize
   PUBLIC :: DryDepStateType


   !> \brief DryDepStateType
   !!
   !! DryDepStateType is the process-specific derived type.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param Scheme Scheme Option
   !! \param DryDepSpeciesIndex Effected Chemical Species from DryDep
   !! \param nSpc # of species
   !! \param SpcIDs CATChem species IDs
   !! \param ScaleFactor Scale Factor
   !! \param Resuspension Activate resuspension  (True/False)
   !!
   !! \ingroup core_modules
   !!!>
   TYPE :: DryDepStateType

      ! Process Specific Parameters

      ! Namelist parameters for specific DryDep goes here as well
      !=================================================================
      ! Module specific variables/arrays/data pointers come below
      !=================================================================

      LOGICAL                         :: Activate              ! Activate Process (True/False)
      LOGICAL                         :: Resuspension          ! Activate resuspension  (True/False)
      INTEGER                         :: SchemeOpt             ! Scheme Option (if there is only one SchemeOpt always = 1)
      real                            :: particleradius        ! Particle radius (m)
      real                            :: particledensity       ! Particle density (kg/m^3)
      real, allocatable               :: drydep_frequency(:)      ! could have one per chem species, revisit later
      real, allocatable               :: drydep_vel(:)            ! could have one per chem species, revisit later


   END TYPE DryDepStateType


CONTAINS

   !>
   !! \brief Initialize the CATChem DryDep module
   !!
   !! \param Config       CATCHem configuration options
   !! \param DryDepState   CATCHem PROCESS state
   !! \param ChemState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !! \ingroup catchem_drydep_process
   !!
   !!!>
   SUBROUTINE CCPR_DryDep_Init( Config, DryDepState, ChemState, RC )
      ! USE


      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)    :: Config    ! Module options
      TYPE(ChemStateType) :: ChemState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(DryDepStateType)          :: DryDepState ! DryDep state
      INTEGER,         INTENT(INOUT) :: RC       ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------


      ! Put any local variables here

      !=================================================================
      ! CCPR_DryDep_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_DryDep_INIT (in process/drydep/ccpr_drydep_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%drydep_activate) then

         ! Activate Process
         !------------------
         DryDepState%Activate = .true.
         allocate(DryDepState%drydep_frequency(ChemState%nSpeciesAeroDryDep), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate DryDepState%drydep_frequency(ChemState%nSpeciesAeroDryDep)'
            CALL CC_Error( ErrMsg, RC, ThisLoc )
         ENDIF
         DryDepState%drydep_frequency(1:ChemState%nSpeciesAeroDryDep)=ZERO

         allocate(DryDepState%drydep_vel(ChemState%nSpeciesAeroDryDep), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate DryDepState%drydep_vel(ChemState%nSpeciesAeroDryDep)'
            CALL CC_Error( ErrMsg, RC, ThisLoc )
         ENDIF
         DryDepState%drydep_vel(1:ChemState%nSpeciesAeroDryDep)=ZERO

         ! Set scheme option
         !------------------
         ! For now, the only option is SchemeOpt = 1
         if (Config%drydep_scheme==1) then
            DryDepState%SchemeOpt = 1
         else
            DryDepState%SchemeOpt = 1
         end if
         DryDepState%Resuspension = Config%drydep_resuspension
      else
         DryDepState%Activate = .false.
      end if

   end subroutine CCPR_DryDep_Init

   !>
   !! \brief Run the DryDep
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] DiagState - The DiagState object
   !! \param [INOUT] DryDepState - The DryDepState object
   !! \param [INOUT] ChemState - The ChemState object
   !! \param [OUT] RC Return code
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   SUBROUTINE CCPr_DryDep_Run( MetState, DiagState, DryDepState, ChemState, RC )

      ! USE
      USE constants, only : Cp, g0, VON_KARMAN
      use CCPr_Scheme_GOCART_DryDep_Mod, only : CCPr_Scheme_GOCART_DryDep

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       !< MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DiagStateType), INTENT(INOUT)      :: DiagState       !< DiagState Instance
      TYPE(DryDepStateType), INTENT(INOUT)    :: DryDepState     !< DryDepState Instance
      TYPE(ChemStateType),  INTENT(INOUT)     :: ChemState       !< ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                 ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      INTEGER :: km
      INTEGER :: i !< counter
      real :: radius
      real :: rhop
      real :: drydepf(1,1)
      REAL(fp) :: dqa                                    ! Change in Species due to drydep
      REAL(fp) :: SpecConc                               ! Temporary Species concentration

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DryDep_Run (in process/drydep/ccpr_DryDep_mod.F90)'

      km = MetState%NLEVS

      ! Run the DryDep Scheme
      !-------------------------
      if (DryDepState%Activate) then
         ! Run the DryDep Scheme
         !-------------------------
         if (DryDepState%SchemeOpt == 1) then
            ! Run the DryDep Scheme - Only Applicable to AEROSOL species
            !-------------------------
            if (ChemState%nSpeciesAeroDryDep > 0) then

               ! loop through aerosol species
               do i = 1, ChemState%nSpeciesAeroDryDep

                  radius = ChemState%chemSpecies(ChemState%DryDepIndex(i))%radius
                  rhop = ChemState%chemSpecies(ChemState%DryDepIndex(i))%density

                  call CCPr_Scheme_GOCART_DryDep( MetState%NLEVS,   &
                     MetState%T,       &
                     MetState%AIRDEN,  &
                     MetState%ZMID,    &
                     MetState%LWI,     &
                     MetState%USTAR,   &
                     MetSTate%PBLH,    &
                     MetState%HFLUX,   &
                     VON_KARMAN,       &
                     Cp,               &
                     g0,               &
                     MetState%Z0H,     &
                     drydepf,          &
                     DryDepState%Resuspension, &
                     radius,           &
                     rhop,             &
                     MetState%U10M,    &
                     MetSTate%V10M,    &
                     MetState%FRLAKE,  &
                     MetState%GWETTOP, &
                     RC)


                  if (RC /= 0) then
                     errMsg = 'Error in GOCART DryDeposition'
                     CALL CC_Error( errMsg, RC, thisLoc )
                  endif  !if (RC /= CC_SUCCESS)

                  ! Fill Diagnostic Variables
                  !--------------------------
                  DryDepState%drydep_frequency(ChemState%DryDepIndex(i)) = drydepf(1,1)
                  DryDepState%drydep_vel(ChemState%DryDepIndex(i)) = MetState%ZMID(1) * drydepf(1,1)
                  DiagState%drydep_frequency(i)= drydepf(1,1)
                  DiagState%drydep_vel(i) = MetState%ZMID(1) * drydepf(1,1)

                  ! apply drydep velocities/freq to chem species
                  dqa = 0.
                  SpecConc = ChemState%chemSpecies(ChemState%DryDepIndex(i))%conc(1)
                  dqa = MAX(0.0_fp, SpecConc * (1.-exp(-1*drydepf(1,1) * MetState%TSTEP)))
                  ChemState%chemSpecies(ChemState%DryDepIndex(i))%conc(1) = SpecConc - dqa

               end do ! do i = 1, ChemState%nSpeciesAeroDryDep

            endif  ! if (ChemState%nSpeciesAeroDryDep > 0)

         endif  ! if (DryDepState%SchemeOpt == 1)

         ! TO DO:  apply dry dep velocities/freq to chem species
         write(*,*) 'TODO: Need to figure out how to add back to the chemical species state '

      endif   !  if (DryDepState%Activate)


   end subroutine CCPr_DryDep_Run

   !>
   !! \brief Finalize the DryDep
   !!
   !! \param [INOUT] DryDepState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_DryDep_Finalize( DryDepState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DryDepStateType), INTENT(INOUT) :: DryDepState  ! DryDepState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DryDep_Finalize (in process/drydep/ccpr_DryDep_mod.F90)'

      DEALLOCATE( DryDepState%drydep_frequency, STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate DryDepState%drydep_frequency'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

      DEALLOCATE( DryDepState%drydep_vel, STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate DryDepState%drydep_vel'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
      ENDIF


   end subroutine CCPr_DryDep_Finalize


END MODULE CCPR_DryDep_Mod
