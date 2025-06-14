!> \file Process_driver_template.F90
!! \brief Template for new process drivers in CATChem
!! \ingroup process_modules
!!
!! \author Barry Baker
!! \date 05/2024
!!
!! This file provides a template for creating new atmospheric chemistry
!! process drivers in the CATChem modeling system.
!!
!! \details
!! To use this template:
!! - Replace <PROCESS> with an identifier for the process group (e.g., 'Dust', 'SeaSalt')
!! - Implement the required initialization, run, and finalization routines
!! - Add process-specific configuration options and state variables
!! - Update documentation to describe the specific process being implemented
!!
!! The template includes standard CATChem interfaces for:
!! - Process initialization with configuration
!! - Process execution with meteorological and chemical states
!! - Process finalization and cleanup
!!
MODULE CCPR_<PROCESS>_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Mod,    Only : ConfigType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_<PROCESS>_Init
   PUBLIC :: CCPR_<PROCESS>_Run
   PUBLIC :: CCPR_<PROCESS>_Final


   ! \brief <PROCESS>StateType
   !!
   !! <PROCESS>StateType is the process-specific derived type. It should hold all module
   !! variables and arrays that are required to compute the emissions.
   !! For instance, if the process relies on an input field read through the
   !! CATChem configuration file (e.g. MY_INPUT_FIELD), the data array pointer
   !! to that field should be listed within the instance and NOT outside of it.
   !! This ensures that the same process can be invoked in various instances,
   !! all of them potentially pointing to different data fields.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param SchemeOpt Scheme Option
   !! \param <Process>SpeciesIndex Effected Chemical Species from <Process>
   !! \param nSpc # of species
   !! \param SpcIDs CATChem species IDs
   !! \param ScaleFactor Scale Factor
   !!!>
   TYPE :: <PROCESS>StateType
      LOGICAL                         :: Activate              ! Activate Process (True/False)
      INTEGER                         :: SchemeOpt             ! Scheme Option (if there is only one SchemeOpt always = 1)
      INTEGER                         :: <Process>SpeciesIndex ! Effected Chemical Species from <Process>
      INTEGER                         :: nSpc                  ! # of species
      INTEGER,  ALLOCATABLE           :: SpcIDs(:)             ! CATChem species IDs

      ! Namelist parameters for specific <PROCESS> goes here as well
      !=================================================================
      ! Module specific variables/arrays/data pointers come below
      !=================================================================
      integer                         :: <Process><Scheme>Opt  ! Scheme Option
      real(fp)                        :: ScaleFactor           ! Scale Factor
   END TYPE <PROCESS>StateType


CONTAINS

   !>
   !! \brief Initialize the CATChem <PROCESS> module
   !!
   !! \param Config_Opt       CATCHem configuration options
   !! \param <PROCESS>State   CATCHem PROCESS state
   !! \param ChmState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !!!>
   SUBROUTINE CCPR_<yourname>_Init( Config, ChemState, <Process>State, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigOptType), POINTER       :: Config    ! Module options
      TYPE(ChemStateType), POINTER       :: ChemState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(<PROCESS>StateType), POINTER :: <Process>State ! <PROCESS> state
      INTEGER,         INTENT(INOUT)    :: RC         ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------

      ! Put any local variables here

      !=================================================================
      ! CCPR_<yourname>_Init begins here!
      !=================================================================
      ThisLoc = ' -> at CCPR_<PROCESS>_INIT (in process/<PROCESS>/ccpr_<PROCESS>_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%<process>activate) then

         ! Activate Process
         !------------------
         <Process>State%Activate = .true.

         ! Set number of species
         !----------------------
         <Process>State%nSpc = 1

         ! Set scheme option
         !------------------
         <Process>State%SchemeOpt = config%<process><Scheme>Opt

         ! Allocate any arrays here for scheme to run
         ALLOCATE(<Process>State%SpcIDs(<Process>State%nSpc))
         CALL CC_CheckVar(<Process>State%SpcIDs(<Process>State%nSpc), 0, RC)  ! Assumed here that the length of SpcIDs is equal to <Process>State%nSpc
         IF (RC /= CC_SUCCESS) RETURN
         <Process>State%SpcIDs(<Process>State%nSpc) = 1

      else

         <Process>State%Activate = .false.

      endif

   end subroutine CCPR_<yourname>_Init

   !>
   !! \brief Run the <Process>
   !!
   !! \param
   !! \param
   !! \param
   !! \param
   !! \param
   !!!>
   SUBROUTINE CCPr_Dust_Run( MetState, DiagState, <PROCESS>State, ChemState, RC )

      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetState_type),  INTENT(IN) :: MetState       ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DiagState_type), INTENT(INOUT)      :: DiagState       ! DiagState Instance
      TYPE(<PROCESS>State_type), INTENT(INOUT) :: <PROCESS>State  ! <PROCESS>State Instance
      TYPE(ChemState_type),  INTENT(INOUT)     :: ChemState       ! ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_<PROCESS>_Run (in process/<PROCESS>/ccpr_<PROCESS>_mod.F90)'

      ! Run the <PROCESS> Scheme
      !-------------------------
      if (<PROCESS>State%Activate) then
         ! Run the <PROCESS> Scheme
         !-------------------------
         if (<PROCESS>State%SchemeOpt == 1) then
            ! Run the <PROCESS> Scheme
            !-------------------------
            CCPr_<Process>_Scheme_<SCHEME_NAME>( MetState, DiagState, <PROCESS>State, ChemState, RC)

         elseif (<PROCESS>State%SchemeOpt == 2) then
            ! Run the <PROCESS> Scheme
            !-------------------------
            CCPr_<Process>_Scheme_<2ndSCHEME_NAME>( MetState, DiagState, <PROCESS>State, ChemState, RC)
         endif

      endif

   end subroutine CCPr_<PROCESS>_Run

   !>
   !! \brief Finalize the <Process>
   !!
   !! \param
   !! \param
   !!!>
   SUBROUTINE CCPr_<PROCESS>_Final( <PROCESS>State, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(<PROCESS>State_type), INTENT(INOUT) :: <PROCESS>State  ! <PROCESS>State Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_<PROCESS>_Final (in process/<PROCESS>/ccpr_<PROCESS>_mod.F90)'

      ! Deallocate any arrays here
      DELLOCATE( <PROCESS>State%SpcIDs, STAT=RC )
      CALL CC_CheckVar('<PROCESS>State%SpcIDs', 0, RC)
      IF (RC /= CC_SUCCESS) RETURN
      <PROCESS>State%SpcIDs => NULL()

   end subroutine CCPr_<PROCESS>_Final

END MODULE CCPR_<yourname>_Mod
