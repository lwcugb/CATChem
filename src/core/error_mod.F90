! \file error_mod.F90
!! \brief Error handling and diagnostic routines for CATChem
!!
!! This module provides standardized error handling, warning messages,
!! and variable checking routines for the CATChem atmospheric chemistry model.
!!
!! \author Barry Baker
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!! \ingroup core_modules
!!
!! \details
!! The error handling module provides consistent error reporting across
!! CATChem with standardized return codes and message formatting.
!! It includes routines for error messages, warnings, and variable validation.
!!
!! \section error_usage Usage Example
!! \code{.f90}
!! use error_mod
!! integer :: rc
!! call CC_Error('Invalid input parameter', rc, 'my_subroutine')
!! if (rc /= CC_SUCCESS) return
!! \endcode
!!
! \brief Error handling and diagnostic routines
!!
!! This module provides error handling, warnings, and validation routines
MODULE Error_Mod
   !
   ! !USES:
   !
   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   !
   PUBLIC :: CC_Error
   PUBLIC :: CC_Warning
   PUBLIC :: CC_CheckVar
   !
   ! !DEFINED PARAMETERS:
   !
   ! \name Return Codes
   !! \brief Standard return codes for CATChem routines
   !! \{
   INTEGER, PUBLIC, PARAMETER :: CC_SUCCESS =  0   ! Routine completed successfully
   INTEGER, PUBLIC, PARAMETER :: CC_FAILURE = -1   ! Routine failed to complete
   ! \}

CONTAINS

   ! \brief Display error message and set failure return code
   !!
   !! This subroutine prints a formatted error message to standard output
   !! and sets the return code to indicate failure. The message includes
   !! optional location and instruction information.
   !!
   !! \param ErrMsg Error message to display
   !! \param RC Return code (set to CC_FAILURE)
   !! \param ThisLoc Optional location where error occurred
   !! \param Instr Optional additional instructions for user
   !!
   !! \par Example:
   !! \code{.f90}
   !! call CC_Error('Invalid temperature value', rc, 'temperature_check', &
   !!               'Check input data file')
   !! \endcode
   SUBROUTINE CC_Error( ErrMsg, RC, ThisLoc, Instr )
      !
      ! !USES:
      !
      USE Charpak_Mod,    ONLY : WordWrapPrint
      !
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN)            :: ErrMsg  ! Message to display
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: ThisLoc ! Location of error
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: Instr   ! Other instructions
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT)            :: RC      ! Error code

      CHARACTER(LEN=1000) :: Message
      !=======================================================================
      ! CC_ERROR begins here
      !=======================================================================

      ! Construct error message

      ! Separator
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Print error message to log
      Message =  'CATChem ERROR: ' // TRIM( ErrMsg )
      CALL WordWrapPrint( Message, 78 )

      ! Print error location to log
      IF ( PRESENT( ThisLoc ) ) THEN
         Message = 'ERROR LOCATION: ' // TRIM( ThisLoc )
         WRITE( 6, '(a)' ) TRIM( ThisLoc )
      ENDIF

      ! Print additional instructions to log
      IF ( PRESENT( Instr ) ) THEN
         WRITE( 6, '(a)' )
         CALL WordWrapPrint( Instr, 78 )
      ENDIF

      ! Separators
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) ''


      ! Force the message to be flushed to the log file
      CALL Flush( 6 )

      ! Return with failure, but preserve existing error code
      IF ( RC == CC_SUCCESS ) THEN
         RC = CC_FAILURE
      ENDIF

   END SUBROUTINE CC_Error

   !>
   !! \brief CC_Warning
   !!
   !! This subroutine prints a warning message and sets RC to CC_SUCCESS.
   !!
   !! \ingroup core_modules
   !!
   !! \param WarnMsg The warning message
   !! \param RC The return code
   !! \param ThisLoc The location of the warning
   !! \param Instr Other instructions
   !!!>
   SUBROUTINE CC_Warning( WarnMsg, RC, ThisLoc, Instr )
      !
      ! !USES:
      !
      USE Charpak_Mod, ONLY : WordWrapPrint
      !!
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN   )            :: WarnMsg
      CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: ThisLoc
      CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Instr
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT)            :: RC

      CHARACTER(LEN=1000) :: Message

      !=======================================================================
      ! CC_ERROR begins here
      !=======================================================================

      ! Separator
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Print error message to log
      Message =  'CATChem WARNING: ' // TRIM( WarnMsg )
      CALL WordWrapPrint( Message, 78 )

      ! Print error location to log
      IF ( PRESENT( ThisLoc ) ) THEN
         Message = 'WARNING LOCATION: ' // TRIM( ThisLoc )
         WRITE( 6, '(a)' ) TRIM( ThisLoc )
      ENDIF

      ! Print additional instructions to log
      IF ( PRESENT( Instr ) ) THEN
         WRITE( 6, '(a)' )
         CALL WordWrapPrint( Instr, 78 )
      ENDIF

      ! Separators
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) ''

      ! Force the message to be flushed to the log file
      CALL Flush( 6 )

      ! Return with success, since this is only a warning message
      RC = CC_SUCCESS

   END SUBROUTINE CC_Warning

   !>
   !! \brief CC_CheckVar
   !!
   !! This subroutine checks if a variable is allocated.
   !!
   !! \param Variable The variable to check
   !! \param Operation 0=Allocate 1=Register 2=Deallocate
   !! \param RC The return code
   !! \ingroup core_modules
   !!!>
   SUBROUTINE CC_CheckVar( Variable, Operation, RC )
      !
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN)    :: Variable   ! Name of variable to check
      INTEGER,          INTENT(IN)    :: Operation  ! 0=Allocate
      ! 1=Register
      ! 2=Deallocate
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT) :: RC         ! Success or failure
      !
      ! !LOCAL VARIABLES:
      !
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=========================================================================
      ! Initialize
      !=========================================================================

      ! Define error message
      SELECT CASE( Operation )
       CASE( 1 )
         ErrMsg = 'Could not register '   // TRIM( Variable ) // '!'
       CASE( 2 )
         ErrMsg = 'Could not deallocate ' // TRIM( Variable ) // '!'
       CASE DEFAULT
         ErrMsg = 'Could not allocate '   // TRIM( Variable ) // '!'
      END SELECT

      ! Define location string
      ThisLoc   = ' -> at CC_CheckVar (in Headers/errcode_mod.F90)'

      !=========================================================================
      ! Display error message if necessary
      !=========================================================================
      IF ( RC /= CC_SUCCESS ) THEN
         CALL CC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

   END SUBROUTINE CC_CheckVar
   !EOC
END MODULE Error_Mod
