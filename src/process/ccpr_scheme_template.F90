!> \file ccpr_scheme_template.F90
!! \brief Template for new process schemes in CATChem
!! \ingroup process_modules
!!
!! \author Barry Baker
!! \date 05/2024
!!
!! This file provides a template for creating new atmospheric chemistry
!! process schemes in the CATChem modeling system.
!!
!! \details
!! To use this template:
!! - Replace SCHEME with a scheme identifier (e.g., 'Fengsha', 'Gong97')
!! - Replace PROCESS with the process group identifier (e.g., 'Dust', 'SeaSalt')
!! - Implement the scheme-specific calculations and parameterizations
!! - Add appropriate scientific references and citations
!! - Update documentation to describe the scheme methodology
!!
!! The template provides the basic module structure for implementing
!! scientific parameterizations within the CATChem process framework.
!!
!! \section scheme_reference Reference
!! Please include appropriate scientific citations here
!!
module CCPr_Scheme_SCHEME_Mod

   implicit none

   private

   public :: CCPr_Scheme_SCHEME

contains

   ! \brief Brief description of the subroutine
   !!
   !! \param MetState     Meteorological Variables
   !! \param DiagState    Diagnostic Variables
   !! \param DustState    Dust Variables
   !! \param RC           Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>
   subroutine CCPr_Scheme_SCHEME(MetState, DiagState, RC)

      ! Uses
      USE Constants,     Only : g0 ! Example to pull in a constant from the CONSTANTS MODULE < Modify as needed >
      use precision_mod, only : fp, ZERO  ! Example to pull in a precision from the PRECISION MODULE < Modify as needed >
      Use MetState_Mod,  Only : MetStateType  ! Needed to access Meteorological Variables
      Use DiagState_Mod, Only : DiagStateType ! Diagnostic Variables are added through DiagState below
      Use Error_Mod,     Only : CC_SUCCESS    ! Error Check Success
      Use CCPr_PROCESS_mod, Only : PROCESSStateType  ! Overall PROCESS State Type - Controlling PROCESS

      IMPLICIT NONE

      ! Arguments
      type(MetStateType),  intent(in) :: MetState     ! Meteorological Variables
      type(DiagStateType), intent(in) :: DiagState    ! Diagnostic Variables
      type(PROCESSStateType), intent(in) :: PROCESSState  ! PROCESS Variables

      integer, intent(out) :: RC                      ! Success or Failure

      ! Local Variables
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      integer :: SCHEMEInt   ! Add Local Scheme Specific Variables < Modify as needed >
      real(fp) :: SCHEMEFlt  ! Add Local Scheme Specific Variables < Modify as needed >

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_SCHEME (in CCPr_Scheme_SCHEME_mod.F90)'
      RC = CC_SUCCESS

      !------------------
      ! Begin Scheme Code
      !------------------

      ! Begin SCHEME Code

      ! SCHEME Options comes in from PROCESSState
      ! Diagnostic Variables are added through DiagState below

      ! End SCHEME Code

   end subroutine CCPr_Scheme_SCHEME

end module CCPr_Scheme_SCHEME_Mod
