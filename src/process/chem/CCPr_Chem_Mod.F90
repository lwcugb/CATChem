!> \file CCPr_Chem_Mod.F90
!! \brief Chemistry process module for CATChem
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2023
!!
!! This module provides interfaces to chemical mechanism solvers
!! for the CATChem atmospheric chemistry model, including integration
!! with the MICM (Multiphase International Chemical Mechanism) solver.
!!
!! \details
!! The chemistry module handles chemical kinetics calculations and
!! provides version information for the underlying chemical mechanism
!! solver implementations.
!!
module CCPr_Chem_mod
   implicit none

   private
   public :: get_micm_version

contains

   !> Get MICM version information
   !!
   !! This function retrieves the version string of the MICM chemical
   !! mechanism solver currently being used.
   !!
   !! @return res Version string of the MICM solver
   function get_micm_version() result(res)
      use musica_util, only: string_t
      use musica_micm, only: get_micm_version_ => get_micm_version
      character(len=256) :: res

      type(string_t) :: micm_version_

      micm_version_ = get_micm_version_()
      res = micm_version_%get_char_array()

   end function get_micm_version

end module CCPr_Chem_mod
