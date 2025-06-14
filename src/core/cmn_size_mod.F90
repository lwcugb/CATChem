!> \file cmn_size_mod.F90
!! \brief Module for common size parameters
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!!
!! This module defines common size parameters and constants used throughout
!! the CATChem atmospheric chemistry model for grid dimensions, surface types,
!! and numerical parameters.
!!
!! \details
!! The common size module provides standardized size parameters for:
!! - Atmospheric pressure levels
!! - Surface and vegetation type classifications
!! - Polynomial fitting coefficients
!!
!! \section cmn_size_usage Usage Example
!! \code{.f90}
!! use cmn_size_mod
!! integer :: surface_types(NSURFTYPE)
!! real(fp) :: top_pressure = PTOP
!! \endcode
!!
MODULE CMN_SIZE_MOD
   !
   ! !USES:
   !
   USE PRECISION_MOD
   IMPLICIT NONE
   PUBLIC
   !
   ! !DEFINED PARAMETERS:
   !
   ! \name Atmospheric Parameters
   !! \brief Parameters related to atmospheric structure
   !! \{
   REAL(fp), PARAMETER :: PTOP = 0.01_fp  !< Model top pressure [hPa]
   ! \}

   ! \name Surface Classification Parameters
   !! \brief Parameters for surface and vegetation type classifications
   !! \{
   INTEGER, PARAMETER :: NSURFTYPE = 73   !< Maximum number of surface types (Olson land use categories)
   INTEGER, PARAMETER :: NTYPE = 20       !< Maximum number of vegetation types in a CTM grid box
   ! \}

   ! \name Numerical Parameters
   !! \brief Parameters for numerical calculations
   !! \{
   INTEGER, PARAMETER :: NPOLY = 20       !< Number of coefficients for polynomial fits
   ! \}

END MODULE CMN_SIZE_MOD
