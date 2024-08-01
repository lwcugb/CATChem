!> \file cmn_size_mod.F90
!! \brief Module for common size parameters
!!
!! \ingroup core_modules
!!!>
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
   REAL(fp), PARAMETER :: PTOP = 0.01_fp

   ! Maximum number of surface types: 73 olson
   INTEGER, PARAMETER :: NSURFTYPE = 73

   ! Maximum number of veg types in a CTM grid box
   INTEGER, PARAMETER :: NTYPE=20

   !Number of coefficients for polynomial fits
   INTEGER, PARAMETER :: NPOLY = 20

END MODULE CMN_SIZE_MOD
