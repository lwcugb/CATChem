! \file precision_mod.F90
!! \brief Floating-point precision control and missing value definitions
!!
!! This module provides compile-time control over floating-point precision
!! throughout CATChem and defines standard missing value constants.
!!
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!! \ingroup core_modules
!!
!! \details
!! The precision module allows switching between single (4-byte) and double
!! (8-byte) precision at compile time using the USE_REAL8 preprocessor flag.
!! It also defines standard missing value constants for different data types.
!!
!! @section precision_usage Usage Example
!! @code{.f90}
!! use precision_mod
!! real(fp) :: temperature = 298.15_fp  ! Uses selected precision
!! if (abs(data - MISSING) < epsilon(data)) then
!!   ! Handle missing data
!! end if
!! @endcode
!!
!! \section precision_compilation Compilation
!! - Default: single precision (fp = f4)
!! - With -DUSE_REAL8: double precision (fp = f8)
!!
! \brief Precision control and missing value constants for CATChem
!!
!! This module provides compile-time precision control and standard constants
module precision_mod
   implicit none

   ! \name Floating-Point Kind Parameters
   !! \brief KIND parameters for different precision levels
   !! \{
   INTEGER, PARAMETER, PUBLIC :: f4 = KIND( 0.0_4 ) ! KIND parameter for single precision (4-byte)
   INTEGER, PARAMETER, PUBLIC :: f8 = KIND( 0.0_8 ) ! KIND parameter for double precision (8-byte)
   ! \}

   ! \name Selected Precision
   !! \brief Main precision kind used throughout CATChem
   !! \details Controlled by USE_REAL8 preprocessor flag
   !! \{
#ifdef USE_REAL8
   INTEGER, PARAMETER, PUBLIC :: fp = f8 ! Selected precision: double precision (8-byte)
#else
   INTEGER, PARAMETER, PUBLIC :: fp = f4 ! Selected precision: single precision (4-byte)
#endif
   ! \}

   ! \name Missing Value Constants
   !! \brief Standard missing value representations for different data types
   !! \{
   LOGICAL,          PARAMETER, PUBLIC :: MISSING_BOOL = .FALSE.   ! Missing boolean value
   INTEGER,          PARAMETER, PUBLIC :: MISSING_INT  = -999      ! Missing integer value
   REAL(fp),         PARAMETER, PUBLIC :: MISSING      = -999.0_fp ! Missing real value (selected precision)
   REAL(f4),         PARAMETER, PUBLIC :: MISSING_REAL = -999.0_f4 ! Missing single precision real value
   REAL(f8),         PARAMETER, PUBLIC :: MISSING_DBLE = -999.0_f8 ! Missing double precision real value
   CHARACTER(LEN=7), PARAMETER, PUBLIC :: MISSING_STR  = "UNKNOWN" ! Missing string value
   ! \}

   ! \name Zero Value Constants
   !! \brief Standard zero value representations
   !! \{
   REAL(fp),         PARAMETER, PUBLIC :: ZERO         =  0.0_fp   ! Zero value (selected precision)
   REAL(f4),         PARAMETER, PUBLIC :: ZERO_REAL    =  0.0_f4   ! Zero value (single precision)
   REAL(f8),         PARAMETER, PUBLIC :: ZERO_DBLE    =  0.0_f8   ! Zero value (kind=f8)
   ! \}
   !=========================================================================
   ! Parameters for very tiny numbers
   !=========================================================================
   ! \name Tiny Value Constants
   !! \brief Standard tiny value representations for numerical stability
   !! \details These values are used to avoid numerical issues in calculations
   !! \{
   REAL(f4),         PARAMETER, PUBLIC :: TINY_REAL    =  1.0e-16_f4 ! A small value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: TINY_DBLE    =  1.0e-31_f8 ! A small value (kind=f8)
#ifdef USE_REAL8
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = TINY_DBLE
#else
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = TINY_REAL
#endif
   ! \}

   !=========================================================================
   ! Parameters for one
   !=========================================================================
   ! \name One Value Constants
   !! \brief Standard one value representations for calculations
   !! \details These values are used in calculations where a unit value is needed
   !! \{
   REAL(fp),         PARAMETER, PUBLIC :: ONE          =  1.0_fp ! One value (kind=fp)
   REAL(f4),         PARAMETER, PUBLIC :: ONE_REAL     =  1.0_f4 ! One value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: ONE_DBLE     =  1.0_f8 ! One value (kind=f8)
   ! \}

   !=========================================================================
   ! Interface for Real Approximately Equal (RAE) functions
   !=========================================================================
   ! \name Real Approximately Equal Interface
   !! \brief Interface for real approximately equal functions
   !! \details Provides functions to check if two real numbers are approximately equal
   !! \{
   interface rae
      module procedure rae_f4, rae_f8
   end interface rae
   ! \}

contains

   ! Real approximately equal: `abs(a - b) < tiny(a)`
   logical function rae_f4(a, b) result(res)
      real(f4), intent(in) :: a, b
      real(f4) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f4

   ! Real approximately equal: `abs(a - b) < tiny(a)`
   logical function rae_f8(a, b) result(res)
      real(f8), intent(in) :: a, b
      real(f8) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f8

end module precision_mod
