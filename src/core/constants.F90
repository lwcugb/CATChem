!> \file constants.F90
!! \brief Physical and mathematical constants for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!!
!! This module provides fundamental physical and mathematical constants
!! used throughout the CATChem atmospheric chemistry model.
!!
!! \details
!! The constants module defines all physical constants, conversion factors,
!! and mathematical constants used in atmospheric chemistry calculations.
!! All values are given in SI units unless otherwise specified.
!!
!! \section constants_usage Usage Example
!! \code{.f90}
!! use constants
!! real(fp) :: air_density
!! air_density = pressure / (Rd * temperature)
!! \endcode
!!
module constants

   USE precision_mod

   implicit none

   public

   ! \name Atmospheric Properties
   !! \brief Constants related to atmospheric composition and properties
   !! \{
   REAL(fp), PARAMETER :: Cp = 1.0046e+3_fp                !< Specific heat of dry air at constant pressure [J/kg/K]
   REAL(fp), PARAMETER :: Cv = 7.1760e+2_fp                !< Specific heat of dry air at constant volume [J/kg/K]
   REAL(fp), PARAMETER :: AIRMW = 28.9644_fp               !< Average molecular weight of dry air [g/mol]
   REAL(fp), PARAMETER :: H2OMW = 18.016_fp                !< Molecular weight of water [g/mol]
   REAL(fp), PARAMETER :: Rd   = 287.0_fp                  !< Gas constant for dry air [J/K/kg]
   REAL(fp), PARAMETER :: Rdg0 = Rd / g0                   !< Gas constant for dry air divided by gravity
   REAL(fp), PARAMETER :: Rv = 461.00_fp                   !< Gas constant for water vapor [J/K/kg]
   REAL(fp), PARAMETER :: SCALE_HEIGHT = 7600.0_fp         !< Atmospheric scale height [m]
   REAL(fp), PARAMETER :: VON_KARMAN = 0.41_fp             !< Von Karman's constant (dimensionless)
   REAL(fp), PARAMETER :: ATM = 1.01325e+5_fp              !< Standard atmospheric pressure [Pa]
   REAL(fp), PARAMETER :: XNUMOLAIR = AVO / ( AIRMW * 1.e-3_fp )  !< Molecules of dry air per kg dry air
   ! \}

   ! \name Fundamental Physical Constants
   !! \brief Universal physical constants
   !! \{
   REAL(fp), PARAMETER :: AVO = 6.022140857e+23_fp         !< Avogadro's number [particles/mol]
   REAL(fp), PARAMETER :: g0     = 9.80665e+0_fp           !< Standard gravity acceleration [m/s^2]
   REAL(fp), PARAMETER :: g0_100 = 100.0_fp / g0           !< 100 divided by standard gravity
   REAL(fp), PARAMETER :: Re = 6.3710072e+6_fp             !< Earth's radius [m]
   REAL(fp), PARAMETER :: RSTARG = 8.3144598_fp            !< Universal gas constant [J/K/mol]
   REAL(fp), PARAMETER :: BOLTZ = 1.38064852e-23_fp        !< Boltzmann's constant [J/K]
   REAL(fp), PARAMETER :: PLANCK = 6.62606957e-34_fp       !< Planck's constant [J⋅s]
   REAL(fp), PARAMETER :: CCONST = 2.99792458e+8_fp        !< Speed of light in vacuum [m/s]
   ! \}

   ! \name Mathematical Constants
   !! \brief Mathematical constants and conversion factors
   !! \{
   REAL(fp), PARAMETER :: PI     = 3.14159265358979323_fp  !< Pi (dimensionless)
   REAL(fp), PARAMETER :: PI_180 = PI / 180.0_fp           !< Radians per degree conversion factor
   REAL(fp), PARAMETER :: E = 2.718281828459045235360287471352_fp  !< Euler's number (dimensionless)
   ! \}

   ! \name Chemistry-Specific Constants
   !! \brief Constants for atmospheric chemistry calculations
   !! \{
   REAL(fp), PARAMETER :: CONSVAP = 6.1078e+03_fp / ( BOLTZ * 1e+7_fp ) !< Condensation vapor pressure factor
   REAL(fp), PARAMETER :: RGASLATM = 8.2057e-2_fp          !< Gas constant in L⋅atm/(K⋅mol)
   REAL(fp), PARAMETER :: MWCARB = 12.01e-3_fp             !< Molecular weight of carbon [kg/mol]
   ! \}

end module constants
