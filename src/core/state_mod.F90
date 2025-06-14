!> \file state_mod.F90
!! \brief Core state management module for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!!
!! This module provides the main state objects that manage the complete
!! state of the CATChem atmospheric chemistry modeling system, including
!! grid, meteorology, chemistry, emissions, and diagnostics.
!!
!! \details
!! The state module centralizes all major state objects used throughout
!! CATChem, providing global access to:
!! - Grid configuration and dimensions (GridState)
!! - Meteorological fields and parameters (MetState)
!! - Chemical species concentrations (ChemState)
!! - Configuration options and settings (Config)
!! - Emission fluxes and mappings (EmisState)
!! - Diagnostic output variables (DiagState)
!!
!! \section state_usage Usage Example
!! \code{.f90}
!! use state_mod
!! ! Grid and meteorology are now globally accessible
!! call initialize_grid(GridState)
!! call setup_meteorology(MetState)
!! \endcode
!!
module state_mod
   use precision_mod
   use Config_Opt_Mod, only : ConfigType
   use GridState_Mod,  only : GridStateType
   use MetState_Mod,   only : MetStateType
   use ChemState_Mod,  only : ChemStateType
   use EmisState_Mod,  only : EmisStateType
   use DiagState_Mod,  only : DiagStateType

   IMPLICIT NONE

   ! PUBLIC
   type(GridStateType), PUBLIC :: GridState  !< Global grid state object
   type(MetStateType),  PUBLIC :: MetState   !< Global meteorological state object
   type(ChemStateType), PUBLIC :: ChemState  !< Global chemical state object
   type(ConfigType),    PUBLIC :: Config     !< Global configuration object
   type(EmisStateType), PUBLIC :: EmisState  !< Global emission state object
   type(DiagStateType), PUBLIC :: DiagState  !< Global diagnostic state object

end module state_mod
