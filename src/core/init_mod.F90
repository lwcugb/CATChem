! \file init_mod.F90
!! \brief Initialization module for the program.
!!
!! This module contains subroutines and functions related to the initialization of the program.
!! It includes subroutines for initializing the grid, the time step, and the solution.
!!
!! \ingroup core_modules
!!!>
module init_mod

   implicit none

   PUBLIC :: Init_Met
   PUBLIC :: Init_Diag
   PUBLIC :: Init_Emis

contains

   ! \brief Initialize the met state
   !!
   !! This subroutine allocates the met state.
   !!
   !! \param GridState The grid state containing information about the grid.
   !! \param MetState The met state to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Met(GridState, MetState, RC)
      !
      use GridState_Mod, Only : GridStateType
      use MetState_Mod
      USE Error_Mod
      implicit none

      ! Arguments
      TYPE(GridStateType), INTENT(IN)  :: GridState
      TYPE(MetStateType), INTENT(INOUT) :: MetState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Init_Met (in core/state_mod.F90)'

      call Met_Allocate(GridState, MetState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating met state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Met

   ! \brief Initialize the emission state
   !!
   !! This subroutine allocates the emission state.
   !!
   !! \param GridState The grid state containing information about the grid
   !! \param EmisState The emission state to be initialized
   !! \param ChemState The chemical state containing species information
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Emis(GridState, EmisState, ChemState, RC)
      !
      use GridState_Mod, Only : GridStateType
      use EmisState_Mod
      use ChemState_Mod, Only : ChemStateType
      USE Error_Mod
      implicit none

      ! Arguments
      TYPE(GridStateType), INTENT(IN)  :: GridState
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      TYPE(ChemStateType), INTENT(INOUT) :: ChemState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Init_Emis (in core/init_mod.F90)'

      call Emis_Allocate(GridState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating emission state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

      !Find mapping index from emission to chemical species
      call Emis_Find_Chem_Map_Index(EmisState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error finding emission to chemical species mapping index'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Emis

   ! \brief Initialize the diag state
   !!
   !! This subroutine allocates the diag state.
   !!
   !! \param Config The configuration options
   !! \param DiagState The diag state to be initialized
   !! \param ChemState The chemical state containing species information
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Diag(Config, DiagState, ChemState, RC)
      use DiagState_Mod
      use Config_Opt_Mod, Only : ConfigType
      ! use GridState_Mod, Only : GridStateType
      use ChemState_Mod, Only : ChemStateType
      use Error_Mod

      implicit none

      ! Arguments
      TYPE(ConfigType),    INTENT(IN)    :: Config
      ! TYPE(GridStateType), INTENT(IN)    :: GridState
      TYPE(DiagStateType), INTENT(INOUT) :: DiagState
      TYPE(ChemStateType), INTENT(INOUT)    :: ChemState
      INTEGER,         INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at Init_Diag (in core/init_mod.F90)'

      call Diag_Allocate(Config, DiagState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating diag state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Diag


end module init_mod
