!> \file
!! \brief CATChem core data types and routines.
!!
!! \defgroup catchem_api
!! CATChem core data types and routines.
!!
!!!>
module CATChem

   !---------------
   ! CATChem States
   !---------------
   use ChemState_Mod,  only: ChemStateType    !< Chemical State
   !
   use ChemState_Mod, only: ChemStateType
   use ChemState_Mod, only: cc_find_species_by_name => FindSpecByName
   use ChemState_Mod, only: cc_get_species_conc => GetSpecConc
   use ChemState_Mod, only: cc_get_species_conc_by_name => GetSpecConcByName
   use ChemState_Mod, only: cc_get_species_conc_by_index => GetSpecConcByIndex
   !
   use Config_Opt_Mod, only: ConfigType
   !
   use Config_Mod, only: cc_read_config => Read_Input_File
   !
   use DiagState_Mod, only: DiagStateType
   !
   use Error_Mod, only: cc_check_var => CC_CheckVar
   use Error_Mod, only: cc_emit_error => CC_Error
   use Error_Mod, only: CC_FAILURE
   use Error_Mod, only: CC_SUCCESS
   use Error_Mod, only: cc_emit_warning => CC_Warning
   !
   use GridState_Mod, only: GridStateType
   !
   use DiagState_Mod,  only: DiagStateType    !< Diagnostic State
   !
   use EmisState_Mod,  only: EmisStateType    !< Emission State
   !
   use GridState_Mod,  only: GridStateType    !< Grid State
   !
   use MetState_Mod,   only: MetStateType     !< Meteorology State
   !
   use species_mod,    only: SpeciesType      !< Species State

   !-------------------
   ! Configuration Read
   !-------------------
   use Config_Mod, only: cc_read_config => Read_Input_File  !< Method for reading the configuration file

   !---------------
   ! Error Handling
   !---------------
   use Error_Mod, only: cc_check_var => CC_CheckVar     !< Method for checking variables
   use Error_Mod, only: cc_emit_error => CC_Error       !< Method for emitting errors
   use Error_Mod, only: CC_FAILURE                      !< CATCHem Failure return code
   use Error_Mod, only: CC_SUCCESS                      !< CATChem Successful return code
   use Error_Mod, only: cc_emit_warning => CC_Warning   !< Method for emitting warnings
   !
   use init_mod, only: cc_init_diag => Init_Diag        !< Method for initializing the diag state
   use init_mod, only: cc_init_met => Init_Met          !< Method for initializing the met state

   !------------------
   ! CATChem Precision
   !------------------
   use precision_mod, only: cc_rk => fp                 !< Real Precision


   !------------------
   ! CATChem Processes
   !------------------
   ! Dust
   use CCPr_Dust_Common_Mod, only: DustStateType                                !< Dust State
   use CCPr_Dust_mod, only: cc_dust_init => CCPr_Dust_Init                      !< Dust Process Initialization Routine
   use CCPr_Dust_mod, only: cc_dust_run => CCPr_Dust_Run                        !< Dust Process Run Routine
   use CCPr_Dust_mod, only: cc_dust_finalize => CCPr_Dust_Finalize              !< Dust Process Finalization Routine
   ! Seasalt
   use CCPr_SeaSalt_Common_Mod, only: SeaSaltStateType                          !< SeaSalt State
   use CCPr_SeaSalt_mod, only: cc_seasalt_init => CCPr_SeaSalt_Init             !< SeaSalt Process Initialization Routine
   use CCPr_SeaSalt_mod, only: cc_seasalt_run => CCPr_SeaSalt_Run               !< SeaSalt Process Run Routine
   use CCPr_SeaSalt_mod, only: cc_seasalt_finalize => CCPr_SeaSalt_Finalize     !< SeaSalt Process Finalization Routine
   ! Megan
   use CCPr_Megan_Common_Mod, only: MeganStateType                              !< Megan State
   use CCPr_Megan_mod, only: cc_megan_init => CCPr_Megan_Init                   !< Megan Process Initialization Routine
   use CCPr_Megan_mod, only: cc_megan_run => CCPr_Megan_Run                     !< Megan Process Run Routine
   use CCPr_Megan_mod, only: cc_megan_finalize => CCPr_Megan_Final              !< Megan Process Finalization Routine

   !
   implicit none

   public

end module CATChem
