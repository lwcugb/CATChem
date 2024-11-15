program test_dust
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use state_mod
   use ccpr_megan_common_mod, only : MeganStateType
   !use EmisState_Mod, only: Emis_Allocate

   implicit none


   type(MeganStateType) :: MeganState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title

   integer :: c ! Loop counter for emission Cats
   integer :: s ! Loop counter for emitted species

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=18), PARAMETER :: configFile ='CATChem_config.yml'

   thisLoc = 'test_megan -> at read CATChem_Config.yml'
   errMsg = ''
   rc = CC_SUCCESS

   write(*,*) '   CCCCC      A     TTTTTTT   CCCCC  H'
   write(*,*) '  C          A A       T     C       H       CCCC   EEEE   M       M'
   write(*,*) '  C         AAAAA      T     C       HHHHH  C      E    E  M M   M M'
   write(*,*) '  C        A     A     T     C       H   H  C      E EE    M   M   M'
   write(*,*) '   CCCCC  A       A    T      CCCCC  H   H   CCCC   EEEEE  M       M'
   write(*,*) ''
   write(*,*) ''

   !----------------------------
   ! Test 1
   !----------------------------

   ! Read input file and initialize grid
   call cc_read_config(Config, GridState, EmisState, ChemState, rc)
   if (rc /= CC_success) then
      errMsg = 'Error reading configuration file: ' // TRIM( configFile )
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   endif
   title = 'Megan Test 1 | Read Config'
   write(*,*) 'title = ', title
   write(*,*) 'Config%megan_activate = ', Config%megan_activate
   write(*,*) 'Config%megan_CO2_Inhib_Opt = ', Config%megan_CO2_Inhib_Opt
   write(*,*) 'Config%megan_CO2_conc_ppm = ', Config%megan_CO2_conc_ppm
   write(*,*) 'EmisState%nCats = ', EmisState%nCats
   !write(*,*) 'EmisState%Cats = ', EmisState%Cats  !cannot write allocatable variables here; need to put in the subroutin at the bottom

   !call Emis_Allocate(GridState, EmisState, RC) !Not sure why cannot call it even if I have made it public and used EmisStateMod in this module
   if (EmisState%nCats > 0) then
      do c = 1, EmisState%nCats
         do s = 1, EmisState%Cats(c)%nSpecies
            ALLOCATE(EmisState%Cats(c)%Species(s)%Flux(GridState%number_of_levels), STAT=RC)
            if (RC /= CC_SUCCESS) then
               ErrMsg = 'Error allocating "EmisState%Cats%Species%Flux"!'
               call cc_emit_error(ErrMsg, RC, ThisLoc)
               stop 1  !!Note here is not 'return'
            endif
         end do
      end do
   end if

   !----------------------------
   ! Test 2
   !----------------------------
   ! Set number of dust species to zero for now (TODO: not sure whether should be used)
   !ChemState%nSpeciesSeaSalt = 0

   ! Meteorological State (get from 20190620 18:00:00 HEMCO log out)
   MetState%LAI= 3.9094312191009521
   allocate(MetState%PFT_16(16))
   MetState%PFT_16= (/0.00, 0.11120668053627014, 0.00, 0.00, 0.00, 0.00, 0.00, 0.35108909010887146, &
      0.00, 0.00, 0.00, 0.00, 0.00,0.18369837105274200,6.9862455129623413E-002,0.26875913143157959/)
   MetState%PMISOLAI= 3.90559316  ! needs to multiply PFTSUM
   MetState%Q_DIR_2= 368.02691650390625
   MetState%Q_DIFF_2= 55.688854217529297
   MetState%PARDR_LASTXDAYS= 69.8014755
   MetState%PARDF_LASTXDAYS= 37.4157791
   MetState%TS= 300.42892456054688
   MetState%T_LASTXDAYS= 294.497833
   MetState%T_LAST24H= 294.919861
   MetState%GWETROOT= 0.76981580257415771
   MetState%SUNCOS= 0.96700188446067026
   MetState%LAT=  38.00
   MetState%DOY= 171
   MetState%LocalHour= 12.00
   MetState%D_BTW_M=  1.00
   MetState%AEF_ISOP= 1.8055753025901623E-009
   MetState%AEF_MBOX= 4.9856540616165277E-013
   MetState%AEF_BPIN= 1.3127712426530056E-011
   MetState%AEF_CARE= 3.6563258621377127E-012
   MetState%AEF_LIMO= 6.3985420662872528E-012
   MetState%AEF_OCIM= 3.0705341449874024E-011
   MetState%AEF_SABI= 1.0054971341792413E-011

   title = "Megan Test 2 | Test isoprene first"
   Config%megan_activate = .TRUE.

   call cc_megan_init(Config, ChemState, EmisState, MeganState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_megan_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_megan_run(MetState, EmisState, DiagState, MeganState, ChemState, RC )
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_megan_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, MeganState, MetState, title)
   call assert(MeganState%TotalEmission > 0.0_fp, "Test Megan isoprene")
   MeganState%TotalEmission = 0.0_fp


contains

   subroutine print_info(Config_, MeganState_, MetState_, title_)
      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(MeganStateType), intent(in) :: MeganState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Megan%activate = ', MeganState_%activate
      write(*,*) 'Megan%CatIndex = ', MeganState_%CatIndex
      write(*,*) 'MeganState%CO2Inhib = ', MeganState_%CO2Inhib
      write(*,*) 'MeganState%CO2conc  = ', MeganState_%CO2conc
      write(*,*) 'MetState%LAI =', MetState_%LAI
      write(*,*) 'MetState%DOY =', MetState_%DOY
      write(*,*) 'MetState%AEF_ISOP =', MetState_%AEF_ISOP
      write(*,*) 'MetState%PFT_16 =', MetState_%PFT_16
      write(*,*) 'MeganState%MeganSpeciesName=', MeganState_%MeganSpeciesName
      write(*,*) 'MeganState%EmissionPerSpecies=', MeganState_%EmissionPerSpecies
      write(*,*) 'MeganState%TotalEmission = ', MeganState_%TotalEmission

   end subroutine print_info

end program test_dust
