!> \file catchem_wrapper_utils.F90
!> \brief CATCHEM-CCPP interface utilities module
!>
!> \details
!> Provides wrapper utilities for interfacing CATCHEM chemistry model with CCPP
!> framework. Handles data transformation and management between host model and
!> CATCHEM chemistry calculations.
!>
!> \author Barry Baker and Wei Li
!>
!> \date 11/2024
!>
!> \ingroup catchem_ccpp_group
!!!>
module catchem_wrapper_utils

    use CATChem, only: ConfigType, MetStateType, ChemStateType, GridStateType, &
                      EmisStateType, DiagStateType, cc_read_config, &
                      cc_allocate_metstate, cc_allocate_chemstate, &
                      cc_allocate_diagstate, cc_allocate_emisstate, &
                      cc_emit_error, cc_check_var, cc_emit_warning, &
                      cc_init_met, cc_init_diag, cc_init_emis, cc_init_process, &
                      cc_run_process, cc_finalize_process, &
                      CC_SUCCESS
    use catchem_types, only: catchem_container_type
    use machine, only: kind_phys


    implicit none

    public :: transform_ccpp_to_catchem
    public :: cc_to_ccpp, ccpp_to_cc
    public :: catchem_init, catchem_run, catchem_finalize

    ! CCPP data structure
    !type(ccpp_t), save, target :: cdata

    !temorary for testing
    !integer, parameter :: kind_phys = 4

contains


    !> Allocates memory for the entire catchem container
    !!
    !! \param[in]    config         CATChem configuration
    !! \param[inout] catchem_states Container for all CATChem states
    !! \param[in]    im            Number of horizontal points
    !! \param[in]    kme           Number of vertical levels
    !! \param[in]    nsoil         Number of soil levels
    !! \param[out]   errflg        Error flag (0=success, non-zero=failure)
    !! \param[out]   errmsg        Error message if errflg is non-zero
    !!
    !! \ingroup catchem_ccpp_group
    !!!>

    subroutine catchem_init(config, catchem_states, im, config_file, errflg, errmsg, &
                            DustState, SeaSaltState, DryDepState)
        use CATChem, only: DustStateType, SeaSaltStateType, DryDepStateType
        implicit none

        type(ConfigType), intent(inout) :: config
        type(catchem_container_type), intent(inout) :: catchem_states
        type(DustStateType), intent(inout) :: DustState
        type(SeaSaltStateType), intent(inout) :: SeaSaltState
        type(DryDepStateType), intent(inout) :: DryDepState
        integer, intent(in) :: im
        !integer, intent(in) :: kme
        !integer, intent(in) :: nsoil
        !integer, intent(in) :: nLandTypes
        character(len=*), intent(in) :: config_file
        integer, intent(out) :: errflg
        character(len=*), intent(out) :: errmsg

        !local variables
        type(EmisStateType) :: EmisState !temporary type to be copied to catchem_states
        type(ChemStateType) :: ChemState
        type(MetStateType)  :: MetState
        type(DiagStateType) :: DiagState
        type(GridStateType) :: GridState
        integer :: i

        ! Initialize
        errflg = 0
        errmsg = ''

        call cc_read_config(config, &
                          GridState,    &
                          EmisState, &
                          ChemState, &
                          errflg, &
                          config_file)
        if (errflg /= 0) then
            errmsg = 'Error reading CATChem configuration file: '//trim(config_file)
        end if

        ! Initialize the States and processes
        call cc_init_met(GridState, MetState, errflg)
        call cc_init_emis(GridState, EmisState, ChemState, errflg)
        !call cc_init_chem(GridState, ChemState, errflg) !TODO: ChemState already allocated in cc_read_config file
        call cc_init_diag( config, DiagState, ChemState, errflg)
        call cc_init_process(config, ChemState, EmisState, DustState, SeaSaltState, DryDepState, errflg)

        ! Allocate State arrays in the container
        if (.not. allocated(catchem_states%MetState)) then
            allocate(catchem_states%MetState(im), stat=errflg)
            if (check_allocation_error('catchem_states%MetState', errflg, errmsg)) return

            ! Initialize each MetState using CC API
            do i = 1, im
                catchem_states%MetState(i) = MetState  !TODO: is there a better way? Duplicated information
            end do
        endif

        if (.not. allocated(catchem_states%EmisState)) then
            allocate(catchem_states%EmisState(im), stat=errflg)
            if (check_allocation_error('catchem_states%EmisState', errflg, errmsg)) return

            ! Initialize each EmisState using CC API
            do i = 1, im
                catchem_states%EmisState(i) = EmisState
            end do
        endif

        if (.not. allocated(catchem_states%ChemState)) then
            allocate(catchem_states%ChemState(im), stat=errflg)
            if (check_allocation_error('catchem_states%ChemState', errflg, errmsg)) return

            ! Initialize each ChemState using CC API
            do i = 1, im
                catchem_states%ChemState(i) = ChemState
            end do
        endif

        if (.not. allocated(catchem_states%DiagState)) then
            allocate(catchem_states%DiagState(im), stat=errflg)
            if (check_allocation_error('catchem_states%DiagState', errflg, errmsg)) return

            ! Initialize each DiagState using CC API
            do i = 1, im
                catchem_states%DiagState(i) = DiagState
            end do
        endif


    !end subroutine allocate_catchem_container
    end subroutine catchem_init

    !> Transforms CCPP meteorological arrays into CATChem meteorological and chemistry states
    !!
    !! \param[in] im         Number of horizontal grid points
    !! \param[in] kme        Number of vertical levels
    !! \param[in] temp       3D temperature field (K)
    !! \param[in] spechum    3D specific humidity field (kg/kg)
    !! \param[in] pfull      3D full level pressure (Pa)
    !! \param[in] phalf      3D half level pressure (Pa)
    !! \param[in] u          3D zonal wind (m/s)
    !! \param[in] v          3D meridional wind (m/s)
    !! \param[in] delp       3D pressure thickness (Pa)
    !! \param[in] zh         3D geopotential height (m)
    !! \param[in] kh         3D vertical diffusivity (m2/s)
    !! \param[in] prsl       3D layer mean pressure (Pa)
    !! \param[in] prslk      3D Exner function
    !! \param[in] u10m       10m u wind (m/s)
    !! \param[in] v10m       10m v wind (m/s)
    !! \param[in] tskin      Surface skin temperature (K)
    !! \param[in] tsfco      Sea surface temperature (K)
    !! \param[in] ps         Surface pressure (Pa)
    !! \param[in] precip     Precipitation rate (kg/m2/s)
    !! \param[in] slmsk      Land-sea mask (1=land,0=sea,2=ice)
    !! \param[in] snowh      Snow depth (m)
    !! \param[in] vegtype    Vegetation type
    !! \param[in] soiltyp    Soil type
    !! \param[in] hf         Sensible heat flux (W/m2)
    !! \param[in] ust        Friction velocity (m/s)
    !! \param[in] zpbl       PBL height (m)
    !! \param[in] coszen     Cosine of solar zenith angle
    !! \param[in] albedo     Surface albedo
    !! \param[in] emis       Surface emissivity
    !! \param[in] ustar      Friction velocity (m/s)
    !! \param[in] shflx      Surface sensible heat flux (W/m2)
    !! \param[in] lhflx      Surface latent heat flux (W/m2)
    !! \param[in] snowc      Snow cover fraction (0-1)
    !! \param[in] vegfrac    Green vegetation fraction (0-1)
    !! \param[in] swdn       Surface downward SW radiation (W/m2)
    !! \param[in] swup       Surface upward SW radiation (W/m2)
    !! \param[in] lwdn       Surface downward LW radiation (W/m2)
    !! \param[in] lwup       Surface upward LW radiation (W/m2)
    !! \param[in] swdnc      Clear-sky surface downward SW radiation (W/m2)
    !! \param[in] swupc      Clear-sky surface upward SW radiation (W/m2)
    !! \param[in] lwdnc      Clear-sky surface downward LW radiation (W/m2)
    !! \param[in] lwupc      Clear-sky surface upward LW radiation (W/m2)
    !! \param[inout] MetState   CATChem meteorology state
    !! \param[inout] ChemState  CATChem chemistry state
    !! \param[out] errmsg     Error message
    !! \param[out] errflg     Error flag
    !!
    !! \ingroup catchem_ccpp_group
    !!!>
    subroutine transform_ccpp_to_catchem(im, kme, kte, nsoil, nlndcat, nsoilcat, lat, lon, &      ! Grid Information
                                        dt, jdate, garea, &  ! Grid Information
                                        lwi, &  ! Model Options
                                        temp, spechum, pfull,pfull_wet, phalf, rh, &  ! Meteorological Variables
                                        u, v, delp, zh, &  ! Meteorological Variables
                                        u10m, v10m, tskin, tsfco, ps, ts, precip, &  ! Meteorological Variables
                                        cldf, airden, delp_dry, &
                                        pfl_lsan, pfl_isan, &  ! precipitation variables
                                        snowh, vegtype, soiltyp, soilmoist, &  ! Surface Variables
                                        zpbl, coszen, &  ! Surface Variables
                                        ustar, shflx, lhflx, &  ! Near-Surface Meteorology
                                        snowc, vegfrac, lai, frlanduse, frsoil, pores, resid, &  ! Surface Variables
                                        z0, landfrac, oceanfrac, lakefrac, seaicefrac, &
                                        swdn, nirbmdi, nirdfdi, visbmdi, visdfdi, &  ! Radiation Fluxes
                                        sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif,&  ! surface albedo
                                        dust_in, &  ! Dust Emission Inputs
                                        MetState, &  ! CATChem MetStates
                                        errmsg, errflg)  ! Error Handling

      use CATChem, only: MetStateType, ChemStateType, DiagStateType, EmisStateType
      implicit none

      !! Transform CCPP meteorological arrays to CATChem states
      integer,  intent(in)    :: im              !> number of horizontal points
      integer,  intent(in)    :: kme             !> number of vertical interfaces
      integer,  intent(in)    :: kte             !> number of vertical levels
      !integer, intent(in)     :: nVert           !> number of vertical levels
      !integer, intent(in)     :: nVertInterface  !> number of vertical interfaces
      integer, intent(in)     :: nsoil           !> number of soil levels
      integer, intent(in)     :: nlndcat         !> number of land categories
      integer, intent(in)     :: nsoilcat        !> number of soil categories
      !real(kind_phys), intent(in) :: rlat(:)     !> latitude (radian)
      !real(kind_phys), intent(in) :: rlon(:)     !> longitude (radian)
      real(kind_phys), intent(in) :: lat(:)      !> latitude (degrees)
      real(kind_phys), intent(in) :: lon(:)      !> longitude (degrees)
      !integer, intent(in)     :: ktau            !> number of timestep
      real(kind_phys), intent(in) :: dt          !> physics timestep (s)
      integer, intent(in)     :: jdate(8)        !> current forecast date and time
      !integer, intent(in)     :: tile_num        !> index of cubed sphere tile
      real(kind_phys), intent(in) :: garea(:)    !> grid cell area
      !integer, intent(in) :: aero_rad_freq_opt   !> catchem aer radiation frequency
      !integer, intent(in) :: aero_feedback_opt   !> catchem aerosol radiation feedback option
      !integer, intent(in) :: plmrise_freq_opt    !> catchem plmrise frequency option
      integer, intent(in) :: lwi(:)               !> sea/land/ice=0/1/2 index
      !integer, intent(in) :: dluse(:)            !> land use index
      real(kind_phys), intent(in) :: pores(30)   !> maximum soil moisture content
      real(kind_phys), intent(in) :: resid(30)   !> minimum soil moisture content

      ! 3D/Layer Variables (dim(:,:))
      real(kind_phys), intent(in) :: temp(:,:)       !> temperature (K)
      real(kind_phys), intent(in) :: spechum(:,:)    !> specific humidity (kg/kg)
      real(kind_phys), intent(in) :: pfull(:,:)      !> full level pressure of dry air (Pa)
      real(kind_phys), intent(in) :: pfull_wet(:,:)  !> full level pressure of moist air (Pa)
      real(kind_phys), intent(in) :: phalf(:,:)      !> half level pressure of dry air(Pa)
      real(kind_phys), intent(in) :: u(:,:)          !> zonal wind (m/s)
      real(kind_phys), intent(in) :: v(:,:)          !> meridional wind (m/s)
      real(kind_phys), intent(in) :: delp(:,:)       !> pressure thickness (Pa)
      real(kind_phys), intent(in) :: delp_dry(:,:)   !> dry air pressure thickness (Pa)
      real(kind_phys), intent(in) :: zh(:,:)         !> geopotential (m2/s2)
      !real(kind_phys), intent(in) :: kh(:,:)         !> vertical diffusivity (m2/s)
      !real(kind_phys), intent(in) :: prsl(:,:)       !> layer mean pressure (Pa)
      !real(kind_phys), intent(in) :: prslk(:,:)      !> Exner function
      real(kind_phys), intent(in) :: rh(:,:)         !> relative humidity [fraction]
      real(kind_phys), intent(in) :: pfl_lsan(:,:)   !>  liquid flux from large scale precipitation (kg/m2/s)
      real(kind_phys), intent(in) :: pfl_isan(:,:)   !>  ice flux from large scale precipitation (kg/m2/s)
      real(kind_phys), intent(in) :: airden(:,:)   !> dry air density [kg/m3]

      ! Surface Variables (dim(:))
      real(kind_phys), intent(in) :: ps(:)           !> surface pressure (Pa)
      real(kind_phys), intent(in) :: tskin(:)        !> skin temperature (K)
      real(kind_phys), intent(in) :: tsfco(:)        !> sea surface temperature (K)
      real(kind_phys), intent(in) :: ts(:)           !> surface temperature (K)
      !real(kind_phys), intent(in) :: slmsk(:)        !> land-sea mask (1=land,0=sea,2=ice)
      real(kind_phys), intent(in) :: snowh(:)        !> snow depth (m)
      integer,         intent(in) :: vegtype(:)      !> vegetation type
      real(kind_phys), intent(in) :: lai(:)          !> leaf area index (m2/m2)
      real(kind_phys), intent(in) :: z0(:)           !> surface roughness length (cm)
      integer,         intent(in) :: soiltyp(:)      !> soil type
      real(kind_phys), intent(in) :: soilmoist(:,:)    !> soil moisture (m3/m3)
      real(kind_phys), intent(in) :: snowc(:)        !> snow cover fraction over land (0-1)
      real(kind_phys), intent(in) :: vegfrac(:)      !> green vegetation fraction (0-1)
      real(kind_phys), intent(in) :: frlanduse(:,:)  !> fraction of each land type
      real(kind_phys), intent(in) :: frsoil(:,:)     !> fraction of each soil type
      real(kind_phys), intent(in) :: landfrac(:)     !> fraction of land
      real(kind_phys), intent(in) :: oceanfrac(:)    !> fraction of ocean
      real(kind_phys), intent(in) :: lakefrac(:)     !> fraction of lake
      !real(kind_phys), intent(in) :: landicefrac(:)  !> fraction of land ice
      real(kind_phys), intent(in) :: seaicefrac(:)   !> fraction of sea ice

      ! Near-Surface Meteorology (dim(:))
      real(kind_phys), intent(in) :: u10m(:)         !> 10m u wind (m/s)
      real(kind_phys), intent(in) :: v10m(:)         !> 10m v wind (m/s)
      real(kind_phys), intent(in) :: ustar(:)        !> friction velocity (m/s)
      real(kind_phys), intent(in) :: zpbl(:)         !> PBL height (m)

      ! Surface Fluxes (dim(:))
      !real(kind_phys), intent(in) :: hf(:)           !> sensible heat flux (W/m2)
      real(kind_phys), intent(in) :: shflx(:)        !> surface sensible heat flux (W/m2)
      real(kind_phys), intent(in) :: lhflx(:)        !> surface latent heat flux (W/m2)
      real(kind_phys), intent(in) :: precip(:)       !> precipitation rate (kg/m2/s)

      ! Radiation Properties (dim(:))
      real(kind_phys), intent(in) :: coszen(:)            !> cosine of solar zenith angle
      real(kind_phys), intent(in) :: cldf(:)              !> fraction of grid box area in which updrafts occur
      real(kind_phys), intent(in) :: sfc_alb_nir_dir(:)   !> surface near-infrared direct albedo
      real(kind_phys), intent(in) :: sfc_alb_nir_dif(:)   !> surface near-infrared diffuse albedo
      real(kind_phys), intent(in) :: sfc_alb_uvvis_dir(:) !> surface visible + uv direct albedo
      real(kind_phys), intent(in) :: sfc_alb_uvvis_dif(:) !> surface visible + uv diffuse albedo
      !real(kind_phys), intent(in) :: emis(:)              !> surface emissivity

      ! Radiation Fluxes (dim(:))
      real(kind_phys), intent(in) :: swdn(:)         !> downward shortwave radiation at surface (W/m2)
      !real(kind_phys), intent(in) :: swup(:)         !> upward shortwave radiation at surface (W/m2)
      !real(kind_phys), intent(in) :: lwdn(:)         !> downward longwave radiation at surface (W/m2)
      !real(kind_phys), intent(in) :: lwup(:)         !> upward longwave radiation at surface (W/m2)
      !real(kind_phys), intent(in) :: swdnc(:)        !> clear-sky downward shortwave radiation (W/m2)
      !real(kind_phys), intent(in) :: swupc(:)        !> clear-sky upward shortwave radiation (W/m2)
      !real(kind_phys), intent(in) :: lwdnc(:)        !> clear-sky downward longwave radiation (W/m2)2)
      !real(kind_phys), intent(in) :: lwupc(:)        !> clear-sky upward longwave radiation (W/m2)
      real(kind_phys), intent(in) :: nirbmdi(:)      !> surface near-infrared beam shortwave radiation (W/m2)
      real(kind_phys), intent(in) :: nirdfdi(:)      !> surface near-infrared diffuse shortwave radiation (W/m2)
      real(kind_phys), intent(in) :: visbmdi(:)      !> surface visible + uv beam shortwave radiation (W/m2)
      real(kind_phys), intent(in) :: visdfdi(:)      !> surface visible + uv diffuse shortwave radiation (W/m2)


      ! Emissions
      !real(kind_phys), intent(in) :: emi_in(:)         !> emissions
      real(kind_phys), dimension(im, 12, 5), intent(in) :: dust_in         !> dust emission inputs
      ! CATChem States
      type(MetStateType),  intent(inout) :: MetState(:)    !> CATChem meteorology state
      !type(ChemStateType), intent(inout) :: ChemState(:)   !> CATChem chemistry state
      !type(EmisStateType), intent(inout) :: EmisState(:)   !> CATChem emission state
      !type(DiagStateType), intent(inout) :: DiagState(:)   !> CATChem diagnostic state

      ! Error handling
      character(len=*), intent(out) :: errmsg    !> error message
      integer,          intent(out) :: errflg    !> error flag

      ! Local variables
      integer :: i, k, k_rev, month_now
      real :: FRLAND_NOSNO_NOICE, FRWATER, FRICE, FRSNO
      real :: frac_temp
      ! Initialize error handling
      errmsg = ''
      errflg = 0

      ! Check dimensions
      if (size(MetState) /= im) then
        errmsg = 'MetState dimension mismatch'
        errflg = 1
        return
      endif
      
      ! Transform data for each horizontal point
      horiz: do i = 1, im
        
        ! Verify vertical dimension
        if (MetState(i)%nLEVS /= kte) then
            errmsg = 'Vertical dimension mismatch'
            errflg = 1
            return
        endif

        !Grid state variables
        MetState(i)%LAT = lat(i)
        MetState(i)%LON = lon(i)
        MetState(i)%AREA_M2 = garea(i)
        MetState(i)%TSTEP = dt
        !TODO: assume the 8 dimension of jdate is something like: (/ 2025, 5, 9, 14, 30, 0, 0, 0 /)
        MetState(i)%YMD = jdate(1) * 10000 + jdate(2) * 100 + jdate(3) !YYYYMMDD
        MetState(i)%HMS = jdate(4) * 10000 + jdate(5) * 100 + jdate(6) !HHMMSS
        month_now = jdate(2)

        ! 3D/Layer center Variables; no need to reverse vertical index for CATChem
        !TODO: make sure layer index starts from 1 not 0
        k_rev = 1 !kte
        vert1: do k = 1, kte
            MetState(i)%T(k_rev)        = temp(i,k)   !tk3d
            MetState(i)%SPHU(k_rev)     = spechum(i,k) /1000.0 !q3d (kg/kg in CC --> g/kg in MetState)
            MetState(i)%PMID_DRY(k_rev) = phalf(i,k) !prl3d_dry
            MetState(i)%U(k_rev)        = u(i,k) !us3d TODO: why use new_state in META???
            MetState(i)%V(k_rev)        = v(i,k)  !vs3d
            MetState(i)%DELP(k_rev)     = delp(i,k)
            MetState(i)%DELP_DRY(k_rev) = delp_dry(i,k)
            MetState(i)%ZMID(k_rev)     = zh(i, k)/9.80665 !phl3d
            !MetState(i)%kh(k_rev)       = kh(i,k) !TODO: MetState does not have kh defined yet.
            !MetState(i)%prsl(k_rev)    = prsl(i,k) !assume the same as phalf
            !MetState(i)%prslk(k_rev)   = prslk(i,k) !TODO: MetState does not have prslk defined yet.
            MetState(i)%AIRDEN(k_rev)  = phalf(i,k)/(287.04*temp(i,k)*(1.+.608*spechum(i,k)))   !< Dry air density [kg/m3]
            MetState(i)%RH(k_rev)      = 50.0 !rh(i,k) * 100 !TODO: not found
            MetState(i)%BXHEIGHT(k_rev)  = delp (i,k) / 9.80665

            k_rev = k_rev + 1
            
        end do vert1

        ! 3D/Layer interface Variables; no need to reverse vertical index for CATChem
        k_rev = 1 !kme
        vert2: do k = 1, kme
            MetState(i)%PEDGE_DRY(k_rev) = pfull(i,k) !pr3d
            MetState(i)%PEDGE(k_rev)     = pfull_wet(i,k) !pfull_wet(i,k) !pr3d TODO: not found
            ! if ( k == kme ) then
            !     MetState(i)%PFLLSAN(k_rev)   = 0.0 !TODO: vertical_layer_dimension not found vertical_interface_dimension
            !     MetState(i)%PFILSAN(k_rev)   = 0.0       ! so giving zero for the first layer for now
            ! else
            !     MetState(i)%PFLLSAN(k_rev)   = 0.0 !pfl_lsan(i,k)  !TODO: not found (need to put after mp_thomson)
            !     MetState(i)%PFILSAN(k_rev)   = 0.0 !pfl_isan(i,k)
            ! end if

            !!TODO: seems no need to use if else above
            MetState(i)%PFLLSAN(k_rev)   = 0.0 !pfl_lsan(i,k)  !TODO: not found (need to put after mp_thomson)
            MetState(i)%PFILSAN(k_rev)   = 0.0 !pfl_isan(i,k)

            k_rev = k_rev + 1
        end do vert2

        ! Surface Variables
        MetState(i)%PS       = ps(i)  !prsfc
        MetState(i)%TS       = ts(i)  !ts
        MetState(i)%TSKIN    = tskin(i) !tskin
        MetState(i)%LWI      = lwi(i) ! TODO: slmsk same as lwi and slmsk can be deleted?
        MetState(i)%SNODP    = snowh(i) !snowdepth
        MetState(i)%DLUSE    = vegtype(i) !vtype
        MetState(i)%DSOILTYPE = soiltyp(i) !stype
        MetState(i)%FRVEG    = vegfrac(i) !gvf
        MetState(i)%LAI      = lai(i) !lai
        MetState(i)%FRLANDUSE = frlanduse(i,:) !< Fraction of each land type
        MetState(i)%FRSOIL = frsoil(i,:) !< Fraction of each soil type
        MetState(i)%FRLAND = landfrac(i)       !< Fraction of land [1]
        !https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html
        if (nlndcat == 20) then !modified NoahMP 20 categories; only choice for now (ivegsrc = 1)
            MetState%LUCNAME = 'NOAH'
        else
            errmsg = 'number of land category is not 20 and not supported'
            errflg = 1
            return
        endif
        MetState(i)%FRLAI = frlanduse(i,:) * lai(i)
        MetState(i)%FRLAI(15:17) = 0.0 !manually give index 15(snow and ice), 16(barren), 17(water) zeros
        MetState(i)%FROCEAN = oceanfrac(i)        !< Fraction of ocean [1]
        MetState(i)%FRSEAICE = seaicefrac(i)  !fice TODO
        !MetState(i)%FRLANDIC = landicefrac(i) !TODO: seems no land ice fraction in CCPP, so give zero for now
        MetState(i)%FRLANDIC = 0.0
        MetState(i)%FRLAKE = lakefrac(i)          !< Fraction of lake [1]
        MetState(i)%FRSNO = snowc(i) !snowfrac
        ! Land without snow or land ice (TODO: here we assume FRLAND is without ice already)
        FRLAND_NOSNO_NOICE = MetState(i)%FRLAND - MetState(i)%FRSNO !- MetState(i)%FRLANDIC
        ! Water without sea ice
        FRWATER = MetState(i)%FRLAKE + MetState(i)%FROCEAN - MetState(i)%FRSEAICE
        ! Land and sea ice
        FRICE = MetState(i)%FRLANDIC + MetState(i)%FRSEAICE
        ! Land snow
        FRSNO = MetState(i)%FRSNO
        ! Set IsLand, IsWater, IsIce, IsSnow based on max fractional area (adoped from GEOS-Chem)
        !MetState(i)%IsLand  = (FRLAND_NOSNO_NOICE > MAX(FRWATER, FRICE, FRSNO))
        MetState(i)%IsWater = (FRWATER > MAX(FRLAND_NOSNO_NOICE, FRICE, FRSNO))
        !MetState(i)%IsIce   = (FRICE > MAX(FRLAND_NOSNO_NOICE, FRWATER, FRSNO))
        !define it based on LWI since we do not have land ice fraction in ccpp
        MetState(i)%IsLand = (lwi(i) == 1) ! land
        MetState(i)%IsIce   = (lwi(i) == 2) ! ice
        MetState(i)%IsSnow  = (FRSNO > MAX(FRLAND_NOSNO_NOICE, FRWATER, FRICE))
        !soil moisture
        MetState(i)%SOILM = soilmoist(i, :)
        MetState(i)%GWETROOT = 0.0
        MetState(i)%GWETTOP  = 0.0
        if ( (nsoilcat == 16) .or. (nsoilcat == 19) ) then ! 16 or 19 STATSGO soil type; only choices for now
            do k =1, nsoilcat
                if ( soiltyp(i) == k .and. (pores(k) - resid(k)) > 0.0 ) then
                    !here we asume nsoil is the lowest layer and index=1 is the top layer
                    MetState(i)%GWETROOT = (soilmoist(i,nsoil) - resid(k)) / (pores(k) - resid(k))
                    MetState(i)%GWETTOP = (soilmoist(i,1) - resid(k)) / (pores(k) - resid(k))
                    MetState(i)%GWETROOT = max(0.0, min(1.0, MetState(i)%GWETROOT))
                    MetState(i)%GWETTOP = max(0.0, min(1.0, MetState(i)%GWETTOP))
                end if
            end do
        else
            errmsg = 'number of soil category is not 16 or 19 and not supported'
            errflg = 1
            return
        end if

        !TODO: need to turn them on for real simulations
        MetState(i)%CLAYFRAC = 0.5 ! dust_in(i, month_now, 1)        !< Fraction of clay [1]
        MetState(i)%RDRAG = 0.8 ! dust_in(i, month_now, 2) !Drag Partition [1]
        MetState(i)%SANDFRAC = 0.5 ! dust_in(i, month_now, 3)       !< Fraction of sand [1]
        MetState(i)%SSM = 0.5 ! dust_in(i, month_now, 4) ! Sediment Supply Map [1]
        MetState(i)%USTAR_THRESHOLD = 0.1 ! dust_in(i, month_now, 5) !< Threshold friction velocity [m/s]

        ! Near-Surface Meteorology
        MetState(i)%U10M     = u10m(i)
        MetState(i)%V10M     = v10m(i)
        MetState(i)%PBLH     = zpbl(i) !pblh
        MetState(i)%USTAR    = ustar(i)
        MetState(i)%CLDFRC   = cldf(i)
        MetState(i)%SST      = tsfco(i) 

        ! Surface Fluxes
        MetState(i)%Z0 = z0(i) *1.0e-2 !cm --> m
        MetState(i)%Z0H = z0(i) *1.0e-2 !< roughness height, for heat (thermal roughness) [m] TODO: give the same value as z0 for now
        MetState(i)%HFLUX    = 0.5 !shflx(i) !hf2d !TODO: need to turn on cplflx in GFS to make it work
        MetState(i)%EFLUX    = 0.0 !lhflx(i) !lf2d
        MetState(i)%PRECLSC  = 0.0 !precip(i) !rain_cplchm
        !calculate the Monin-Obukhov length following a simplified method as used in GOCART
        !in the ESMF/GOCART_GridComp/O3_GridComp/O3_GridCompMod.F90 file
        ! check MAPL/Python/MAPL/constants.py for MAPL_CP and MAPL_GRAV values
        if (abs(MetState(i)%HFLUX) > 1.00e-32) then
            MetState(i)%OBK  = - MetState(i)%AIRDEN(1)  * 1004.68307 * temp(i,kte) * ustar(i)**3. / (0.40 * 9.80665 * MetState(i)%HFLUX)
        else
            MetState(i)%OBK = 1.00e+05
        end if

        ! Radiation Properties
        MetState(i)%SUNCOS   = coszen(i) !xcosz
        frac_temp = visbmdi(i) / (visbmdi(i) + visdfdi(i))
        MetState(i)%ALBD_VIS   = frac_temp * sfc_alb_uvvis_dir(i) + (1.0 - frac_temp) * sfc_alb_uvvis_dif(i) !albedo_vis TODO: should be vis and uv
        frac_temp = nirbmdi(i) / (nirbmdi(i) + nirdfdi(i))
        MetState(i)%ALBD_NIR   = frac_temp * sfc_alb_nir_dir(i) + (1.0 - frac_temp) * sfc_alb_nir_dif(i) !albedo_nir

        ! Radiation Fluxes
        MetState(i)%PARDR = 0.47 * ( nirbmdi(i) + visbmdi(i) ) !TODO: empirical factor from papers
        MetState(i)%PARDF = 0.57 * ( nirdfdi(i) + visdfdi(i) )
        MetState(i)%SWGDN     = swdn(i) !dswsfc

        !MetState(i)%swup     = swup(i) !not definded in MetState
        !MetState(i)%lwdn     = lwdn(i)
        !MetState(i)%lwup     = lwup(i)
        !MetState(i)%swdnc    = swdnc(i)
        !MetState(i)%swupc    = swupc(i)
        !MetState(i)%lwdnc    = lwdnc(i)
        !MetState(i)%lwupc    = lwupc(i)

      end do horiz

    end subroutine transform_ccpp_to_catchem



    subroutine catchem_run (im, catchem_states, DustState, SeaSaltState, DryDepState, errflg)
        use CATChem, only: DustStateType, SeaSaltStateType, DryDepStateType
        implicit none

        integer, intent(in) :: im
        type(catchem_container_type), intent(inout) :: catchem_states
        type(DustStateType), intent(inout) :: DustState
        type(SeaSaltStateType), intent(inout) :: SeaSaltState
        type(DryDepStateType), intent(inout) :: DryDepState
        integer, intent(inout) :: errflg

        integer:: i ! loop index

        ! Error handling
        !---------------
        CHARACTER(LEN=255)    :: ErrMsg
        CHARACTER(LEN=255)    :: ThisLoc
        ThisLoc = ' -> at catchem_process_run (in drivers/ccpp/catchem_wrapper_utils.F90)'

        do i = 1, im
            call cc_run_process(catchem_states%MetState(i), catchem_states%DiagState(i),   &
                                catchem_states%ChemState(i), catchem_states%EmisState(i),  &
                                DustState, SeaSaltState, DryDepState, errflg)
            if (errflg /= CC_SUCCESS) then
                ErrMsg = 'Error in catchem_process_run'
                call cc_emit_error(ErrMsg, errflg, ThisLoc ) !TODO: consider rename the subroutine
                return
            end if
        end do

    end subroutine catchem_run


    subroutine catchem_finalize (catchem_states, DustState, SeaSaltState, DryDepState, errflg)
        use CATChem, only: DustStateType, SeaSaltStateType, DryDepStateType
        implicit none

        type(catchem_container_type), intent(inout) :: catchem_states
        type(DustStateType), intent(inout) :: DustState
        type(SeaSaltStateType), intent(inout) :: SeaSaltState
        type(DryDepStateType), intent(inout) :: DryDepState
        integer, intent(inout) :: errflg

        ! Error handling
        !---------------
        CHARACTER(LEN=255)    :: ErrMsg
        CHARACTER(LEN=255)    :: ThisLoc
        ThisLoc = ' -> at catchem_finalize (in drivers/ccpp/catchem_wrapper_utils.F90)'

        call cc_finalize_process(DustState, SeaSaltState, DryDepState, errflg)
        if (errflg /= CC_SUCCESS) then
            ErrMsg = 'Error in cc_finalize_process'
            call cc_emit_error(ErrMsg, errflg, ThisLoc )
            return
        end if

        ! Finalize CATChem states
        if (allocated(catchem_states%MetState)) then
            deallocate(catchem_states%MetState, stat=errflg)
            if (errflg /= CC_SUCCESS) then
                ErrMsg = 'Error in deallocate catchem_states%MetState'
                call cc_emit_error(ErrMsg, errflg, ThisLoc )
                return
            end if
        end if

        if (allocated(catchem_states%EmisState)) then
            deallocate(catchem_states%EmisState, stat=errflg)
            if (errflg /= CC_SUCCESS) then
                ErrMsg = 'Error in deallocate catchem_states%EmisState'
                call cc_emit_error(ErrMsg, errflg, ThisLoc )
                return
            end if
        end if

        if (allocated(catchem_states%ChemState)) then
            deallocate(catchem_states%ChemState, stat=errflg)
            if (errflg /= CC_SUCCESS) then
                ErrMsg = 'Error in deallocate catchem_states%ChemState'
                call cc_emit_error(ErrMsg, errflg, ThisLoc )
                return
            end if
        end if

        if (allocated(catchem_states%DiagState)) then
            deallocate(catchem_states%DiagState, stat=errflg)
            if (errflg /= CC_SUCCESS) then
                ErrMsg = 'Error in deallocate catchem_states%DiagState'
                call cc_emit_error(ErrMsg, errflg, ThisLoc )
                return
            end if
        end if

    end subroutine catchem_finalize


    subroutine cc_to_ccpp(im, kte, ind_start, ntchm, ChemState, chemarr)
        use CATChem, only: ChemStateType
        implicit none

        integer, intent(in) :: im
        integer, intent(in) :: kte
        integer, intent(in) :: ind_start
        integer, intent(in) :: ntchm
        type(ChemStateType), intent(in) :: ChemState(:)
        real(kind_phys), intent(inout) :: chemarr(:,:,:)

        !local variables
        integer :: i, n
        real(kind_phys) :: con(kte), unit_conv !, con_rev(kte)

        !TODO: here we assume the order of the chemical tracers in ChemState is the same as in chemarr
        do i = 1, im
            do n = 1, ntchm
                con = ChemState(i)%ChemSpecies(n)%conc
                !seems non eed to reverse the vertical index for CATChem
                !con_rev = con(size(con):1:-1)
                !unite conversion
                if(ChemState(i)%ChemSpecies(n)%is_gas) then
                    unit_conv = 1.0e3 * ChemState(i)%ChemSpecies(n)%mw_g / 28.9644  ! convert from ppm to ug/kg for gases
                else
                    unit_conv = 1.00 ! convert from ug/kg to ug/kg for aerosols
                end if
                chemarr(i,:,n+ind_start-1) = con * unit_conv
            end do
            !write(*,*) 'cc_to_ccpp i=:', i, n, any(con_rev > 0.0), unit_conv
        end do

    end subroutine cc_to_ccpp

    subroutine ccpp_to_cc(im, kte, ind_start, ntchm, ChemState, chemarr)
        use CATChem, only: ChemStateType
        implicit none

        integer, intent(in) :: im
        integer, intent(in) :: kte
        integer, intent(in) :: ind_start
        integer, intent(in) :: ntchm
        type(ChemStateType), intent(inout) :: ChemState(:)
        real(kind_phys), intent(in) :: chemarr(:,:,:)

        !local variables
        integer :: i, n
        real(kind_phys) :: con(kte),  unit_conv !, con_rev(kte)

        !TODO: here we assume the order of the chemical tracers in ChemState is the same as in chemarr
        do i = 1, im
            do n = 1, ntchm
                con = chemarr(i,:,n+ind_start-1)
                !seems non eed to reverse the vertical index for CATChem
                !con_rev = con(size(con):1:-1)
                !unite conversion
                if(ChemState(i)%ChemSpecies(n)%is_gas) then
                    unit_conv = 28.9644  / ChemState(i)%ChemSpecies(n)%mw_g * 1.0e-3  ! convert from ug/kg to ppm for gases
                else
                    unit_conv = 1.00 ! convert from ug/kg to ug/kg for aerosols
                end if
                ChemState(i)%ChemSpecies(n)%conc = con * unit_conv
            end do
        end do

    end subroutine ccpp_to_cc


    !> Checks for allocation errors and sets error message
    !!
    !! \param[in] state_name Name of the state being allocated
    !! \param[in] errflg     Error flag from allocation
    !! \param[out] errmsg    Output error message if allocation failed
    !! \return has_error     True if allocation error occurred
    !! \ingroup catchem_ccpp_group
    !!!>
    function check_allocation_error(state_name, errflg, errmsg) result(has_error)
        implicit none

        character(len=*), intent(in) :: state_name
        integer, intent(in) :: errflg
        character(len=*), intent(out) :: errmsg
        logical :: has_error

        has_error = (errflg /= 0)
        if (has_error) then
            errmsg = 'Error allocating '//trim(state_name)//' - catchem_wrapper_init'
        end if
    end function check_allocation_error

end module catchem_wrapper_utils
