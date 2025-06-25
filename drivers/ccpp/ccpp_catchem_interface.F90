!> \file ccpp_catchem_interface.F90
!! \brief CATCHEM-CCPP interface utilities module
!!
!! \details
!! This is the CCPP-Compliant wrapper for interfacing CATCHEM chemistry model with CCPP
!! framework. Handles data transformation and management between host model and
!! CATCHEM chemistry calculations. This module contains the init run and finalize functions and
!! subroutines that facilitate the data exchange and coordinate transformations
!! required for CATCHEM integration within CCPP-compliant models.
!!
!! \author Barry Baker and Wei Li, NOAA/OAR/ARL
!!
!! \date November 2024
!! \defgroup catchem_ccpp_group CATChem CCPP Interface
!! \ingroup catchem_ccpp_group
!!
!! \note This is part of the CATCHEM-CCPP interface layer that enables
!!       chemistry calculations within the CCPP framework
!!!>
module ccpp_catchem_interface

  use catchem_types, only: catchem_container_type !> CATChem container type

  !use physcons, only: g => con_g, pi => con_pi
  use machine, only: kind_phys
  !use catchem_config
  use CATChem
  use catchem_wrapper_utils

implicit none

private

public :: ccpp_catchem_interface_init, ccpp_catchem_interface_run, ccpp_catchem_interface_finalize

type(ConfigType) :: Config                          !> CATChem configuration object
type(DustStateType) :: DustState
type(SeaSaltStateType) :: SeaSaltState
type(DryDepStateType) :: DryDepState
type(catchem_container_type) :: CATChemStates       !> Container for all CATChem states

!   integer :: im    !> Number of horizontal points
!   integer :: kme   !> Number of vertical levels
!   integer :: nsoil !> Number of soil layers

contains



   !> \section arg_table_ccpp_catchem_interface_init Argument Table
   !! \htmlinclude ccpp_catchem_interface_init.html
   !!
   !! \brief Initialize the CATChem container
   !! \param[in] catchem_configfile_in Name of the CATChem configuration file
   !! \param[in] do_catchem_in Flag to enable CATChem calculations
   !! \param[in] export_catchem_diags_in Flag to export CATChem diagnostics
   !! \param[in] n_dbg_lines_in Number of debug lines
   !! \param[in] im Number of horizontal points
   !! \param[in] kme Number of vertical levels
   !! \param[in] nsoil Number of soil layers
   !! \param[out] errmsg Error message
   !! \param[out] errflg Error flag
   !!
   !! \note This subroutine initializes the CATChem container and reads the configuration
   !!       file to set up the chemistry model
   !!
   !! \ingroup catchem_ccpp_group
   !!!>
   subroutine ccpp_catchem_interface_init

   end subroutine ccpp_catchem_interface_init

  !> \brief Brief description of the subroutine
  !!
  !! \section arg_table_ccpp_catchem_interface_finalize Argument Table
  !! \htmlinclude ccpp_catchem_interface_finalize.html
  !!
  !>
  subroutine ccpp_catchem_interface_finalize(do_catchem, errmsg, errflg)

   implicit none

   ! Input parameters
   logical, intent(in) :: do_catchem
   ! Output parameters
   character(len=*), intent(out) :: errmsg
   integer, intent(out) :: errflg
   ! Local variables
   CHARACTER(LEN=255)    :: ThisLoc

   errmsg = ''
   errflg = 0
   ThisLoc = ' -> at catchem_interface_finalize (in drivers/ccpp/ccpp_catchem_inerface.F90)'

   if (.not. do_catchem) return

   call catchem_finalize (CATChemStates, DustState, SeaSaltState, DryDepState, errflg)
   if (errflg /= CC_SUCCESS) then
         ErrMsg = 'Error in catchem_finalize'
         call cc_emit_error(ErrMsg, errflg, ThisLoc )
         return
   end if

  end subroutine ccpp_catchem_interface_finalize

  !>
  !! This is the Configurable ATmospheric Chemistry (CATChem)
  !! This is the CATChem interface Module
  !! \section arg_table_ccpp_catchem_interface_run Argument Table
  !! \htmlinclude ccpp_catchem_interface_run.html
  !!
  !>
  subroutine ccpp_catchem_interface_run(im, kte, kme, garea, nsoil, nlndcat, nsoilcat, &
     lat, lon, & ! Grid information
     do_catchem, catchem_configfile_in, & ! CATChem Flag on
     dt, jdate, & ! Time information
     xcosz, &
     lwi, frlanduse, gvf, seaicefrac, oceanfrac, lakefrac, landfrac, & ! land water specific variables
     stype, vtype, snowdepth, frsnow, lai, frsoil, pores, resid, & ! land surface variables
     ustar, u10m, v10m, tskin, tsfco, ts, hf2d, lf2d, znt, prsfc, pblh, & !surface variables
     dswsfc, nirbmdi, nirdfdi, visbmdi, visdfdi, &  ! Radiation Fluxes
     sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, & ! surface albedo
     soilmoist, pr3d, phl3d, prl3d, tk3d, q3d, us3d, vs3d, rh, & ! 3D Variables
     delp, airden, pfl_lsan, pfl_isan, & ! 3D Variables
     rain_cpl, cldf, & ! cloud variables
     dust_in, & ! Emissions
     ntrac, ntchs, ntchm, chemarr_phys, chemarr, & ! Chemistry Variables
     errmsg, errflg)

     implicit none
     ! ARGUMENTS
     !----------

     ! CCPP Interface Variables
     !-------------------------

     ! Grid information
     integer, intent(in) :: im       ! horizontal loop extent
     integer, intent(in) :: kte      ! vertical layer dimension
     integer, intent(in) :: kme      ! vertical interface dimension = number of vertical layers + 1
     !integer, intent(in) :: ktau     ! current forecast iteration
     !integer, intent(in) :: tile_num ! index of cube sphere tile
     integer, intent(in) :: nsoil    !> vertical_dimension_of_soil
     integer, intent(in) :: nlndcat !> number_of_vegetation_categories
     integer, intent(in) :: nsoilcat !> number of soil categories
     real(kind_phys), dimension(im), intent(in) :: garea !> grid area (m^2) of each grid cell
     !real(kind_phys), dimension(im), intent(in) :: rlat  !> radian latitude
     !real(kind_phys), dimension(im), intent(in) :: rlon  !> radian longitude
     real(kind_phys), dimension(im), intent(in) :: lat   !> latitude
     real(kind_phys), dimension(im), intent(in) :: lon   !> longitude
     real(kind_phys), dimension(im), intent(in) :: xcosz !> cosine of solar zenith angle

     ! Time information
     !integer, intent(in) :: idat(8)  ! initialization date and time (in iso order)
     integer, intent(in) :: jdate(8)    !> julian date
     real(kind_phys), intent(in) :: dt      ! physics timestep (s)
     !real(kind_phys), intent(in) :: julian  ! forecast julian day (days)

     ! MPI information
     !type(MPI_comm), intent(in) :: mpicomm
     !integer, intent(in) :: mpirank
     !integer, intent(in) :: mpiroot

     ! Logicals
     logical, intent(in) :: do_catchem
     character(len=*), intent(in) :: catchem_configfile_in

     ! tracer information
     integer, intent(in) :: ntrac    ! total number of tracers
     integer, intent(in) :: ntchs      ! index of first chemical tracer in tracer concentration array
     integer, intent(in) :: ntchm      ! number of chemical tracers
     !integer, intent(in) :: ntaero      ! number of aerosol tracers
     real(kind_phys), dimension(im, kte, ntrac), intent(inout) :: chemarr_phys
     real(kind_phys), dimension(im, kte, ntrac), intent(inout) :: chemarr

     ! Emissions
     !real(kind_phys), dimension(im, kte, 3), intent(in) :: emi_in ! 3 should be temporary... need to replace with either an input namelist option or have it be a pointer with variable dimension
     !real(kind_phys), dimension(im, kte, ntrac), intent(in) :: kemit_in !> emission inputs
     real(kind_phys), dimension(im, 12, 5), intent(in) :: dust_in         !> dust emission inputs
     !real(kind_phys), dimension(im), intent(in) :: pert_scale_anthro !> anthropogenic emission perturbation scale factor
     !real(kind_phys), dimension(im), intent(in) :: pert_scale_plume !> plume emission perturbation scale factor
     !real(kind_phys), dimension(im), intent(in) :: pert_scale_dust !> dust emission perturbation scale factor
     !real(kind_phys), dimension(im), intent(in) :: emis_amp_anthro !> anthropogenic emission amplitude
     !real(kind_phys), dimension(im), intent(in) :: emis_amp_plume !> plume emission amplitude
     !real(kind_phys), dimension(im), intent(in) :: emis_amp_seas !> seasonal emission amplitude
     !real(kind_phys), dimension(im), intent(in) :: emis_amp_dust !> dust emission amplitude
     !integer(kind_phys), intent(in) :: biomass_burn_opt_in !> biomass burning option
     !integer(kind_phys), intent(in) :: plumerise_flag_in !> plume rise flag
     !integer(kind_phys), intent(in) :: plumerisefire_freq_in !> plume rise fire frequency
     !real(kind_phys), dimension(im, kte), intent(inout) :: drydep !> dry deposition
     !real(kind_phys), dimension(im, kte), intent(inout) :: wetdpl !> wet deposition
     !real(kind_phys), dimension(im, kte), intent(inout) :: abem !> aerosol backscatter extinction mass
     !real(kind_phys), dimension(im, kte), intent(in) :: ca_emis_plume !> plume emission
     !logical, dimension(im, kte), intent(in) :: ca_sgs_emis !> SGS emission
     !logical, dimension(im, kte), intent(in) :: ca_sgs !> SGS
     !real(kind_phys), dimension(im, kte), intent(in) :: ca_sgs_gbbepx_frp !> SGS GBBEPx FRP
     !real(kind_phys), dimension(im, kte), intent(in) :: fire_GBBEPx !> GBBEPx fire emissions
     !real(kind_phys), dimension(im, kte), intent(in) :: fire_MODIS !> MODIS fire emissions
     !real(kind_phys), dimension(im, kte), intent(in) :: ebu !> EBU
     !real(kind_phys), dimension(im, kte), intent(in) :: sppt_wts !> SPPT weights
     !logical, intent(in) :: do_sppt_emis !> SPPT emission flag

     ! land surface information
     integer, dimension(im), intent(in)                :: lwi         !> sea land ice mask (sea = 0, land = 1, ice = 2)
     integer, dimension(im), intent(in)                :: stype       !> soil type
     integer, dimension(im), intent(in)                :: vtype       !> vegatation type

     real(kind_phys), dimension(im, nlndcat), intent(in) :: frlanduse     !> fraction of each land surface category
     real(kind_phys), dimension(im, nsoilcat), intent(in) :: frsoil       !> fraction of each soil type
     real(kind_phys), dimension(30), intent(in)        :: pores           !> porosity of each soil type
     real(kind_phys), dimension(30), intent(in)        :: resid           !> residual water content of each soil type
     real(kind_phys), dimension(im), intent(in)        :: seaicefrac      !> fractin of ice cover over ocean
     !real(kind_phys), dimension(im), intent(in)        :: landicefrac     !> fractin of ice cover over land
     real(kind_phys), dimension(im), intent(in)        :: oceanfrac       !> fraction of ocean cover
     real(kind_phys), dimension(im), intent(in)        :: frsnow        !> fraction of snow cover over land
     real(kind_phys), dimension(im), intent(in)        :: lakefrac        !> fraction of lake cover
     real(kind_phys), dimension(im), intent(in)        :: landfrac        !> fraction of land cover
     real(kind_phys), dimension(im), intent(in)        :: gvf             !> green vegetative fraction
     real(kind_phys), dimension(im), intent(in)        :: lai             !> leaf area index

     real(kind_phys), dimension(im, nsoil), intent(in) :: soilmoist     !> volumetric fraction of soil moisture for lsm
     !real(kind_phys), dimension(im, nsoil), intent(in) :: soiltemp      !> soil temperature (K)
     real(kind_phys), dimension(im), intent(in)        :: snowdepth     !> water equivalent snow depth (mm)
     real(kind_phys), dimension(im), intent(in)        :: prsfc          !> pressure at the surface (Pa)
     !real(kind_phys), dimension(im), intent(in)        :: prslp         !> sea level pressure (Pa)
     real(kind_phys), dimension(im), intent(in)        :: pblh          !> PBL Thickness determined by the PBL scheme (m)
     !real(kind_phys), dimension(im), intent(in)        :: kpbl          !> PBL level
     !real(kind_phys), dimension(im), intent(in)        :: hpbl_thetav   !> PBL Height based on modified parcel method (m)
     real(kind_phys), dimension(im), intent(in)        :: u10m          !> 10 m wind speed (m/s)
     real(kind_phys), dimension(im), intent(in)        :: v10m          !> 10 m wind speed (m/s)
     real(kind_phys), dimension(im), intent(in)        :: ustar         !> friction velocity (m/s)
     !real(kind_phys), dimension(im), intent(in)        :: psim          !> Monin-Obukhov similarity parameter for momentum at 10m
     !real(kind_phys), dimension(im), intent(in)        :: psih          !> Monin-Obukhov similarity parameter for heat at 10m
     real(kind_phys), dimension(im), intent(in)        :: tskin         !> skin temperature (K)
     real(kind_phys), dimension(im), intent(in)        :: tsfco         !> sea surface temperature (K)
     !real(kind_phys), dimension(im), intent(in)        :: t2m           !> 2 m temperature (K)
     real(kind_phys), dimension(im), intent(in)        :: ts            !> surface temperature (K)
     !real(kind_phys), dimension(im), intent(in)        :: dpt2m         !> 2 m dew point temperature (K)
     real(kind_phys), dimension(im), intent(in)        :: hf2d          !> Sensible heat flux (W m-2)
     real(kind_phys), dimension(im), intent(in)        :: lf2d          !> Latent heat flux (W m-2)
     real(kind_phys), dimension(im), intent(in)        :: znt           !> surface roughness length in (cm)
     real(kind_phys), dimension(im), intent(in)        :: dswsfc        !> downward short wave flux (W m-2)
     !real(kind_phys), dimension(im), intent(in)        :: recmol        !> one over obukhov length (m-1)
     real(kind_phys), dimension(im), intent(in)        :: sfc_alb_nir_dir    !> surface near-infrared direct albedo
     real(kind_phys), dimension(im), intent(in)        :: sfc_alb_nir_dif    !> surface near-infrared diffuse albedo
     real(kind_phys), dimension(im), intent(in)        :: sfc_alb_uvvis_dir    !> surface ultraviolet-visible direct albedo
     real(kind_phys), dimension(im), intent(in)        :: sfc_alb_uvvis_dif    !> surface ultraviolet-visible diffuse albedo
     real(kind_phys), dimension(im), intent(in)        :: nirbmdi      !> surface near-infrared beam shortwave radiation (W/m2)
     real(kind_phys), dimension(im), intent(in)        :: nirdfdi      !> surface near-infrared diffuse shortwave radiation (W/m2)
     real(kind_phys), dimension(im), intent(in)        :: visbmdi      !> surface visible beam shortwave radiation (W/m2)
     real(kind_phys), dimension(im), intent(in)        :: visdfdi      !> surface visible diffuse shortwave radiation (W/m2)

     real(kind_phys), dimension(im, kme), intent(in) :: pr3d            !> air pressure at model layer interfaces (Pa)
     real(kind_phys), dimension(im, kte), intent(in) :: prl3d           !> pressure at the model level (Pa)
     !real(kind_phys), dimension(im, kme), intent(in) :: pr3d_dry        !> dry air pressure at model layer interfaces (Pa)
     !real(kind_phys), dimension(im, kte), intent(in) :: prl3d_dry       !> dry air pressure at the model level (Pa)
     real(kind_phys), dimension(im, kte), intent(in) :: delp            !> air pressure thickness at the model level (Pa)
     !real(kind_phys), dimension(im, kte), intent(in) :: delp_dry        !> dry air pressure thickness at the model level (Pa)
     !real(kind_phys), dimension(im, kme), intent(in) :: ph3d            !> geopotential at the model level interfaces (m2 s-2)
     real(kind_phys), dimension(im, kte), intent(in) :: phl3d           !> geopotential at the model layer (m2 s-2)
     real(kind_phys), dimension(im, kte), intent(in) :: tk3d            !> temperature at the model level (K)
     real(kind_phys), dimension(im, kte), intent(in) :: us3d            !> zonal wind at the model level (m/s)
     real(kind_phys), dimension(im, kte), intent(in) :: vs3d            !> meridional wind at the model level (m/s)
     real(kind_phys), dimension(im, kte), intent(in) :: q3d             !> specific humidity at the model level (kg/kg)
     real(kind_phys), dimension(im, kte), intent(in) :: airden         !> dry air density (kg/m3)
     !real(kind_phys), dimension(im, kte), intent(in) :: w               !> lagrangian_tendency_of_air_pressure
     !real(kind_phys), dimension(im, kte), intent(in) :: exch            !> atmospheric heat diffusivity
     real(kind_phys), dimension(im, kte), intent(in) :: rh              !> relative humidity
     real(kind_phys), dimension(im, kte), intent(in) :: pfl_lsan   !>  liquid flux from large scale precipitation (kg/m2/s)
     real(kind_phys), dimension(im, kte), intent(in) :: pfl_isan   !>  ice flux from large scale precipitation (kg/m2/s)

     ! precipitation information
     real(kind_phys), dimension(im), intent(in)        :: rain_cpl        !> total rain at this time step (m)
     !real(kind_phys), dimension(im), intent(in)        :: rainc_cpl       !> convective rain at this time step (m)
     real(kind_phys), dimension(im), intent(in)        :: cldf            !> total cloud fraction
     !real(kind_phys), dimension(im, kte), intent(in)   :: dqdt            !> instantaneous_water_vapor_specific_humidity_tendency_due_to_convection
     !integer, intent(in)                               :: chem_conv_tr_in !> catchem convective transport option

     ! Radiation
     !integer, intent(in) :: aer_ra_feedback_in  !> catchem aer radiation feedback option
     !integer, intent(in) :: aer_ra_frq_in       !> catchem_aer_ra_frq

     ! Output
     !-------
     character(len=*), intent(out) :: errmsg
     integer, intent(out) :: errflg
     ! Local
     !------
     CHARACTER(LEN=255)    :: ThisLoc
     logical, save :: firsttime = .true. !> first time flag

     errmsg = ''
     errflg = 0
     ThisLoc = ' -> at catchem_interface_run (in drivers/ccpp/ccpp_catchem_interface.F90)'


     if (.not. do_catchem) return

     if (firsttime) then
        firsttime = .false.
        !Initialize CATChem container
        !It is not put in ccpp_catchem_interface_init because horizontal_loop_extent cannot be used in interface init
        call catchem_init(Config, CATChemStates, im, catchem_configfile_in, errflg, errmsg, &
        DustState, SeaSaltState, DryDepState)
        if (errflg /= CC_SUCCESS) then
          ErrMsg = 'Error in catchem_init'
          call cc_emit_error(ErrMsg, errflg, ThisLoc )
          return
        end if
     end if

     !give tracer from CCPP to CATChem
     call ccpp_to_cc(im, kte, ntchs, ntchm, CATChemStates%ChemState, chemarr_phys)

     ! Fill MetState Arrays
     call transform_ccpp_to_catchem(im, kme, kte, nsoil, nlndcat, nsoilcat, lat, lon, &      ! Grid Information
                                        dt, jdate, garea, &  ! Grid Information
                                        lwi, &  ! Model Options
                                        tk3d, q3d, pr3d, pr3d, prl3d, rh, &  ! Meteorological Variables; TODO: not using dry air pressre and delp for now
                                        us3d, vs3d, delp, phl3d, &  ! Meteorological Variables
                                        u10m, v10m, tskin, tsfco, prsfc, ts, rain_cpl, &  ! Meteorological Variables
                                        cldf, airden, delp, &
                                        pfl_lsan, pfl_isan, &  ! precipitation variables
                                        snowdepth, vtype, stype, soilmoist, &  ! Surface Variables
                                        pblh, xcosz, &  ! Surface Variables
                                        ustar, hf2d, lf2d, &  ! Near-Surface Meteorology
                                        frsnow, gvf, lai, frlanduse, frsoil, pores, resid, &  ! Surface Variables
                                        znt, landfrac, oceanfrac, lakefrac, seaicefrac, &
                                        dswsfc, nirbmdi, nirdfdi, visbmdi, visdfdi, &  ! Radiation Fluxes
                                        sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif,&  ! surface albedo
                                        dust_in, & !Emissions
                                        CATChemStates%MetState, & ! CATChem States
                                        errmsg, errflg)
      if (errflg /= CC_SUCCESS) then
         call cc_emit_error(ErrMsg, errflg, ThisLoc )
         return
      end if

     ! Run CATChem
     call catchem_run(im, CATChemStates, DustState, SeaSaltState, DryDepState, errflg)
     if (errflg /= CC_SUCCESS) then
         ErrMsg = 'Error in catchem_run'
         call cc_emit_error(ErrMsg, errflg, ThisLoc ) !TODO: consider rename the subroutine
         return
     end if

     ! give chemical tracer back to CCPP
     call cc_to_ccpp(im, kte, ntchs, ntchm, CATChemStates%ChemState, chemarr)

   end subroutine ccpp_catchem_interface_run


end module ccpp_catchem_interface
