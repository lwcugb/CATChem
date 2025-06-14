!> \file ccpr_scheme_ginoux_mod.F90
!! \brief Ginoux windblown dust emission scheme
!! \ingroup catchem_dust_process
!!
!! \author Barry Baker
!! \date 05/2024
!!
!! This module implements the Ginoux dust emission scheme for calculating
!! windblown dust emissions in the CATChem atmospheric chemistry model.
!!
!! \details
!! The Ginoux scheme calculates dust emissions based on surface wind speed,
!! soil moisture, and surface characteristics. It provides size-resolved
!! dust emission fluxes for integration with atmospheric transport models.
!!
!! \section ginoux_reference Reference
!! Ginoux, P., "Sources and distributions of dust aerosols simulated with the GOCART model",
!! Journal of Geophysical Research, vol. 106, no. D17, pp. 20,255–20,273, 2001.
!! doi:10.1029/2000JD000053.
!!!! \brief Contains the Ginoux windblown dust emission scheme
!!
!! Ginoux, P., “Sources and distributions of dust aerosols simulated with the GOCART model”,
!! <Journal of Geophysical Research, vol. 106, no. D17, pp. 20, 255–20, 273, 2001.
!! doi:10.1029/2000JD000053.
!!
!! \author Barry baker
!! \date 05/2024
!! \ingroup catchem_dust_process
!!!>

module CCPr_Scheme_Ginoux_Mod

   implicit none

   private

   public :: CCPr_Scheme_Ginoux

contains

   !>
   !! \brief Calculates the dust emission flux in ug m-2 s-1
   !!
   !! \param DSOILTYPE Dominant soil type index
   !! \param SSM Sediment supply map (0-1)
   !! \param TSKIN Skin temperature [K]
   !! \param FROCEAN Fraction of ocean
   !! \param FRSNO Fraction of snow
   !! \param airden Air density [kg/m3]
   !! \param u10m 10m zonal wind speed [m/s]
   !! \param v10m 10m meridional wind speed [m/s]
   !! \param gwettop Top level volumetric soil moisture [m3/m3]
   !! \param AlphaScaleFactor Alpha scaling parameter
   !! \param EffectiveRadius Effective particle radius array
   !! \param DustDensity Dust particle density array
   !! \param TotalEmission Total emission rate [ug/m2/s]
   !! \param EmissionPerSpecies Emission rate per species [ug/m2/s]
   !! \param RC Return code
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine CCPr_Scheme_Ginoux(DSOILTYPE,          &
      SSM,                &
      TSKIN,              &
      FROCEAN,            &
      FRSNO,              &
      airden,             &
      u10m,               &
      v10m,               &
      gwettop,            &
      AlphaScaleFactor,   &
      EffectiveRadius,    &
      DustDensity,        &
      TotalEmission,      &
      EmissionPerSpecies, &
      RC)

      ! Uses
      Use CCPr_Dust_Common_Mod
      use precision_mod, only : fp, ZERO

      implicit none

      ! Arguments
      integer, intent(in) :: DSOILTYPE         ! Dominant Soil Type
      real(fp), intent(in) :: SSM              ! Sediment Supply Map (0-1)
      real(fp), intent(in) :: TSKIN            ! Skin Temperature [K]
      real(fp), intent(in) :: FROCEAN          ! Fraction of Ocean
      real(fp), intent(in) :: FRSNO            ! Fraction of Snow
      real(fp), dimension(:), intent(in) :: airden           ! Air density [kg/m3]
      real(fp), intent(in) :: u10m             ! 10m wind speed [m/s]
      real(fp), intent(in) :: v10m             ! 10m wind speed [m/s]
      real(fp), intent(in) :: gwettop          ! Top level volumetric soil moisture [m3/m3]
      real(fp), intent(in) :: AlphaScaleFactor ! Alpha Scaling Parameter
      real(fp), dimension(:), intent(in) :: EffectiveRadius  ! Effective Particle Radius
      real(fp), dimension(:), intent(in) :: DustDensity      ! Dust Particle Density

      real(fp), intent(inout)             :: TotalEmission ! Total Emission Rate [ug/m2/s]
      real(fp), dimension(:), intent(inout) :: EmissionPerSpecies ! Emission Rate per Species [ug/m2/s]
      integer, intent(out) :: RC

      ! Local Variables
      logical :: do_dust                               ! Enable Dust Calculation Flag
      integer :: n                                     ! Bin index
      integer :: nbins                                 ! number of dust bins
      real(fp) :: ginoux_scaling                       ! Ginoux scaling
      real(fp) :: u_thresh0                            ! Dry threshold wind speed [m/s]
      real(fp) :: u_thresh                             ! Moisture Corrected threshold wind speed [m/s]
      real(fp) :: w10m                                 ! 10m wind speed [m/s]

      ! Initialize
      RC = 0
      EmissionPerSpecies = ZERO
      TotalEmission = ZERO

      nbins = size(EffectiveRadius)

      !--------------------------------------------------------------------
      ! Don't do dust over certain criteria
      !--------------------------------------------------------------------
      do_dust = .true. ! Default value for all cases

      ! Don't do dust over bedrock, lava, or Permanent Ice (15, 16, 18)
      !----------------------------------------------------------------
      if (DSOILTYPE == 15 .or. DSOILTYPE == 16 .or. DSOILTYPE == 18) then
         do_dust = .false.
      endif

      ! Check for valid inputs from SSM and RDRAG
      !------------------------------------------
      if (SSM < 0.15_fp .or. SSM > 1.0_fp) then
         do_dust = .false.
      endif

      ! Don't do dust over frozen soil
      !--------------------------------
      if (TSKIN <= 273.15_fp) then
         do_dust = .false.
      endif

      if (do_dust) then

         ! get the scaling factor following Ginoux et al. (2001)
         ginoux_scaling = (1 - FROCEAN) * (1 - FRSNO) * SSM * AlphaScaleFactor

         do n = 1, nbins

            ! get threshold friction velocity following MB97
            call MB97_threshold_velocity(DustDensity(n), AIRDEN(1), EffectiveRadius(n), u_thresh0)

            ! get 10m mean wind speed
            w10m = sqrt(U10M ** 2 + V10M ** 2)

            if (GWETTOP .lt. 0.5_fp) then

               ! add the moisture correction following Ginoux et al. (2001)
               u_thresh = amax1(0., u_thresh0 * (1.2_fp + 0.2_fp*alog10(max(1.e-3_fp, GWETTOP))) )

               if (w10m .gt. u_thresh) then
                  EmissionPerSpecies(n) = ginoux_scaling * w10m ** 2 * max(0.,(w10m - u_thresh) )
               endif

            endif ! GWETTOP < 0.5
         enddo ! nbins

         TotalEmission = sum(EmissionPerSpecies)
         ! DiagState%dust_total_flux = TotalEmission

      endif ! do_dust

   end subroutine CCPr_Scheme_Ginoux

end module CCPr_Scheme_Ginoux_Mod
