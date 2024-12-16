!> \brief This module contains the routines for the SOFIEV plumerise scheme
!! \ingroup catchem_plumerise_process
!!!>
module ccpr_scheme_sofiev_mod

   use precision_mod, only: rae
   implicit none

   PUBLIC :: CCPr_Sofiev_Plmrise

contains

   subroutine find_height_index(GEOHGT, hgt, index)
      ! simple function to find the height index given the target height in
      ! a profile given the geopotential height.

      ! Arguments
      real, intent(in) :: GEOHGT (:)
      real, intent(in) :: hgt
      integer, intent(out) :: index

      ! Local Variables:
      integer :: i
      real, dimension(:), allocatable :: diffs

      allocate(diffs(size(geohgt)))

      do i = 1, size(GEOHGT)
         diffs(i) = abs(geohgt(i) - hgt)
      end do

      index = minloc(diffs, 1)
      ! function
      ! index = FINDLOC(geohgt, hgt)
      !index = FINDLOC(geohgt, hgt)

      return

   end subroutine find_height_index

   subroutine distribute_conc_linear(ZF, plmHGT, base_emiss, emis)

      ! Arguments
      real, intent(in) :: ZF(:)      ! Height of full layer (m)

      real, intent(in) :: plmHGT  ! plume rise height
      real, intent(in) :: base_emiss ! emission from file
      real, intent(out) :: emis(:) ! output column emission

      ! Local Variables:
      real :: hgt_prev     ! place holder for previous height index
      real :: dz           ! layer thickness
      real :: total_height ! height of of plume
      real :: column_frac  ! fraction of layer in total plume height
      integer :: plmHGT_index ! plume rise height index in Z
      integer :: z            ! loop index


      ! initialize emission array to zero
      emis = 0.

      ! find the plume height index
      call find_height_index(ZF, plmHGT, plmHGT_index)

      ! get the total height
      total_height = ZF(plmHGT_index)

      hgt_prev = 0.
      do z=1,  plmHGT_index
         dz = ZF(z) - hgt_prev
         column_frac = dz / total_height
         emis(z) = column_frac * base_emiss
         hgt_prev = zf(z)
      end do

      return

   end subroutine distribute_conc_linear

   !> \brief Distribute emissions around a given height
   !!
   !! This function will distribute the emission around a given height
   !! with a total plumePct, % of emission, in the given height and the rest split between the layer obove and below
   !!
   !! \param[in] ZF Height of full layer (m)
   !! \param[in] plmHGT plume rise height
   !! \param[in] plmPct_opt % of emission
   !! \param[in] underplumePct % of emission
   !! \param[in] base_emiss emission from file
   !! \param[out] emis output column emission
   !!
   !! \ingroup catchem_plumerise_process
   !!
   !!!>
   subroutine CCPr_Sofiev_Distribute(ZF, plmHGT, EFrac, RC, plmPct_opt, SfcPlmPct_opt)

      ! Arguments
      real, intent(in) :: ZF(:)      ! Height of full layer (m)
      real, intent(in) :: plmHGT     ! plume rise height
      real, intent(out) :: EFrac(:)   ! Fractional Column Emission

      integer, intent(inout) :: RC

      real, optional, intent(in) :: plmPct_opt
      real, optional, intent(in) :: SfcPlmPct_opt

      ! Local Variables:
      ! real :: hgt_prev                         ! place holder for previous height index
      ! real :: dz                               ! layer thickness
      real :: total_height                     ! height of of plume
      ! integer :: lev0                          ! bottom level index index
      ! integer :: lev1                          ! upper level index
      ! real :: column_frac                      ! fraction of layer in total plume height
      ! real, allocatable :: LayPlmPct(:)  ! fraction of layer in total plume height
      real :: plmPct
      real :: SfcPlmPct
      integer :: plmHGT_index

      ! initialize emission array to zero
      EFrac = 0.
      plmPct = 0.
      SfcPlmPct = 0.

      ! find the plume height index
      call find_height_index(ZF, plmHGT, plmHGT_index)

      ! get the total height
      total_height = ZF(plmHGT_index)

      if (present(plmPct_opt)) then
         plmPct = plmPct_opt
      else
         plmPct = .75
      endif

      if (present(SfcPlmPct_opt)) then
         SfcPlmPct = SfcPlmPct_opt
      else
         SfcPlmPct = .05
      endif

      if (plmPct + SfcPlmPct > 1) then
         RC = 0
         return
      endif

      ! Emissions are emitted at plmHGT_index +- 1 level and at the sfc
      if (plmPct < 1. .and. SfcPlmPct > 0.) then
         if (plmHGT_index == 1) then
            EFrac(1) = plmPct
            EFrac(2) = (1. - plmPct)
         else if (plmHGT_index == 2) then
            EFrac(plmHGT_index) = plmPct
            EFrac(plmHGT_index-1) = (1. - plmPct/2)
            EFrac(plmHGT_index+1) = (1. - plmPct/2)
         else
            EFrac(plmHGT_index) = plmPct
            EFrac(plmHGT_index-1) = (1. - plmPct/2. - SfcPlmPct)
            EFrac(plmHGT_index+1) =  (1. - plmPct/2. - SfcPlmPct)
            EFrac(1) = SfcPlmPct
         endif
      else if (rae(plmPct, 1.) .and. rae(SfcPlmPct, 0.)) then
         EFrac(plmHGT_index) = 1.
      endif

      return

   end subroutine CCPr_Sofiev_Distribute

   !> \brief Find the plume height index
   !!
   !! \ingroup catchem_plumerise_process
   !!
   !! \param Z
   !! \param T
   !! \param P
   !! \param PBLH
   !! \param psfc
   !! \param frp
   !! \param plmHGT
   !! \param EmisFrac
   !! \param RC
   !!
   !! \ingroup catchem_plumerise_process
   !! \author B. Baker
   !! \date 2024
   !!!>
   subroutine CCPr_Sofiev_Plmrise(Z,        &
      T,        &
      P,        &
      PBLH,     &
      psfc,     &
      frp,      &
      plmHGT,   &
      EmisFrac, &
      RC)
      real, intent(in) :: z(:)
      real, intent(in) :: T(:)
      real, intent(in) :: P(:)
      real, intent(in) :: PBLH
      real, intent(in) :: psfc
      real, intent(in) :: frp
      real, intent(out) :: plmHGT
      real, intent(out) :: EmisFrac(:)
      integer, intent(out) :: RC

      integer :: pblx2_index
      real :: PT1, PT2, LayerDepth

      !find the index of 2x the pbl
      call find_height_index(z, pblh * 2, pblx2_index)
      pblx2_index = MAX(pblx2_index, 2)
      ! now get the plume height
      PT1 = T(pblx2_index - 1) * (psfc / 100. / p(pblx2_index - 1))**(2./7.)
      PT2 = T(pblx2_index) * (psfc / 100. / p(pblx2_index))**(2./7.)
      LayerDepth = z(pblx2_index) - z(pblx2_index - 1)

      call plumeRiseSofiev(PT1, PT2, LayerDepth, frp, PBLH, plmHGT)

      call CCPr_Sofiev_Distribute(Z, plmHGT, EmisFrac, RC)
      if (RC /= 0) then
         write(*,*) 'Error in CCPr_Sofiev_PlmriseHgt'
         stop
      endif

      return
   end subroutine CCPr_Sofiev_Plmrise


   !> \brief Sofiev Plume Rise Algorithm
   !!
   !! Ref: M. Sofiev et al., Evaluation of the smoke-injection
   !! height from wild-land fires using remote sensing data.
   !!    Atmos. Chem. Phys., 12, 1995-2006, 2012.
   !!
   !! \param[in] PT1 Potential Temperature right below PBL height
   !! \param[in] PT2 Potential Temperature right above PBL height
   !! \param[in] laydepth Layer depth
   !! \param[in] frp Fire radiative power
   !! \param[in] pblh PBL height
   !! \param[out] Hp Plume height
   !!
   !! \ingroup catchem_plumerise_process
   !! \author Barry Baker
   !! \date 07/2024
   !!!>
   subroutine plumeRiseSofiev(PT1, PT2,laydepth,frp,pblh,Hp)

      !  This subroutine implements the Sofiev plume rise algorithm
      !  History: 09/16/2019: Prototype by Daniel Tong (DT)
      !           10/15/2019: bug fix based on feedback from M. Sofiev, DT
      !	         11/2020: parameterization options, Yunyao Li (YL)
      !           07/2024: implemented into CATChem, Barry Baker (BB)

      real Hp       !< plume height (m)
      real pblh     !< PBL height (m)
      real frp      !< fire radiative power (W)
      real NFT_sq   !< N square in Free Troposphere (@ z = 2pblh)
      real PT1, PT2 !< Potential Temperature right below and above PBL height
      real laydepth !< depth of the layer at the PBL height
      real grav     !< gravity

      real Pf0      !< reference fire power (W)
      real N0_sq    !< Brunt-Vaisala frequency (s-2)
      real alpha    !< part of ABL passed freely
      real beta     !< weights contribution of fire intensity
      real gama     !< power-law dependence on FRP
      real delta    !< dependence on stability in the FT

      ! ... Initial values.
      ! ... predefined values parameter set 3 to estimate whether hp higher
      ! than abl
      alpha     = 0.15
      beta      = 102
      gama      = 0.49
      delta     = 0

      Pf0       = 1000000.0
      N0_sq     = 0.00025
      grav      = 9.8

      ! ! ... calculate PT from T and P
      !       PT1 = T1 * (1000/P1)**0.286
      !       PT2 = T2 * (1000/P2)**0.286

      ! ... calculate Brunt-Vaisala frequency
      NFT_sq = grav/PT1*abs(PT1-PT2)/laydepth

      ! ... calculate first guess plume rise top height
      Hp = alpha*pblh + beta*(frp/Pf0)**gama * exp(-delta*NFT_sq/N0_sq)
      ! ... compare Hp with ABL
      if (Hp .lt. pblh) then
         alpha     = 0.24
         beta      = 170
         gama      = 0.35
         delta     = 0.6
         Hp = alpha*pblh + beta*(frp/Pf0)**gama*exp(-delta*NFT_sq/N0_sq)
      else
         alpha     = 0.93
         beta      = 298
         gama      = 0.13
         delta     = 0.7
         Hp = alpha*pblh + beta*(frp/Pf0)**gama*exp(-delta*NFT_sq/N0_sq)
      end if

      ! print *, "The height of fire plume (m) is: ", Hp

   end subroutine

end module ccpr_scheme_sofiev_mod
