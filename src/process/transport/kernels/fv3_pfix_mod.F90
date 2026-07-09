!***********************************************************************
!*                 GEOS-Chem / GMI pressure-fixer provenance
!*
!* This file is a FAITHFUL PORT (numerics unchanged) of the Philip
!* Cameron-Smith / LLNL "PJC" pressure fixer from the GEOS-Chem Global
!* Chemical Transport Model:
!*   https://github.com/geoschem/geos-chem  GeosCore/pjc_pfix_mod.F90
!*   (routines Do_Pjc_Pfix, Adjust_Press, Init_Press_Fix, Set_Press_Terms,
!*    Calc_Horiz_Mass_Flux, Calc_Divergence, Do_Divergence_Pole_Sum,
!*    Do_Press_Fix_Llnl, Average_Press_Poles, Calc_Advection_Factors,
!*    Init_Pjc_Pfix, Xpavg)
!*   Original code: Shian-Jiann Lin (DAO); Philip Cameron-Smith and
!*   John Tannahill (GMI project @ LLNL, 2003); Brendan Field and
!*   Bob Yantosca (Harvard, 5/8/03).  GEOS-Chem is distributed under a
!*   permissive (BSD-style) license; see the GEOS-Chem repository.
!*
! CATChem adaptation notes
! ------------------------
! The published GEOS-Chem module is tightly bound to State_Grid, Pressure_Mod
! (GET_AP/GET_BP), Error_Mod and module-level SAVEd geometry.  Here the SAME
! numerics are repackaged as a single self-contained driver, `do_pjc_pfix`,
! that receives everything it needs as plain array arguments.  No GEOS-Chem
! framework, no module state, no I/O.
!
! What the fixer does (unchanged from GEOS-Chem):
!   Given the archived A-grid winds (U,V) and the "before"/"after" true
!   surface pressures (P1,P2) over a dynamic step, it builds the horizontal
!   pressure mass-fluxes and adds a *barotropic* (column, distributed in the
!   vertical by the hybrid dSigma weight dbk) correction so that the
!   vertically-integrated horizontal mass divergence exactly equals the
!   prescribed surface-pressure tendency (P2-P1).  This is precisely the
!   consistency an offline transport model needs: after advection the
!   Lagrangian surface pressure closes onto the met surface pressure.
!
! Grid assumptions (identical to GEOS-Chem's pfix):
!   * regular global lon-lat grid, i = 1..nx PERIODIC in longitude,
!   * j = 1..ny running SOUTH pole (j=1) -> NORTH pole (j=ny);
!     caller supplies grid-box centre latitudes `clat(1:ny)` [rad],
!     monotonically increasing.
!   * enlarged polar cap with j1p = 3 (the two southern-most and two
!     northern-most latitude rows are treated as the polar cap; their
!     pressures are area-averaged, exactly as GEOS-Chem does).
!
! Vertical index convention:
!   Everything vertical here (dap,dbk,delpm,u,v,xmass,ymass) is carried in the
!   host's own layer order.  The only vertically-resolved operations are (a)
!   the per-level flux build and (b) the dbk-weighted distribution of the
!   column correction; the divergence is summed over all layers (order-
!   invariant).  So CATChem's surface-first order is used directly: pass
!   dap(k)=ap(k)-ap(k+1) >0 and dbk(k)=bp(k)-bp(k+1) >0 for surface-first
!   hybrid edges ap(1..nz+1),bp(1..nz+1).
!
! Precision: compiled inside CATChem_process_transport with
! -fdefault-real-8 -fdefault-double-8, so bare `real` == REAL(8), matching
! the GEOS-Chem REAL(fp)=REAL8 build.
!***********************************************************************

module fv3_pfix_mod

   implicit none
   private
   public :: do_pjc_pfix
   public :: pfix_correction

   ! Physical constants (GEOS-Chem PhysConstants: Re, PI).
   real, parameter :: PJC_PI = 3.14159265358979323846_8
   real, parameter :: PJC_RE = 6.371e+6_8    !< Earth radius [m]

contains

   !> \brief Philip Cameron-Smith / LLNL pressure fixer (driver).
   !!
   !! Faithful port of GEOS-Chem `Do_Pjc_Pfix` + `Adjust_Press` +
   !! `Init_Press_Fix` + `Do_Press_Fix_Llnl`.  Produces the pressure mass
   !! fluxes `xmass` (E-W) and `ymass` (N-S) AFTER the barotropic fix, so that
   !! the vertically-integrated divergence of (xmass,ymass) equals the met
   !! surface-pressure tendency (P2-P1) to machine precision (up to the global
   !! mean removed for whole-atmosphere mass conservation, and the polar-cap
   !! averaging).
   !!
   !! All pressures in Pa (any consistent unit works; the fixer is linear).
   !! Winds in m/s on the A-grid.  Areas in m^2.  `dt` in seconds.
   subroutine do_pjc_pfix(nx, ny, nz, area, clat, dap, dbk, dt, &
                          p1_in, p2_in, u, v, xmass, ymass, dps_ctm)
      integer, intent(in)  :: nx, ny, nz
      real,    intent(in)  :: area(nx, ny)     !< grid-box area [m^2]
      real,    intent(in)  :: clat(ny)         !< centre latitude [rad], S->N
      real,    intent(in)  :: dap(nz)          !< ap(k)-ap(k+1) [Pa] (>0)
      real,    intent(in)  :: dbk(nz)          !< bp(k)-bp(k+1) [-]  (>0)
      real,    intent(in)  :: dt               !< dynamic step [s]
      real,    intent(in)  :: p1_in(nx, ny)    !< true sfc P before step [Pa]
      real,    intent(in)  :: p2_in(nx, ny)    !< true sfc P after  step [Pa]
      real,    intent(in)  :: u(nx, ny, nz)    !< A-grid zonal wind [m/s]
      real,    intent(in)  :: v(nx, ny, nz)    !< A-grid meridional wind [m/s]
      real,    intent(out) :: xmass(nx, ny, nz)!< fixed E-W pressure flux
      real,    intent(out) :: ymass(nx, ny, nz)!< fixed N-S pressure flux
      !> Optional: vertically-integrated divergence of the FIXED fluxes
      !! [Pa]; on return this matches (P2-P1) after global-mean/pole averaging
      !! and is provided for verification/diagnostics.
      real, optional, intent(out) :: dps_ctm(nx, ny)

      ! Geometry (built each call; cost is O(ny), negligible vs advection).
      real :: rel_area(nx, ny)
      real :: elat(ny+1), sine(ny+1), cose(ny+1)
      real :: cosp(ny), gw(ny)
      real :: geofac(ny), geofac_pc
      real :: dp, dl, ri2, total_area

      ! Working pressures / fluxes.
      real :: p1(nx, ny), p2(nx, ny)
      real :: dps(nx, ny), dps_c(nx, ny)
      real :: delpm(nx, ny, nz)
      real :: dpi(nx, ny, nz)
      real :: dgpress

      integer :: i, j, k, j1p, j2p

      j1p = 3
      j2p = ny - 2
      ri2 = real(nx)
      dl  = 2.0_8 * PJC_PI / ri2
      dp  = PJC_PI / real(ny - 1)

      !---------------------------------------------------------------------
      ! Geometry  (GEOS-Chem Init_Pjc_Pfix + Calc_Advection_Factors)
      !---------------------------------------------------------------------
      total_area = sum(area)
      rel_area(:,:) = area(:,:) / total_area

      ! Latitude edges (S pole -> N pole) with sines/cosines.
      elat(1) = -0.5_8 * PJC_PI
      sine(1) = -1.0_8
      cose(1) =  0.0_8
      do j = 2, ny
         elat(j) = 0.5_8 * (clat(j-1) + clat(j))
         sine(j) = sin(elat(j))
         cose(j) = cos(elat(j))
      end do
      elat(ny+1) =  0.5_8 * PJC_PI
      sine(ny+1) =  1.0_8
      cose(ny+1) =  0.0_8

      ! gw = d(sin lat); cosp = gw/dlat (area-consistent cos at centres).
      gw(1)  = sine(2) - sine(1)
      cosp(1) = gw(1) / (2.0_8 * (elat(2) - elat(1)))
      do j = 2, ny-1
         gw(j)   = sine(j+1) - sine(j)
         cosp(j) = gw(j) / (elat(j+1) - elat(j))
      end do
      gw(ny)  = sine(ny+1) - sine(ny)
      cosp(ny) = gw(ny) / (2.0_8 * (elat(ny+1) - elat(ny)))

      ! Meridional geometrical factors (Calc_Advection_Factors).
      do j = 1, ny
         geofac(j) = dp / (2.0_8 * rel_area(1,j) * ri2)
      end do
      geofac_pc = dp / (2.0_8 * sum(rel_area(1,1:2)) * ri2)

      !---------------------------------------------------------------------
      ! Adjust_Press: remove global-mean pressure discrepancy (whole-
      ! atmosphere mass conservation), then average the polar caps.
      !---------------------------------------------------------------------
      p1(:,:) = p1_in(:,:)
      p2(:,:) = p2_in(:,:)

      dgpress = sum((p2(:,:) - p1(:,:)) * rel_area(:,:))
      p2(:,:) = p2(:,:) - dgpress

      call average_press_poles(nx, ny, rel_area, p1)
      call average_press_poles(nx, ny, rel_area, p2)

      ! Target surface-pressure tendency.
      dps(:,:) = p2(:,:) - p1(:,:)

      !---------------------------------------------------------------------
      ! Set_Press_Terms: layer thickness at the mid-step (used for fluxes).
      !---------------------------------------------------------------------
      do k = 1, nz
         delpm(:,:,k) = dap(k) + dbk(k) * 0.5_8 * (p1(:,:) + p2(:,:))
      end do

      !---------------------------------------------------------------------
      ! Calc_Horiz_Mass_Flux: raw pressure mass fluxes from the winds.
      !---------------------------------------------------------------------
      call calc_horiz_mass_flux(nx, ny, nz, dt, dl, dp, cosp, cose, &
                                delpm, u, v, xmass, ymass)

      !---------------------------------------------------------------------
      ! Calc_Divergence + vertical sum: CTM surface-pressure tendency from
      ! the raw winds.
      !---------------------------------------------------------------------
      call calc_divergence(nx, ny, nz, j1p, j2p, geofac_pc, geofac, &
                           xmass, ymass, dpi)
      do j = 1, ny
         do i = 1, nx
            dps_c(i,j) = sum(dpi(i,j,:))
         end do
      end do

      !---------------------------------------------------------------------
      ! Do_Press_Fix_Llnl: barotropic correction so the fixed-flux divergence
      ! matches dps (= P2-P1).  Overwrites xmass,ymass with the fixed fluxes.
      !---------------------------------------------------------------------
      call do_press_fix_llnl(nx, ny, nz, j1p, j2p, geofac_pc, geofac, &
                             dbk, dps, dps_c, rel_area, xmass, ymass)

      !---------------------------------------------------------------------
      ! Optional verification diagnostic: divergence of the FIXED fluxes.
      !---------------------------------------------------------------------
      if (present(dps_ctm)) then
         call calc_divergence(nx, ny, nz, j1p, j2p, geofac_pc, geofac, &
                              xmass, ymass, dpi)
         do j = 1, ny
            do i = 1, nx
               dps_ctm(i,j) = sum(dpi(i,j,:))
            end do
         end do
      end if

   end subroutine do_pjc_pfix

   !> \brief Barotropic pressure-fix TERMS for coupling to an external
   !!        (e.g. FV3) advection core.
   !!
   !! Same LLNL solver as `do_pjc_pfix`, but instead of computing the model's
   !! surface-pressure tendency from its own winds, the caller supplies the
   !! advection core's ALREADY vertically-integrated horizontal mass divergence
   !! `dps_ctm` [Pa] (the surface-pressure change the core's fluxes would
   !! produce over the step).  The routine returns the per-column correction
   !! terms that, when added to the core's mass fluxes (weighted by the hybrid
   !! dSigma dbk in the vertical), drive the total divergence to the prescribed
   !! `p2-p1`:
   !!    dmfx(i,j,k) = xcolmass_fix(i,j) * dbk(k) * area(i,j)
   !!    dmfy(i,j,k) = mmf(j)           * dbk(k) * cn
   !! where the FV3 layer-thickness tendency is (dmfx(i)-dmfx(i+1) +
   !! dmfy(j)-dmfy(j+1)) / area(i,j).  On a regular lon-lat grid this reproduces
   !! the native pfix correction exactly (cn/area(j) == geofac(j)).
   subroutine pfix_correction(nx, ny, area, p1_in, p2_in, dps_ctm, &
                              xcolmass_fix, mmf, cn)
      integer, intent(in)  :: nx, ny
      real,    intent(in)  :: area(nx, ny)      !< grid-box area [m^2]
      real,    intent(in)  :: p1_in(nx, ny)     !< sfc P before step [Pa]
      real,    intent(in)  :: p2_in(nx, ny)     !< sfc P after  step [Pa]
      real,    intent(in)  :: dps_ctm(nx, ny)   !< core's integrated divergence [Pa]
      real,    intent(out) :: xcolmass_fix(nx, ny) !< zonal column correction
      real,    intent(out) :: mmf(ny)           !< meridional column correction
      real,    intent(out) :: cn                !< meridional flux scale [m^2]

      real :: rel_area(nx, ny)
      real :: geofac(ny), geofac_pc
      real :: p1(nx, ny), p2(nx, ny)
      real :: dps(nx, ny), dpsc(nx, ny), ddps(nx, ny)
      real :: mmfd(ny), fxintegral(nx+1)
      real :: dp, ri2, total_area, dgpress, fxmean
      integer :: i, j, j1p, j2p

      j1p = 3
      j2p = ny - 2
      ri2 = real(nx)
      dp  = PJC_PI / real(ny - 1)

      !--- Geometry (meridional factors) ---
      total_area = sum(area)
      rel_area(:,:) = area(:,:) / total_area

      do j = 1, ny
         geofac(j) = dp / (2.0_8 * rel_area(1,j) * ri2)
      end do
      geofac_pc = dp / (2.0_8 * sum(rel_area(1,1:2)) * ri2)

      ! Meridional flux scale: cn/area(j) == geofac(j) on a regular lon-lat grid.
      cn = dp * total_area / (2.0_8 * ri2)

      !--- Prescribed target tendency (Adjust_Press + pole averaging) ---
      p1(:,:) = p1_in(:,:)
      p2(:,:) = p2_in(:,:)
      dgpress = sum((p2(:,:) - p1(:,:)) * rel_area(:,:))
      p2(:,:) = p2(:,:) - dgpress
      call average_press_poles(nx, ny, rel_area, p1)
      call average_press_poles(nx, ny, rel_area, p2)
      dps(:,:) = p2(:,:) - p1(:,:)

      !--- Core's divergence, polar caps area-averaged (pfix pole semantics) ---
      dpsc(:,:) = dps_ctm(:,:)
      call average_press_poles(nx, ny, rel_area, dpsc)

      !--- LLNL direct fix (terms only; cf. do_press_fix_llnl) ---
      ddps(:,:) = dps(:,:) - dpsc(:,:)
      dgpress   = sum(ddps(:,:) * rel_area(:,:))

      mmfd(:) = 0.0_8
      do j = j1p, j2p
         mmfd(j) = -(sum(ddps(:,j)) / ri2 - dgpress)
      end do
      mmfd(1)    = -(ddps(1,1)    - dgpress)
      mmfd(2)    = -(ddps(1,2)    - dgpress)
      mmfd(ny-1) = -(ddps(1,ny-1) - dgpress)
      mmfd(ny)   = -(ddps(1,ny)   - dgpress)

      mmf(:)   = 0.0_8
      mmf(j1p) = mmfd(1) / geofac_pc
      do j = j1p, j2p
         mmf(j+1) = mmf(j) + mmfd(j) / geofac(j)
      end do

      xcolmass_fix(:,:) = 0.0_8
      do j = j1p, j2p
         fxintegral(:) = 0.0_8
         do i = 1, nx
            fxintegral(i+1) = fxintegral(i) - (ddps(i,j) - dgpress) - mmfd(j)
         end do
         fxmean = sum(fxintegral(2:nx+1)) / ri2
         do i = 1, nx
            xcolmass_fix(i,j) = fxintegral(i) - fxmean
         end do
      end do
   end subroutine pfix_correction


   !! two southern-most and two northern-most rows equal (area-weighted).
   !! GEOS-Chem `Average_Press_Poles`.
   subroutine average_press_poles(nx, ny, rel_area, press)
      integer, intent(in)    :: nx, ny
      real,    intent(in)    :: rel_area(nx, ny)
      real,    intent(inout) :: press(nx, ny)
      real :: meanp

      meanp = sum(rel_area(:,1:2) * press(:,1:2)) / sum(rel_area(:,1:2))
      press(:,1:2) = meanp

      meanp = sum(rel_area(:,ny-1:ny) * press(:,ny-1:ny)) / sum(rel_area(:,ny-1:ny))
      press(:,ny-1:ny) = meanp
   end subroutine average_press_poles

   !> Horizontal pressure mass fluxes from A-grid winds.
   !! GEOS-Chem `Calc_Horiz_Mass_Flux`.  xmass at the west face of cell i;
   !! ymass at the south face of cell j (cose = cos of south edge).
   subroutine calc_horiz_mass_flux(nx, ny, nz, dt, dl, dp, cosp, cose, &
                                   delpm, u, v, xmass, ymass)
      integer, intent(in)  :: nx, ny, nz
      real,    intent(in)  :: dt, dl, dp
      real,    intent(in)  :: cosp(ny), cose(ny+1)
      real,    intent(in)  :: delpm(nx, ny, nz)
      real,    intent(in)  :: u(nx, ny, nz), v(nx, ny, nz)
      real,    intent(out) :: xmass(nx, ny, nz), ymass(nx, ny, nz)

      integer :: i, j, k
      real    :: factx, facty

      facty = 0.5_8 * dt / (PJC_RE * dp)

      ! E-W flux (periodic in i).
      do k = 1, nz
         do j = 1, ny
            factx = 0.5_8 * dt / (dl * PJC_RE * cosp(j))
            xmass(1,j,k) = factx * (u(1,j,k)*delpm(1,j,k) + u(nx,j,k)*delpm(nx,j,k))
            do i = 2, nx
               xmass(i,j,k) = factx * (u(i,j,k)*delpm(i,j,k) + u(i-1,j,k)*delpm(i-1,j,k))
            end do
         end do
      end do

      ! N-S flux.  ymass(:,1) = south-pole edge (cose(1)=0 -> zero).
      do k = 1, nz
         do i = 1, nx
            ymass(i,1,k) = facty * cose(1) * (v(i,1,k)*delpm(i,1,k))
         end do
         do j = 2, ny
            do i = 1, nx
               ymass(i,j,k) = facty * cose(j) * &
                    (v(i,j,k)*delpm(i,j,k) + v(i,j-1,k)*delpm(i,j-1,k))
            end do
         end do
      end do
   end subroutine calc_horiz_mass_flux

   !> Horizontal divergence of the pressure mass fluxes (per level).
   !! GEOS-Chem `Calc_Divergence` + `Do_Divergence_Pole_Sum` (poles).
   subroutine calc_divergence(nx, ny, nz, j1p, j2p, geofac_pc, geofac, &
                              xmass, ymass, dpi)
      integer, intent(in)  :: nx, ny, nz, j1p, j2p
      real,    intent(in)  :: geofac_pc, geofac(ny)
      real,    intent(in)  :: xmass(nx, ny, nz), ymass(nx, ny, nz)
      real,    intent(out) :: dpi(nx, ny, nz)

      integer :: i, j, k
      real    :: sumsp, sumnp

      dpi(:,:,:) = 0.0_8

      ! N-S then E-W divergence over the non-polar band.
      do k = 1, nz
         do j = j1p, j2p
            do i = 1, nx
               dpi(i,j,k) = (ymass(i,j,k) - ymass(i,j+1,k)) * geofac(j)
            end do
            do i = 1, nx-1
               dpi(i,j,k) = dpi(i,j,k) + xmass(i,j,k) - xmass(i+1,j,k)
            end do
            dpi(nx,j,k) = dpi(nx,j,k) + xmass(nx,j,k) - xmass(1,j,k)
         end do
      end do

      ! Poles (Do_Divergence_Pole_Sum): south = j=1, north = j=ny.
      do k = 1, nz
         sumsp = 0.0_8
         sumnp = 0.0_8
         do i = 1, nx
            sumsp = sumsp + ymass(i,j1p,k)
            sumnp = sumnp + ymass(i,j2p+1,k)
         end do
         dpi(:,1, k) = -sumsp / real(nx) * geofac_pc
         dpi(:,ny,k) =  sumnp / real(nx) * geofac_pc
      end do

      ! Polar-cap enlargement: copy the pole divergence into the inner ring.
      dpi(:,2,   :) = dpi(:,1, :)
      dpi(:,ny-1,:) = dpi(:,ny,:)
   end subroutine calc_divergence

   !> The LLNL direct pressure fixer.  GEOS-Chem `Do_Press_Fix_Llnl`.
   !! Computes the barotropic column correction (meridional `mmf`, zonal
   !! `xcolmass_fix`) that drives the flux divergence to `dps`, and adds it
   !! (weighted by dbk) to xmass,ymass.
   subroutine do_press_fix_llnl(nx, ny, nz, j1p, j2p, geofac_pc, geofac, &
                                dbk, dps, dps_c, rel_area, xmass, ymass)
      integer, intent(in)    :: nx, ny, nz, j1p, j2p
      real,    intent(in)    :: geofac_pc, geofac(ny), dbk(nz)
      real,    intent(in)    :: dps(nx, ny), dps_c(nx, ny), rel_area(nx, ny)
      real,    intent(inout) :: xmass(nx, ny, nz), ymass(nx, ny, nz)

      real    :: ddps(nx, ny)
      real    :: mmfd(ny), mmf(ny)
      real    :: fxintegral(nx+1)
      real    :: xcolmass_fix(nx, ny)
      real    :: dgpress, fxmean, ri2
      integer :: i, j, k

      ri2 = real(nx)

      ! Discrepancy between prescribed and wind-derived pressure tendency.
      ddps(:,:) = dps(:,:) - dps_c(:,:)
      dgpress   = sum(ddps(:,:) * rel_area(:,:))

      !--- Mean meridional flux divergence (zonal-mean pressure change) ---
      mmfd(:) = 0.0_8
      do j = j1p, j2p
         mmfd(j) = -(sum(ddps(:,j)) / ri2 - dgpress)
      end do
      ! Poles (rows already area-averaged, so element (1,j) is representative).
      mmfd(1)    = -(ddps(1,1)    - dgpress)
      mmfd(2)    = -(ddps(1,2)    - dgpress)
      mmfd(ny-1) = -(ddps(1,ny-1) - dgpress)
      mmfd(ny)   = -(ddps(1,ny)   - dgpress)

      !--- Mean meridional fluxes cos(e)*fy  (meridional integration) ---
      mmf(:)   = 0.0_8
      mmf(j1p) = mmfd(1) / geofac_pc
      do j = j1p, j2p
         mmf(j+1) = mmf(j) + mmfd(j) / geofac(j)
      end do

      !--- Zonal integration for the E-W column correction ---
      xcolmass_fix(:,:) = 0.0_8
      do j = j1p, j2p
         fxintegral(:) = 0.0_8
         do i = 1, nx
            fxintegral(i+1) = fxintegral(i) - (ddps(i,j) - dgpress) - mmfd(j)
         end do
         fxmean = sum(fxintegral(2:nx+1)) / ri2
         do i = 1, nx
            xcolmass_fix(i,j) = fxintegral(i) - fxmean
         end do
      end do

      !--- Distribute the column corrections in the vertical (weight dbk) ---
      do k = 1, nz
         do j = j1p, j2p
            do i = 1, nx
               xmass(i,j,k) = xmass(i,j,k) + xcolmass_fix(i,j) * dbk(k)
            end do
         end do
         do j = j1p, j2p+1
            do i = 1, nx
               ymass(i,j,k) = ymass(i,j,k) + mmf(j) * dbk(k)
            end do
         end do
      end do
   end subroutine do_press_fix_llnl

end module fv3_pfix_mod
