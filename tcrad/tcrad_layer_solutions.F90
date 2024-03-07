! tcrad_layer_solutions.F90 - Two-stream and related layer solutions for TCRAD package
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module tcrad_layer_solutions

  use parkind1, only : jpim, jprb

  implicit none
  public

  ! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined
#ifdef NG_LW
  integer, parameter, private :: ng = NG_LW
#else
#define ng ng_lw_in
#endif

  ! Two-stream scheme can be Elsasser: radiance field assumed to be
  ! two pencil beams travelling at cosine-zenith angles (mu) of
  ! +/-1/1.66; or Eddington: radiance field assumed to be
  ! L(mu)=L0+mu*L1.
  enum, bind(c)
    enumerator ITwoStreamElsasser, ITwoStreamEddington
  end enum

  ! Two stream scheme currently in use 
  integer(jpim) :: i_two_stream_scheme = ITwoStreamElsasser
  
  ! The effective factor by which the zenith optical depth needs to be
  ! multiplied to account for longwave transmission at all angles
  ! through the atmosphere.  Alternatively think of
  ! acos(1/lw_diffusivity) to be the effective zenith angle of
  ! longwave radiation.
  real(jprb) :: lw_diffusivity = 1.66_jprb ! Elsasser default

  ! To avoid division by near-zero values use simpler formulae in the
  ! low optical depth regime
  real(jprb), parameter :: OD_THRESH_2STREAM = 1.0e-3_jprb
  real(jprb), parameter :: OD_THRESH = 1.0e-3_jprb

  integer(jpim), parameter :: MAX_GAUSS_LEGENDRE_POINTS = 8

contains

  !---------------------------------------------------------------------
  ! Set the two-stream scheme; this overwrites global module variables
  ! so should normally be called from outside a parallel block
  subroutine set_two_stream_scheme(i_scheme)
    
    integer(jpim), intent(in) :: i_scheme

    ! Only overwrite global module variables if they need changing
    if (i_scheme == ITwoStreamEddington &
         .and. i_two_stream_scheme /= ITwoStreamEddington) then
      ! Toon et al. (1989), Table 1
      i_two_stream_scheme = ITwoStreamEddington
      lw_diffusivity = 2.0_jprb
    else if (i_two_stream_scheme /= ITwoStreamElsasser) then
      ! Elsasser (1942)
      i_two_stream_scheme = ITwoStreamElsasser
      lw_diffusivity = 1.66_jprb
    end if
    
  end subroutine set_two_stream_scheme
  
  
  !---------------------------------------------------------------------
  ! Return Gauss-Legendre quadrature points, or "optimized" values for
  ! radiation if npoint == -1 or -2
  subroutine gauss_legendre(npoint, xpoint, weight)

    integer(jpim), intent(in)  :: npoint
    real(jprb),    intent(out) :: xpoint(:), weight(:)

    ! Set default values
    xpoint = 0.5_jprb
    weight = 0.0_jprb

    if (npoint == 1 .or. npoint == 0) then
      xpoint(1) = 0.5_jprb
      weight(1) = 1.0_jprb
    else if (npoint == -1) then
      ! Optimized values minimizing error in transmission over the
      ! full range of optical depths; this is likely to be updated in
      ! future
      xpoint(1) = 0.6158_jprb;
      weight(1) = 0.5_jprb / xpoint(1)
    else if (npoint == 2) then
      xpoint(1:2) = [0.211324865405187_jprb, 0.788675134594813_jprb]
      weight(1:2) = [0.5_jprb, 0.5_jprb]
    else if (npoint == -2) then
      ! Optimized values minimizing error in transmission over the
      ! full range of optical depths; this is likely to be updated in
      ! future
      xpoint(1:2) = [0.2704_jprb, 0.8018_jprb];
      weight(1:2) = 1.0_jprb / sum(xpoint(1:2))
    else if (npoint == 3) then
      xpoint(1:3) = [0.112701665379258_jprb, 0.5_jprb, &
           &         0.887298334620742_jprb]
      weight(1:3) = [0.277777777777777_jprb, 0.444444444444444_jprb, &
           &         0.277777777777777_jprb]
    else if (npoint == 4) then
      xpoint(1:4) = [0.0694318442029737_jprb, 0.330009478207572_jprb, &
           &         0.669990521792428_jprb,  0.930568155797026_jprb]
      weight(1:4) = [0.173927422568727_jprb, 0.326072577431273_jprb, &
           &         0.326072577431273_jprb, 0.173927422568727_jprb]
    else if (npoint == 5) then
      xpoint(1:5) = [0.9530899230_jprb, 0.7692346551_jprb, &
           &  0.5000000000_jprb, 0.2307653449_jprb, 0.0469100770_jprb]
      weight(1:5) = [0.1184634425_jprb, 0.2393143352_jprb, &
           &  0.2844444444_jprb, 0.2393143352_jprb, 0.1184634425_jprb]
    else if (npoint == 6) then
      xpoint(1:6) = [0.9662347571_jprb, 0.8306046932_jprb, 0.6193095930_jprb, &
           &         0.3806904070_jprb, 0.1693953068_jprb, 0.0337652429_jprb]
      weight(1:6) = [0.0856622462_jprb, 0.1803807865_jprb, 0.2339569673_jprb, &
           &         0.2339569673_jprb, 0.1803807865_jprb, 0.0856622462_jprb]
      
    else if (npoint == 7) then
      xpoint(1:7) = [0.9745539562_jprb, 0.8707655928_jprb, 0.7029225757_jprb, &
           &         0.5000000000_jprb, 0.2970774243_jprb, 0.1292344072_jprb, &
           &         0.0254460438_jprb]
      weight(1:7) = [0.0647424831_jprb, 0.1398526957_jprb, 0.1909150253_jprb, &
           &         0.2089795918_jprb, 0.1909150253_jprb, 0.1398526957_jprb, &
           &         0.0647424831_jprb]
    else
      xpoint(1:8) = [0.9801449282_jprb, 0.8983332387_jprb, 0.7627662050_jprb, &
           &         0.5917173212_jprb, 0.4082826788_jprb, 0.2372337950_jprb, &
           &         0.1016667613_jprb, 0.0198550718_jprb]
      weight(1:8) = [0.0506142681_jprb, 0.1111905172_jprb, 0.1568533229_jprb, &
           &         0.1813418917_jprb, 0.1813418917_jprb, 0.1568533229_jprb,&
           &         0.1111905172_jprb, 0.0506142681_jprb]
    end if

  end subroutine gauss_legendre


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver (1980) two-stream formulas, as
  ! well as the upward flux at the top and the downward flux at the
  ! base of the layer due to emission from within the layer assuming a
  ! linear variation of Planck function within the layer.
  subroutine calc_reflectance_transmittance(ng_lw_in, nlev, nreg, &
       &  region_fracs, planck_hl, od, ssa, asymmetry, &
       &  reflectance, transmittance, source_up, source_dn)
        
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: ng_lw_in, nlev, nreg

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(ng,nlev+1) :: planck_hl

    ! Optical depth in each region and layer
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: od

    ! Single scattering albedo in each cloudy region
    real(jprb), intent(in), dimension(ng,2:nreg,nlev) :: ssa

    ! Asymmetry factor of clouds
    real(jprb), intent(in), dimension(ng,nlev) :: asymmetry

    ! Outputs

    ! Layer reflectance and transmittance
    real(jprb), intent(out), dimension(ng,nreg,nlev) :: reflectance, transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer, in
    ! units of Watts of power per square metre of the entire gridbox,
    ! so emission is proportional to the size of each region
    real(jprb), intent(out), dimension(ng,nreg,nlev) :: source_up, source_dn

    ! Two-stream exchange coefficients
#define OPT_REFTRANS 
#ifdef OPT_REFTRANS 
    real(jprb), dimension(ng) :: gamma1, gamma2, od_loc, k_exponent

    ! Working variables
    real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
    real(jprb) :: factor, exponential2, reftrans_factor
    real(jprb), parameter :: Half= 0.5_jprb
    real(jprb), parameter :: One = 1.0_jprb
    real(jprb), parameter :: Two = 2.0_jprb
    ! For optical depths lower than this value, use the linear
    ! approximation
    real(jprb), parameter :: OdThresholdLw = 1.0e-3_jprb
#ifdef PARKIND1_SINGLE
  ! real(jprb), parameter :: KMin = 1.e4_jprb * epsilon(1._jprb)
  real(jprb), parameter :: KMinSw = 1.e-4_jprb
  real(jprb), parameter :: KMinLw = 1.e-4_jprb
#else
  real(jprb), parameter :: KMinSw = 1.e-12_jprb
  real(jprb), parameter :: KMinLw = 1.e-12_jprb
#endif
#else
    real(jprb) :: gamma1, gamma2

    ! Working variables
    real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
    real(jprb) :: factor, exponential, exponential2, k_exponent, reftrans_factor
#endif
    ! Loop indices
    integer(jpim) :: jg, jreg, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance',0,hook_handle)

    ! Set cloudy regions to default values
    reflectance(:,2:,:)   = 0.0_jprb
    transmittance(:,2:,:) = 1.0_jprb
    source_up(:,2:,:)     = 0.0_jprb
    source_dn(:,2:,:)     = 0.0_jprb

    do jlev = 1,nlev
      ! No-scattering solution in clear-sky region: compute upward and
      ! downward emission assuming the Planck function to vary
      ! linearly with optical depth within the layer (e.g. Wiscombe,
      ! JQSRT 1976).
#ifndef OPT_REFTRANS 
      do jg = 1,ng
        reflectance(jg,1,jlev) = 0.0_jprb
        if (od(jg,1,jlev) > OD_THRESH_2STREAM) then
          coeff = lw_diffusivity*od(jg,1,jlev)
          transmittance(jg,1,jlev) = exp(-coeff)
          coeff = (planck_hl(jg,jlev+1)-planck_hl(jg,jlev)) / coeff
          coeff_up_top  =  coeff + planck_hl(jg,jlev)
          coeff_up_base =  coeff + planck_hl(jg,jlev+1)
          coeff_dn_top  = -coeff + planck_hl(jg,jlev)
          coeff_dn_base = -coeff + planck_hl(jg,jlev+1)
          source_up(jg,1,jlev) = coeff_up_top &
               &  - transmittance(jg,1,jlev) * coeff_up_base
          source_dn(jg,1,jlev) = coeff_dn_base &
               &  - transmittance(jg,1,jlev) * coeff_dn_top
        else
          ! Linear limit at low optical depth
          coeff = lw_diffusivity*od(jg,1,jlev)
          transmittance(jg,1,jlev) = 1.0_jprb - coeff
          source_up(jg,1,jlev) = coeff * 0.5_jprb &
               &  * (planck_hl(jg,jlev)+planck_hl(jg,jlev+1))
          source_dn(jg,1,jlev) = source_up(jg,1,jlev)
        end if
        ! Scale the sources by the area fraction of the region
        source_up(jg,1,jlev) = region_fracs(1,jlev)*source_up(jg,1,jlev)
        source_dn(jg,1,jlev) = region_fracs(1,jlev)*source_dn(jg,1,jlev)
      end do

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! We have some cloud
        do jreg = 2,nreg
          ! Scattering solution
          do jg = 1,ng
            if (i_two_stream_scheme == ITwoStreamElsasser) then
              ! See Fu et al. (1997), Eqs. 2.9 and 2.10
              factor = (lw_diffusivity * 0.5_jprb) * ssa(jg,jreg,jlev)
              gamma1 = lw_diffusivity - factor*(1.0_jprb + asymmetry(jg,jlev))
              gamma2 = factor * (1.0_jprb - asymmetry(jg,jlev))
            else
              ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
              gamma1 = 1.75_jprb - ssa(jg,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jg,jlev))
              gamma2 = ssa(jg,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jg,jlev)) - 0.25_jprb
            end if

            k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
                 1.E-12_jprb)) ! Eq 18 of Meador & Weaver (1980)
            if (od(jg,jreg,jlev) > OD_THRESH_2STREAM) then
              exponential = exp(-k_exponent*od(jg,jreg,jlev))
              exponential2 = exponential*exponential
              reftrans_factor = 1.0 / (k_exponent + gamma1&
                   &  + (k_exponent - gamma1)*exponential2)
              ! Meador & Weaver (1980) Eq. 25
              reflectance(jg,jreg,jlev) = gamma2 &
                   &  * (1.0_jprb - exponential2) * reftrans_factor
              ! Meador & Weaver (1980) Eq. 26
              transmittance(jg,jreg,jlev) = 2.0_jprb * k_exponent &
                   &  * exponential * reftrans_factor
      
              ! Compute upward and downward emission assuming the
              ! Planck function to vary linearly with optical depth
              ! within the layer (e.g. Wiscombe , JQSRT 1976).

              ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
              coeff = (planck_hl(jg,jlev+1)-planck_hl(jg,jlev)) &
                   / (od(jg,jreg,jlev)*(gamma1+gamma2))
              coeff_up_top  =  coeff + planck_hl(jg,jlev)
              coeff_up_base =  coeff + planck_hl(jg,jlev+1)
              coeff_dn_top  = -coeff + planck_hl(jg,jlev)
              coeff_dn_base = -coeff + planck_hl(jg,jlev+1)
              source_up(jg,jreg,jlev) = coeff_up_top &
                   &  - reflectance(jg,jreg,jlev) * coeff_dn_top &
                   &  - transmittance(jg,jreg,jlev) * coeff_up_base
              source_dn(jg,jreg,jlev) = coeff_dn_base &
                   &  - reflectance(jg,jreg,jlev) * coeff_up_base &
                   &  - transmittance(jg,jreg,jlev) * coeff_dn_top
            else
              ! Low optical depth approximation
              reflectance(jg,jreg,jlev) = gamma2 * od(jg,jreg,jlev)
              transmittance(jg,jreg,jlev) = (1.0_jprb - k_exponent*od(jg,jreg,jlev)) &
                   &  / (1.0_jprb + od(jg,jreg,jlev)*(gamma1-k_exponent))
              source_up(jg,jreg,jlev) = (1.0_jprb - reflectance(jg,jreg,jlev) &
                   &  - transmittance(jg,jreg,jlev)) &
                   &       * 0.5 * (planck_hl(jg,jlev) + planck_hl(jg,jlev+1))
              source_dn(jg,jreg,jlev) = source_up(jg,jreg,jlev)
            end if
            ! Scale the sources by the area fraction of the region
            source_up(jg,jreg,jlev) = region_fracs(jreg,jlev)*source_up(jg,jreg,jlev)
            source_dn(jg,jreg,jlev) = region_fracs(jreg,jlev)*source_dn(jg,jreg,jlev)
          end do
        end do
      end if

#else

      associate(planck_top=>planck_hl(:,jlev), planck_bot=>planck_hl(:,jlev+1), exponential=>od_loc, source_up_low_od=>gamma1)

      od_loc = lw_diffusivity*od(:,1,jlev)
      transmittance(:,1,jlev) = exp(-od_loc)
      do jg = 1, ng
        reflectance(jg,1,jlev) = 0.0_jprb
        source_up_low_od(jg) = od_loc(jg)*0.5_jprb*(planck_top(jg)+planck_bot(jg))
      ! Compute upward and downward emission assuming the Planck
      ! function to vary linearly with optical depth within the layer
      ! (e.g. Wiscombe , JQSRT 1976).
        od_loc(jg) = (planck_bot(jg)-planck_top(jg)) / od_loc(jg)
        coeff_up_top  =  od_loc(jg) + planck_top(jg)
        coeff_up_base  =  od_loc(jg) + planck_bot(jg)
        coeff_dn_top  = -od_loc(jg) + planck_top(jg)
        coeff_dn_base  = -od_loc(jg) + planck_bot(jg)
        source_up(jg,1,jlev) =  coeff_up_top - transmittance(jg,1,jlev) * coeff_up_base
        source_dn(jg,1,jlev) =  coeff_dn_base - transmittance(jg,1,jlev) * coeff_dn_top
      end do
        
      do jg = 1, ng
        if (od(jg,1,jlev) < 1.0e-3) then
          ! Linear limit at low optical depth
          ! coeff = LwDiffusivity*od(jg)
          ! transmittance(jg) = 1.0_jprb - coeff
          ! source_up(jg) = coeff * 0.5_jprb * (planck_top(jg)+planck_bot(jg))
          source_up(jg,1,jlev) = source_up_low_od(jg)
          source_dn(jg,1,jlev) = source_up(jg,1,jlev)
        end if
      end do

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! We have some cloud
        do jreg = 2,nreg
          ! Scattering solution
          if (i_two_stream_scheme == ITwoStreamElsasser) then
            do jg = 1, ng
              factor = (lw_diffusivity * 0.5_jprb) * ssa(jg,jreg,jlev)
              gamma1(jg) = lw_diffusivity - factor*(1.0_jprb + asymmetry(jg,jlev))
              gamma2(jg) = factor * (1.0_jprb - asymmetry(jg,jlev))
              k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
                KMinLw)) ! Eq 18 of Meador & Weaver (1980)
            end do
          else
            do jg = 1, ng
              ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
              gamma1 = 1.75_jprb - ssa(jg,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jg,jlev))
              gamma2 = ssa(jg,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jg,jlev)) - 0.25_jprb
              k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
                KMinLw)) ! Eq 18 of Meador & Weaver (1980)
            end do 
          end if
          exponential = exp(-k_exponent*od(:,jreg,jlev))


          do jg = 1, ng
            exponential2 = exponential(jg)*exponential(jg)
            reftrans_factor = 1.0_jprb / (k_exponent(jg)  + gamma1(jg) + (k_exponent(jg) - gamma1(jg))*exponential2)
            ! Meador & Weaver (1980) Eq. 25
            reflectance(jg,jreg,jlev) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
            ! Meador & Weaver (1980) Eq. 26
            transmittance(jg,jreg,jlev) = 2.0_jprb * k_exponent(jg) * exponential(jg) * reftrans_factor
            !
            ! Toon et al. (JGR 1989) Eqs 26-27
            !
            !   if (od(jg,jreg,jlev) > 1.0e-8_jprb) then
            coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg,jreg,jlev)*(gamma1(jg)+gamma2(jg)))
            coeff_up_top  =  coeff + planck_top(jg)
            coeff_up_base  =  coeff + planck_bot(jg)
            coeff_dn_top  = -coeff + planck_top(jg)
            coeff_dn_base  = -coeff + planck_bot(jg)
            source_up(jg,jreg,jlev) = (coeff_up_top - reflectance(jg,jreg,jlev) * coeff_dn_top - &
              & transmittance(jg,jreg,jlev) * coeff_up_base)
            source_dn(jg,jreg,jlev) = (coeff_dn_base - reflectance(jg,jreg,jlev) * coeff_up_base - &
              & transmittance(jg,jreg,jlev) * coeff_dn_top)
          !   else
          !     source_up(jg) = 0._jprb
          !     source_dn(jg) = 0._jprb
          !   end if
          end do

          do jg = 1, ng
            if (od(jg,jreg,jlev) < OdThresholdLw) then
              if (od(jg,jreg,jlev) < 1.0e-8_jprb) then 
                source_up(jg,jreg,jlev) = 0._jprb
                source_dn(jg,jreg,jlev) = 0._jprb
              else 
                ! In between e-3..e-8, must use Robin's original code, otherwise not secure 
                ! for instance when LW aerosol scattering is turned on 
                source_up(jg,jreg,jlev) = (One - reflectance(jg,jreg,jlev) - transmittance(jg,jreg,jlev)) &
                &       * Half * (planck_top(jg) + planck_bot(jg))
                source_dn(jg,jreg,jlev) = source_up(jg,jreg,jlev)
              end if 
            end if 
          end do
          
        end do
      end if

      ! Scale the sources by the area fraction of the region
      do jreg = 1, nreg
        do jg = 1, ng
          source_up(jg,jreg,jlev) = region_fracs(jreg,jlev)*source_up(jg,jreg,jlev)
          source_dn(jg,jreg,jlev) = region_fracs(jreg,jlev)*source_dn(jg,jreg,jlev)
        end do 
      end do 
      end associate 
#endif
    end do

    ! print *, "R min max", minval(reflectance), maxval(reflectance), "T", minval(transmittance), maxval(transmittance)

    if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance',1,hook_handle)

  end subroutine calc_reflectance_transmittance


  !---------------------------------------------------------------------
  ! Compute the rates of emission and scattering into a particular
  ! zenith angle cosine (mu) at the top and base of each layer and
  ! region
  subroutine calc_radiance_rates(ng_lw_in, nlev, nreg, &
       &  mu, region_fracs, planck_hl, ssa, asymmetry, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
       &  rate_up_top, rate_up_base, rate_dn_top, rate_dn_base)
    
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: ng_lw_in, nlev, nreg

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(ng,nlev+1) :: planck_hl

    ! Single scattering albedo in each cloudy region
    real(jprb), intent(in), dimension(ng,2:nreg,nlev) :: ssa

    ! Asymmetry factor of clouds
    real(jprb), intent(in), dimension(ng,nlev) :: asymmetry

    ! Upward and downward fluxes at the top and base of each layer and
    ! region, in Watts of power per square metre of the entire
    ! gridbox, so the energy is scaled by the size of each region
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: flux_up_base, flux_dn_base
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: flux_up_top, flux_dn_top
  
    ! Outputs

    ! Rate of emission/scattering upwards or downwards at the top or
    ! base of each layer and region, in Watts of power per square
    ! metre of the entire gridbox. These terms are available instead
    ! of source_up and source_dn for the 3D radiance calculation. Note
    ! that this is the rate of emission/scattering along the radiance
    ! direction per unit optical depth in that layer region, but since
    ! optical depth is dimensionless, these rates still have units of
    ! W m-2.
    real(jprb), intent(out), dimension(ng,nreg,nlev), optional &
         &  :: rate_up_top, rate_up_base, rate_dn_top, rate_dn_base

    ! Local variables

    ! Working variables in W m-2
    real(jprb), dimension(ng,nreg) :: planck_top, planck_base
    real(jprb), dimension(ng,nreg) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, factor, coeff, scaling, coeffplus, coeffmin

    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jg, jreg
    integer(jpim), dimension(2) :: jmax

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_rates',0,hook_handle)

    secant = 1.0_jprb / mu

    ! 0.5: half the scattering goes up and half down
    factor = 0.5_jprb * 3.0_jprb * mu / lw_diffusivity

    if (present(rate_up_top) .and. present(rate_up_base) &
      & .and. present(rate_dn_top) .and. present(rate_dn_base)) then
        ! Optimized code
        do jlev = 1,nlev
          ! CHECK / FIX: is the source that per unit optical depth in the
          ! vertical or along the direction to the sensor?
  
          ! Clear-sky region first
          rate_up_top(:,1,jlev)  = planck_hl(:,jlev) * region_fracs(1,jlev)
          rate_up_base(:,1,jlev) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
          rate_dn_top(:,1,jlev)  = planck_hl(:,jlev) * region_fracs(1,jlev)
          rate_dn_base(:,1,jlev) = planck_hl(:,jlev+1) * region_fracs(1,jlev)

          if (region_fracs(1,jlev) < 1.0_jprb) then
            ! Cloudy layer: scale the Planck terms by the region fraction
            ! and also by the single-scattering co-albedo
            do jreg = 2, nreg 
              do jg = 1, ng
                scaling = (1.0_jprb - ssa(jg,jreg,jlev)) * region_fracs(jreg,jlev)
                coeffplus =  0.5_jprb + factor*asymmetry(jg,jlev)
                coeffmin  =  0.5_jprb - factor*asymmetry(jg,jlev)

                ! Compute the rate of energy emitted or scattered in the
                ! upward direction mu at the top and base of the layer
                rate_up_top(jg,jreg,jlev) = planck_hl(jg,jlev)*scaling + ssa(jg,jreg,jlev) &
                    &  * (flux_up_top(jg,jreg,jlev) * coeffplus + flux_dn_top(jg,jreg,jlev) * coeffmin)
                rate_up_base(jg,jreg,jlev) =  planck_hl(jg,jlev+1)*scaling  + ssa(jg,jreg,jlev) &
                    &  * (flux_up_base(jg,jreg,jlev) * coeffplus + flux_dn_base(jg,jreg,jlev) * coeffmin)        

                ! Compute the rate of energy emitted or scattered in the
                ! downward direction mu at the top and base of the layer
                rate_dn_top(jg,jreg,jlev) = planck_hl(jg,jlev)*scaling + ssa(jg,jreg,jlev) &
                    &  * (flux_up_top(jg,jreg,jlev) * coeffmin + flux_dn_top(jg,jreg,jlev) * coeffplus)
                rate_dn_base(jg,jreg,jlev) = planck_hl(jg,jlev+1)*scaling + ssa(jg,jreg,jlev) &
                    &  * (flux_up_base(jg,jreg,jlev) * coeffmin + flux_dn_base(jg,jreg,jlev) * coeffplus)     
              end do 
            end do 
          end if
  
        end do

    else 

      do jlev = 1,nlev

        ! CHECK / FIX: is the source that per unit optical depth in the
        ! vertical or along the direction to the sensor?

        if (region_fracs(1,jlev) < 1.0_jprb) then
          ! Cloudy layer: scale the Planck terms by the region fraction
          ! and also by the single-scattering co-albedo
          max_reg = nreg
          planck_top(:,1) = planck_hl(:,jlev) * region_fracs(1,jlev)
          planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
          do jreg = 2, nreg 
            do jg = 1, ng
              coeff = (1.0_jprb - ssa(jg,jreg,jlev)) * region_fracs(jreg,jlev)
              planck_top(jg,jreg) = planck_hl(jg,jlev) * coeff 
              planck_base(jg,jreg) = planck_hl(jg,jlev+1) * coeff
              end do 
          end do 
        else
          ! Clear layer
          max_reg = 1
          planck_top(:,1)  = planck_hl(:,jlev)
          planck_base(:,1) = planck_hl(:,jlev+1)
          planck_top(:,2:) = 0.0_jprb 
          planck_base(:,2:) = 0.0_jprb
        end if

        if (present(rate_up_top) .and. present(rate_up_base)) then
          ! Compute the rate of energy emitted or scattered in the
          ! upward direction mu at the top and base of the layer: first
          ! the Planck emission for all regions...
          rate_up_top(:,1,jlev)  = planck_top(:,1)
          rate_up_base(:,1,jlev) = planck_base(:,1)
          ! ...then scattering from the scattering source function, but
          ! only in cloudy regions
          if (max_reg > 1) then
            do jreg = 2, nreg 
              do jg = 1, ng 
                coeffplus =  0.5_jprb + factor*asymmetry(jg,jlev)
                coeffmin =  0.5_jprb - factor*asymmetry(jg,jlev)

                rate_up_top(jg,jreg,jlev) = planck_top(jg,jreg) + ssa(jg,jreg,jlev) &
                    &  * (flux_up_top(jg,jreg,jlev) * coeffplus + flux_dn_top(jg,jreg,jlev) * coeffmin)
                rate_up_base(jg,jreg,jlev) = planck_base(jg,jreg) + ssa(jg,jreg,jlev) &
                    &  * (flux_up_base(jg,jreg,jlev) * coeffplus + flux_dn_base(jg,jreg,jlev) * coeffmin)    
              end do 
            end do 
          end if
        end if
        if (present(rate_dn_top) .and. present(rate_dn_base)) then
          ! Compute the rate of energy emitted or scattered in the
          ! downward direction mu at the top and base of the layer:
          ! first the Planck emission for all regions...
          rate_dn_top(:,1,jlev)  = planck_top(:,1)
          rate_dn_base(:,1,jlev) = planck_base(:,1)
          ! ...then scattering from the scattering source function, but
          ! only in cloudy regions
          if (max_reg > 1) then
            do jreg = 2, nreg 
              do jg = 1, ng 
                coeffmin =  0.5_jprb - factor*asymmetry(jg,jlev)
                coeffplus =  0.5_jprb + factor*asymmetry(jg,jlev)

                rate_dn_top(jg,jreg,jlev) = planck_top(jg,jreg) + ssa(jg,jreg,jlev) &
                    &  * (flux_up_top(jg,jreg,jlev) * coeffmin + flux_dn_top(jg,jreg,jlev) * coeffplus)
                rate_dn_base(jg,jreg,jlev) = planck_base(jg,jreg) + ssa(jg,jreg,jlev) &
                    &  * (flux_up_base(jg,jreg,jlev) * coeffmin + flux_dn_base(jg,jreg,jlev) * coeffplus)
              end do 
            end do 

          end if
        end if
      end do
    end if 
    if (lhook) call dr_hook('tcrad:calc_radiance_rates',1,hook_handle)

  end subroutine calc_radiance_rates


  !---------------------------------------------------------------------
  ! Compute the transmittance to a beam of radiation at a particular
  ! zenith angle cosine (mu), as well as optionally the source from
  ! the layer in that direction up and/or down. The latter includes
  ! only emission, so is suitable to be used in a no-scattering
  ! radiance calculation.
  subroutine calc_no_scattering_radiance_source(ng_lw_in, nlev, nreg, &
       &  mu, region_fracs, planck_hl, od, &
       &  transmittance, source_up, source_dn)
      
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: ng_lw_in, nlev, nreg

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(ng,nlev+1) :: planck_hl

    ! Optical depth in each region and layer
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: od
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(ng,nreg,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(ng,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables in W m-2
    real(jprb) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, coeff
   
    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jg, jreg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance_source',0,hook_handle)

    secant = 1.0_jprb / mu

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! Cloudy layer: scale the Planck terms by the region fraction
        ! and also by the single-scattering co-albedo
        max_reg = nreg
        transmittance(:,1:max_reg,jlev) = exp(-od(:,1:max_reg,jlev)*secant)
      else
        ! Clear layer
        max_reg = 1
        transmittance(:,1,jlev) = exp(-od(:,1,jlev)*secant)
        transmittance(:,2:nreg,jlev) = 1.0_jprb
      end if

      if (present(source_up)) then
        ! Compute upward source from layer top due to Planck emission
        ! within the layer
        do jreg = 1,max_reg
          do jg = 1,ng
            source_top  = planck_hl(jg,jlev)   * region_fracs(jreg,jlev)
            source_base = planck_hl(jg,jlev+1) * region_fracs(jreg,jlev)
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (source_base - source_top) &
                   &   * mu / od(jg,jreg,jlev)
              source_up(jg,jreg,jlev) = coeff + source_top &
                   - transmittance(jg,jreg,jlev) * (coeff + source_base)
            else
              source_up(jg,jreg,jlev) = od(jg,jreg,jlev) &
                   &  * 0.5_jprb * (source_top+source_base) / mu
            end if
          end do
        end do
        if (max_reg == 1) then
          ! In a clear layer, set cloudy values to zero
          source_up(:,2:nreg,jlev) = 0.0_jprb
        end if
      end if

      if (present(source_dn)) then
        ! Compute downward source from layer base due to Planck emission
        ! within the layer
        do jreg = 1,max_reg
          do jg = 1,ng
            source_top  = planck_hl(jg,jlev)   * region_fracs(jreg,jlev)
            source_base = planck_hl(jg,jlev+1) * region_fracs(jreg,jlev)
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (source_top - source_base) &
                   &   * mu / od(jg,jreg,jlev)
              source_dn(jg,jreg,jlev) = coeff + source_base &
                   - transmittance(jg,jreg,jlev) * (coeff + source_top)
            else
              source_dn(jg,jreg,jlev) = od(jg,jreg,jlev) &
                   &  * 0.5_jprb * (source_top+source_base) / mu
            end if
          end do
        end do
        if (max_reg == 1) then
          ! In a clear layer, set cloudy values to zero
          source_dn(:,2:nreg,jlev) = 0.0_jprb
        end if
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance_source',1,hook_handle)

  end subroutine calc_no_scattering_radiance_source


  !---------------------------------------------------------------------
  ! Compute the clear-sky transmittance to a beam of radiation at a
  ! particular zenith angle cosine (mu), as well as optionally the
  ! source from the layer in that direction up and/or down. The latter
  ! includes only emission, so is suitable to be used in a clear-sky
  ! no-scattering radiance calculation.
  subroutine calc_clear_sky_trans_source(ng_lw_in, nlev, &
       &  mu, planck_hl, od, &
       &  transmittance, source_up, source_dn)
      
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: ng_lw_in, nlev

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(ng,nlev+1) :: planck_hl

    ! Optical depth in each layer
    real(jprb), intent(in), dimension(ng,nlev) :: od
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(ng,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(ng,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables
    real(jprb) :: secant, coeff, coeff2, mu_over_od
   
    ! Loop indices for level and spectral interval
    integer(jpim) :: jlev, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source',0,hook_handle)

    secant = 1.0_jprb / mu

    transmittance = exp(-od*secant)

    do jlev = 1,nlev

      if (present(source_up) .and. present(source_dn)) then
        ! optimized version with only one loop having conditional; which is further implemented as a
        ! postprocessing loop to vectorize the main computations
        do jg = 1,ng
          mu_over_od =  mu / od(jg,jlev)
          coeff = (planck_hl(jg,jlev+1) - planck_hl(jg,jlev)) * mu_over_od
          coeff2 = (planck_hl(jg,jlev) - planck_hl(jg,jlev+1)) * mu_over_od
          ! Compute upward source from layer top and downward source from layer base
          ! due to Planck emission within the layer
          source_up(jg,jlev) = coeff + planck_hl(jg,jlev) &
              - transmittance(jg,jlev) * (coeff + planck_hl(jg,jlev+1))
          source_dn(jg,jlev) = coeff2 + planck_hl(jg,jlev+1) &
              - transmittance(jg,jlev) * (coeff2 + planck_hl(jg,jlev))
        end do 
        
        do jg = 1,ng
          if (od(jg,jlev) < OD_THRESH) then
            source_up(jg,jlev) = od(jg,jlev) &
                &  * 0.5_jprb * (planck_hl(jg,jlev)+planck_hl(jg,jlev+1)) / mu
            source_dn(jg,jlev) = source_up(jg,jlev)
          end if 
        end do

      else

        if (present(source_up)) then
          ! Compute upward source from layer top due to Planck emission
          ! within the layer
          do jg = 1,ng
            if (od(jg,jlev) > OD_THRESH) then
              coeff = (planck_hl(jg,jlev+1) - planck_hl(jg,jlev)) &
                  &   * mu / od(jg,jlev)
              source_up(jg,jlev) = coeff + planck_hl(jg,jlev) &
                  - transmittance(jg,jlev) * (coeff + planck_hl(jg,jlev+1))
            else
              source_up(jg,jlev) = od(jg,jlev) &
                  &  * 0.5_jprb * (planck_hl(jg,jlev)+planck_hl(jg,jlev+1)) / mu
            end if
          end do
        end if

        if (present(source_dn)) then
          ! Compute downward source from layer base due to Planck emission
          ! within the layer
          do jg = 1,ng
            if (od(jg,jlev) > OD_THRESH) then
              coeff = (planck_hl(jg,jlev) - planck_hl(jg,jlev+1)) &
                  &   * mu / od(jg,jlev)
              source_dn(jg,jlev) = coeff + planck_hl(jg,jlev+1) &
                  - transmittance(jg,jlev) * (coeff + planck_hl(jg,jlev))
            else
              source_dn(jg,jlev) = od(jg,jlev) &
                  &  * 0.5_jprb * (planck_hl(jg,jlev)+planck_hl(jg,jlev+1)) / mu
            end if
          end do
        end if

      end if 

    end do

    if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source',1,hook_handle)

  end subroutine calc_clear_sky_trans_source

  ! subroutine calc_clear_sky_trans_source_collapse(ng_lw_in, nlev, &
  !      &  mu, planck_top, planck_bot, od, &
  !      &  transmittance, source_up, source_dn)
      
  !   use yomhook,  only           : lhook, dr_hook, jphook

  !   ! Inputs

  !   ! Number of spectral intervals and levels
  !   integer(jpim), intent(in) :: ng_lw_in, nlev

  !   ! Cosine of the zenith angle (positive)
  !   real(jprb) :: mu

  !   ! Planck function integrated over each spectral interval at each
  !   ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
  !   ! black-body surface)
  !   real(jprb), intent(in), dimension(ng*nlev) :: planck_top, planck_bot

  !   ! Optical depth in each layer
  !   real(jprb), intent(in), dimension(ng*nlev) :: od
  
  !   ! Outputs

  !   ! Layer transmittance at the requested zenith angle
  !   real(jprb), intent(out), dimension(ng*nlev) :: transmittance

  !   ! Source term up from the top of the layer or down from its base,
  !   ! in Watts of power per square metre of the entire gridbox, so the
  !   ! energy is scaled by the size of each region. Since the user may
  !   ! only require a radiance up or down, these output arguments are
  !   ! optional.
  !   real(jprb), intent(out), dimension(ng*nlev), optional &
  !        &  :: source_up, source_dn

  !   ! Working variables
  !   real(jprb) :: secant, coeff, coeff2, mu_over_od
   
  !   ! Loop indices for level and spectral interval
  !   integer(jpim) :: jlev, jg

  !   real(jphook) :: hook_handle

  !   if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source',0,hook_handle)

  !   secant = 1.0_jprb / mu

  !   transmittance = exp(-od*secant)

  !   if (present(source_up) .and. present(source_dn)) then
  !     ! optimized version with only one loop having conditional; which is further implemented as a
  !     ! postprocessing loop to vectorize the main computations
  !     do jg = 1,ng*nlev
  !       mu_over_od =  mu / od(jg)
  !       coeff = (planck_bot(jg) - planck_top(jg)) * mu_over_od
  !       coeff2 = (planck_top(jg) - planck_bot(jg)) * mu_over_od
  !       ! Compute upward source from layer top and downward source from layer base
  !       ! due to Planck emission within the layer
  !       source_up(jg) = coeff + planck_top(jg) &
  !           - transmittance(jg) * (coeff + planck_bot(jg))
  !       source_dn(jg) = coeff2 + planck_bot(jg) &
  !           - transmittance(jg) * (coeff2 + planck_top(jg))
  !     end do 
      
  !     do jg = 1,ng*nlev
  !       if (od(jg) < OD_THRESH) then
  !         source_up(jg) = od(jg) &
  !             &  * 0.5_jprb * (planck_top(jg)+planck_bot(jg)) / mu
  !         source_dn(jg) = source_up(jg)
  !       end if 
  !     end do

  !   else

  !     if (present(source_up)) then
  !       ! Compute upward source from layer top due to Planck emission
  !       ! within the layer
  !       do jg = 1,ng*nlev
  !         ! if (od(jg) > OD_THRESH) then
  !           coeff = (planck_bot(jg) - planck_top(jg)) &
  !               &   * mu / od(jg)
  !           source_up(jg) = coeff + planck_top(jg) &
  !               - transmittance(jg) * (coeff + planck_bot(jg))
  !         ! else
  !         !   source_up(jg) = od(jg) &
  !         !       &  * 0.5_jprb * (planck_top(jg)+planck_bot(jg)) / mu
  !         ! end if
  !       end do
  !       do jg = 1,ng*nlev
  !         if (od(jg) < OD_THRESH) then
  !           source_up(jg) = od(jg) &
  !               &  * 0.5_jprb * (planck_top(jg)+planck_bot(jg)) / mu
  !         end if
  !       end do
  !     end if

  !     if (present(source_dn)) then
  !       ! Compute downward source from layer base due to Planck emission
  !       ! within the layer
  !       do jg = 1,ng*nlev
  !         ! if (od(jg) > OD_THRESH) then
  !           coeff = (planck_top(jg) - planck_bot(jg)) &
  !               &   * mu / od(jg)
  !           source_dn(jg) = coeff + planck_bot(jg) &
  !               - transmittance(jg) * (coeff + planck_top(jg))
  !         ! else
  !         !   source_dn(jg) = od(jg) &
  !         !       &  * 0.5_jprb * (planck_top(jg)+planck_bot(jg)) / mu
  !         ! end if
  !       end do
  !       do jg = 1,ng*nlev
  !         if (od(jg) < OD_THRESH) then
  !           source_dn(jg) = od(jg) &
  !               &  * 0.5_jprb * (planck_top(jg)+planck_bot(jg)) / mu
  !         end if
  !       end do
  !     end if

  !   end if 


  !   if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source_collapse',1,hook_handle)

  ! end subroutine calc_clear_sky_trans_source_collapse
  
  !---------------------------------------------------------------------
  ! Calculate the transmittance of each layer and region along a path
  ! with consine of zenith angle "mu", as well as (optionally) the
  ! emission up from the top of the layer and down through its base,
  ! computed neglecting 3D effects
  subroutine calc_radiance_trans_source(ng_lw_in, nlev, nreg, &
             &  mu, region_fracs, od, transmittance, &
             &  rate_up_top, rate_up_base, rate_dn_top, rate_dn_base, &
             &  source_up, source_dn)

    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: ng_lw_in, nlev, nreg

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Optical depth in each region and layer
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: od
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(ng,nreg,nlev) :: transmittance

    ! Optional inputs

    ! Rate of emission/scattering upwards or downwards at the top or
    ! base of each layer and region, in Watts of power per square
    ! metre of the entire gridbox. These terms are available instead
    ! of source_up and source_dn for the 3D radiance calculation. Note
    ! that this is the rate of emission/scattering along the radiance
    ! direction per unit optical depth in that layer region, but since
    ! optical depth is dimensionless, these rates still have units of
    ! W m-2.
    real(jprb), intent(in), dimension(ng,nreg,nlev), optional &
         &  :: rate_up_top, rate_up_base, rate_dn_top, rate_dn_base

    ! Optional outputs

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(ng,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Local variables

    ! Working variables
    real(jprb) :: secant, coeff

    ! Maximum number of active regions in layer
    integer(jpim) :: max_reg

    ! Loop indices
    integer(jpim) :: jlev, jreg, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',0,hook_handle)

    secant = 1.0_jprb / mu

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        max_reg = nreg;
      else
        max_reg = 1;
      end if

      do jreg = 1,max_reg
        transmittance(:,jreg,jlev) = exp(-od(:,jreg,jlev)*secant)

        if (present(source_dn)) then
          do jg = 1,ng
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (rate_dn_top(jg,jreg,jlev) - rate_dn_base(jg,jreg,jlev)) &
                   &   * mu / od(jg,jreg,jlev)
              source_dn(jg,jreg,jlev) = coeff + rate_dn_base(jg,jreg,jlev) &
                   - transmittance(jg,jreg,jlev) * (coeff + rate_dn_top(jg,jreg,jlev))
            else
              source_dn(jg,jreg,jlev) = od(jg,jreg,jlev) * 0.5_jprb &
                   &  * (rate_dn_top(jg,jreg,jlev) + rate_dn_base(jg,jreg,jlev)) / mu
            end if
          end do
        end if
        if (present(source_up)) then
          do jg = 1,ng
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (rate_up_base(jg,jreg,jlev) - rate_up_top(jg,jreg,jlev)) &
                   &   * mu / od(jg,jreg,jlev)
              source_up(jg,jreg,jlev) = coeff + rate_up_top(jg,jreg,jlev) &
                   - transmittance(jg,jreg,jlev) * (coeff + rate_up_base(jg,jreg,jlev))
            else
              source_up(jg,jreg,jlev) = od(jg,jreg,jlev) * 0.5_jprb &
                   &  * (rate_up_top(jg,jreg,jlev) + rate_up_base(jg,jreg,jlev)) / mu
            end if
          end do
        end if
      end do

      if (max_reg == 1) then
        transmittance(:,2:,jlev) = 1.0_jprb
        if (present(source_dn)) then
          source_dn(:,2:,jlev) = 0.0_jprb
        end if
        if (present(source_up)) then
          source_up(:,2:,jlev) = 0.0_jprb
        end if
      end if
      
    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',1,hook_handle)

  end subroutine calc_radiance_trans_source


!#warning "calc_radiance_source is deprecated"

  !---------------------------------------------------------------------
  ! Compute the transmittance to a beam of radiation at a particular
  ! zenith angle cosine (mu), as well as optionally the source from
  ! the layer in that direction up and/or down. The latter includes
  ! both emission, and the source term from scattering using the
  ! incoming flux terms from a previous two-stream calculation.
  subroutine calc_radiance_source(ng_lw_in, nlev, nreg, &
       &  mu, region_fracs, planck_hl, od, ssa, asymmetry, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
       &  transmittance, source_up, source_dn)
    
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: ng_lw_in, nlev, nreg

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(ng,nlev+1) :: planck_hl

    ! Optical depth in each region and layer
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: od

    ! Single scattering albedo in each cloudy region
    real(jprb), intent(in), dimension(ng,2:nreg,nlev) :: ssa

    ! Asymmetry factor of clouds
    real(jprb), intent(in), dimension(ng,nlev) :: asymmetry

    ! Upward and downward fluxes at the top and base of each layer and
    ! region, in Watts of power per square metre of the entire
    ! gridbox, so the energy is scaled by the size of each region
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: flux_up_base, flux_dn_base
    real(jprb), intent(in), dimension(ng,nreg,nlev) :: flux_up_top, flux_dn_top
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(ng,nreg,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(ng,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables in W m-2
    real(jprb), dimension(ng,nreg) :: planck_top, planck_base
    real(jprb), dimension(ng,nreg) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, factor, coeff

    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jg, jreg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_source',0,hook_handle)

    secant = 1.0_jprb / mu

    ! 0.5: half the scattering goes up and half down
    factor = 0.5_jprb * 3.0_jprb * mu / lw_diffusivity

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! Cloudy layer: scale the Planck terms by the region fraction
        ! and also by the single-scattering co-albedo
        max_reg = nreg
        planck_top(:,1) = planck_hl(:,jlev) * region_fracs(1,jlev)
        planck_top(:,2:nreg) = spread(planck_hl(:,jlev),2,nreg-1) &
             &  * (1.0_jprb - ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,ng)
        planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
        planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2,nreg-1) &
             &  * (1.0_jprb - ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,ng)
        ! Compute transmittance in all regions
        transmittance(:,1:nreg,jlev) = exp(-od(:,1:nreg,jlev)*secant)
      else
        ! Clear layer
        max_reg = 1
        planck_top(:,1)  = planck_hl(:,jlev)
        planck_base(:,1) = planck_hl(:,jlev+1)
        ! Compute transmittance in active region
        transmittance(:,1,jlev) = exp(-od(:,1,jlev)*secant)
        transmittance(:,2:nreg,jlev) = 1.0_jprb
      end if

      if (present(source_up)) then
        ! Compute the rate of energy emitted or scattered in the
        ! upward direction mu at the top and base of the layer: first
        ! the Planck emission for all regions...
        source_top  = planck_top
        source_base = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          source_top(:,2:nreg) = source_top(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
          source_base(:,2:nreg) = source_base(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
        else
          source_up(:,2:nreg,jlev) = 0.0_jprb
        end if

        ! Compute the energy making it to the top of the layer
        ! accounting for attenuation within it
        do jreg = 1,max_reg
          do jg = 1,ng
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (source_base(jg,jreg)-source_top(jg,jreg)) &
                   &   * mu / od(jg,jreg,jlev)
              source_up(jg,jreg,jlev) = coeff + source_top(jg,jreg) &
                   - transmittance(jg,jreg,jlev) * (coeff + source_base(jg,jreg))
            else
              source_up(jg,jreg,jlev) = od(jg,jreg,jlev) * 0.5_jprb &
                   &  * (source_base(jg,jreg)+source_top(jg,jreg)) / mu
            end if
          end do
        end do
      end if

      if (present(source_dn)) then
        ! Compute the rate of energy emitted or scattered in the
        ! downward direction mu at the top and base of the layer:
        ! first the Planck emission for all regions...
        source_top  = planck_top
        source_base = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          source_top(:,2:nreg) = source_top(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
          source_base(:,2:nreg) = source_base(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
        else
          source_dn(:,2:nreg,jlev) = 0.0_jprb
        end if
        ! Compute the energy making it to the top of the layer
        ! accounting for attenuation within it
        do jreg = 1,max_reg
          do jg = 1,ng
            if (od(jg,jreg,jlev) > OD_THRESH) then
              coeff = (source_top(jg,jreg)-source_base(jg,jreg)) &
                   &   * mu / od(jg,jreg,jlev)
              source_dn(jg,jreg,jlev) = coeff + source_base(jg,jreg) &
                   - transmittance(jg,jreg,jlev) * (coeff + source_top(jg,jreg))
            else
              source_dn(jg,jreg,jlev) = od(jg,jreg,jlev) * 0.5_jprb &
                   &  * (source_base(jg,jreg)+source_top(jg,jreg)) / mu
            end if
          end do
        end do
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_source',1,hook_handle)

  end subroutine calc_radiance_source


end module tcrad_layer_solutions
