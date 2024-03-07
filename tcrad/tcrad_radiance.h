! tcrad_radiance.h - Calculate radiances in TCRAD -*- f90 -*-
!
! (C) Copyright 2021- ECMWF.
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
! This file is included in the modules specifying the NREGION
! parameter (typically 2 or 3) which makes this routine use a
! "doubleclouds" or "tripleclouds" assumption.
!

!---------------------------------------------------------------------
! Compute upward radiance profile by solving the Schwarzschild
! radiative transfer equation in each region of each layer assuming
! the source term to vary linearly with optical depth in each
! layer. The overlap matrix is used to translate the radiances exiting
! the top of one layer into radiances in the base of the layer above,
! consistent with the Tripleclouds approximation. This routine adds to
! any existing radiance profile, useful if used as part of a
! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
! et al. (1997).
subroutine calc_radiance_up(ng_lw_in, nlev, &
     &  weight, surf_up, &
     &  transmittance, source_up, u_overlap, radiance_up)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: ng_lw_in, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Surface upwelling flux in W m-2
  real(jprb), intent(in),  dimension(ng,NREGION) :: surf_up

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: transmittance

  ! Upward source from the top of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: source_up

  ! Upward overlap matrix - see Hogan et al. (JGR 2016) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap

  ! Output

  ! Upward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(ng,nlev+1) :: radiance_up

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(ng,NREGION) :: radiance_base, radiance_top

  ! Loop index for level
  integer(jpim) :: jlev

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_up',0,hook_handle)

  ! Set surface upward radiance
  radiance_base = weight * surf_up

  ! Save radiance profile averaged over regions, adding to existing
  ! values
  radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + sum(radiance_base,2)

  do jlev = nlev,1,-1
    ! Solution to Schwarzschild equation
    radiance_top = transmittance(:,:,jlev)*radiance_base + weight * source_up(:,:,jlev)
    ! Overlap rules to obtain radiances at base of the layer above
    radiance_base = singlemat_x_vec(ng,u_overlap(:,:,jlev),radiance_top)
    ! Save radiances
    radiance_up(:,jlev) = radiance_up(:,jlev) + sum(radiance_base,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_up',1,hook_handle)

end subroutine calc_radiance_up


!---------------------------------------------------------------------
! Compute downward radiance profile by solving the Schwarzschild
! radiative transfer equation in each region of each layer assuming
! the source term to vary linearly with optical depth in each
! layer. The overlap matrix is used to translate the radiances exiting
! the base of one layer into radiances in the top of the layer below,
! consistent with the Tripleclouds approximation. This routine adds to
! any existing radiance profile, useful if used as part of a
! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
! et al. (1997).
subroutine calc_radiance_dn(ng_lw_in, nlev, &
     &  weight, transmittance, source_dn, v_overlap, radiance_dn)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: ng_lw_in, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: transmittance

  ! Down source from the base of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: source_dn

  ! Downward overlap matrix - see Shonk and Hogan (2008) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap

  ! Output

  ! Downward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(ng,nlev+1) :: radiance_dn

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(ng,NREGION) :: radiance_base, radiance_top

  ! Loop index for level
  integer(jpim) :: jlev

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_dn',0,hook_handle)

  ! Start with zero at TOA
  radiance_top = 0.0_jprb

  do jlev = 1,nlev
    ! Solution to Schwarzschild equation
    radiance_base = transmittance(:,:,jlev)*radiance_top + weight * source_dn(:,:,jlev)
    ! Overlap rules to obtain radiances at base of the layer above
    radiance_top = singlemat_x_vec(ng,v_overlap(:,:,jlev+1),radiance_base)
    ! Save radiances
    radiance_dn(:,jlev+1) = radiance_dn(:,jlev+1) + sum(radiance_top,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_dn',1,hook_handle)

end subroutine calc_radiance_dn


!---------------------------------------------------------------------
! As calc_radiance_up but with a matrix for the transmission,
! representing lateral exchange of radiation between regions
subroutine calc_radiance_up_3d(ng_lw_in, nlev, &
     &  weight, surf_up, &
     &  transmittance, source_up, u_overlap, radiance_up)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: ng_lw_in, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Surface upwelling flux in W m-2
  real(jprb), intent(in),  dimension(ng,NREGION) :: surf_up

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(ng,NREGION,NREGION,nlev) :: transmittance

  ! Upward source from the top of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: source_up

  ! Upward overlap matrix - see Hogan et al. (JGR 2016) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap

  ! Output

  ! Upward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(ng,nlev+1) :: radiance_up

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(ng,NREGION) :: radiance_base, radiance_top

  ! Loop index
  integer(jpim) :: jlev, j

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_up_3d',0,hook_handle)

  ! Set surface upward radiance
  radiance_base = weight * surf_up

  ! Save radiance profile averaged over regions, adding to existing
  ! values
  radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + sum(radiance_base,2)

  do jlev = nlev,1,-1
    if (NREGION==3) then 
      associate(A=>transmittance(:,:,:,jlev), rt=>radiance_top, rb=>radiance_base, C=>u_overlap(:,:,jlev))
        do j = 1, ng
          ! Transmittance including exchange between regions
          ! both inner and outer loop of the matrix loops j1 and j2 unrolled
          rt(j,1) = A(j,1,1)*rb(j,1) + A(j,1,2)*rb(j,2) + A(j,1,3)*rb(j,3) + weight * source_up(j,1,jlev)
          rt(j,2) = A(j,2,1)*rb(j,1) + A(j,2,2)*rb(j,2) + A(j,2,3)*rb(j,3) + weight * source_up(j,2,jlev)
          rt(j,3) = A(j,3,1)*rb(j,1) + A(j,3,2)*rb(j,2) + A(j,3,3)*rb(j,3) + weight * source_up(j,3,jlev)
           
          ! Overlap rules to obtain radiances at base of the layer above
          rb(j,1) = C(1,1)*rt(j,1) + C(1,2)*rt(j,2) + C(1,3)*rt(j,3) ! j1=1
          rb(j,2) = C(2,1)*rt(j,1) + C(2,2)*rt(j,2) + C(2,3)*rt(j,3) ! j1=2
          rb(j,3) = C(3,1)*rt(j,1) + C(3,2)*rt(j,2) + C(3,3)*rt(j,3) ! j1=3

          ! Save radiances 
          radiance_up(j,jlev) = radiance_up(j,jlev) + rb(j,1) + rb(j,2) + rb(j,3)
        end do
      end associate
    else 
      ! Transmittance including exchange between regions
      radiance_top = mat_x_vec(ng,transmittance(:,:,:,jlev),radiance_base) &
          &       + weight * source_up(:,:,jlev)
      ! Overlap rules to obtain radiances at base of the layer above
      radiance_base = singlemat_x_vec(ng,u_overlap(:,:,jlev),radiance_top)
      ! Save radiances
      radiance_up(:,jlev) = radiance_up(:,jlev) + sum(radiance_base,2)
    end if 
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_up_3d',1,hook_handle)

end subroutine calc_radiance_up_3d


!---------------------------------------------------------------------
! As calc_radiance_dn but with a matrix for the transmission,
! representing lateral exchange of radiation between regions
subroutine calc_radiance_dn_3d(ng_lw_in, nlev, &
     &  weight, transmittance, source_dn, v_overlap, radiance_dn)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: ng_lw_in, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(ng,NREGION,NREGION,nlev) :: transmittance

  ! Down source from the base of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(ng,NREGION,nlev) :: source_dn

  ! Downward overlap matrix - see Shonk and Hogan (2008) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap

  ! Output

  ! Downward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(ng,nlev+1) :: radiance_dn

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(ng,NREGION) :: radiance_base, radiance_top

  ! Loop index
  integer(jpim) :: jlev, j

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_3d',0,hook_handle)

  ! Start with zero at TOA
  radiance_top = 0.0_jprb

  do jlev = 1,nlev
    if (NREGION==3) then 
      associate(A=>transmittance(:,:,:,jlev), rt=>radiance_top, rb=>radiance_base, C=>v_overlap(:,:,jlev+1))
        do j = 1, ng 
          ! Solution to Schwarzschild equation
          ! both inner and outer loop of the matrix loops j1 and j2 unrolled
          ! inner loop:   j2=1               j2=2                  j2=3 
          rb(j,1) = A(j,1,1)*rt(j,1) + A(j,1,2)*rt(j,2) + A(j,1,3)*rt(j,3) + weight * source_dn(j,1,jlev)
          rb(j,2) = A(j,2,1)*rt(j,1) + A(j,2,2)*rt(j,2) + A(j,2,3)*rt(j,3) + weight * source_dn(j,2,jlev)
          rb(j,3) = A(j,3,1)*rt(j,1) + A(j,3,2)*rt(j,2) + A(j,3,3)*rt(j,3) + weight * source_dn(j,3,jlev)
           
          ! Overlap rules to obtain radiances at base of the layer above
          rt(j,1) = C(1,1)*rb(j,1) + C(1,2)*rb(j,2) + C(1,3)*rb(j,3) ! j1=1
          rt(j,2) = C(2,1)*rb(j,1) + C(2,2)*rb(j,2) + C(2,3)*rb(j,3) ! j1=2
          rt(j,3) = C(3,1)*rb(j,1) + C(3,2)*rb(j,2) + C(3,3)*rb(j,3) ! j1=3

          ! Save radiances 
          radiance_dn(j,jlev+1) = radiance_dn(j,jlev+1) + rt(j,1) + rt(j,2) + rt(j,3)
        end do
      end associate
    else 
      ! Solution to Schwarzschild equation
      radiance_base = mat_x_vec(ng,transmittance(:,:,:,jlev),radiance_top) &
          &        + weight * source_dn(:,:,jlev)
      ! Overlap rules to obtain radiances at base of the layer above
      radiance_top = singlemat_x_vec(ng,v_overlap(:,:,jlev+1),radiance_base)

      ! Save radiances
      radiance_dn(:,jlev+1) = radiance_dn(:,jlev+1) + sum(radiance_top,2)
    end if 

  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_3d',1,hook_handle)

end subroutine calc_radiance_dn_3d




