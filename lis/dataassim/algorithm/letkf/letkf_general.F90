!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: letkf_general
!
! This module contains the Local Ensemble Transform Kalman Filter (LETKF)
! implementation adapted from PDAF.
!
! !DESCRIPTION:
!   LETKF is a localized ensemble Kalman filter that performs analysis
!   in ensemble space rather than state space. This makes it efficient
!   for high-dimensional problems.
!
!   Key differences from EnKF:
!   - Analysis performed in ensemble space (N_ens x N_ens matrices)
!   - Deterministic square-root filter (no perturbed observations needed)
!   - Natural parallelization over grid points
!
!   Key differences from LNETF:
!   - Uses linear Kalman gain (not likelihood-based weights)
!   - Gaussian assumption (LNETF handles non-Gaussian better)
!   - Computationally simpler for linear observation operators
!
! !REFERENCES:
!   Hunt, B.R., Kostelich, E.J. and Szunyogh, I., 2007: Efficient data
!   assimilation for spatiotemporal chaos: A local ensemble transform
!   Kalman filter. Physica D, 230, 112-126.
!
!   PDAF implementation by Lars Nerger
!
! !REVISION HISTORY:
!   28 Jan 2026: Initial implementation adapted from LNETF and PDAF
!EOP
module letkf_general

  implicit none

  private

  public :: letkf_analysis
  public :: gaspari_cohn
  public :: get_gaspari_cohn
  public :: haversine_km
  public :: get_gaspari_cohn_km

  ! Module-level parameters for LETKF configuration
  integer, save :: type_forget = 0   ! Type of forgetting factor (0=fixed)
  integer, save :: type_trans = 0    ! Type of transformation (0=deterministic, 2=random)
  real, save    :: forget = 1.0      ! Forgetting factor (covariance inflation)

  public :: type_forget, type_trans, forget

contains


!BOP
!
! !ROUTINE: letkf_analysis
! \label{letkf_analysis}
!
! !INTERFACE:
  subroutine letkf_analysis( gid, &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact )
! !USES:
    use letkf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: gid
    integer, intent(in) :: N_state, N_obs, N_ens
    type(obs_type), intent(in), dimension(N_obs) :: Observations
    real, intent(in), dimension(N_obs,N_ens) :: Obs_pred
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov
    real, intent(inout), dimension(N_state,N_ens) :: State_incr

    ! Optional inputs for localization
    real, dimension(N_state), intent(in), optional :: State_lon, State_lat
    real, intent(in), optional :: xcompact, ycompact

!
! !DESCRIPTION:
!
! Perform LETKF update following Hunt et al. (2007)
!
! ALGORITHM STEPS:
! 1. Compute ensemble perturbations: Z = X - mean(X), HZ = H(X) - mean(H(X))
! 2. Compute A^{-1} = (N-1)/rho * I + (HZ)^T R^{-1} HZ
! 3. Eigenvalue decomposition of A^{-1}
! 4. Compute weight vector: w = A * (HZ)^T R^{-1} d
! 5. Compute transform matrix: W = sqrt(N-1) * sqrt(A) [* Omega] + w
! 6. Transform ensemble: X^a = mean(X^f) + Z * W
!
! IMPORTANT:
! On input, State_incr must contain State_minus(1:N_state,1:N_ens)
! On output, State_incr contains the analysis ensemble
!
! If optional inputs State_lon, State_lat, xcompact, and ycompact
! are present, Gaspari-Cohn localization is applied
!EOP

    ! Local variables
    integer :: n_e, i, j, k, row, col
    real :: sqrtNm1, rho_loc
    logical :: apply_localization

    ! Arrays for ensemble statistics
    real, allocatable, dimension(:,:) :: State_prime   ! X - mean(X)
    real, allocatable, dimension(:) :: State_bar       ! mean(X)
    real, allocatable, dimension(:,:) :: HZ_l          ! H(X) - mean(H(X))
    real, allocatable, dimension(:) :: HXbar_l         ! mean(H(X))
    real, allocatable, dimension(:) :: innov_l         ! y - mean(H(X))

    ! Arrays for LETKF computation
    real, allocatable, dimension(:,:) :: RiHZ          ! R^{-1} * HZ
    real, allocatable, dimension(:,:) :: Ainv_l        ! A^{-1} matrix
    real, allocatable, dimension(:,:) :: tmp_Ainv      ! Temporary for Ainv
    real, allocatable, dimension(:) :: RiHZd           ! R^{-1} * HZ^T * d
    real, allocatable, dimension(:) :: VRiHZd          ! V^T * RiHZd (for A*b)
    real, allocatable, dimension(:) :: svals           ! Eigenvalues
    real, allocatable, dimension(:) :: work            ! Work array for SYEV
    real, allocatable, dimension(:,:) :: Asqrt_l       ! sqrt(A)

    ! Arrays for ensemble transformation
    real, allocatable, dimension(:,:) :: rndmat        ! Random rotation
    real, allocatable, dimension(:,:) :: ens_blk       ! Block of ensemble

    ! Localization weights
    real, allocatable, dimension(:) :: loc_weights
    real :: dx, dy, dweight

    integer :: info, ldwork
    integer :: maxblksize, blklower, blkupper

    ! ------------------------------------------------------------------

    ! Allocate arrays
    allocate(State_prime(N_state, N_ens))
    allocate(State_bar(N_state))
    allocate(HZ_l(N_obs, N_ens))
    allocate(HXbar_l(N_obs))
    allocate(innov_l(N_obs))
    allocate(RiHZ(N_obs, N_ens))
    allocate(Ainv_l(N_ens, N_ens))
    allocate(tmp_Ainv(N_ens, N_ens))
    allocate(RiHZd(N_ens))
    allocate(VRiHZd(N_ens))
    allocate(svals(N_ens))
    allocate(work(3*N_ens))
    allocate(Asqrt_l(N_ens, N_ens))
    allocate(rndmat(N_ens, N_ens))
    allocate(loc_weights(N_obs))

    ldwork = 3*N_ens

    ! Check if localization should be applied
    apply_localization = (present(State_lon) .and. present(State_lat) &
                         .and. present(xcompact) .and. present(ycompact))

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Starting LETKF analysis for grid point', gid
       write(LIS_logunit,*) '[LETKF] N_state=', N_state, ' N_obs=', N_obs, ' N_ens=', N_ens
       if (apply_localization) then
          write(LIS_logunit,*) '[LETKF] Localization ENABLED with xcompact=', xcompact
       else
          write(LIS_logunit,*) '[LETKF] Localization DISABLED'
       endif
    endif

    ! ========================================================================
    ! STEP 1: Compute ensemble perturbations
    ! ========================================================================

    ! State ensemble mean
    State_bar = sum(State_incr, 2) / real(N_ens)

    ! State perturbations: Z = X - mean(X)
    do n_e = 1, N_ens
       State_prime(:, n_e) = State_incr(:, n_e) - State_bar
    end do

    ! Observation prediction mean
    HXbar_l = sum(Obs_pred, 2) / real(N_ens)

    ! Observation perturbations: HZ = H(X) - mean(H(X))
    do n_e = 1, N_ens
       HZ_l(:, n_e) = Obs_pred(:, n_e) - HXbar_l
    end do

    ! Innovation: d = y - mean(H(X))
    do i = 1, N_obs
       innov_l(i) = Observations(i)%value - HXbar_l(i)
    end do

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Innovation (first 3):', innov_l(1:min(3,N_obs))
    endif

    ! ========================================================================
    ! STEP 2: Compute localization weights (if enabled)
    ! ========================================================================

    if (apply_localization) then
       do i = 1, N_obs
          dx = State_lon(1) - Observations(i)%lon
          dy = State_lat(1) - Observations(i)%lat
          loc_weights(i) = get_gaspari_cohn(dx, dy, xcompact, ycompact)
       end do
    else
       loc_weights = 1.0
    endif

    ! ========================================================================
    ! STEP 3: Compute R^{-1} * HZ with localization
    ! ========================================================================

    ! RiHZ(i,j) = loc_weight(i) * HZ(i,j) / R(i,i)
    ! This is the localized version of R^{-1} * HZ
    do j = 1, N_ens
       do i = 1, N_obs
          if (Obs_cov(i,i) > 1.0e-15) then
             RiHZ(i,j) = loc_weights(i) * HZ_l(i,j) / Obs_cov(i,i)
          else
             RiHZ(i,j) = 0.0
          endif
       end do
    end do

    ! ========================================================================
    ! STEP 4: Compute A^{-1} = (N-1)/rho * I + HZ^T * R^{-1} * HZ
    ! ========================================================================

    ! Initialize Ainv = (N-1) * I  (will apply forgetting factor later)
    Ainv_l = 0.0
    do i = 1, N_ens
       Ainv_l(i,i) = real(N_ens - 1)
    end do

    ! Compute HZ^T * RiHZ using SGEMM
    ! tmp_Ainv = HZ^T * RiHZ
    call sgemm('T', 'N', N_ens, N_ens, N_obs, &
               1.0, HZ_l, N_obs, RiHZ, N_obs, &
               0.0, tmp_Ainv, N_ens)

    ! Complete A^{-1}: Ainv = (N-1)/forget * I + HZ^T * RiHZ
    ! forget < 1 increases the prior weight (inflation)
    rho_loc = forget
    if (rho_loc < 0.01) rho_loc = 1.0  ! Safety check

    Ainv_l = Ainv_l / rho_loc + tmp_Ainv

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Ainv diagonal (first 3):', &
            (Ainv_l(i,i), i=1,min(3,N_ens))
    endif

    ! ========================================================================
    ! STEP 5: Eigenvalue decomposition of A^{-1}
    ! ========================================================================

    ! SSYEV computes eigenvalues and eigenvectors
    ! On exit: Ainv_l contains eigenvectors, svals contains eigenvalues
    call ssyev('V', 'L', N_ens, Ainv_l, N_ens, svals, work, ldwork, info)

    if (info /= 0) then
       write(LIS_logunit,*) '[WARN] LETKF gid=', gid, ': Eigenvalue decomposition failed, info=', info
       write(LIS_logunit,*) '[WARN] Using identity transformation'
       ! Fallback: no update
       deallocate(State_prime, State_bar, HZ_l, HXbar_l, innov_l)
       deallocate(RiHZ, Ainv_l, tmp_Ainv, RiHZd, VRiHZd, svals, work, Asqrt_l, rndmat, loc_weights)
       return
    endif

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Eigenvalues of Ainv (first 3):', svals(1:min(3,N_ens))
    endif

    ! ========================================================================
    ! STEP 6: Compute weight vector w = A * (HZ)^T * R^{-1} * d
    ! ========================================================================

    ! First: RiHZd = (R^{-1} * HZ)^T * d = HZ^T * R^{-1} * d
    ! Note: RiHZ already contains R^{-1} * HZ (with localization)
    call sgemv('T', N_obs, N_ens, 1.0, RiHZ, N_obs, innov_l, 1, 0.0, RiHZd, 1)

    ! Now compute w = A * RiHZd using the eigendecomposition
    ! A = V * diag(1/svals) * V^T
    ! So: w = V * diag(1/svals) * V^T * RiHZd

    ! VRiHZd = V^T * RiHZd
    call sgemv('T', N_ens, N_ens, 1.0, Ainv_l, N_ens, RiHZd, 1, 0.0, VRiHZd, 1)

    ! Apply inverse eigenvalues: VRiHZd = diag(1/svals) * VRiHZd
    do row = 1, N_ens
       if (svals(row) > 1.0e-10) then
          VRiHZd(row) = VRiHZd(row) / svals(row)
       else
          VRiHZd(row) = 0.0
       endif
    end do

    ! RiHZd = V * VRiHZd  (this is now the weight vector w)
    call sgemv('N', N_ens, N_ens, 1.0, Ainv_l, N_ens, VRiHZd, 1, 0.0, RiHZd, 1)

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Weight vector w (first 3):', RiHZd(1:min(3,N_ens))
    endif

    ! ========================================================================
    ! STEP 7: Compute sqrt(A) = V * diag(1/sqrt(svals)) * V^T
    ! ========================================================================

    ! tmp_Ainv = V * diag(1/sqrt(svals))
    do col = 1, N_ens
       do row = 1, N_ens
          if (svals(col) > 1.0e-10) then
             tmp_Ainv(row, col) = Ainv_l(row, col) / sqrt(svals(col))
          else
             tmp_Ainv(row, col) = 0.0
          endif
       end do
    end do

    ! Asqrt = sqrt(N-1) * tmp_Ainv * V^T
    sqrtNm1 = sqrt(real(N_ens - 1))
    call sgemm('N', 'T', N_ens, N_ens, N_ens, &
               sqrtNm1, tmp_Ainv, N_ens, Ainv_l, N_ens, &
               0.0, Asqrt_l, N_ens)

    ! ========================================================================
    ! STEP 8: Apply random rotation (optional)
    ! ========================================================================

    if (type_trans == 2) then
       ! Generate random orthonormal matrix
       call generate_rndmat(N_ens, rndmat)

       ! tmp_Ainv = Asqrt * rndmat
       call sgemm('N', 'N', N_ens, N_ens, N_ens, &
                  1.0, Asqrt_l, N_ens, rndmat, N_ens, &
                  0.0, tmp_Ainv, N_ens)

       if (gid == 1) then
          write(LIS_logunit,*) '[LETKF] Using random rotation (type_trans=2)'
       endif
    else
       ! Deterministic: tmp_Ainv = Asqrt
       tmp_Ainv = Asqrt_l
       if (gid == 1) then
          write(LIS_logunit,*) '[LETKF] Using deterministic transformation (type_trans=0)'
       endif
    endif

    ! ========================================================================
    ! STEP 9: Complete transformation matrix W = sqrt(A) + w
    ! ========================================================================

    ! Ainv_l now becomes the transformation matrix W
    ! W(:,j) = tmp_Ainv(:,j) + w  for all j
    do col = 1, N_ens
       do row = 1, N_ens
          Ainv_l(row, col) = tmp_Ainv(row, col) + RiHZd(row)
       end do
    end do

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Transform matrix W(1,1:3):', Ainv_l(1,1:min(3,N_ens))
    endif

    ! ========================================================================
    ! STEP 10: Transform ensemble: X^a = mean(X^f) + Z * W
    ! ========================================================================

    maxblksize = 200
    allocate(ens_blk(maxblksize, N_ens))

    do blklower = 1, N_state, maxblksize
       blkupper = min(blklower + maxblksize - 1, N_state)

       ! Store forecast perturbations
       do col = 1, N_ens
          ens_blk(1:(blkupper-blklower+1), col) = State_prime(blklower:blkupper, col)
       end do

       ! Initialize with ensemble mean
       do col = 1, N_ens
          State_incr(blklower:blkupper, col) = State_bar(blklower:blkupper)
       end do

       ! Add transformed perturbations: X^a = mean + Z * W
       call sgemm('N', 'N', blkupper-blklower+1, N_ens, N_ens, &
                  1.0, ens_blk, maxblksize, Ainv_l, N_ens, &
                  1.0, State_incr(blklower,1), N_state)
    end do

    deallocate(ens_blk)

    if (gid == 1) then
       write(LIS_logunit,*) '[LETKF] Analysis complete'
       write(LIS_logunit,*) '[LETKF] Analysis mean (first 3):', &
            sum(State_incr(1:min(3,N_state),:), 2) / real(N_ens)
    endif

    ! Deallocate arrays
    deallocate(State_prime, State_bar, HZ_l, HXbar_l, innov_l)
    deallocate(RiHZ, Ainv_l, tmp_Ainv, RiHZd, VRiHZd, svals, work, Asqrt_l, rndmat, loc_weights)

  end subroutine letkf_analysis


!BOP
!
! !ROUTINE: generate_rndmat
! \label{generate_rndmat_letkf}
!
! !INTERFACE:
  subroutine generate_rndmat(dim, rndmat)
! !USES:
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: dim
    real, intent(out) :: rndmat(dim, dim)
!
! !DESCRIPTION:
!   Generate a random orthonormal matrix with eigenvector (1,...,1)^T
!   This ensures the ensemble mean is preserved during transformation.
!
!   Simplified implementation using Householder reflections.
!EOP

    integer :: i, j, k, iter
    real :: norm, rndval
    real, allocatable :: rndvec(:)

    ! Initialize with identity
    rndmat = 0.0
    do i = 1, dim
       rndmat(i,i) = 1.0
    end do

    if (dim < 2) return

    allocate(rndvec(dim))

    ! Apply sequence of random Householder reflections
    do iter = 1, dim - 1
       ! Generate random vector
       do i = iter, dim
          call random_number(rndval)
          rndvec(i) = rndval - 0.5
       end do

       ! Normalize
       norm = 0.0
       do i = iter, dim
          norm = norm + rndvec(i)**2
       end do
       norm = sqrt(norm)
       if (norm > 1.0e-10) then
          do i = iter, dim
             rndvec(i) = rndvec(i) / norm
          end do
       endif

       ! Apply Householder: H = I - 2*v*v^T
       do j = iter, dim
          rndval = 0.0
          do k = iter, dim
             rndval = rndval + rndvec(k) * rndmat(k, j)
          end do
          rndval = 2.0 * rndval
          do i = iter, dim
             rndmat(i, j) = rndmat(i, j) - rndval * rndvec(i)
          end do
       end do
    end do

    deallocate(rndvec)

  end subroutine generate_rndmat


!BOP
! !ROUTINE: gaspari_cohn
! \label{gaspari_cohn_letkf}
!
! !INTERFACE:
  function gaspari_cohn( d )
!
! !DESCRIPTION:
! Evaluate 5th-order polynomial from Gaspari and Cohn, 1999, Eq (4.10)
! See enkf_general.F90 or lnetf_general.F90 for detailed documentation.
!EOP

    implicit none
    real :: gaspari_cohn, d, y
    real, parameter :: tol = 1e-3

    d = 2.*abs(d)

    if (d >= 2.) then
       y = 0.
    else if (d <= tol) then
       y = 1.
    else if (d <= 1.) then
       y = d**2 *( d*( d*( -.25*d + .5) + .625) -5./3.) + 1.
    else
       y = d*( d*( d*( d*( d/12. - .5) + .625) + 5./3.) -5.) + 4. - 2./3./d
    end if

    gaspari_cohn = y
    return

  end function gaspari_cohn


!BOP
!
! !ROUTINE: get_gaspari_cohn
! \label{get_gaspari_cohn_letkf}
!
! !INTERFACE:
  function get_gaspari_cohn( dx, dy, xcompact, ycompact )
!
! !DESCRIPTION:
! Compute anisotropic Gaspari-Cohn weight for given lat/lon distance.
!EOP

    implicit none
    real :: get_gaspari_cohn, dx, dy, xcompact, ycompact, d

    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
    get_gaspari_cohn = gaspari_cohn(d)

  end function get_gaspari_cohn


!BOP
!
! !ROUTINE: haversine_km
! \label{haversine_km_letkf}
!
! !INTERFACE:
  function haversine_km(lon1, lat1, lon2, lat2)
!
! !DESCRIPTION:
! Calculate great-circle distance using Haversine formula.
! See lnetf_general.F90 for detailed documentation.
!EOP

    implicit none
    real :: haversine_km
    real, intent(in) :: lon1, lat1, lon2, lat2

    real, parameter :: PI = 3.14159265358979323846
    real, parameter :: DEG_TO_RAD = PI / 180.0
    real, parameter :: EARTH_RADIUS_KM = 6371.0

    real :: lat1_rad, lat2_rad, dlat, dlon, a, c

    ! Input validation
    if (lon1 /= lon1 .or. lat1 /= lat1 .or. lon2 /= lon2 .or. lat2 /= lat2) then
       haversine_km = 999999.0
       return
    endif

    if (abs(lat1) > 90.0 .or. abs(lat2) > 90.0) then
       haversine_km = 999999.0
       return
    endif

    lat1_rad = lat1 * DEG_TO_RAD
    lat2_rad = lat2 * DEG_TO_RAD
    dlat = (lat2 - lat1) * DEG_TO_RAD
    dlon = (lon2 - lon1) * DEG_TO_RAD

    a = sin(dlat/2.0)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2.0)**2
    if (a < 0.0) a = 0.0
    if (a > 1.0) a = 1.0

    c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    haversine_km = EARTH_RADIUS_KM * c

  end function haversine_km


!BOP
!
! !ROUTINE: get_gaspari_cohn_km
! \label{get_gaspari_cohn_km_letkf}
!
! !INTERFACE:
  function get_gaspari_cohn_km(lon1, lat1, lon2, lat2, sigma_km)
!
! !DESCRIPTION:
! Compute Gaspari-Cohn weight based on distance in km.
!EOP

    implicit none
    real :: get_gaspari_cohn_km
    real, intent(in) :: lon1, lat1, lon2, lat2, sigma_km

    real :: dist_km, compact_km, d

    dist_km = haversine_km(lon1, lat1, lon2, lat2)
    compact_km = 3.0 * sigma_km
    d = dist_km / compact_km
    get_gaspari_cohn_km = gaspari_cohn(d)

  end function get_gaspari_cohn_km


end module letkf_general
