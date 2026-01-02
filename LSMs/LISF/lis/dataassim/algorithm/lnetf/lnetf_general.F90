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
! !MODULE: lnetf_general
!
! This module contains the Local Nonlinear Ensemble Transform Filter (LNETF)
! implementation adapted from PDAF.
!
! !DESCRIPTION:
!   LNETF is a nonlinear ensemble data assimilation method that uses
!   particle weights based on observation likelihoods instead of the
!   Kalman gain. It is particularly suitable for non-Gaussian error
!   distributions and nonlinear observation operators.
!
!   This implementation now supports spatial localization using the
!   Gaspari-Cohn (1999) 5th-order polynomial function, consistent
!   with the EnKF implementation.
!
! !REFERENCES:
!   Toedter, J. and Ahrens, B., 2015: A Second-Order Exact Ensemble Square
!   Root Filter for Nonlinear Data Assimilation. Mon. Wea. Rev., 143, 1347-1367.
!
!   Gaspari, G. and Cohn, S.E., 1999: Construction of correlation functions
!   in two and three dimensions. Q.J.R. Meteorol. Soc., 125, 723-757.
!
!   PDAF implementation by Lars Nerger and Paul Kirchgessner
!
! !REVISION HISTORY:
!   27 Dec 2024: Initial implementation adapted from PDAF
!   01 Jan 2026: Added Gaspari-Cohn localization support
!EOP
module lnetf_general

  implicit none

  private

  public :: lnetf_analysis
  public :: compute_likelihood_gaussian
  public :: compute_likelihood_localized
  public :: compute_effective_sample_size
  public :: inflate_weights
  public :: generate_rndmat
  public :: gaspari_cohn
  public :: get_gaspari_cohn
  public :: haversine_km
  public :: get_gaspari_cohn_km

  ! Module-level parameters for LNETF configuration
  ! These can be set via namelist or configuration file
  integer, save :: type_winf = 1       ! Type of weights inflation (0=off, 1=on)
  real, save    :: limit_winf = 0.0    ! Limit for N_eff/N below which to inflate
  integer, save :: type_trans = 2      ! Type of transformation (1=random, 2=deterministic)
  real, save    :: forget = 1.0        ! Forgetting factor for covariance inflation

  public :: type_winf, limit_winf, type_trans, forget

contains


!BOP
!
! !ROUTINE: lnetf_analysis
! \label{lnetf_analysis}
!
! !INTERFACE:
  subroutine lnetf_analysis( gid, &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_err, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact, &
       n_eff_out, mean_loc_weight_out, max_loc_weight_out, min_loc_weight_out )
! !USES:
    use lnetf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: gid
    integer, intent(in) :: N_state, N_obs, N_ens
    type(obs_type), intent(in), dimension(N_obs) :: Observations
    real, intent(in), dimension(N_obs,N_ens) :: Obs_pred
    real, intent(in), dimension(N_obs,N_ens) :: Obs_err
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov
    real, intent(inout), dimension(N_state,N_ens) :: State_incr

    ! optional inputs for localization
    real, dimension(N_state), intent(in), optional :: State_lon, State_lat
    real, intent(in), optional :: xcompact       ! [deg] longitude
    real, intent(in), optional :: ycompact       ! [deg] latitude

    ! optional outputs for diagnostics
    real, intent(out), optional :: n_eff_out           ! Effective sample size
    real, intent(out), optional :: mean_loc_weight_out ! Mean Gaspari-Cohn weight
    real, intent(out), optional :: max_loc_weight_out  ! Max Gaspari-Cohn weight
    real, intent(out), optional :: min_loc_weight_out  ! Min Gaspari-Cohn weight

!
! !DESCRIPTION:
!
! Perform LNETF update using particle weights and ensemble transformation
!
! ALGORITHM STEPS:
! 1. Compute particle weights (likelihood) for each ensemble member
!    - If localization is enabled, apply Gaspari-Cohn distance weighting
! 2. Normalize weights and handle zero-weight cases
! 3. Compute transform matrix A = diag(w) - w*w^T
! 4. Calculate square root of A via eigenvalue decomposition
! 5. Apply random rotation (or deterministic transformation)
! 6. Transform ensemble: X^a = X^f * W
!
! IMPORTANT:
! On input, State_incr must contain State_minus(1:N_state,1:N_ens)
! On output, State_incr contains the analysis ensemble
!
! If optional inputs State_lon, State_lat, xcompact, and ycompact
! are present, Gaspari-Cohn localization is applied to the likelihood
!EOP

    ! Local variables
    integer :: n_e, i, j, k
    real :: total_weight, weight, n_eff
    real :: fac, dx, dy, dweight
    logical :: apply_localization

    ! Log-sum-exp variables for numerical stability
    real :: log_weight, max_log_weight
    real, dimension(N_ens) :: log_weights

    ! Large arrays - allocate on heap to avoid stack overflow
    real, allocatable, dimension(:,:) :: State_prime, State_bar_mat
    real, allocatable, dimension(:) :: State_bar

    real, allocatable, dimension(:) :: innov_i, obs_values
    real, allocatable, dimension(:) :: loc_weights  ! Localization weights

    real, dimension(N_ens) :: weights
    real, dimension(N_ens,N_ens) :: A_matrix, T_matrix, T_tmp, rndmat
    real, dimension(N_ens) :: eigenvalues
    real, dimension(3*N_ens) :: work
    integer :: info, lwork

    integer :: maxblksize, blklower, blkupper
    real, allocatable :: ens_blk(:,:)

    ! ------------------------------------------------------------------

    ! Allocate large arrays on heap
    allocate(State_prime(N_state, N_ens))
    allocate(State_bar_mat(N_state, N_ens))
    allocate(State_bar(N_state))
    allocate(innov_i(N_obs))
    allocate(obs_values(N_obs))
    allocate(loc_weights(N_obs))

    ! Check if localization should be applied
    apply_localization = (present(State_lon) .and. present(State_lat) &
                         .and. present(xcompact) .and. present(ycompact))

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Starting LNETF analysis for grid point', gid
       write(LIS_logunit,*) '[LNETF] N_state=', N_state, ' N_obs=', N_obs, ' N_ens=', N_ens
       if (apply_localization) then
          write(LIS_logunit,*) '[LNETF] Localization ENABLED with xcompact=', xcompact, ' ycompact=', ycompact
       else
          write(LIS_logunit,*) '[LNETF] Localization DISABLED (1D filter mode)'
       endif
    endif

    ! ========================================================================
    ! STEP 0: Calculate localization weights (if enabled)
    ! ========================================================================

    if (apply_localization) then
       ! Calculate Gaspari-Cohn weight for each observation based on distance
       ! from the state location (first element, as all state variables are
       ! at the same grid point)
       do i = 1, N_obs
          dx = State_lon(1) - Observations(i)%lon
          dy = State_lat(1) - Observations(i)%lat
          loc_weights(i) = get_gaspari_cohn(dx, dy, xcompact, ycompact)
       end do

       if (gid == 1) then
          write(LIS_logunit,*) '[LNETF] Localization weights (first 3):', &
               loc_weights(1:min(3,N_obs))
       endif
    else
       ! No localization - all weights are 1.0
       loc_weights = 1.0
    endif

    ! ========================================================================
    ! STEP 1: Compute particle weights (LOG-likelihoods for numerical stability)
    ! ========================================================================

    ! Extract observation values
    do i = 1, N_obs
       obs_values(i) = Observations(i)%value
    end do

    ! Compute LOG-likelihood for each ensemble member
    do n_e = 1, N_ens
       ! Calculate innovation: obs - H(x_i)
       innov_i = obs_values - Obs_pred(:, n_e)

       ! Compute localized Gaussian LOG-likelihood (returns log value)
       call compute_likelihood_localized(N_obs, innov_i, Obs_cov, loc_weights, log_weight)

       log_weights(n_e) = log_weight
    end do

    ! Debug: Show observation predictions and innovations for first few members
    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Observations:', obs_values(1:min(3,N_obs))
       write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Obs_pred for first 3 members:'
       do n_e = 1, min(3, N_ens)
          write(LIS_logunit,*) '  Member', n_e, ':', Obs_pred(1:min(3,N_obs), n_e)
       end do
       ! Show ensemble spread in predictions
       if (N_obs >= 1) then
          write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Obs_pred stats for first obs:'
          write(LIS_logunit,*) '  Min:', minval(Obs_pred(1,:)), ' Max:', maxval(Obs_pred(1,:))
          write(LIS_logunit,*) '  Mean:', sum(Obs_pred(1,:))/real(N_ens)
       end if
    end if

    if (gid == 1) write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Log weights:', log_weights(1:min(5,N_ens))

    ! ========================================================================
    ! STEP 2: Normalize weights using log-sum-exp trick for numerical stability
    ! This avoids underflow when likelihoods are very small
    ! ========================================================================

    ! Find maximum log-weight to prevent overflow in exp()
    max_log_weight = maxval(log_weights)

    if (gid == 1) write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Max log weight:', max_log_weight

    ! Check if all log-weights are -infinity (all innovations are huge)
    if (max_log_weight < -700.0) then
       ! All likelihoods underflowed - use equal weights
       write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': All log-weights < -700, using equal weights'
       weights = 1.0 / real(N_ens)
    else
       ! Log-sum-exp trick: w_i = exp(log_w_i - max_log_w) / sum(exp(log_w_j - max_log_w))
       ! This is numerically stable because the largest exp() argument is 0
       do n_e = 1, N_ens
          weights(n_e) = exp(log_weights(n_e) - max_log_weight)
       end do

       ! Normalize
       total_weight = sum(weights)
       if (total_weight > 1.0e-15) then
          weights = weights / total_weight
       else
          write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': Zero total weight after exp - using equal weights'
          weights = 1.0 / real(N_ens)
       endif
    endif

    if (gid == 1) write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Normalized weights:', weights(1:min(5,N_ens))

    ! ========================================================================
    ! STEP 2.5: Compute effective sample size and apply weights inflation
    ! ========================================================================

    call compute_effective_sample_size(N_ens, weights, n_eff)

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Effective sample size: N_eff =', n_eff, &
                            ' (', 100.0*n_eff/real(N_ens), '%)'
    endif

    ! Apply weights inflation if enabled and N_eff/N < limit_winf
    ! This prevents filter degeneracy by spreading out the weights
    if (type_winf == 1 .and. limit_winf > 0.0) then
       if (n_eff / real(N_ens) < limit_winf) then
          if (gid == 1) then
             write(LIS_logunit,*) '[LNETF] Applying weights inflation (N_eff/N < limit_winf)'
             write(LIS_logunit,*) '[LNETF]   N_eff/N =', n_eff/real(N_ens), ' limit =', limit_winf
          endif
          call inflate_weights(N_ens, limit_winf, weights)
          ! Recompute effective sample size after inflation
          call compute_effective_sample_size(N_ens, weights, n_eff)
          if (gid == 1) then
             write(LIS_logunit,*) '[LNETF] After inflation: N_eff =', n_eff, &
                                  ' (', 100.0*n_eff/real(N_ens), '%)'
          endif
       endif
    endif

    ! ========================================================================
    ! STEP 3: Compute transform matrix A = diag(w) - w*w^T
    ! ========================================================================

    if (gid == 1) write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Building transform matrix A'

    ! Initialize A = -w*w^T
    do j = 1, N_ens
       do i = 1, N_ens
          A_matrix(i,j) = -weights(i) * weights(j)
       end do
    end do

    ! Add diag(w) to get A = diag(w) - w*w^T
    do i = 1, N_ens
       A_matrix(i,i) = A_matrix(i,i) + weights(i)
    end do

    ! Check for NaN/Inf in A_matrix
    if (gid == 1) then
       do i = 1, N_ens
          if (A_matrix(i,i) /= A_matrix(i,i)) then
             write(LIS_logunit,*) '[ERR] LNETF: A_matrix diagonal has NaN at', i
          endif
       end do
       write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' A diagonal:', (A_matrix(i,i), i=1,min(3,N_ens))
       write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Sum of A:', sum(A_matrix)
    endif

    ! ========================================================================
    ! STEP 4: Compute sqrt(A) via eigenvalue decomposition
    ! ========================================================================

    ! Eigenvalue decomposition of A
    ! A_matrix will be overwritten with eigenvectors
    lwork = 3*N_ens

    ! Call LAPACK's symmetric eigenvalue routine
    ! Note: 'V' means compute eigenvectors, 'L' means lower triangle
    ! Initialize eigenvalues to zero in case ssyev fails
    eigenvalues = 0.0

    ! Use SSYEV (single precision) to match Fortran REAL type
    call ssyev('V', 'L', N_ens, A_matrix, N_ens, eigenvalues, work, lwork, info)

    if (info /= 0) then
       write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': Eigenvalue decomposition info=', info
       write(LIS_logunit,*) '[WARN] Note: A is rank-deficient (rank N-1) by design'
       write(LIS_logunit,*) '[WARN] Using fallback: identity transform'

       ! Fallback: Use identity transformation (no update)
       ! Set eigenvalues to safe values
       eigenvalues = 0.0
       eigenvalues(1) = 1.0  ! At least one non-zero eigenvalue

       ! Set A_matrix to identity
       A_matrix = 0.0
       do i = 1, N_ens
          A_matrix(i,i) = 1.0
       end do
    else
       if (gid == 1) write(LIS_logunit,*) '[LNETF] gid=', gid, ' Eigenvalue decomposition successful'
    endif

    ! Print first 3 eigenvalues safely
    if (gid == 1) then
       do i = 1, min(3, N_ens)
          write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Eigenvalue(', i, ')=', eigenvalues(i)
       end do
    endif

    ! Check for small/negative eigenvalues
    ! Note: A should have one zero eigenvalue by design (rank = N_ens - 1)
    do i = 1, N_ens
       if (eigenvalues(i) < 1.0e-15 .or. eigenvalues(i) /= eigenvalues(i)) then
          ! Zero or NaN
          eigenvalues(i) = 0.0
       else if (eigenvalues(i) < 0.0) then
          ! Negative eigenvalue - should not happen for symmetric matrix
          if (gid == 1) write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': Negative eigenvalue', i, eigenvalues(i)
          eigenvalues(i) = 0.0
       else
          eigenvalues(i) = sqrt(eigenvalues(i))
       endif
    end do

    if (gid == 1) then
       do i = 1, min(3, N_ens)
          write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' sqrt(Eigenvalue(', i, '))=', eigenvalues(i)
       end do
    endif

    ! Compute T = eigenvectors * diag(sqrt(eigenvalues))
    do j = 1, N_ens
       do i = 1, N_ens
          T_matrix(i,j) = A_matrix(i,j) * eigenvalues(j)
       end do
    end do

    ! Compute sqrt(A) = T * T^T
    ! Use SGEMM (single precision) to match Fortran REAL type
    call sgemm('N', 'T', N_ens, N_ens, N_ens, 1.0, &
               T_matrix, N_ens, A_matrix, N_ens, 0.0, T_tmp, N_ens)

    ! ========================================================================
    ! STEP 5: Generate random rotation matrix and apply scaling
    ! ========================================================================

    ! Generate orthonormal matrix based on type_trans setting
    ! type_trans = 1: Random orthonormal matrix (stochastic transformation)
    ! type_trans = 2: Identity matrix (deterministic transformation)
    if (type_trans == 1) then
       ! Generate random orthonormal matrix with eigenvector (1,...,1)^T
       ! This preserves ensemble mean while adding random rotations
       call generate_rndmat(N_ens, rndmat)
       if (gid == 1) then
          write(LIS_logunit,*) '[LNETF] Using random rotation matrix (type_trans=1)'
       endif
    else
       ! Deterministic transformation (identity matrix)
       rndmat = 0.0
       do i = 1, N_ens
          rndmat(i,i) = 1.0
       end do
       if (gid == 1) then
          write(LIS_logunit,*) '[LNETF] Using deterministic transformation (type_trans=2)'
       endif
    endif

    ! Multiply T by sqrt(N_ens) to get unbiased ensemble
    fac = sqrt(real(N_ens))

    ! Apply forgetting factor if needed (covariance inflation)
    ! This inflates the analysis ensemble spread
    if (forget < 1.0 .and. forget > 0.0) then
       fac = fac / sqrt(forget)
       if (gid == 1) then
          write(LIS_logunit,*) '[LNETF] Applying forgetting factor:', forget
       endif
    endif

    ! Compute final transformation: T = sqrt(N_ens) * sqrt(A) * rndmat
    ! Use SGEMM (single precision) to match Fortran REAL type
    call sgemm('N', 'N', N_ens, N_ens, N_ens, fac, &
               T_tmp, N_ens, rndmat, N_ens, 0.0, T_matrix, N_ens)

    ! Add weights to each column: W = T + w
    do j = 1, N_ens
       do i = 1, N_ens
          T_matrix(i,j) = T_matrix(i,j) + weights(i)
       end do
    end do

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Transformation matrix W(1,1:5):', T_matrix(1,1:min(5,N_ens))
    endif

    ! ========================================================================
    ! STEP 6: Transform ensemble: X^a = X^f * W
    ! ========================================================================

    ! Use block formulation for transformation to save memory
    maxblksize = 200

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Transforming ensemble with block size', maxblksize
    endif

    allocate(ens_blk(maxblksize, N_ens))

    do blklower = 1, N_state, maxblksize
       blkupper = min(blklower + maxblksize - 1, N_state)

       ! Store forecast ensemble block
       do j = 1, N_ens
          ens_blk(1:(blkupper-blklower+1), j) = State_incr(blklower:blkupper, j)
       end do

       ! Transform: X^a = X^f * W
       ! Use SGEMM (single precision) to match Fortran REAL type
       call sgemm('N', 'N', blkupper-blklower+1, N_ens, N_ens, 1.0, &
                  ens_blk, maxblksize, T_matrix, N_ens, 0.0, &
                  State_incr(blklower:blkupper,1), N_state)
    end do

    deallocate(ens_blk)

    ! Compute ensemble mean and convert to increments
    State_bar = sum(State_incr, 2) / real(N_ens)

    ! Create matrix with State_bar replicated
    do n_e = 1, N_ens
       State_bar_mat(:,n_e) = State_bar
    end do

    ! Compute increments: X^a - X^f
    ! Note: On input, State_incr contained X^f, now it contains X^a
    ! We need to compute the actual increments
    ! This is handled in the calling routine

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Analysis complete'
       write(LIS_logunit,*) '[LNETF] Analysis mean (first 3):', State_bar(1:min(3,N_state))
    endif

    ! ========================================================================
    ! Return diagnostic outputs if requested
    ! ========================================================================
    if (present(n_eff_out)) then
       n_eff_out = n_eff
    endif

    if (present(mean_loc_weight_out)) then
       ! Compute mean of localization weights for this tile
       if (N_obs > 0) then
          mean_loc_weight_out = sum(loc_weights) / real(N_obs)
       else
          mean_loc_weight_out = 0.0
       endif
    endif

    if (present(max_loc_weight_out)) then
       ! Compute max of localization weights (closest observation influence)
       if (N_obs > 0) then
          max_loc_weight_out = maxval(loc_weights)
       else
          max_loc_weight_out = 0.0
       endif
    endif

    if (present(min_loc_weight_out)) then
       ! Compute min of localization weights (farthest observation influence)
       if (N_obs > 0) then
          min_loc_weight_out = minval(loc_weights)
       else
          min_loc_weight_out = 0.0
       endif
    endif

    ! Deallocate arrays
    deallocate(State_prime)
    deallocate(State_bar_mat)
    deallocate(State_bar)
    deallocate(innov_i)
    deallocate(obs_values)
    deallocate(loc_weights)

  end subroutine lnetf_analysis


!BOP
!
! !ROUTINE: compute_likelihood_gaussian
! \label{compute_likelihood_gaussian}
!
! !INTERFACE:
  subroutine compute_likelihood_gaussian(N_obs, innov, Obs_cov, likelihood)
! !USES:
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: N_obs
    real, intent(in), dimension(N_obs) :: innov
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov
    real, intent(out) :: likelihood
!
! !DESCRIPTION:
!   Compute Gaussian likelihood: exp(-0.5 * innov^T * R^{-1} * innov)
!
!   For simplicity, assumes diagonal R (observation errors are uncorrelated)
!   This can be extended to full R matrix using Cholesky decomposition
!EOP

    integer :: i
    real :: maha_dist  ! Mahalanobis distance

    ! Simplified version assuming diagonal Obs_cov
    maha_dist = 0.0

    do i = 1, N_obs
       if (Obs_cov(i,i) > 1.0e-15) then
          maha_dist = maha_dist + (innov(i)**2) / Obs_cov(i,i)
       endif
    end do

    ! Gaussian likelihood
    likelihood = exp(-0.5 * maha_dist)

  end subroutine compute_likelihood_gaussian


!BOP
!
! !ROUTINE: compute_likelihood_localized
! \label{compute_likelihood_localized}
!
! !INTERFACE:
  subroutine compute_likelihood_localized(N_obs, innov, Obs_cov, loc_weights, log_likelihood)
! !USES:
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: N_obs
    real, intent(in), dimension(N_obs) :: innov
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov
    real, intent(in), dimension(N_obs) :: loc_weights  ! Gaspari-Cohn weights
    real, intent(out) :: log_likelihood  ! Changed to LOG-likelihood to avoid underflow
!
! !DESCRIPTION:
!   Compute localized Gaussian LOG-likelihood with distance-based weighting
!
!   log_likelihood = -0.5 * sum_i( w_i * innov_i^2 / R_ii )
!
!   where w_i is the Gaspari-Cohn localization weight for observation i
!
!   IMPORTANT: Returns LOG-likelihood to avoid numerical underflow.
!   The calling routine should use the log-sum-exp trick to normalize weights.
!
!   This follows the PDAF approach where localization weights are applied
!   to the observation errors in the likelihood calculation.
!
!   References:
!   - PDAF likelihood_l_pdaf.F90 template
!   - Gaspari and Cohn (1999), Q.J.R. Meteorol. Soc.
!EOP

    integer :: i
    real :: maha_dist  ! Localized Mahalanobis distance

    ! Compute localized Mahalanobis distance
    maha_dist = 0.0

    do i = 1, N_obs
       if (Obs_cov(i,i) > 1.0e-15 .and. loc_weights(i) > 1.0e-10) then
          ! Apply localization weight to the likelihood contribution
          ! Weight acts on the observation error variance: R_eff = R / w
          ! This makes distant observations have less influence
          maha_dist = maha_dist + loc_weights(i) * (innov(i)**2) / Obs_cov(i,i)
       endif
    end do

    ! Return LOG-likelihood to avoid underflow
    ! (exp(-0.5 * maha_dist) can underflow for large maha_dist)
    log_likelihood = -0.5 * maha_dist

  end subroutine compute_likelihood_localized


!BOP
!
! !ROUTINE: compute_effective_sample_size
! \label{compute_effective_sample_size}
!
! !INTERFACE:
  subroutine compute_effective_sample_size(N_ens, weights, n_eff)
! !USES:
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: N_ens
    real, intent(in), dimension(N_ens) :: weights
    real, intent(out) :: n_eff
!
! !DESCRIPTION:
!   Compute effective sample size: N_eff = 1 / sum(w_i^2)
!
!   This diagnostic indicates how many "effective" ensemble members
!   contribute to the analysis. Low N_eff suggests weight collapse.
!EOP

    real :: sum_w_squared
    integer :: i

    sum_w_squared = 0.0
    do i = 1, N_ens
       sum_w_squared = sum_w_squared + weights(i)**2
    end do

    if (sum_w_squared > 1.0e-15) then
       n_eff = 1.0 / sum_w_squared
    else
       n_eff = 0.0
    endif

  end subroutine compute_effective_sample_size


!BOP
! !ROUTINE: gaspari_cohn
! \label{gaspari_cohn_lnetf}
!
! !INTERFACE:
  function gaspari_cohn( d )
!
! !DESCRIPTION:
! Evaluate 5th-order polynomial from Gaspari and Cohn, 1999, Eq (4.10)
!
! On input, d = separation distance relative to the distance at which
! all correlations vanish. In the isotropic case, Gaspari and Cohn, 1999,
! Eq. (4.10)
!
! \begin{verbatim}
!    d = sqrt(dx**2 + dy**2) / (2*c) = |z| / (2*c)
! \end{verbatim}
!
! or in the anisotropic case
!
! \begin{verbatim}
!    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
! \end{verbatim}
!
! \begin{verbatim}
! *** Use |z|/c = 2*d. All correlations vanish for d > 1. ***
! \end{verbatim}
!
! This function is identical to the one in enkf_general.F90
!EOP

    implicit none

    real :: gaspari_cohn, d, y

    real, parameter :: tol = 1e-3

    ! Get rid of possibly negative distances.
    ! Multiply with 2. to return to the Gaspari and Cohn, 1999, notation.

    d = 2.*abs(d)

    if (d >= 2.) then

       y = 0.

    else if (d <= tol) then

       y = 1.

    else if (d <= 1.) then

       ! y = -.25*d**5 + .5*d**4 + .625*d**3 - 5./3.*d**2 + 1.

       y = d**2 *( d*( d*( -.25*d + .5) + .625) -5./3.) + 1.

    else

       ! y = d**5/12. - .5*d**4 + .625*d**3 + 5./3.*d**2 - 5.*d + 4. - 2./3./d

       y = d*( d*( d*( d*( d/12. - .5) + .625) + 5./3.) -5.) + 4. - 2./3./d

    end if

    gaspari_cohn = y

    return

  end function gaspari_cohn


!BOP
!
! !ROUTINE: get_gaspari_cohn
!  \label{get_gaspari_cohn_lnetf}
!
! !INTERFACE:
  function get_gaspari_cohn( dx, dy, xcompact, ycompact )

! !DESCRIPTION:
! For a given lat/lon distance, compute the
! anisotropic compact support (Gaspari and Cohn) weights
!
! Input coordinates must be in degrees latitude/longitude
!
!  dx = longitude separation of two points \newline
!  dy = latitude separation of two points
!
!  xcompact = longitude scale of compact support \newline
!  ycompact = latitude scale of compact support
!
! For the isotropic Gaspari and Cohn function, the relative
! distance is in Gaspari and Cohn, 1999, notation (Eq (4.10))
!
!      d = sqrt(dx**2 + dy**2) / 2*c
!
! get\_gaspari\_cohn() uses a generalized anisotropic Gaspari and Cohn
! approach (essentially coordinate stretching)
!
!  d = sqrt((dx/xcompact)**2 + (dy/ycompact)**2  )
!
! All correlations vanish outside of an ellipse with semi-axes
! xcompact and ycompact, ie Gaspari and Cohn weights vanish
! for d > 1 (note the factor 2!)
!
! When the anisotropic case is reduced back to the isotropic case,
! (ie if xcompact==ycompact) then c = xcompact/2 = ycompact/2.
!
! This function is identical to the one in enkf_general.F90
!EOP

    implicit none

    real :: get_gaspari_cohn, dx, dy, xcompact, ycompact, d

    ! Compute (anisotropic) distance relative to compact support

    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
    get_gaspari_cohn = gaspari_cohn(d)

  end function get_gaspari_cohn


!BOP
!
! !ROUTINE: inflate_weights
! \label{inflate_weights}
!
! !INTERFACE:
  subroutine inflate_weights(dim_ens, alpha, weights)
! !USES:
    use LIS_logMod, only : LIS_logunit
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: dim_ens           ! Ensemble size
    real, intent(in) :: alpha                ! Minimum limit of N_eff / N
    real, intent(inout) :: weights(dim_ens)  ! Weights (before and after inflation)
!
! !DESCRIPTION:
!   Inflate particle weights to achieve N_eff/N > alpha
!
!   This routine is adapted from PDAF_inflate_weights (PDAF by Lars Nerger).
!   It uses an iterative approach to find the optimal inflation factor
!   that brings N_eff/N above the specified threshold alpha.
!
!   The inflation is applied as: w_new = w^(1-a) * (1/N)^a
!   where a is the inflation factor (0 <= a <= 1)
!
!   Algorithm:
!   1. Start with a=0 (no inflation)
!   2. Incrementally increase a until N_eff/N >= alpha
!   3. At a=1, weights become uniform (1/N)
!
! !REFERENCES:
!   PDAF - Parallel Data Assimilation Framework
!   https://pdaf.awi.de
!EOP

    ! Local variables
    integer :: i
    real :: alpha_iter
    real :: a_step, tot_weight
    real :: alpha_lim, n_eff
    real, allocatable :: logw(:)      ! Logarithm of particle weights
    real, allocatable :: aweights(:)  ! Temporary weights

    ! Allocate local arrays
    allocate(logw(dim_ens))
    allocate(aweights(dim_ens))

    ! Get logarithm of weights
    do i = 1, dim_ens
       if (weights(i) > 1.0e-30) then
          logw(i) = log(weights(i))
       else
          logw(i) = -69.0  ! log(1e-30)
       endif
    end do

    ! Store initial weights
    aweights = weights

    write(LIS_logunit,*) '[LNETF] Inflating weights to achieve N_eff/N >', alpha

    ! Set limit: target N_eff
    alpha_lim = alpha * real(dim_ens)

    ! Iteratively find inflation factor
    alpha_iter = 0.0
    a_step = 0.05

    inflation_loop: do
       if (alpha_iter >= 1.0) then
          alpha_iter = 1.0
          exit inflation_loop
       endif

       ! Apply inflation: w_new = w^(1-a) * (1/N)^a
       ! In log space: log(w_new) = (1-a)*log(w) + a*log(1/N)
       do i = 1, dim_ens
          aweights(i) = exp((1.0 - alpha_iter) * logw(i) &
                           + alpha_iter * log(1.0/real(dim_ens)))
       end do

       ! Normalize weights
       tot_weight = sum(aweights)
       if (tot_weight > 1.0e-15) then
          aweights = aweights / tot_weight
       else
          aweights = 1.0 / real(dim_ens)
       endif

       ! Compute effective sample size
       n_eff = 0.0
       do i = 1, dim_ens
          n_eff = n_eff + aweights(i)**2
       end do
       if (n_eff > 1.0e-15) then
          n_eff = 1.0 / n_eff
       else
          n_eff = 0.0
       endif

       ! Check if target is reached
       if (n_eff >= alpha_lim) then
          exit inflation_loop
       endif

       ! Increase inflation factor
       alpha_iter = alpha_iter + a_step

    end do inflation_loop

    ! Apply final inflated weights
    weights = aweights

    write(LIS_logunit,*) '[LNETF] Inflation factor used:', alpha_iter
    write(LIS_logunit,*) '[LNETF] Final N_eff:', n_eff, ' (target:', alpha_lim, ')'

    deallocate(logw, aweights)

  end subroutine inflate_weights


!BOP
!
! !ROUTINE: generate_rndmat
! \label{generate_rndmat}
!
! !INTERFACE:
  subroutine generate_rndmat(dim, rndmat)
! !USES:
    use LIS_logMod, only : LIS_logunit
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: dim              ! Size of matrix
    real, intent(out) :: rndmat(dim, dim)   ! Random orthonormal matrix
!
! !DESCRIPTION:
!   Generate a random orthonormal matrix with eigenvector (1,...,1)^T
!
!   This routine is adapted from PDAF_generate_rndmat (PDAF by Lars Nerger).
!   It uses Householder reflections to construct an orthonormal matrix
!   that has (1,...,1)^T as an eigenvector with eigenvalue 1.
!
!   This property ensures that the ensemble mean is preserved during
!   the transformation: if x_mean = sum(x_i)/N, then
!   after transformation y = x * rndmat, y_mean = x_mean
!
!   The algorithm:
!   1. Generate (dim-1) random orthonormal matrix using Householder
!   2. Embed it into dim x dim matrix with (1,...,1)^T as first eigenvector
!
! !REFERENCES:
!   PDAF - Parallel Data Assimilation Framework
!   https://pdaf.awi.de
!EOP

    ! Local variables
    integer :: i, j, k, iter
    integer :: dimrnd
    real :: norm, rndval
    real, allocatable :: rndvec(:)
    real, allocatable :: house(:,:)
    real, allocatable :: temp1(:,:), temp2(:,:)
    real, allocatable :: matU(:,:), matB(:,:)
    real, allocatable :: matUB(:,:), matUBB(:,:)
    integer, save :: iseed(4) = (/1, 5, 7, 9/)
    integer, save :: first = 1

    ! Initialize random seed on first call
    if (first == 1) then
       ! Use system-dependent seed initialization
       call random_seed()
       first = 0
    endif

    ! For matrix with eigenvector (1,...,1)^T, we build (dim-1) random part
    dimrnd = dim - 1

    if (dimrnd < 1) then
       ! For dim=1, just return identity
       rndmat(1,1) = 1.0
       return
    endif

    ! Allocate arrays
    allocate(rndvec(dim))
    allocate(house(dim+1, dim))
    allocate(temp1(dim, dim), temp2(dim, dim))
    allocate(matU(dimrnd, dimrnd))
    allocate(matB(dim, dimrnd))
    allocate(matUB(dim, dimrnd))
    allocate(matUBB(dim, dim))

    ! ================================================================
    ! Step 1: Generate random orthonormal matrix of size (dim-1)
    ! Using Householder reflections
    ! ================================================================

    ! Initialize with identity
    matU = 0.0
    do i = 1, dimrnd
       matU(i,i) = 1.0
    end do

    ! Apply sequence of Householder reflections
    do iter = 1, dimrnd - 1
       ! Generate random vector
       do i = 1, dimrnd - iter + 1
          call random_number(rndval)
          rndvec(i) = rndval - 0.5
       end do

       ! Normalize
       norm = 0.0
       do i = 1, dimrnd - iter + 1
          norm = norm + rndvec(i)**2
       end do
       norm = sqrt(norm)
       if (norm > 1.0e-10) then
          do i = 1, dimrnd - iter + 1
             rndvec(i) = rndvec(i) / norm
          end do
       endif

       ! Build Householder matrix H = I - 2*v*v^T
       ! and apply to matU
       do j = iter, dimrnd
          rndval = 0.0
          do k = iter, dimrnd
             rndval = rndval + rndvec(k-iter+1) * matU(k, j)
          end do
          rndval = 2.0 * rndval
          do i = iter, dimrnd
             matU(i, j) = matU(i, j) - rndval * rndvec(i-iter+1)
          end do
       end do
    end do

    ! ================================================================
    ! Step 2: Embed into dim x dim matrix with (1,...,1)^T eigenvector
    ! ================================================================

    ! Create matrix B that maps (dim-1) space to orthogonal complement of (1,...,1)
    ! B is dim x (dim-1), with columns orthogonal to (1,...,1)^T
    matB = 0.0
    do j = 1, dimrnd
       do i = 1, dim
          if (i == j) then
             matB(i, j) = 1.0 - 1.0/real(dim)
          else if (i == dim) then
             matB(i, j) = -1.0/real(dim)
          else if (i > j) then
             matB(i, j) = -1.0/real(dim)
          endif
       end do
       ! Normalize column
       norm = 0.0
       do i = 1, dim
          norm = norm + matB(i,j)**2
       end do
       norm = sqrt(norm)
       if (norm > 1.0e-10) then
          do i = 1, dim
             matB(i,j) = matB(i,j) / norm
          end do
       endif
    end do

    ! Apply Gram-Schmidt to ensure orthonormality
    do j = 2, dimrnd
       do k = 1, j-1
          rndval = 0.0
          do i = 1, dim
             rndval = rndval + matB(i,j) * matB(i,k)
          end do
          do i = 1, dim
             matB(i,j) = matB(i,j) - rndval * matB(i,k)
          end do
       end do
       ! Normalize
       norm = 0.0
       do i = 1, dim
          norm = norm + matB(i,j)**2
       end do
       norm = sqrt(norm)
       if (norm > 1.0e-10) then
          do i = 1, dim
             matB(i,j) = matB(i,j) / norm
          end do
       endif
    end do

    ! Compute matUB = matB * matU (dim x dimrnd * dimrnd x dimrnd = dim x dimrnd)
    matUB = 0.0
    do j = 1, dimrnd
       do i = 1, dim
          do k = 1, dimrnd
             matUB(i,j) = matUB(i,j) + matB(i,k) * matU(k,j)
          end do
       end do
    end do

    ! Compute rndmat = matUB * matB^T + (1/dim) * ones
    ! This gives orthonormal matrix with (1,...,1)^T as eigenvector
    rndmat = 0.0
    do j = 1, dim
       do i = 1, dim
          do k = 1, dimrnd
             rndmat(i,j) = rndmat(i,j) + matUB(i,k) * matB(j,k)
          end do
          ! Add projection onto (1,...,1)^T
          rndmat(i,j) = rndmat(i,j) + 1.0/real(dim)
       end do
    end do

    deallocate(rndvec, house, temp1, temp2)
    deallocate(matU, matB, matUB, matUBB)

  end subroutine generate_rndmat


!BOP
!
! !ROUTINE: haversine_km
! \label{haversine_km}
!
! !INTERFACE:
  function haversine_km(lon1, lat1, lon2, lat2)
!
! !DESCRIPTION:
!   Calculate the great-circle distance between two points on Earth
!   using the Haversine formula.
!
!   This is used for accurate km-based localization following
!   Seo et al. (2021) and Tak et al. (2025) LETKF implementations.
!
!   Input: longitudes and latitudes in degrees
!   Output: distance in kilometers
!
! !REFERENCES:
!   Seo, E., et al. (2021): Assimilation of SMAP and ASCAT soil moisture
!   retrievals into the JULES land surface model using the LETKF.
!   Tak, Y.-J., et al. (2025): Multi-Sensor Data Assimilation for
!   Global Soil Moisture Estimation Using the LETKF.
!EOP

    implicit none

    real :: haversine_km
    real, intent(in) :: lon1, lat1, lon2, lat2

    real, parameter :: PI = 3.14159265358979323846
    real, parameter :: DEG_TO_RAD = PI / 180.0
    real, parameter :: EARTH_RADIUS_KM = 6371.0

    real :: lat1_rad, lat2_rad, dlat, dlon
    real :: a, c

    ! Convert to radians
    lat1_rad = lat1 * DEG_TO_RAD
    lat2_rad = lat2 * DEG_TO_RAD
    dlat = (lat2 - lat1) * DEG_TO_RAD
    dlon = (lon2 - lon1) * DEG_TO_RAD

    ! Haversine formula
    a = sin(dlat/2.0)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))

    haversine_km = EARTH_RADIUS_KM * c

  end function haversine_km


!BOP
!
! !ROUTINE: get_gaspari_cohn_km
! \label{get_gaspari_cohn_km}
!
! !INTERFACE:
  function get_gaspari_cohn_km(lon1, lat1, lon2, lat2, sigma_km)
!
! !DESCRIPTION:
!   Compute Gaspari-Cohn localization weight based on distance in km.
!
!   This is the km-based version for LETKF-style localization
!   following Seo et al. (2021) and Tak et al. (2025).
!
!   The localization weight is computed as:
!   1. Calculate haversine distance between points (in km)
!   2. Normalize by compact support radius (3*sigma_km)
!   3. Apply Gaspari-Cohn 5th-order polynomial
!
!   For sigma_km = 30 km (paper default):
!   - At 30 km: weight ≈ 0.7
!   - At 60 km: weight ≈ 0.2
!   - At 90 km (compact support): weight = 0
!
!   This results in "almost 1D filter" behavior where distant
!   observations have very little influence.
!
! !REFERENCES:
!   Gaspari, G. and Cohn, S.E. (1999): Construction of correlation
!   functions in two and three dimensions. Q.J.R. Meteorol. Soc.
!EOP

    implicit none

    real :: get_gaspari_cohn_km
    real, intent(in) :: lon1, lat1, lon2, lat2
    real, intent(in) :: sigma_km

    real :: dist_km, compact_km, d

    ! Calculate actual distance in km
    dist_km = haversine_km(lon1, lat1, lon2, lat2)

    ! Compact support cutoff at 3*sigma
    ! (weight = 0 at this distance and beyond)
    compact_km = 3.0 * sigma_km

    ! Normalized distance
    d = dist_km / compact_km

    ! Apply Gaspari-Cohn function
    get_gaspari_cohn_km = gaspari_cohn(d)

  end function get_gaspari_cohn_km


end module lnetf_general
