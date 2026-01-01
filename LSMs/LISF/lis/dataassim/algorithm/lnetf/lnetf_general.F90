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
! !REFERENCES:
!   Toedter, J. and Ahrens, B., 2015: A Second-Order Exact Ensemble Square
!   Root Filter for Nonlinear Data Assimilation. Mon. Wea. Rev., 143, 1347-1367.
!
!   PDAF implementation by Lars Nerger and Paul Kirchgessner
!
! !REVISION HISTORY:
!   27 Dec 2024: Initial implementation adapted from PDAF
!EOP
module lnetf_general

  implicit none

  private

  public :: lnetf_analysis
  public :: compute_likelihood_gaussian
  public :: compute_effective_sample_size

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
       State_lon, State_lat, xcompact, ycompact )
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

!
! !DESCRIPTION:
!
! Perform LNETF update using particle weights and ensemble transformation
!
! ALGORITHM STEPS:
! 1. Compute particle weights (likelihood) for each ensemble member
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
! are present, localization is applied
!EOP

    ! Local variables
    integer :: n_e, i, j, k
    real :: total_weight, weight, n_eff
    real :: fac, dx, dy, dweight
    logical :: apply_localization

    ! Large arrays - allocate on heap to avoid stack overflow
    real, allocatable, dimension(:,:) :: State_prime, State_bar_mat
    real, allocatable, dimension(:) :: State_bar

    real, allocatable, dimension(:) :: innov_i, obs_values

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

    ! Check if localization should be applied
    apply_localization = (present(State_lon) .and. present(State_lat) &
                         .and. present(xcompact) .and. present(ycompact))

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Starting LNETF analysis for grid point', gid
       write(LIS_logunit,*) '[LNETF] N_state=', N_state, ' N_obs=', N_obs, ' N_ens=', N_ens
    endif

    ! ========================================================================
    ! STEP 1: Compute particle weights (likelihoods)
    ! ========================================================================

    ! Extract observation values
    do i = 1, N_obs
       obs_values(i) = Observations(i)%value
    end do

    ! Compute likelihood for each ensemble member
    do n_e = 1, N_ens
       ! Calculate innovation: obs - H(x_i)
       innov_i = obs_values - Obs_pred(:, n_e)

       ! Compute Gaussian likelihood
       call compute_likelihood_gaussian(N_obs, innov_i, Obs_cov, weight)

       weights(n_e) = weight
    end do

    ! Debug: Show observation predictions and innovations for first few members
    if (gid == 16) then
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

    write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Raw weights:', weights(1:min(3,N_ens))

    ! ========================================================================
    ! STEP 2: Normalize weights
    ! ========================================================================

    total_weight = sum(weights)
    write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Total weight:', total_weight

    if (total_weight > 1.0e-15) then
       weights = weights / total_weight
    else
       ! ERROR: All weights are zero - use equal weights
       write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': Zero total weight - using equal weights'
       weights = 1.0 / real(N_ens)
    endif

    write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Normalized weights:', weights(1:min(3,N_ens))

    ! ========================================================================
    ! STEP 2.5: Compute effective sample size
    ! ========================================================================

    call compute_effective_sample_size(N_ens, weights, n_eff)

    if (gid == 1) then
       write(LIS_logunit,*) '[LNETF] Effective sample size: N_eff =', n_eff, &
                            ' (', 100.0*n_eff/real(N_ens), '%)'
    endif

    ! TODO: Implement weights inflation if N_eff/N_ens > limit_winf
    ! Currently skipped for initial implementation

    ! ========================================================================
    ! STEP 3: Compute transform matrix A = diag(w) - w*w^T
    ! ========================================================================

    if (gid == 16) write(LIS_logunit,*) '[LNETF-DEBUG] gid=', gid, ' Building transform matrix A'

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
    if (gid == 16) then
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
       write(LIS_logunit,*) '[LNETF] gid=', gid, ' Eigenvalue decomposition successful'
    endif

    ! Print first 3 eigenvalues safely
    if (gid == 16) then
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
          if (gid == 16) write(LIS_logunit,*) '[WARN] LNETF gid=', gid, ': Negative eigenvalue', i, eigenvalues(i)
          eigenvalues(i) = 0.0
       else
          eigenvalues(i) = sqrt(eigenvalues(i))
       endif
    end do

    if (gid == 16) then
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

    ! Generate random orthonormal matrix with eigenvector (1,...,1)^T
    ! For now, use identity matrix (deterministic transformation)
    ! TODO: Implement proper random rotation matrix generation
    rndmat = 0.0
    do i = 1, N_ens
       rndmat(i,i) = 1.0
    end do

    ! Multiply T by sqrt(N_ens) to get unbiased ensemble
    fac = sqrt(real(N_ens))

    ! Apply forgetting factor if needed (covariance inflation)
    ! fac = fac / sqrt(forget) for analysis inflation
    ! Currently using default forget = 1.0

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

    ! Deallocate arrays
    deallocate(State_prime)
    deallocate(State_bar_mat)
    deallocate(State_bar)
    deallocate(innov_i)
    deallocate(obs_values)

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

end module lnetf_general
