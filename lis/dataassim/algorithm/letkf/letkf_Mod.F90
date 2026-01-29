!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module letkf_Mod
!BOP
!
! !MODULE: letkf_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines that control
!   the incorporation of a data set using the Local Ensemble Transform
!   Kalman Filter (LETKF) method, into a land surface model.
!
!  The LETKF algorithm is based on the work of:
!  Hunt, B.R., Kostelich, E.J. and Szunyogh, I., 2007: Efficient data
!  assimilation for spatiotemporal chaos: A local ensemble transform
!  Kalman filter. Physica D, 230, 112-126.
!
!  Implementation adapted from PDAF (Parallel Data Assimilation Framework)
!  by Lars Nerger, and LIS LNETF implementation.
!
!  Key differences from LNETF:
!  - Uses linear Kalman gain (not likelihood-based weights)
!  - Gaussian assumption (LNETF handles non-Gaussian better)
!  - Deterministic square-root filter (no perturbed observations)
!
!  NOTES: Data assimilation is currently only supported for land surface
!  models (and not across different surface model types)
!
! !REVISION HISTORY:
!   28 Jan 2026: Initial implementation adapted from LNETF and PDAF
!
! !USES:
  use ESMF
  use letkf_types
  use letkf_general
  use my_matrix_functions
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_surfaceModelMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: letkf_init  ! Initialization for LETKF
  public :: letkf_setup
  public :: letkf_increments ! compute analysis increments
  public :: letkf_update ! apply the analysis increments
  public :: letkf_diagnostics ! write LETKF related diagnostics
  public :: letkf_final ! Finalization for LETKF
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: letkf_struc ! data structure containing LETKF diagnostics
!EOP

  type, public ::  letkf_dec
     logical     :: fileOpen
     real, allocatable :: innov(:)
     real, allocatable :: forecast_var(:) !HPHt
     real, allocatable :: anlys_res(:)
     real, allocatable :: anlys_incr(:,:)
     real, allocatable :: norm_innov(:)
     real, allocatable :: transform(:,:)  ! ensemble transformation matrix
     real :: localization_factor = 5.0    ! for localization (grid-relative)

     ! Localization parameters (km-based, following Seo et al. 2021 / Tak et al. 2025)
     ! Reference: Seo, E., et al. (2021) - JULES + LETKF with σ=30km, patch=150km
     !            Tak, Y.-J., et al. (2025) - Multi-sensor SMAP+ASCAT LETKF
     real :: localization_scale_km = 30.0   ! Gaspari-Cohn scale σ (km)
     logical :: use_km_localization = .false.

     ! LETKF specific parameters (Hunt et al. 2007)
     integer :: type_forget = 0           ! Type of forgetting factor
     integer :: type_trans = 0            ! Type of ensemble transformation (0=deterministic, 2=random)
     real :: forget = 1.0                 ! Forgetting factor (covariance inflation)

     !--------------------------------------------------------------------------
     ! Localization diagnostic arrays (for output to increment file)
     !--------------------------------------------------------------------------
     real, allocatable :: n_local_obs(:)
     real, allocatable :: mean_loc_weight(:)
     real, allocatable :: max_loc_weight(:)
     real, allocatable :: min_loc_weight(:)

  end type letkf_dec
!EOP

  type(letkf_dec), allocatable :: letkf_struc(:,:)

contains

!BOP
!
! !ROUTINE: letkf_init
! \label{letkf_init}
!
! !INTERFACE:
  subroutine letkf_init()
! !USES:

!
! !DESCRIPTION:
!  This method performs the required initializations for the
!  LETKF method.
!
!EOP
    allocate(letkf_struc(LIS_rc%nnest, LIS_rc%ndas))

  end subroutine letkf_init

!BOP
!
! !ROUTINE: letkf_setup
! \label{letkf_setup}
!
! !INTERFACE:
  subroutine letkf_setup(k)
! !USES:

!
! !DESCRIPTION:
!  This method performs the required setup for the
!  LETKF method. The method reads the runtime settings from
!  the LIS configuration file.
!
!EOP
    integer                      :: n
    integer                      :: k
    integer                      :: status
    integer                      :: Nobjs
    integer                      :: N_obs_size

    do n=1,LIS_rc%nnest
       if(LIS_rc%nensem(n).le.1) then
          write(LIS_logunit,*) '[ERR] Please set number of ensembles '
          write(LIS_logunit,*) '[ERR] to greater than 1 for LETKF',LIS_rc%nensem(n)
          call LIS_endrun
       endif
       letkf_struc(n,k)%fileOpen = .false.
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed in letkf_Mod')

       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet failed in letkf_Mod')

       if(LIS_rc%winnov(k).eq.1) then
          allocate(letkf_struc(n,k)%norm_innov(Nobjs*N_obs_size))
          allocate(letkf_struc(n,k)%innov(Nobjs*N_obs_size))
          allocate(letkf_struc(n,k)%anlys_res(Nobjs*N_obs_size))
          allocate(letkf_struc(n,k)%forecast_var(Nobjs*N_obs_size))

          ! Initialize innovation arrays to undefined
          letkf_struc(n,k)%norm_innov = LIS_rc%udef
          letkf_struc(n,k)%innov = LIS_rc%udef
          letkf_struc(n,k)%anlys_res = LIS_rc%udef
          letkf_struc(n,k)%forecast_var = LIS_rc%udef
       endif
       allocate(letkf_struc(n,k)%anlys_incr(LIS_rc%nstvars(k),&
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)))

       !----------------------------------------------------------------------
       ! Allocate localization diagnostic arrays
       !----------------------------------------------------------------------
       allocate(letkf_struc(n,k)%n_local_obs( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(letkf_struc(n,k)%mean_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(letkf_struc(n,k)%max_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(letkf_struc(n,k)%min_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))

       ! Initialize to undefined
       letkf_struc(n,k)%n_local_obs = LIS_rc%udef
       letkf_struc(n,k)%mean_loc_weight = LIS_rc%udef
       letkf_struc(n,k)%max_loc_weight = LIS_rc%udef
       letkf_struc(n,k)%min_loc_weight = LIS_rc%udef
    enddo

!----------------------------------------------------------------------------
! Read localization parameters from config file
!----------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       ! First try km-based localization (preferred, following papers)
       call ESMF_ConfigGetAttribute(LIS_config, &
            letkf_struc(n,k)%localization_scale_km, &
            label="LETKF localization scale (km):", rc=status)
       if(status.eq.0) then
          letkf_struc(n,k)%use_km_localization = .true.
          write(LIS_logunit,*) '[INFO] =============================================='
          write(LIS_logunit,*) '[INFO] LETKF using km-based localization'
          write(LIS_logunit,*) '[INFO]   Gaspari-Cohn scale (σ): ', &
               letkf_struc(n,k)%localization_scale_km, ' km'
          write(LIS_logunit,*) '[INFO]   Compact support cutoff (3σ): ', &
               3.0 * letkf_struc(n,k)%localization_scale_km, ' km'
          write(LIS_logunit,*) '[INFO] Reference: Hunt et al. (2007), Seo et al. (2021)'
          write(LIS_logunit,*) '[INFO] =============================================='
       else
          ! Fall back to factor-based localization
          call ESMF_ConfigGetAttribute(LIS_config, &
               letkf_struc(n,k)%localization_factor, &
               label="LETKF localization radius factor:", rc=status)
          if(status.ne.0) then
             letkf_struc(n,k)%localization_factor = 5.0
             write(LIS_logunit,*) '[INFO] LETKF localization radius factor not found in config'
             write(LIS_logunit,*) '[INFO] Using default factor: 5.0'
          else
             write(LIS_logunit,*) '[INFO] LETKF localization radius factor set to:', &
                  letkf_struc(n,k)%localization_factor
          endif
          letkf_struc(n,k)%use_km_localization = .false.
       endif

       ! Read forgetting factor (covariance inflation)
       call ESMF_ConfigGetAttribute(LIS_config, &
            letkf_struc(n,k)%forget, &
            label="LETKF forgetting factor:", rc=status)
       if(status.ne.0) then
          letkf_struc(n,k)%forget = 1.0
       else
          write(LIS_logunit,*) '[INFO] LETKF forgetting factor set to:', &
               letkf_struc(n,k)%forget
       endif

       ! Read transformation type
       call ESMF_ConfigGetAttribute(LIS_config, &
            letkf_struc(n,k)%type_trans, &
            label="LETKF transformation type:", rc=status)
       if(status.ne.0) then
          letkf_struc(n,k)%type_trans = 0  ! deterministic
       else
          write(LIS_logunit,*) '[INFO] LETKF transformation type set to:', &
               letkf_struc(n,k)%type_trans
       endif
    enddo

  end subroutine letkf_setup

!BOP
!
! !ROUTINE: letkf_increments
! \label{letkf_increments}
!
! !INTERFACE:
  subroutine letkf_increments(n,k)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
    integer, intent(IN)    :: k
!
! !DESCRIPTION:
!  This routine computes the analysis increments for LETKF using
!  PDAF-style 3D spatial localization.
!
!  LETKF Algorithm (Hunt et al. 2007):
!  1. Compute ensemble perturbations
!  2. For each tile, find observations within localization radius
!  3. Compute A^{-1} = (N-1)/rho * I + (HZ)^T R^{-1} HZ
!  4. Eigenvalue decomposition for sqrt(A)
!  5. Compute weight vector w and transform matrix W
!  6. Transform ensemble: X^a = mean(X^f) + Z * W
!
!EOP
    logical                           :: data_status
    integer                           :: status
    integer                           :: Nobjs
    integer                           :: state_size
    integer                           :: Nobs
    integer                           :: N_obs_size
    integer                           :: N_obs_actual
    integer                           :: N_selected_obs
    integer                           :: N_local_obs
    integer                           :: N_ens
    integer                           :: N_state
    type(obs_type), allocatable       :: Observations(:)
    type(obs_type), allocatable       :: Observations_filtered(:)
    type(obs_type), allocatable       :: obs_da(:)
    type(obs_param_type), allocatable :: obs_param(:)
    real,         allocatable         :: Obs_pred(:,:)
    real,         allocatable         :: Obs_pred_filtered(:,:)
    real,         allocatable         :: obspred_da(:,:)
    real,         allocatable         :: Obs_pert(:,:)
    real,         allocatable         :: Obs_pert_filtered(:,:)
    real,         allocatable         :: obspert_da(:,:)
    real,         allocatable         :: Obs_cov(:,:)
    integer                           :: i,j,v,tileid
    integer                           :: st_id, en_id, sid,eid
    real,         allocatable         :: stvar(:,:)
    real,         pointer             :: stdata(:)
    real,         pointer             :: stincrdata(:)
    real,         allocatable         :: state_incr(:,:)
    real                              :: innov,std_innov(1)
    integer                           :: kk,m,p
    logical                           :: assim, obspred_flag
    integer                           :: gid, t

    integer                           :: jj, mm
    real,         allocatable         :: state_tmp(:,:)
    real                              :: dx,dy,xcompact,ycompact
    real,         allocatable         :: lons(:), lats(:)
    real,         allocatable         :: state_lat(:), state_lon(:)

    ! PDAF-style Local LETKF variables
    real,         allocatable         :: obs_lon_arr(:), obs_lat_arr(:)
    real                              :: max_dist, dist, tile_lon, tile_lat
    integer,      allocatable         :: local_obs_idx(:)
    integer                           :: local_count

    ! km-based localization variables
    real                              :: mean_lat, km_per_deg_lon, km_per_deg_lat
    real                              :: sigma_km, xcompact_km, ycompact_km
    real, parameter                   :: PI = 3.14159265358979323846
    real, parameter                   :: DEG_TO_RAD = PI / 180.0
    real, parameter                   :: EARTH_RADIUS_KM = 6371.0
    logical                           :: use_km_loc

    ! Diagnostic counters
    integer                           :: tiles_with_obs, tiles_no_obs
    integer                           :: total_local_obs_count
    integer                           :: min_local_obs, max_local_obs
    real                              :: avg_local_obs

    integer                           :: n_tiles


!----------------------------------------------------------------------------
!  Check if the observation state is updated or not.
!----------------------------------------------------------------------------

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Update Status",&
         value=data_status,rc=status)
    call LIS_verify(status, &
         'ESMF_AttributeGet: Data Update Status failed in letkf_increments')

    call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.false.)

    letkf_struc(n,k)%anlys_incr = 0.0

    ! Initialize localization diagnostic arrays to undefined
    letkf_struc(n,k)%n_local_obs = LIS_rc%udef
    letkf_struc(n,k)%mean_loc_weight = LIS_rc%udef
    letkf_struc(n,k)%max_loc_weight = LIS_rc%udef
    letkf_struc(n,k)%min_loc_weight = LIS_rc%udef

    if(data_status) then
       write(LIS_logunit,*) &
            '[INFO] Assimilating Observations using LETKF for DA instance',k

       ! Get observation parameters
       allocate(obs_param(LIS_rc%nobtypes(k)))
       call LIS_surfaceModel_DAGetObsParams(n,k,obs_param)

       ! Get number of observations
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed in letkf_increments')

       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet failed in letkf_increments')

       Nobs = Nobjs*N_obs_size
       N_ens = LIS_rc%nensem(n)

       ! Allocate observation arrays
       allocate(Observations(Nobs))
       allocate(Obs_pred(Nobs, N_ens))
       allocate(Obs_pert(Nobs, N_ens))

       ! Get observations
       call generateObservations(n,k,Nobjs,Nobs,LIS_OBS_State(n,k), &
            LIS_OBS_Pert_State(n,k),Observations)

       ! Get observation predictions from model
       call LIS_surfaceModel_DAGetObsPred(n,k,Obs_pred)

       ! Get observation perturbations
       call getObsPert(LIS_OBS_Pert_State(n,k),N_obs_size,Nobs,N_ens,Obs_pert)

       ! Count observations to assimilate
       N_obs_actual = 0
       do i=1,Nobs
          if(Observations(i)%assim) N_obs_actual = N_obs_actual + 1
       enddo

       write(LIS_logunit,*) '[INFO] Total observations:', Nobs
       write(LIS_logunit,*) '[INFO] Observations with assim=true:', N_obs_actual

       if(N_obs_actual == 0) then
          write(LIS_logunit,*) '[INFO] No observations to assimilate, skipping LETKF'
          deallocate(obs_param)
          deallocate(Observations)
          deallocate(Obs_pred)
          deallocate(Obs_pert)
          return
       endif

       ! Filter to keep only assim=true observations
       allocate(Observations_filtered(N_obs_actual))
       allocate(Obs_pred_filtered(N_obs_actual, N_ens))
       allocate(Obs_pert_filtered(N_obs_actual, N_ens))

       j = 0
       do i=1,Nobs
          if(Observations(i)%assim) then
             j = j + 1
             Observations_filtered(j) = Observations(i)
             Obs_pred_filtered(j,:) = Obs_pred(i,:)
             Obs_pert_filtered(j,:) = Obs_pert(i,:)
          endif
       enddo

       ! Compute innovations for diagnostics
       if(LIS_rc%winnov(k).eq.1) then
          do i=1,Nobs
             if(Observations(i)%assim) then
                innov = Observations(i)%value - sum(Obs_pred(i,:))/real(N_ens)
                call row_variance(1,N_ens,Obs_pred(i,:),std_innov(1))
                letkf_struc(n,k)%forecast_var(i) = std_innov(1)
                std_innov = std_innov+(Observations(i)%std)**2
                std_innov = sqrt(std_innov)
                letkf_struc(n,k)%norm_innov(i) = innov/std_innov(1)
                letkf_struc(n,k)%innov(i) = innov
             else
                letkf_struc(n,k)%norm_innov(i) = LIS_rc%udef
                letkf_struc(n,k)%innov(i) = LIS_rc%udef
                letkf_struc(n,k)%forecast_var(i) = LIS_rc%udef
             endif
          enddo
       endif

       ! Get state variables
       state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)
       N_state = LIS_rc%nstvars(k)

       allocate(stvar(N_state, state_size))
       allocate(state_incr(N_state, state_size))
       allocate(state_tmp(N_state, state_size))

       call LIS_surfaceModel_DAScaleStateVar(n,k)
       call LIS_surfaceModel_DAgetStateVar(n,k,stvar)

       state_incr = stvar
       state_tmp = stvar

       ! Get tile coordinates
       allocate(lons(state_size))
       allocate(lats(state_size))
       call LIS_surfaceModel_DAgetLonLat(n,k,lons,lats)

       allocate(state_lat(state_size))
       allocate(state_lon(state_size))
       state_lat = lats
       state_lon = lons

       ! Get observation coordinates
       allocate(obs_lon_arr(N_obs_actual))
       allocate(obs_lat_arr(N_obs_actual))
       do i=1,N_obs_actual
          obs_lon_arr(i) = Observations_filtered(i)%lon
          obs_lat_arr(i) = Observations_filtered(i)%lat
       enddo

       ! Setup localization parameters
       use_km_loc = letkf_struc(n,k)%use_km_localization

       if(use_km_loc) then
          sigma_km = letkf_struc(n,k)%localization_scale_km
          max_dist = 6.0 * sigma_km  ! Search distance in km

          ! Compute km to degree conversion at mean latitude
          mean_lat = sum(state_lat) / real(state_size)
          km_per_deg_lat = EARTH_RADIUS_KM * DEG_TO_RAD
          km_per_deg_lon = EARTH_RADIUS_KM * cos(mean_lat * DEG_TO_RAD) * DEG_TO_RAD

          xcompact_km = 3.0 * sigma_km
          ycompact_km = 3.0 * sigma_km
          xcompact = xcompact_km / km_per_deg_lon
          ycompact = ycompact_km / km_per_deg_lat

          write(LIS_logunit,*) '[INFO] LETKF km-based localization:'
          write(LIS_logunit,*) '[INFO]   Sigma (km):', sigma_km
          write(LIS_logunit,*) '[INFO]   Compact support (km):', xcompact_km
          write(LIS_logunit,*) '[INFO]   At mean_lat=', mean_lat
          write(LIS_logunit,*) '[INFO]   xcompact (deg):', xcompact
          write(LIS_logunit,*) '[INFO]   ycompact (deg):', ycompact
       else
          xcompact = letkf_struc(n,k)%localization_factor * LIS_rc%gridDesc(n,9)
          ycompact = letkf_struc(n,k)%localization_factor * LIS_rc%gridDesc(n,10)
          max_dist = 2.0 * max(xcompact, ycompact) * 111.0  ! Approximate km

          write(LIS_logunit,*) '[INFO] LETKF factor-based localization:'
          write(LIS_logunit,*) '[INFO]   Factor:', letkf_struc(n,k)%localization_factor
          write(LIS_logunit,*) '[INFO]   xcompact (deg):', xcompact
          write(LIS_logunit,*) '[INFO]   ycompact (deg):', ycompact
       endif

       ! Set LETKF module parameters
       forget = letkf_struc(n,k)%forget
       type_trans = letkf_struc(n,k)%type_trans

       ! Initialize diagnostic counters
       tiles_with_obs = 0
       tiles_no_obs = 0
       total_local_obs_count = 0
       min_local_obs = N_obs_actual
       max_local_obs = 0

       n_tiles = state_size / N_ens

       !------------------------------------------------------------------------
       ! Loop over tiles for local LETKF analysis
       !------------------------------------------------------------------------
       write(LIS_logunit,*) '[INFO] Starting local LETKF loop over', n_tiles, 'tiles'

       do i = 1, n_tiles
          tile_lon = state_lon((i-1)*N_ens+1)
          tile_lat = state_lat((i-1)*N_ens+1)

          ! Find local observations
          allocate(local_obs_idx(N_obs_actual))
          local_count = 0

          do j = 1, N_obs_actual
             if(use_km_loc) then
                dist = haversine_km(tile_lon, tile_lat, &
                     obs_lon_arr(j), obs_lat_arr(j))
                if(dist <= max_dist) then
                   local_count = local_count + 1
                   local_obs_idx(local_count) = j
                endif
             else
                dx = abs(tile_lon - obs_lon_arr(j))
                dy = abs(tile_lat - obs_lat_arr(j))
                if(dx <= 2.0*xcompact .and. dy <= 2.0*ycompact) then
                   local_count = local_count + 1
                   local_obs_idx(local_count) = j
                endif
             endif
          enddo

          N_local_obs = local_count

          ! Store localization diagnostics
          letkf_struc(n,k)%n_local_obs(i) = real(N_local_obs)

          if(N_local_obs > 0) then
             tiles_with_obs = tiles_with_obs + 1
             total_local_obs_count = total_local_obs_count + N_local_obs
             if(N_local_obs < min_local_obs) min_local_obs = N_local_obs
             if(N_local_obs > max_local_obs) max_local_obs = N_local_obs

             ! Extract local observations
             allocate(obs_da(N_local_obs))
             allocate(obspred_da(N_local_obs, N_ens))
             allocate(obspert_da(N_local_obs, N_ens))
             allocate(obs_cov(N_local_obs, N_local_obs))

             do j = 1, N_local_obs
                obs_da(j) = Observations_filtered(local_obs_idx(j))
                obspred_da(j,:) = Obs_pred_filtered(local_obs_idx(j),:)
                obspert_da(j,:) = Obs_pert_filtered(local_obs_idx(j),:)
             enddo

             ! Assemble local observation error covariance
             call assemble_obs_cov(LIS_rc%nobtypes(k), N_local_obs, &
                  obs_param, obs_da, obs_cov)

             ! Copy forecast state to temporary
             state_tmp(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
                  stvar(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens)

             ! Call LETKF analysis
             call letkf_analysis(i, &
                  N_state, N_local_obs, N_ens, &
                  obs_da, obspred_da, obs_cov, &
                  state_tmp(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens), &
                  state_lon((i-1)*N_ens+1:(i-1)*N_ens+N_ens), &
                  state_lat((i-1)*N_ens+1:(i-1)*N_ens+N_ens), &
                  xcompact, ycompact )

             ! Compute increment
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
                  state_tmp(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) - &
                  stvar(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens)

             deallocate(obs_da)
             deallocate(obspred_da)
             deallocate(obspert_da)
             deallocate(obs_cov)
          else
             tiles_no_obs = tiles_no_obs + 1
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0

             letkf_struc(n,k)%mean_loc_weight(i) = 0.0
             letkf_struc(n,k)%max_loc_weight(i) = 0.0
             letkf_struc(n,k)%min_loc_weight(i) = 0.0
          endif

          letkf_struc(n,k)%anlys_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
               state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens)

          deallocate(local_obs_idx)
       enddo

       !------------------------------------------------------------------------
       ! Summary Statistics
       !------------------------------------------------------------------------
       if(tiles_with_obs > 0) then
          avg_local_obs = real(total_local_obs_count) / real(tiles_with_obs)
       else
          avg_local_obs = 0.0
          min_local_obs = 0
       endif

       write(LIS_logunit,*) '[INFO] =========================================='
       write(LIS_logunit,*) '[INFO] LETKF Summary:'
       write(LIS_logunit,*) '[INFO]   Total tiles processed:', n_tiles
       write(LIS_logunit,*) '[INFO]   Tiles with local observations:', tiles_with_obs
       write(LIS_logunit,*) '[INFO]   Tiles without observations:', tiles_no_obs
       if(tiles_with_obs > 0) then
          write(LIS_logunit,*) '[INFO]   Min local obs per tile:', min_local_obs
          write(LIS_logunit,*) '[INFO]   Max local obs per tile:', max_local_obs
          write(LIS_logunit,*) '[INFO]   Avg local obs per tile:', avg_local_obs
       endif
       write(LIS_logunit,*) '[INFO] =========================================='

       call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.true.)

       ! Update state
       call LIS_surfaceModel_DAsetAnlysisUpdates(n,k,N_state,state_size,&
            stvar,state_incr)

       call LIS_surfaceModel_DADescaleStatevar(n,k)

       ! Cleanup
       deallocate(obs_param)
       deallocate(stvar)
       deallocate(State_incr)
       deallocate(state_tmp)
       deallocate(Observations)
       deallocate(Obs_pred)
       deallocate(Obs_pert)
       deallocate(Observations_filtered)
       deallocate(Obs_pred_filtered)
       deallocate(Obs_pert_filtered)
       deallocate(state_lat)
       deallocate(state_lon)
       deallocate(lats)
       deallocate(lons)
       if(allocated(obs_lon_arr)) deallocate(obs_lon_arr)
       if(allocated(obs_lat_arr)) deallocate(obs_lat_arr)
    end if

  end subroutine letkf_increments

!BOP
!
! !ROUTINE: letkf_update
! \label{letkf_update}
!
! !INTERFACE:
  subroutine letkf_update(n,k)
! !USES:

    implicit none

    integer,       intent(in)    :: n
    integer,       intent(in)    :: k
!
! !DESCRIPTION:
!  This routine updates the model prognostics using the analysis
!  increments computed earlier.
!
!EOP
    integer                       :: status
    logical                       :: fresh_incr
    type(obs_type), allocatable   :: Observations(:)
    integer                       :: i
    integer                       :: N_ens
    integer                       :: Nobjs, Nobs, N_obs_size
    real,     allocatable         :: Obs_pred(:,:)

    call LIS_surfaceModel_DAGetFreshIncrementsStatus(n,k,fresh_incr)

    if(fresh_incr) then
       if(LIS_rc%incroption(k).eq.0) then

          call LIS_surfaceModel_DAUpdateState(n,k)
          call LIS_surfaceModel_DAQCstate(n,k)
          call LIS_surfaceModel_DASetStateVar(n,k)

          ! Compute analysis residuals
          if(LIS_rc%winnov(k).eq.1) then
             N_ens = LIS_rc%nensem(n)
             call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
             call LIS_verify(status, 'ESMF_StateGet failed in letkf_update')

             call ESMF_AttributeGet(LIS_OBS_State(n,k),&
                  name="Number Of Observations",&
                  value=N_obs_size,rc=status)
             call LIS_verify(status, &
                  'ESMF_AttributeGet: Number of Observations failed in letkf_update')

             Nobs = Nobjs*N_obs_size
             allocate(Observations(Nobs))

             call generateObservations(n, k, Nobjs, Nobs, LIS_OBS_State(n,k), &
                  LIS_OBS_Pert_State(n,k), Observations)

             allocate(Obs_pred(Nobs,N_ens))
             call LIS_surfaceModel_DAGetObsPred(n,k,Obs_pred)

             do i=1,Nobs
                if(Observations(i)%assim) then
                   letkf_struc(n,k)%anlys_res(i) = Observations(i)%value - &
                        sum(Obs_pred(i,:))/real(LIS_rc%nensem(n))
                else
                   letkf_struc(n,k)%anlys_res(i) = LIS_rc%udef
                endif
             enddo

             deallocate(Observations)
             deallocate(Obs_pred)
          endif
       endif

       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
            .true., rc=status)
       call LIS_verify(status, 'ESMF_AttributeSet failed in letkf_update')

       write(LIS_logunit,*) '[INFO] Finished assimilating Observations using LETKF'

       call print_increment_statistics(n, k)
    else
       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status, 'ESMF_AttributeSet failed in letkf_update')
    endif

  end subroutine letkf_update


!BOP
!
! !ROUTINE: letkf_diagnostics
!  \label{letkf_diagnostics}
!
! !INTERFACE:
  subroutine letkf_diagnostics(n,k)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
    integer, intent(IN)    :: k
!
! !DESCRIPTION:
!  This subroutine generates the LETKF related diagnostics and outputs
!  it to a file.
!
!EOP

    integer                :: v
    integer                :: status
    logical                :: assim_status
    logical                :: alarmCheck
    character*3            :: fda

    write(fda,'(i3.3)') k
    alarmCheck = LIS_isAlarmRinging(LIS_rc,"LIS DA output "//trim(fda))

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Assimilate Status",&
         value=assim_status,rc=status)

    if(assim_status) then
       if(.not.letkf_struc(n,k)%fileopen.and.LIS_masterproc) then
          call LIS_create_output_directory("LETKF")
       endif

       call writeInnovationOutput(n,k)
       call writeAnalysisIncr(n,k)
    endif

    if(alarmCheck) then
       if(.not.letkf_struc(n,k)%fileopen.and.LIS_masterproc) then
          call LIS_create_output_directory("LETKF")
       endif
       call writeEnsembleSpread(n,k)
    endif

#endif
  end subroutine letkf_diagnostics

!BOP
!
! !ROUTINE: letkf_final
! \label{letkf_final}
!
! !INTERFACE:
  subroutine letkf_final
!
! !DESCRIPTION:
!  This method performs the finalization for all LETKF
!  related structures and subroutines.
!
!EOP
    deallocate(letkf_struc)
  end subroutine letkf_final


!BOP
!
! !ROUTINE: assemble_obs_cov
! \label{assemble_obs_cov_letkf}
!
! !INTERFACE:
  subroutine assemble_obs_cov(nob, N_obs, obs_param, &
       Observations, Obs_cov)

    implicit none
! !ARGUMENTS:
    integer, intent(in)                              :: nob, N_obs
    type(obs_param_type), dimension(nob), intent(in) :: obs_param
    type(obs_type), dimension(N_obs), intent(in)     :: Observations
    real, intent(out), dimension(N_obs,N_obs)        :: Obs_cov
!
! !DESCRIPTION:
!  Assemble measurements error covariance.
!EOP

    integer :: i, j, i_species, j_species
    real :: fac, xcorr_tmp, ycorr_tmp

    ! Initialize
    Obs_cov = 0.

    ! Diagonal elements
    do i=1,N_obs
       Obs_cov(i,i) = Observations(i)%std**2
    end do

    ! Off-diagonal elements
    do i=1,N_obs
       do j=(i+1),N_obs
          i_species = Observations(i)%species
          j_species = Observations(j)%species

          ! Non-zero correlation only between observations of same type
          if (i_species == j_species) then
             xcorr_tmp = obs_param(i_species)%xcorr
             ycorr_tmp = obs_param(i_species)%ycorr

             if (xcorr_tmp>0. .and. ycorr_tmp>0.) then
                ! Gaussian correlation
                fac =  &
                     ((Observations(i)%lon-Observations(j)%lon)**2 &
                     /xcorr_tmp**2 )                                         &
                     +                                                       &
                     ((Observations(i)%lat-Observations(j)%lat)**2 &
                     /ycorr_tmp**2 )

                fac = exp(-.5*fac)

                Obs_cov(i,j) = Observations(i)%std * Observations(j)%std * fac
                Obs_cov(j,i) = Obs_cov(i,j)
             end if
          end if
       end do
    end do

  end subroutine assemble_obs_cov


!BOP
!
! !ROUTINE: getObsPert
! \label{getObsPert_letkf}
!
! !INTERFACE:
  subroutine getObsPert(OBS_Pert_State, gsize, dim1, dim2, pert)
! !USES:

! !ARGUMENTS:
    type(ESMF_State)      :: OBS_Pert_State
    integer               :: gsize
    integer               :: dim1, dim2
    real                  :: pert(dim1,dim2)
!
! !DESCRIPTION:
!  This routine retrieves the perturbations to be applied to the
!  observations from the Observation perturbation state.
!EOP
    real, pointer         :: obs_temp(:,:)
    integer               :: i
    integer               :: obs_state_count
    integer               :: status
    character*100,allocatable     :: obs_state_objs(:)
    type(ESMF_Field), allocatable :: obs_field(:)

    call ESMF_StateGet(OBS_Pert_State,itemCount=obs_state_count,rc=status)
    call LIS_verify(status, 'ESMF_StateGet failed in getObsPert')

    allocate(obs_state_objs(obs_state_count))
    allocate(obs_field(obs_state_count))

    call ESMF_StateGet(OBS_Pert_State,itemNameList=obs_state_objs,rc=status)
    call LIS_verify(status,'ESMF_StateGet failed in getObsPert')

    do i=1,obs_state_count
       call ESMF_StateGet(OBS_Pert_State,obs_state_objs(i),&
            obs_field(i),rc=status)
       call LIS_verify(status,'ESMF_StateGet failed in getObsPert')
       call ESMF_FieldGet(obs_field(i),localDE=0,farrayPtr=obs_temp,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed in getObsPert')
       pert((i-1)*gsize + 1:i*gsize,:) = obs_temp(:,:)
    enddo
    deallocate(obs_state_objs)
    deallocate(obs_field)

  end subroutine getObsPert

!BOP
!
! !ROUTINE: generateObservations
! \label{generateObservations_letkf}
!
! !INTERFACE:
  subroutine generateObservations(n, k, Nobjs, N_obs_size, LIS_OBS_State, &
       LIS_OBS_Pert_State, Observations )
! !USES:

! !ARGUMENTS:
    integer,     intent(IN)        :: n
    integer,     intent(IN)        :: k
    integer,     intent(IN)        :: Nobjs
    integer,     intent(IN)        :: N_obs_size
    type(ESMF_State)               :: LIS_OBS_State
    type(ESMF_State)               :: LIS_OBS_Pert_State
    type(obs_type)                 :: Observations(N_obs_size)
!
! !DESCRIPTION:
!   This subroutine unpacks the Observation state and packages them into the
!   'obs_type' datastructure.
!EOP

    integer                        :: typ
    integer                        :: gval,gsize
    integer                        :: count1
    character*100                  :: temp
    integer                        :: i, t,g
    character*1                    :: vid(2)
    integer                        :: status
    type(ESMF_Field)               :: valfield
    type(ESMF_Field)               :: pertField
    real, pointer                  :: value1(:)
    integer, allocatable           :: gid1(:)
    real,    allocatable           :: obsstd(:)
    integer                        :: pert_type
    integer, allocatable           :: obsassimflag(:)
    integer                        :: counts, row, col
    integer                        :: pertyp(N_obs_size)
    real                           :: value(N_obs_size)
    integer                        :: gid(N_obs_size)
    real                           :: std(N_obs_size)
    integer                        :: species(N_obs_size)
    integer                        :: assimflag(N_obs_size)

    count1 = 1
    typ    = 1
    do i=1,Nobjs

       write(unit=temp,fmt='(i2.2)') i
       read(unit=temp,fmt='(2a1)') vid

       call ESMF_StateGet(LIS_OBS_State,"Observation"//vid(1)//vid(2),&
            valfield,rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed in generateObservations')

       call ESMF_StateGet(LIS_OBS_Pert_State,&
            "Observation"//vid(1)//vid(2),&
            pertField,rc=status)
       call LIS_verify(status,'ESMF_StateGet failed in generateObservations')

       call ESMF_AttributeGet(pertField, "Perturbation Type",&
            pert_type, rc=status)
       call LIS_verify(status, "Perturbation Type attribute not found in letkf_Mod")

       call ESMF_AttributeGet(valfield,"Grid Number of Observations",&
            gval, rc=status)
       call LIS_verify(status, 'ESMF_FieldGet: Obs gval failed in generateObservations')

       call ESMF_FieldGet(valfield,localDE=0,farrayPtr=value1,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet in letkf')

       allocate(gid1(gval))
       allocate(obsstd(gval))
       allocate(obsassimflag(gval))

       call ESMF_AttributeGet(valfield,"Grid Observation Id",&
            itemCount = counts,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Obs gid1 count failed')

       call ESMF_AttributeGet(valfield,"Grid Observation Id",&
            gid1, itemCount = counts, rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Obs gid failed')

       call ESMF_AttributeGet(valfield,"Standard Deviation of Observations",&
            itemCount = counts,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Obs std count failed')

       call ESMF_AttributeGet(valfield,"Standard Deviation of Observations",&
            obsstd, itemCount = counts, rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Obs std failed')

       call ESMF_AttributeGet(valfield,"Assimilation Flag",&
            itemCount = counts,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Assimflag count failed')

       call ESMF_AttributeGet(valfield,"Assimilation Flag",&
            obsassimflag, itemCount = counts, rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet: Assimflag failed')

       do t=1, gval
          pertyp(count1) = pert_type
          species(count1) = typ
          value(count1) = value1(t)
          gid(count1) = gid1(t)
          std(count1) = obsstd(t)
          assimflag(count1) = obsassimflag(t)
          count1 = count1 + 1
       enddo
       typ = typ + 1
       deallocate(gid1)
       deallocate(obsstd)
       deallocate(obsassimflag)
    enddo

    do i=1,N_obs_size
       Observations(i)%value = value(i)
       Observations(i)%std   = std(i)
       Observations(i)%species = species(i)
       Observations(i)%pert_type = pertyp(i)
       if(assimflag(i).ne.0) then
          Observations(i)%assim = .true.
       else
          Observations(i)%assim = .false.
       endif

       ! Get lon/lat from gid
       row = (gid(i)-1)/LIS_rc%obs_gnc(k)+1
       col = gid(i) - (row-1)*LIS_rc%obs_gnc(k)

       Observations(i)%lon = LIS_rc%obs_gridDesc(k,4) + (col-1)*LIS_rc%obs_gridDesc(k,9)
       Observations(i)%lat = LIS_rc%obs_gridDesc(k,5) + (row-1)*LIS_rc%obs_gridDesc(k,10)
    enddo

  end subroutine generateObservations


!BOP
!
! !ROUTINE: writeInnovationOutput
! \label{writeInnovationOutput_letkf}
!
! !INTERFACE:
  subroutine writeInnovationOutput(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n
    integer,  intent(in)    :: k
!
! !DESCRIPTION:
!  This routine writes the innovation values to an external file.
!EOP
    integer                :: ftn
    character(len=LIS_CONST_PATH_LEN) :: innovfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3), ares_Id, ninnov_Id, innov_id
    integer                :: forecast_sigma_id, aincr_Id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status

    if(LIS_rc%winnov(k).eq.1) then

       shuffle = 1
       deflate = 1
       deflate_level = 9

       if(LIS_masterproc) then
          call LIS_create_innov_filename(n,k, innovfile, 'LETKF')

#if (defined USE_NETCDF4)
          status = nf90_create(path=innovfile,cmode=nf90_hdf5, ncid = ftn)
          call LIS_verify(status,'creating netcdf file failed in letkf_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=innovfile,cmode=nf90_clobber, ncid = ftn)
          call LIS_verify(status,'creating netcdf file failed in letkf_Mod')
#endif

          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%obs_gnc(k),&
               dimID(1)),'nf90_def_dim for east_west failed in letkf_mod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%obs_gnr(k),&
               dimID(2)),'nf90_def_dim for north_south failed in letkf_mod')

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att failed for missing_value in letkf_mod')

          ! Normalized innovations
          write(unit=finst, fmt='(i2.2)') k
          varname = "ninnov_"//trim(finst)
          vardimname = "ninnov_"//trim(finst)//"_levels"
          standard_name = "Normalized_innovations_for_DA_instance_"//trim(finst)

          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for ninnov')
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float, dimids = dimID(1:2), varID=ninnov_Id),&
               'nf90_def_var failed for ninnov')

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,ninnov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for ninnov')
#endif
          call LIS_verify(nf90_put_att(ftn,ninnov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for ninnov')

          ! Innovations
          varname = "innov_"//trim(finst)
          vardimname = "innov_"//trim(finst)//"_levels"
          standard_name = "Innovations_for_DA_instance_"//trim(finst)

          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for innov')
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float, dimids = dimID(1:2), varID=innov_Id),&
               'nf90_def_var failed for innov')

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,innov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for innov failed')
#endif
          call LIS_verify(nf90_put_att(ftn,innov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for innov')

          ! Analysis residuals
          varname = "analysis_residual_"//trim(finst)
          vardimname = "analysis_residual_"//trim(finst)//"_levels"
          standard_name = "Analysis_residuals_for_DA_instance_"//trim(finst)

          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for analysis_residual')
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float, dimids = dimID(1:2), varID=ares_Id),&
               'nf90_def_var failed for analysis_residual')

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,ares_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for analysis_residual')
#endif
          call LIS_verify(nf90_put_att(ftn,ares_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for analysis_residual')

          ! Forecast variance
          varname = "forecast_variance_"//trim(finst)
          vardimname = "forecast_variance_"//trim(finst)//"_levels"
          standard_name = "Forecast_variance_for_DA_instance_"//trim(finst)

          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for forecast_variance')
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float, dimids = dimID(1:2), varID=forecast_sigma_id),&
               'nf90_def_var failed for forecast_variance')

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,forecast_sigma_id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for forecast_variance')
#endif
          call LIS_verify(nf90_put_att(ftn,forecast_sigma_id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for forecast_variance')

          call LIS_verify(nf90_enddef(ftn),&
               'nf90_enddef failed in letkf_mod')

          ! Write data
          call LIS_verify(nf90_put_var(ftn,ninnov_Id,&
               letkf_struc(n,k)%norm_innov),&
               'nf90_put_var failed for norm_innov')
          call LIS_verify(nf90_put_var(ftn,innov_Id,&
               letkf_struc(n,k)%innov),&
               'nf90_put_var failed for innov')
          call LIS_verify(nf90_put_var(ftn,ares_Id,&
               letkf_struc(n,k)%anlys_res),&
               'nf90_put_var failed for anlys_res')
          call LIS_verify(nf90_put_var(ftn,forecast_sigma_id,&
               letkf_struc(n,k)%forecast_var),&
               'nf90_put_var failed for forecast_var')

          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in letkf_mod')
       endif
    endif

  end subroutine writeInnovationOutput


!BOP
!
! !ROUTINE: writeAnalysisIncr
! \label{writeAnalysisIncr_letkf}
!
! !INTERFACE:
  subroutine writeAnalysisIncr(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n
    integer,  intent(in)    :: k
!
! !DESCRIPTION:
!  This routine writes the analysis increments to an external file.
!EOP
    integer                :: ftn
    character(len=LIS_CONST_PATH_LEN) :: incrfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3), aincr_Id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status

    shuffle = 1
    deflate = 1
    deflate_level = 9

    if(LIS_masterproc) then
       call LIS_create_incr_filename(n,k, incrfile, 'LETKF')

#if (defined USE_NETCDF4)
       status = nf90_create(path=incrfile,cmode=nf90_hdf5, ncid = ftn)
       call LIS_verify(status,'creating netcdf file failed in letkf_Mod')
#endif
#if (defined USE_NETCDF3)
       status = nf90_create(path=incrfile,cmode=nf90_clobber, ncid = ftn)
       call LIS_verify(status,'creating netcdf file failed in letkf_Mod')
#endif

       call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
            dimID(1)),'nf90_def_dim for east_west failed in letkf_mod')
       call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
            dimID(2)),'nf90_def_dim for north_south failed in letkf_mod')
       call LIS_verify(nf90_def_dim(ftn,'ensemble',LIS_rc%nensem(n),&
            dimID(3)),'nf90_def_dim for ensemble failed in letkf_mod')

       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att failed for missing_value in letkf_mod')

       ! Analysis increment
       write(unit=finst, fmt='(i2.2)') k
       varname = "analysis_increment_"//trim(finst)
       standard_name = "Analysis_increment_for_DA_instance_"//trim(finst)

       call LIS_verify(nf90_def_var(ftn,varname,&
            nf90_float, dimids = dimID(1:3), varID=aincr_Id),&
            'nf90_def_var failed for analysis_increment')

#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,aincr_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for analysis_increment')
#endif
       call LIS_verify(nf90_put_att(ftn,aincr_Id,&
            "standard_name",standard_name),&
            'nf90_put_att failed for analysis_increment')

       call LIS_verify(nf90_enddef(ftn),'nf90_enddef failed in letkf_mod')

       ! Write data
       call LIS_verify(nf90_put_var(ftn,aincr_Id,&
            letkf_struc(n,k)%anlys_incr),&
            'nf90_put_var failed for anlys_incr')

       call LIS_verify(nf90_close(ftn),'nf90_close failed in letkf_mod')
    endif

  end subroutine writeAnalysisIncr


!BOP
!
! !ROUTINE: writeEnsembleSpread
! \label{writeEnsembleSpread_letkf}
!
! !INTERFACE:
  subroutine writeEnsembleSpread(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n
    integer,  intent(in)    :: k
!
! !DESCRIPTION:
!  This routine writes the ensemble spread diagnostics.
!EOP
    ! Simplified placeholder - full implementation can be added later
    write(LIS_logunit,*) '[INFO] LETKF ensemble spread diagnostics: not yet implemented'

  end subroutine writeEnsembleSpread


!BOP
!
! !ROUTINE: print_increment_statistics
! \label{print_increment_statistics_letkf}
!
! !INTERFACE:
  subroutine print_increment_statistics(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n
    integer,  intent(in)    :: k
!
! !DESCRIPTION:
!  This routine prints statistics about the analysis increments.
!EOP
    integer :: i, j, state_size, N_state, N_ens
    integer :: n_nonzero, n_total
    real :: incr_sum, incr_sum_sq, incr_mean, incr_std
    real :: incr_min, incr_max
    real :: tmp_val

    N_state = LIS_rc%nstvars(k)
    N_ens = LIS_rc%nensem(n)
    state_size = size(letkf_struc(n,k)%anlys_incr, 2)

    write(LIS_logunit,*) '[INFO] =========================================='
    write(LIS_logunit,*) '[INFO] LETKF Analysis Increment Statistics:'

    do j = 1, N_state
       n_nonzero = 0
       n_total = 0
       incr_sum = 0.0
       incr_sum_sq = 0.0
       incr_min = 1.0e30
       incr_max = -1.0e30

       do i = 1, state_size
          tmp_val = letkf_struc(n,k)%anlys_incr(j, i)
          if(abs(tmp_val) > 1.0e-15 .and. tmp_val /= LIS_rc%udef) then
             n_nonzero = n_nonzero + 1
             incr_sum = incr_sum + tmp_val
             incr_sum_sq = incr_sum_sq + tmp_val**2
             if(tmp_val < incr_min) incr_min = tmp_val
             if(tmp_val > incr_max) incr_max = tmp_val
          endif
          n_total = n_total + 1
       enddo

       if(n_nonzero > 0) then
          incr_mean = incr_sum / real(n_nonzero)
          incr_std = sqrt(incr_sum_sq/real(n_nonzero) - incr_mean**2)
       else
          incr_mean = 0.0
          incr_std = 0.0
          incr_min = 0.0
          incr_max = 0.0
       endif

       write(LIS_logunit,*) '[INFO]   State var', j, ':'
       write(LIS_logunit,*) '[INFO]     Non-zero increments:', n_nonzero, '/', n_total
       write(LIS_logunit,*) '[INFO]     Mean increment:', incr_mean
       write(LIS_logunit,*) '[INFO]     Std deviation:', incr_std
       write(LIS_logunit,*) '[INFO]     Min:', incr_min, ' Max:', incr_max
    enddo

    write(LIS_logunit,*) '[INFO] =========================================='

  end subroutine print_increment_statistics

end module letkf_Mod
