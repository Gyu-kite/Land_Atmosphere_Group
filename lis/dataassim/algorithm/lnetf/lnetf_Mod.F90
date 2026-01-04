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
module lnetf_Mod
!BOP
!
! !MODULE: lnetf_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines that control
!   the incorporation of a data set using the Local Nonlinear
!   Ensemble Transform Filter (LNETF) method, into a land surface model.
!
!  The LNETF algorithm is based on the work of:
!  J. Toedter and B. Ahrens, "A Second-Order Exact Ensemble Square
!  Root Filter for Nonlinear Data Assimilation", Mon. Wea. Rev. 143 (2015) 1347
!  Implementation adapted from PDAF (Parallel Data Assimilation Framework)
!  by Lars Nerger and Paul Kirchgessner.
!
!  NOTES: Data assimilation is currently only supported for land surface
!  models (and not across different surface model types)
!
! !REVISION HISTORY:
!   27 Feb 2005: Sujay Kumar; Initial Specification (LNETF)
!   28 Dec 2024: Adapted for LNETF from PDAF implementation
!
! !USES:
  use ESMF
  use lnetf_types
  use lnetf_general
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
  public :: lnetf_init  ! Initialization for LNETF
  public :: lnetf_setup
  public :: lnetf_increments ! compute analysis increments
  public :: lnetf_update ! apply the analysis increments
  public :: lnetf_diagnostics ! write LNETF related diagnostics
  public :: lnetf_final ! Finalization for LNETF
  ! Multi-source DA functions (simultaneous assimilation)
  public :: lnetf_multisource_increments  ! compute increments from combined obs
  public :: lnetf_multisource_update      ! apply combined increments
  public :: lnetf_multisource_diagnostics ! combined diagnostics output
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: lnetf_struc ! data structure containing LNETF diagnostics
!EOP

  type, public ::  lnetf_dec
     logical     :: fileOpen
     real, allocatable :: innov(:)
     real, allocatable :: forecast_var(:) !HPHt
     real, allocatable :: anlys_res(:)
     real, allocatable :: anlys_incr(:,:)
     real, allocatable :: norm_innov(:)
     real, allocatable :: weights(:)      ! particle weights
     real, allocatable :: transform(:,:)  ! ensemble transformation matrix
     real :: localization_factor = 5.0    ! for localization (grid-relative)

     ! Localization parameters (km-based, following Seo et al. 2021 / Tak et al. 2025)
     ! Reference: Seo, E., et al. (2021) - JULES + LETKF with σ=30km, patch=150km
     !            Tak, Y.-J., et al. (2025) - Multi-sensor SMAP+ASCAT LETKF
     real :: localization_scale_km = 30.0   ! Gaspari-Cohn scale σ (km) - default 30km like papers
     logical :: use_km_localization = .false.  ! If true, use km-based; if false, use factor-based

     ! LNETF specific parameters
     integer :: type_forget = 0           ! Type of forgetting factor
     integer :: type_trans = 0            ! Type of ensemble transformation
     integer :: type_winf = 0             ! Type of weights inflation
     real :: forget = 1.0                 ! Forgetting factor
     real :: limit_winf = 0.0             ! Limit for weights inflation
     real :: eff_sample_size = 0.0        ! Effective sample size

     !--------------------------------------------------------------------------
     ! Localization diagnostic arrays (for output to increment file)
     !--------------------------------------------------------------------------
     ! n_local_obs: Number of observations within localization radius for each tile
     !   - Higher values indicate more observation influence
     !   - Zero values mean no observations within localization radius
     real, allocatable :: n_local_obs(:)

     ! mean_loc_weight: Mean Gaspari-Cohn localization weight for each tile
     !   - Ranges from 0.0 to 1.0
     !   - 1.0 = observations very close to tile (strongest influence)
     !   - Near 0.0 = observations at edge of localization radius (weak influence)
     !   - Undefined for tiles with no local observations
     real, allocatable :: mean_loc_weight(:)

     ! max_loc_weight: Maximum Gaspari-Cohn localization weight for each tile
     !   - Shows the closest observation's influence
     !   - Useful for identifying swath boundaries vs centers
     real, allocatable :: max_loc_weight(:)

     ! min_loc_weight: Minimum Gaspari-Cohn localization weight for each tile
     !   - Shows the farthest observation's influence within localization radius
     !   - Large gap between max and min indicates observations at varying distances
     real, allocatable :: min_loc_weight(:)

     ! eff_sample_size_tile: Effective sample size (N_eff) for each tile
     !   - Diagnostic for filter degeneracy: N_eff = 1 / sum(w_i^2)
     !   - Higher values (close to N_ens) = healthy filter with well-spread weights
     !   - Low values (close to 1) = weight collapse, filter degeneracy risk
     !   - Used to trigger weight inflation in LNETF algorithm
     real, allocatable :: eff_sample_size_tile(:)

  end type lnetf_dec
!EOP

  type(lnetf_dec), allocatable :: lnetf_struc(:,:)

contains

!BOP
!
! !ROUTINE: lnetf_init
! \label{lnetf_init}
!
! !INTERFACE:
  subroutine lnetf_init()
! !USES:


!
! !DESCRIPTION:
!  This method performs the required initializations for the
!  LNETF method. The method reads the runtime settings from
!  the LIS configuration file.
!
!EOP
    allocate(lnetf_struc(LIS_rc%nnest, LIS_rc%ndas))
   
  end subroutine lnetf_init
!BOP
! 
! !ROUTINE: lnetf_setup
! \label{lnetf_setup}
!  
! !INTERFACE: 
  subroutine lnetf_setup(k)
! !USES: 

!
! !DESCRIPTION: 
!  This method performs the required initializations for the 
!  GMAO LNETF method. The method reads the runtime settings from 
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
          write(LIS_logunit,*) '[ERR] to greater than 1 for LNETF',LIS_rc%nensem(n)
          call LIS_endrun
       endif
       lnetf_struc(n,k)%fileOpen = .false.
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed in lnetf_Mod')

       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet failed in lnetf_Mod')

       if(LIS_rc%winnov(k).eq.1) then
          allocate(lnetf_struc(n,k)%norm_innov(Nobjs*N_obs_size))
          allocate(lnetf_struc(n,k)%innov(Nobjs*N_obs_size))
          allocate(lnetf_struc(n,k)%anlys_res(Nobjs*N_obs_size))
          allocate(lnetf_struc(n,k)%forecast_var(Nobjs*N_obs_size))

          ! Initialize innovation arrays to undefined (will be set during DA if assim=true)
          lnetf_struc(n,k)%norm_innov = LIS_rc%udef
          lnetf_struc(n,k)%innov = LIS_rc%udef
          lnetf_struc(n,k)%anlys_res = LIS_rc%udef
          lnetf_struc(n,k)%forecast_var = LIS_rc%udef
!          allocate(lnetf_struc(n,k)%k_gain(&
!               LIS_surfaceModel_DAgetStateSpaceSize(n), &
!               LIS_rc%nstvars(k)))
       endif
       allocate(lnetf_struc(n,k)%anlys_incr(LIS_rc%nstvars(k),&
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)))

       !----------------------------------------------------------------------
       ! Allocate localization diagnostic arrays (one value per tile)
       ! Note: state_space_size includes all ensemble members, so we divide
       ! by nensem to get the number of unique tiles
       !----------------------------------------------------------------------
       allocate(lnetf_struc(n,k)%n_local_obs( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(lnetf_struc(n,k)%mean_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(lnetf_struc(n,k)%max_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(lnetf_struc(n,k)%min_loc_weight( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))
       allocate(lnetf_struc(n,k)%eff_sample_size_tile( &
            LIS_surfaceModel_DAgetStateSpaceSize(n,k)/LIS_rc%nensem(n)))

       ! Initialize to undefined
       lnetf_struc(n,k)%n_local_obs = LIS_rc%udef
       lnetf_struc(n,k)%mean_loc_weight = LIS_rc%udef
       lnetf_struc(n,k)%max_loc_weight = LIS_rc%udef
       lnetf_struc(n,k)%min_loc_weight = LIS_rc%udef
       lnetf_struc(n,k)%eff_sample_size_tile = LIS_rc%udef
    enddo
!----------------------------------------------------------------------------
! Read localization parameters from config file
! Following Seo et al. (2021) and Tak et al. (2025) LETKF setup:
!   - σ = 30 km (Gaspari-Cohn scale)
!   - Local patch = 150 km (compact support radius = 2σ at patch edge)
!   - Result: "almost 1D filter" with weights 10^-2 to 10^-3 at edge
!----------------------------------------------------------------------------
     do n=1,LIS_rc%nnest
        ! First try km-based localization (preferred, like in papers)
        call ESMF_ConfigGetAttribute(LIS_config, &
             lnetf_struc(n,k)%localization_scale_km, &
             label="LNETF localization scale (km):", rc=status)
        if(status.eq.0) then
           lnetf_struc(n,k)%use_km_localization = .true.
           write(LIS_logunit,*) '[INFO] =============================================='
           write(LIS_logunit,*) '[INFO] LNETF using km-based localization (like LETKF papers)'
           write(LIS_logunit,*) '[INFO]   Gaspari-Cohn scale (σ): ', &
                lnetf_struc(n,k)%localization_scale_km, ' km'
           write(LIS_logunit,*) '[INFO]   Weight ~0.5 at distance: ', &
                lnetf_struc(n,k)%localization_scale_km, ' km'
           write(LIS_logunit,*) '[INFO]   Compact support cutoff (3σ): ', &
                3.0 * lnetf_struc(n,k)%localization_scale_km, ' km'
           write(LIS_logunit,*) '[INFO]   Max search distance (6σ): ', &
                6.0 * lnetf_struc(n,k)%localization_scale_km, ' km'
           write(LIS_logunit,*) '[INFO] Reference: Seo et al. (2021), Tak et al. (2025)'
           write(LIS_logunit,*) '[INFO] =============================================='
        else
           ! Fall back to factor-based localization
           call ESMF_ConfigGetAttribute(LIS_config, &
                lnetf_struc(n,k)%localization_factor, &
                label="LNETF localization radius factor:", rc=status)
           if(status.ne.0) then
              lnetf_struc(n,k)%localization_factor = 5.0  ! 기본값
              write(LIS_logunit,*) '[INFO] LNETF localization radius factor not found in config'
              write(LIS_logunit,*) '[INFO] Using default factor: 5.0'
              write(LIS_logunit,*) '[INFO] Tip: For LETKF-style localization, add to config:'
              write(LIS_logunit,*) '[INFO]   LNETF localization scale (km): 30.0'
           else
              write(LIS_logunit,*) '[INFO] LNETF localization radius factor set to:', &
                   lnetf_struc(n,k)%localization_factor
           endif
           lnetf_struc(n,k)%use_km_localization = .false.
        endif
     enddo
  
  end subroutine lnetf_setup

!BOP
!
! !ROUTINE: lnetf_increments
! \label{lnetf_increments}
!
! !INTERFACE:
  subroutine lnetf_increments(n,k)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
    integer, intent(IN)    :: k
!
! !DESCRIPTION:
!  This routine computes the analysis increments for LNETF using
!  PDAF-style 3D spatial localization. For each tile, observations
!  within the localization radius are selected and used for the
!  local analysis.
!
!  PDAF-style Local LNETF Algorithm:
!  1. Filter observations to keep only assim=.true.
!  2. For each tile, find observations within localization radius
!  3. Apply LNETF analysis with Gaspari-Cohn weighted likelihood
!  4. Tiles without nearby observations get zero increment
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

    ! PDAF-style Local LNETF variables
    real,         allocatable         :: obs_lon_arr(:), obs_lat_arr(:)
    real                              :: max_dist, dist, tile_lon, tile_lat
    integer,      allocatable         :: local_obs_idx(:)
    integer                           :: local_count

    ! km-based localization variables (Seo et al. 2021 / Tak et al. 2025 style)
    real                              :: mean_lat, km_per_deg_lon, km_per_deg_lat
    real                              :: sigma_km, xcompact_km, ycompact_km
    real, parameter                   :: PI = 3.14159265358979323846
    real, parameter                   :: DEG_TO_RAD = PI / 180.0
    real, parameter                   :: EARTH_RADIUS_KM = 6371.0
    logical                           :: use_km_loc

    ! Diagnostic counters for Local LNETF
    integer                           :: tiles_with_obs, tiles_no_obs
    integer                           :: total_local_obs_count
    integer                           :: min_local_obs, max_local_obs
    real                              :: avg_local_obs

    ! Output diagnostics from lnetf_analysis
    real                              :: tile_n_eff
    real                              :: tile_mean_loc_weight
    real                              :: tile_max_loc_weight
    real                              :: tile_min_loc_weight
    integer                           :: n_tiles


!----------------------------------------------------------------------------
!  Check if the observation state is updated or not. If it is updated,
!  the data is then assimilated.
!----------------------------------------------------------------------------

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Update Status",&
         value=data_status,rc=status)
    call LIS_verify(status, &
         'ESMF_AttributeGet: Data Update Status failed in lnetf_increments')

    call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.false.)

    lnetf_struc(n,k)%anlys_incr = 0.0

    ! Initialize localization diagnostic arrays to undefined
    lnetf_struc(n,k)%n_local_obs = LIS_rc%udef
    lnetf_struc(n,k)%mean_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%max_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%min_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%eff_sample_size_tile = LIS_rc%udef

    if(data_status) then
       write(LIS_logunit,*) &
            '[INFO] Assimilating Observations using PDAF-style Local LNETF for DA instance',k

       call LIS_getDomainResolutions(n,dx,dy)
       state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)

       ! Determine localization method
       use_km_loc = lnetf_struc(n,k)%use_km_localization

       if (use_km_loc) then
          !-------------------------------------------------------------------
          ! km-based localization (Seo et al. 2021 / Tak et al. 2025 style)
          !-------------------------------------------------------------------
          sigma_km = lnetf_struc(n,k)%localization_scale_km

          ! Following LETKF papers:
          ! σ = 30 km → at distance σ, weight ≈ 0.5-0.6
          ! Compact support cutoff = 2 * xcompact in Gaspari-Cohn notation
          ! To match paper's behavior (weight ≈ 0.01 at 75 km with σ=30 km):
          !   xcompact_km = σ (so weight=0 at 2σ=60km? No, that's too short)
          ! Actually, the code uses: d = distance/xcompact, then gc(d)=0 when d>=1
          ! So xcompact IS the cutoff distance where weight becomes 0
          !
          ! For σ=30km paper setup:
          !   - Weight at σ: ~0.5 (half-maximum)
          !   - Compact support: ~90-120 km (weight=0)
          ! We set: xcompact_km = 3*σ (compact support cutoff at 3σ ≈ 90 km)
          ! This gives: at σ, d=1/3, gc(1/3)≈0.7; at 2σ, d=2/3, gc(2/3)≈0.2

          xcompact_km = 3.0 * sigma_km  ! Compact support cutoff at 3σ
          ycompact_km = 3.0 * sigma_km

          ! Convert km to degrees (approximate)
          ! This is used for the main filtering loop
          km_per_deg_lat = 111.0  ! km per degree latitude (constant)

          ! For xcompact (longitude), we need cos(lat) correction
          ! Use domain center latitude for approximation
          ! Will be refined per-tile in the loop
          mean_lat = (LIS_rc%gridDesc(n,4) + LIS_rc%gridDesc(n,7)) / 2.0
          km_per_deg_lon = 111.0 * cos(mean_lat * DEG_TO_RAD)

          xcompact = xcompact_km / km_per_deg_lon  ! degrees
          ycompact = ycompact_km / km_per_deg_lat  ! degrees

          write(LIS_logunit,*) '[INFO] =============================================='
          write(LIS_logunit,*) '[INFO] LNETF km-based localization (LETKF paper style):'
          write(LIS_logunit,*) '[INFO]   σ (Gaspari-Cohn scale): ', sigma_km, ' km'
          write(LIS_logunit,*) '[INFO]   Compact support cutoff (3σ): ', 3.0*sigma_km, ' km'
          write(LIS_logunit,*) '[INFO]   Max search distance (6σ): ', 6.0*sigma_km, ' km'
          write(LIS_logunit,*) '[INFO]   Approximate xcompact: ', xcompact, ' deg'
          write(LIS_logunit,*) '[INFO]   Approximate ycompact: ', ycompact, ' deg'
          write(LIS_logunit,*) '[INFO]   (at mean lat ', mean_lat, ' deg)'
          write(LIS_logunit,*) '[INFO] =============================================='

       else
          !-------------------------------------------------------------------
          ! Factor-based localization (original method)
          !-------------------------------------------------------------------
          xcompact = dx * lnetf_struc(n,k)%localization_factor
          ycompact = dy * lnetf_struc(n,k)%localization_factor

          write(LIS_logunit,*) '[INFO] LNETF factor-based localization:'
          write(LIS_logunit,*) '[INFO]   Localization radius (degrees): ', &
               sqrt(xcompact**2 + ycompact**2)
       endif

       N_state = LIS_rc%nstvars(k)
       N_ens = LIS_rc%nensem(n)

!----------------------------------------------------------------------------
! Get observation state dimensions
!----------------------------------------------------------------------------
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status, &
            'ESMF_StateGet failed in lnetf_increments')

       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status, &
            'ESMF_AttributeGet: Number of Observations failed in lnetf_increments')

       Nobs = Nobjs*N_obs_size
       allocate(Observations(Nobs))

       call generateObservations(n, k, Nobjs, Nobs, LIS_OBS_State(n,k), &
            LIS_OBS_Pert_State(n,k),Observations)

!----------------------------------------------------------------------------
!  Retrieve Obs_pred : model's estimate of the observations (full grid)
!----------------------------------------------------------------------------
       allocate(Obs_pred(Nobs,N_ens))
       call LIS_surfaceModel_DAGetObsPred(n,k,Obs_pred)

!----------------------------------------------------------------------------
!  Retrieve Obs_pert : observation perturbations (full grid)
!----------------------------------------------------------------------------
       allocate(Obs_pert(Nobs,N_ens))
       call getObsPert(LIS_OBS_Pert_State(n,k),N_obs_size,&
            Nobs, N_ens, Obs_pert)

!----------------------------------------------------------------------------
!  Compute and store innovations BEFORE filtering (to preserve spatial order)
!----------------------------------------------------------------------------
       if(LIS_rc%winnov(k).eq.1) then
          do i=1,Nobs
             if(Observations(i)%assim) then
                innov = Observations(i)%value - &
                     sum(Obs_pred(i,:))/real(LIS_rc%nensem(n))
                call row_variance(1,LIS_rc%nensem(n),Obs_pred(i,:),std_innov(1))
                lnetf_struc(n,k)%forecast_var(i) = std_innov(1)
                std_innov = std_innov+(Observations(i)%std)**2
                std_innov = sqrt(std_innov)
                lnetf_struc(n,k)%norm_innov(i) = innov/std_innov(1)
                lnetf_struc(n,k)%innov(i) = innov
             else
                lnetf_struc(n,k)%norm_innov(i) = LIS_rc%udef
                lnetf_struc(n,k)%innov(i) = LIS_rc%udef
                lnetf_struc(n,k)%forecast_var(i) = LIS_rc%udef
             endif
          enddo
       endif

! ========================================================================
! PDAF-style Local LNETF: Filter observations to keep only assim=.true.
! This dramatically reduces computational cost
! ========================================================================
       N_obs_actual = 0
       do i = 1, Nobs
          if(Observations(i)%assim) then
             N_obs_actual = N_obs_actual + 1
          endif
       enddo

       write(LIS_logunit,*) '[INFO] Total obs in grid: ', Nobs
       write(LIS_logunit,*) '[INFO] Actual obs with assimflag=1: ', N_obs_actual

       if(N_obs_actual .eq. 0) then
          write(LIS_logunit,*) '[INFO] No observations with assimflag=1 in this domain'
          write(LIS_logunit,*) '[INFO] Will set all increments to 0'
          deallocate(Observations)
          deallocate(Obs_pred)
          deallocate(Obs_pert)
          N_obs_size = 0
          Nobs = 0
          allocate(Observations(0))
          allocate(Obs_pred(0, N_ens))
          allocate(Obs_pert(0, N_ens))
          allocate(obs_lon_arr(0))
          allocate(obs_lat_arr(0))
          max_dist = 2.0 * sqrt(xcompact**2 + ycompact**2)
       else
          ! Filter observations: keep only assim=.true.
          allocate(Observations_filtered(N_obs_actual))
          allocate(Obs_pred_filtered(N_obs_actual, N_ens))
          allocate(Obs_pert_filtered(N_obs_actual, N_ens))

          j = 1
          do i = 1, Nobs
             if(Observations(i)%assim) then
                Observations_filtered(j) = Observations(i)
                Obs_pred_filtered(j, :) = Obs_pred(i, :)
                Obs_pert_filtered(j, :) = Obs_pert(i, :)
                j = j + 1
             endif
          enddo

          ! Replace with filtered versions
          deallocate(Observations)
          deallocate(Obs_pred)
          deallocate(Obs_pert)

          allocate(Observations(N_obs_actual))
          allocate(Obs_pred(N_obs_actual, N_ens))
          allocate(Obs_pert(N_obs_actual, N_ens))

          Observations = Observations_filtered
          Obs_pred = Obs_pred_filtered
          Obs_pert = Obs_pert_filtered

          deallocate(Observations_filtered)
          deallocate(Obs_pred_filtered)
          deallocate(Obs_pert_filtered)

          ! Update counts
          N_obs_size = N_obs_actual
          Nobs = N_obs_actual

          write(LIS_logunit,*) '[INFO] LNETF using filtered observations: N_obs_size=', &
               N_obs_size

!----------------------------------------------------------------------------
! PDAF-style: Extract observation lat/lon for distance calculation
!----------------------------------------------------------------------------
          allocate(obs_lon_arr(N_obs_size))
          allocate(obs_lat_arr(N_obs_size))
          do jj=1,N_obs_size
             obs_lon_arr(jj) = Observations(jj)%lon
             obs_lat_arr(jj) = Observations(jj)%lat
          enddo

          ! Gaspari-Cohn compact support: weight=0 when d >= 2
          max_dist = 2.0 * sqrt(xcompact**2 + ycompact**2)

          write(LIS_logunit,*) '[INFO] Local LNETF: max observation distance (deg):', max_dist
          write(LIS_logunit,*) '[INFO] Local LNETF: total filtered observations:', N_obs_size

          ! DEBUG: Print observation location range
          write(LIS_logunit,*) '[DEBUG-LLNETF] Observation location range:'
          write(LIS_logunit,*) '[DEBUG-LLNETF]   lon: min=', minval(obs_lon_arr), &
               ' max=', maxval(obs_lon_arr)
          write(LIS_logunit,*) '[DEBUG-LLNETF]   lat: min=', minval(obs_lat_arr), &
               ' max=', maxval(obs_lat_arr)

          ! For km-based localization, calculate max_dist in km
          if (use_km_loc) then
             ! Search radius = 2 * compact_support = 2 * 3σ = 6σ
             ! This ensures we find all observations that might have non-zero weight
             max_dist = 6.0 * lnetf_struc(n,k)%localization_scale_km  ! km
             write(LIS_logunit,*) '[INFO] Local LNETF: max search distance (km):', max_dist
          endif
       endif

!----------------------------------------------------------------------------
!  Assemble observation covariances.
!----------------------------------------------------------------------------
       allocate(obs_param(LIS_rc%nobtypes(k)))
       call generateObsparam(Nobjs, LIS_OBS_Pert_State(n,k),obs_param)

!----------------------------------------------------------------------------
! Retrieve the state variables
!----------------------------------------------------------------------------
       allocate(stvar(N_state,state_size))
       allocate(state_incr(N_state,state_size))
       allocate(state_tmp(N_state,state_size))

       allocate(state_lat(N_state))
       allocate(state_lon(N_state))
       allocate(lats(state_size))
       allocate(lons(state_size))

       call LIS_surfaceModel_DAGetStateVar(n,k)
       call LIS_surfaceModel_DAScaleStateVar(n,k)

       call LIS_surfaceModel_DAextractStateVector(n,k,N_state,state_size,stvar)

       call LIS_surfaceModel_getlatlons(n,k,state_size,lats,lons)

       ! DEBUG: Print tile location range
       write(LIS_logunit,*) '[DEBUG-LLNETF] Tile location range:'
       write(LIS_logunit,*) '[DEBUG-LLNETF]   lon: min=', minval(lons), &
            ' max=', maxval(lons)
       write(LIS_logunit,*) '[DEBUG-LLNETF]   lat: min=', minval(lats), &
            ' max=', maxval(lats)

       state_incr = stvar
       state_tmp  = stvar

       ! Initialize Local LNETF diagnostic counters
       tiles_with_obs = 0
       tiles_no_obs = 0
       total_local_obs_count = 0
       min_local_obs = 999999
       max_local_obs = 0

! ========================================================================
! PDAF-style Local LNETF: Main tile loop
! For each tile, find observations within localization radius
! ========================================================================
       do i=1,state_size/LIS_rc%nensem(n)

          obspred_flag = .true.
          tileid = (i-1)*LIS_rc%nensem(n)+1

          ! Get tile location
          tile_lon = lons(tileid)
          tile_lat = lats(tileid)
          state_lat(:) = tile_lat
          state_lon(:) = tile_lon

          ! Count observations within localization radius
          N_local_obs = 0
          do jj=1,N_obs_size
             if (use_km_loc) then
                ! km-based distance using Haversine formula
                dist = haversine_km(tile_lon, tile_lat, obs_lon_arr(jj), obs_lat_arr(jj))
             else
                ! degree-based distance (original method)
                dist = sqrt((obs_lon_arr(jj)-tile_lon)**2 + (obs_lat_arr(jj)-tile_lat)**2)
             endif
             if(dist < max_dist .and. Observations(jj)%assim) then
                N_local_obs = N_local_obs + 1
             endif
          enddo

          ! Skip if no nearby observations
          if(N_local_obs == 0) then
             tiles_no_obs = tiles_no_obs + 1
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0
             lnetf_struc(n,k)%anlys_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0

             ! Store localization diagnostics for tile with no observations
             lnetf_struc(n,k)%n_local_obs(i) = 0.0
             lnetf_struc(n,k)%mean_loc_weight(i) = 0.0
             lnetf_struc(n,k)%max_loc_weight(i) = 0.0
             lnetf_struc(n,k)%min_loc_weight(i) = 0.0
             lnetf_struc(n,k)%eff_sample_size_tile(i) = 0.0
             cycle
          endif

          ! Update diagnostic counters
          tiles_with_obs = tiles_with_obs + 1
          total_local_obs_count = total_local_obs_count + N_local_obs
          if(N_local_obs < min_local_obs) min_local_obs = N_local_obs
          if(N_local_obs > max_local_obs) max_local_obs = N_local_obs

          ! Allocate local observation index array
          allocate(local_obs_idx(N_local_obs))

          ! Build index of nearby observations
          local_count = 0
          do jj=1,N_obs_size
             if (use_km_loc) then
                ! km-based distance using Haversine formula
                dist = haversine_km(tile_lon, tile_lat, obs_lon_arr(jj), obs_lat_arr(jj))
             else
                ! degree-based distance (original method)
                dist = sqrt((obs_lon_arr(jj)-tile_lon)**2 + (obs_lat_arr(jj)-tile_lat)**2)
             endif
             if(dist < max_dist .and. Observations(jj)%assim) then
                local_count = local_count + 1
                local_obs_idx(local_count) = jj
             endif
          enddo

          ! Use only local observations
          N_selected_obs = N_local_obs
          assim = .true.

          allocate(obs_da(N_selected_obs))
          allocate(obspred_da(N_selected_obs,N_ens))
          allocate(obspert_da(N_selected_obs,N_ens))
          allocate(obs_cov(N_selected_obs, N_selected_obs))

          ! Fill local observation arrays using index
          do kk=1,N_local_obs
             jj = local_obs_idx(kk)
             obs_da(kk) = Observations(jj)
             obspred_da(kk,:) = Obs_pred(jj,:)
             obspert_da(kk,:) = Obs_pert(jj,:)

             ! Check for undefined predictions
             do mm=1,N_ens
                if(obspred_da(kk,mm).eq.LIS_rc%udef) then
                   obspred_flag = .false.
                endif
             enddo
          enddo

          call assemble_obs_cov(LIS_rc%nobtypes(k), N_selected_obs, &
               obs_param,obs_da,Obs_cov)

          if(assim.and.obspred_flag) then

             ! DEBUG: Print observation statistics (for first 5 tiles WITH observations)
             if(tiles_with_obs .le. 5) then
                write(LIS_logunit,*) '[DEBUG-LLNETF] =========================================='
                write(LIS_logunit,*) '[DEBUG-LLNETF] Local LNETF Statistics (Tile ', i, '):'
                write(LIS_logunit,*) '[DEBUG-LLNETF]   Tile location: lon=', tile_lon, ' lat=', tile_lat
                write(LIS_logunit,*) '[DEBUG-LLNETF]   Local observations within radius: ', N_local_obs
             endif

! ========================================================================
! Local LNETF: Call lnetf_analysis with only nearby observations
! Localization (Gaspari-Cohn weighting) is applied inside lnetf_analysis
! ========================================================================
             call lnetf_analysis(i,N_state,N_selected_obs, N_ens, &
                  obs_da,                                        &
                  obspred_da,                       &
                  obspert_da,                       &
                  Obs_cov,              &
                  state_incr(:, ((i-1)*N_ens+1):((i-1)*N_ens+N_ens)),&
                  state_lon, state_lat,xcompact,ycompact, &
                  tile_n_eff, tile_mean_loc_weight, &
                  tile_max_loc_weight, tile_min_loc_weight)

             ! Store localization diagnostics for this tile
             lnetf_struc(n,k)%n_local_obs(i) = real(N_local_obs)
             lnetf_struc(n,k)%mean_loc_weight(i) = tile_mean_loc_weight
             lnetf_struc(n,k)%max_loc_weight(i) = tile_max_loc_weight
             lnetf_struc(n,k)%min_loc_weight(i) = tile_min_loc_weight
             lnetf_struc(n,k)%eff_sample_size_tile(i) = tile_n_eff

             !--------------------------------------------------------------
             ! CRITICAL FIX: Compute actual increment = X^a - X^f
             ! After lnetf_analysis, state_incr contains X^a (analysis state)
             ! state_tmp contains X^f (forecast state)
             ! The increment is the difference: X^a - X^f
             !--------------------------------------------------------------
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
                  state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) - &
                  state_tmp(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens)

          else
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0

             ! Store localization diagnostics for skipped tile
             lnetf_struc(n,k)%n_local_obs(i) = real(N_local_obs)
             lnetf_struc(n,k)%mean_loc_weight(i) = 0.0
             lnetf_struc(n,k)%max_loc_weight(i) = 0.0
             lnetf_struc(n,k)%min_loc_weight(i) = 0.0
             lnetf_struc(n,k)%eff_sample_size_tile(i) = 0.0
          endif

          lnetf_struc(n,k)%anlys_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
               state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens)

          deallocate(obs_da)
          deallocate(obspred_da)
          deallocate(obspert_da)
          deallocate(obs_cov)
          deallocate(local_obs_idx)

       enddo

!----------------------------------------------------------------------------
! Local LNETF Summary Statistics
!----------------------------------------------------------------------------
       if(tiles_with_obs > 0) then
          avg_local_obs = real(total_local_obs_count) / real(tiles_with_obs)
       else
          avg_local_obs = 0.0
          min_local_obs = 0
       endif

       write(LIS_logunit,*) '[INFO] =========================================='
       write(LIS_logunit,*) '[INFO] PDAF-style Local LNETF Summary:'
       write(LIS_logunit,*) '[INFO]   Total tiles processed:', &
            state_size/LIS_rc%nensem(n)
       write(LIS_logunit,*) '[INFO]   Tiles with local observations:', tiles_with_obs
       write(LIS_logunit,*) '[INFO]   Tiles without observations:', tiles_no_obs
       if(tiles_with_obs > 0) then
          write(LIS_logunit,*) '[INFO]   Min local obs per tile:', min_local_obs
          write(LIS_logunit,*) '[INFO]   Max local obs per tile:', max_local_obs
          write(LIS_logunit,*) '[INFO]   Avg local obs per tile:', avg_local_obs
       endif
       write(LIS_logunit,*) '[INFO] =========================================='

       call LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,.true.)

!----------------------------------------------------------------------------
! Updating State vector and increments state
!----------------------------------------------------------------------------
       call LIS_surfaceModel_DAsetAnlysisUpdates(n,k,N_state,state_size,&
            stvar,state_incr)

       call LIS_surfaceModel_DADescaleStatevar(n,k)

!----------------------------------------------------------------------------
! Cleanup
!----------------------------------------------------------------------------
       deallocate(obs_param)
       deallocate(stvar)
       deallocate(State_incr)
       deallocate(state_tmp)
       deallocate(Observations)
       deallocate(Obs_pred)
       deallocate(Obs_pert)
       deallocate(state_lat)
       deallocate(state_lon)
       deallocate(lats)
       deallocate(lons)
       if(allocated(obs_lon_arr)) deallocate(obs_lon_arr)
       if(allocated(obs_lat_arr)) deallocate(obs_lat_arr)
    end if

  end subroutine lnetf_increments

!BOP
! 
! !ROUTINE: lnetf_update
! \label{lnetf_update}
!
! !INTERFACE: 
  subroutine lnetf_update(n,k)
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
       if(LIS_rc%incroption(k).eq.0) then !include analysis increments
         
          call LIS_surfaceModel_DAUpdateState(n,k)
!----------------------------------------------------------------------
!  Update the state variables
!----------------------------------------------------------------------       
          call LIS_surfaceModel_DAQCstate(n,k)

          call LIS_surfaceModel_DASetStateVar(n,k)
!----------------------------------------------------------------------
!  compute analysis residuals
!---------------------------------------------------------------------
          if(LIS_rc%winnov(k).eq.1) then 
             N_ens = LIS_rc%nensem(n)
             call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
             call LIS_verify(status, &
                  'ESMF_StateGet failed in lnetf_update')
             
             call ESMF_AttributeGet(LIS_OBS_State(n,k),&
                  name="Number Of Observations",&
                  value=N_obs_size,rc=status)
             call LIS_verify(status, &
                  'ESMF_AttributeGet: Number of Observations failed in lnetf_update')
             
             Nobs = Nobjs*N_obs_size
             allocate(Observations(Nobs))
             
             call generateObservations(n, k,Nobjs, Nobs, LIS_OBS_State(n,k), &
                  LIS_OBS_Pert_State(n,k),Observations)
             
             allocate(Obs_pred(Nobs,N_ens))               
             
             call LIS_surfaceModel_DAGetObsPred(n,k,Obs_pred)

            
             do i=1,Nobs
                if(Observations(i)%assim) then
                   lnetf_struc(n,k)%anlys_res(i) = Observations(i)%value - &
                        sum(Obs_pred(i,:))/real(LIS_rc%nensem(n))
                else
                   lnetf_struc(n,k)%anlys_res(i) = LIS_rc%udef
                endif
             enddo
             
             deallocate(Observations)
             deallocate(Obs_pred)
          endif
       endif
!----------------------------------------------------------------------
!  Cleanup
!---------------------------------------------------------------------
       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
            .true., rc=status)
       call LIS_verify(status, &
            'ESMF_AttributeSet failed in lnetf_update')
       
       write(LIS_logunit,*) '[INFO] Finished assimilating Observations using LNETF'

       !----------------------------------------------------------------------
       ! Print increment statistics for quality check
       !----------------------------------------------------------------------
       call print_increment_statistics(n, k)
    else
       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status, &
            'ESMF_AttributeSet failed in lnetf_update')
    endif
    
  end subroutine lnetf_update


!BOP
! 
! !ROUTINE: lnetf_diagnostics
!  \label{lnetf_diagnostics}
! 
! !INTERFACE:
  subroutine lnetf_diagnostics(n,k)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n 
    integer, intent(IN)    :: k
!  
! !DESCRIPTION:  
!  This subroutine generates the LNETF related diagnostics and outputs
!  it to a file. This includes a text output of selected ensemble
!  members, their mean, standard deviation, observations, and normalized
!  innovations. The frequency of diagnostic outputs can be specified
!  in the LIS configuration file.
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!    call to create output directory for DA statistics
!   \item[pruneVarname](\ref{pruneVarname_lnetf}) \newline
!    trims the variable name, eliminating white spaces
!   \item[LIS\_create\_stats\_filename](\ref{LIS_create_stats_filename}) \newline
!    creates the filename for statistics 
!   \item[LIS\_create\_innov\_filename](\ref{LIS_create_innov_filename}) \newline
!    creates the name of the innovations file
!   \item[LIS\_surfaceModel_DAgetstatevar]\ref{LIS_surfaceModel_DAgetstatevar}
!    obtain the specified prognostic variables. 
!   \item[LIS\_surfaceModel\_DAextractStateVector](\ref{LIS_surfaceModel_DAextractStateVector}) \newline
!    unpack the state and retrive the data
!  \end{description}
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
!--------------------------------------------------------------------------
! write innovations file with the following entries for all observation 
! types in each data assimilation instance. 
! 
! * Normalized innovations
! * Raw innovations
! * Ensemble spread
! * Analysis increments
! * Forecast variance
!--------------------------------------------------------------------------
       if(.not.lnetf_struc(n,k)%fileopen.and.LIS_masterproc) then
          
          call LIS_create_output_directory("LNETF")
          
       endif
       
       call writeInnovationOutput(n,k)

       call writeAnalysisIncr(n,k)

    endif

    if(alarmCheck) then 
       if(.not.lnetf_struc(n,k)%fileopen.and.LIS_masterproc) then
          
          call LIS_create_output_directory("LNETF")
          
       endif

       call writeEnsembleSpread(n,k)

    endif


#endif
  end subroutine lnetf_diagnostics

!BOP
! 
! !ROUTINE: writeInnovationOutput
! \label{writeInnovationOutput}
!
! !INTERFACE: 
  subroutine writeInnovationOutput(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n 
    integer,  intent(in)    :: k 
!
! !DESCRIPTION: 
! 
!  This routine writes the innovation values (observation minus the model 
!  forecast) to an external file. 
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!   \item[k]    index of the observation datastream
!  \end{description}
! 
!
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
       deflate_level =9

       if(LIS_masterproc) then 
          call LIS_create_innov_filename(n,k, innovfile,&
               'LNETF')
          
#if (defined USE_NETCDF4)
          status = nf90_create(path=innovfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(innovfile)//&
               ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=innovfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(innovfile)//&
               ' failed in lnetf_Mod')
#endif
          
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%obs_gnc(k),&
               dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%obs_gnr(k),&
               dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att failed for missing_value in lnetf_mod')
          
!--------------------------------------------------------------------------
!  Normalized innovations -meta data
!--------------------------------------------------------------------------
          write(unit=finst, fmt='(i2.2)') k
          varname = "ninnov_"//trim(finst)
          vardimname = "ninnov_"//trim(finst)//"_levels"
          standard_name = "Normalized_innovations_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for ninnov_'//trim(finst))
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=ninnov_Id),&
               'nf90_def_var failed for '//trim(varname)//' in lnetf_mod')
          
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               ninnov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for ninnov in lnetf_mod')      
#endif
          call LIS_verify(nf90_put_att(ftn,ninnov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for '//trim(standard_name)//' in lnetf_mod')
          
!--------------------------------------------------------------------------
!  Innovations -meta data
!--------------------------------------------------------------------------
          varname = "innov_"//trim(finst)
          vardimname = "innov_"//trim(finst)//"_levels"
          standard_name = "Innovations_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for innov_'//trim(finst))

          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=innov_Id),&
               'nf90_def_var failed for innov')
          
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               innov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for innov failed in lnetf_mod')             
#endif
          call LIS_verify(nf90_put_att(ftn,innov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for innov in lnetf_mod')
          
!--------------------------------------------------------------------------
!  analysis residuals -meta data
!--------------------------------------------------------------------------
          write(unit=finst, fmt='(i2.2)') k
          varname = "analysis_residual_"//trim(finst)
          vardimname = "analysis_residual_"//trim(finst)//"_levels"
          standard_name = "Analysis_residuals_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for analysis_residual_'//trim(finst))

          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=ares_Id),&
               'nf90_def_var failed for '//trim(varname)//' in lnetf_mod')
          
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               ares_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for analysis_residual in lnetf_mod')      
#endif
          call LIS_verify(nf90_put_att(ftn,ares_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for '//trim(standard_name)//' in lnetf_mod')
                    
!--------------------------------------------------------------------------
!  Forecast-variance -meta data
!--------------------------------------------------------------------------
          varname = "forecast_sigma_"//trim(finst)
          vardimname = "forecast_sigma_"//trim(finst)//"_levels"
          standard_name = "Forecast_variance_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for forecast_sigma_'//trim(finst))

          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=forecast_sigma_Id),&
               'nf90_def_var for forecast_sigma failed in lnetf_mod')
             
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               forecast_sigma_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var for forecast_sigma failed in lnetf_mod')             
#endif
          call LIS_verify(nf90_put_att(ftn,forecast_sigma_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for forecast_sigma failed in lnetf_mod')

          call LIS_verify(nf90_enddef(ftn),&
               'nf90_enddef failed in lnetf_mod')
       endif
       
       call LIS_writevar_innov(ftn,n, k, ninnov_id, &
            lnetf_struc(n,k)%norm_innov)
       call LIS_writevar_innov(ftn,n, k, innov_id, &
            lnetf_struc(n,k)%innov)
       call LIS_writevar_innov(ftn,n, k, ares_id, &
            lnetf_struc(n,k)%anlys_res)

       call LIS_writevar_innov(ftn,n, k, forecast_sigma_id, &
            lnetf_struc(n,k)%forecast_var)
       
       if(LIS_masterproc) then 
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in lnetf_mod')
       endif

    endif
  end subroutine writeInnovationOutput


!BOP
! 
! !ROUTINE: writeEnsembleSpread
! \label{writeEnsembleSpread}
!
! !INTERFACE: 
  subroutine writeEnsembleSpread(n,k)
!
! !DESCRIPTION: 
!  This routine writes the ensemble spread (standard deviation) 
!  of the model state vector. 
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!   \item[k]    index of the observation datastream
!  \end{description}

!EOP
    integer,  intent(in)    :: n 
    integer,  intent(in)    :: k 

    integer                :: ftn 
    integer                :: v
    character(len=LIS_CONST_PATH_LEN) :: spreadfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: ensspread_id(LIS_rc%nstvars(k))
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status
    integer                :: state_size
    real, allocatable      :: stvar(:,:)
    character*100,    allocatable     :: state_objs(:)

    state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)

    if(LIS_rc%wensems(k).eq.1) then 

       shuffle = 1
       deflate = 1
       deflate_level =9
       
       if(LIS_masterproc) then 
          call LIS_create_daspread_filename(n,k,spreadfile,&
               'LNETF')
          
#if (defined USE_NETCDF4)
          status = nf90_create(path=spreadfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=spreadfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in lnetf_Mod')
#endif
          
          if(LIS_rc%wopt.eq."1d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'ngrid',&
                  LIS_rc%glbngrid_red(n),&
                  dimID(1)),'nf90_def_dim for ngrid failed in lnetf_mod')
          elseif(LIS_rc%wopt.eq."2d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
                  dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
             call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
                  dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')
          endif

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att for missing_value failed in lnetf_mod')
          
!--------------------------------------------------------------------------
!  Ensemble spread -meta data
!--------------------------------------------------------------------------
          allocate(state_objs(LIS_rc%nstvars(k)))          
          call LIS_surfaceModel_DAgetStateVarNames(n,k,state_objs)

          do v = 1, LIS_rc%nstvars(k)
             write(unit=finst, fmt='(i2.2)') k
             varname = "ensspread_"//trim(state_objs(v))//"_"//trim(finst)
             vardimname = "ensspread_"//trim(state_objs(v))//&
                  "_"//trim(finst)//"_levels"
             standard_name = "Ensemble_spread_for_DA_instance_"//&
                  trim(state_objs(v))//"_"//&
                  trim(finst)

             if(LIS_rc%wopt.eq."1d gridspace") then           
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1), varID=ensspread_Id(v)),&
                     'nf90_def_var for ensspread failed in lnetf_mod')
                
             elseif(LIS_rc%wopt.eq."2d gridspace") then 
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1:2), varID=ensspread_Id(v)),&
                     'nf90_def_var for ensspread failed in lnetf_mod') 
             endif
          
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  ensspread_Id(v),&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate for ensspread failed in lnetf_mod')             
#endif
             call LIS_verify(nf90_put_att(ftn,ensspread_Id(v),&
                  "standard_name",standard_name),&
                  'nf90_put_att for ensspread failed in lnetf_mod')
             call LIS_verify(nf90_enddef(ftn),&
                  'nf90_enddef failed in lnetf_mod')
          end do
          deallocate(state_objs)          
       endif
       
       allocate(stvar(LIS_rc%nstvars(k),&
            state_size))
              
       call LIS_surfaceModel_DAGetStateVar(n,k)

       call LIS_surfaceModel_DAextractStateVector(n,k,&
            LIS_rc%nstvars(k),state_size, stvar)

       do v=1,LIS_rc%nstvars(k)
          call LIS_writevar_spread(ftn,n,k,ensspread_id(v), &
               stvar(v,:),v)
       enddo
       
       deallocate(stvar)

       if(LIS_masterproc) then 
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in lnetf_mod')
       endif
    endif

  end subroutine writeEnsembleSpread


!BOP
! 
! !ROUTINE: writeAnalysisIncr
! \label{writeAnalysisIncr}
!
! !INTERFACE: 
  subroutine writeAnalysisIncr(n,k)
!
! !DESCRIPTION: 
!  This routine writes the analysis increments from DA
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!   \item[k]    index of the observation datastream
!  \end{description}

!EOP
    integer,  intent(in)    :: n 
    integer,  intent(in)    :: k 

    integer                :: ftn
    integer                :: v
    character(len=LIS_CONST_PATH_LEN) :: incrfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: incr_id(LIS_rc%nstvars(k))
    integer                :: n_local_obs_id, mean_loc_weight_id
    integer                :: max_loc_weight_id, min_loc_weight_id
    integer                :: eff_sample_size_id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status
    real, allocatable      :: stvar(:,:)
    character*100,    allocatable     :: state_objs(:)

    ! Temporary arrays for expanding tile-level diagnostics to patch-level
    real, allocatable      :: diag_expanded(:)
    integer                :: ntiles, npatch, i, j, nensem

    if(LIS_rc%wensems(k).eq.1) then 

       shuffle = 1
       deflate = 1
       deflate_level =9
       
       if(LIS_masterproc) then 
          call LIS_create_incr_filename(n,k,incrfile,&
               'LNETF')
          
#if (defined USE_NETCDF4)
          status = nf90_create(path=incrfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(incrfile)//&
               ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=incrfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(incrfile)//&
               ' failed in lnetf_Mod')
#endif
          
          if(LIS_rc%wopt.eq."1d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'ngrid',&
                  LIS_rc%glbngrid_red(n),&
                  dimID(1)),'nf90_def_dim for ngrid failed in lnetf_mod')
          elseif(LIS_rc%wopt.eq."2d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
                  dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
             call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
                  dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')
          endif

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att for missing_value failed in lnetf_mod')
          
!--------------------------------------------------------------------------
!  Ensemble incr -meta data
!--------------------------------------------------------------------------
          allocate(state_objs(LIS_rc%nstvars(k)))          
          call LIS_surfaceModel_DAgetStateVarNames(n,k,state_objs)

          do v = 1, LIS_rc%nstvars(k)
             write(unit=finst, fmt='(i2.2)') k
             varname = "anlys_incr_"//trim(state_objs(v))//"_"//trim(finst)
             vardimname = "anlys_incr_"//trim(state_objs(v))//&
                  "_"//trim(finst)//"_levels"
             standard_name = "Analysis_incr_for_DA_instance_"//&
                  trim(state_objs(v))//"_"//&
                  trim(finst)

             if(LIS_rc%wopt.eq."1d gridspace") then           
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1), varID=incr_Id(v)),&
                     'nf90_def_var for incr failed in lnetf_mod')
                
             elseif(LIS_rc%wopt.eq."2d gridspace") then 
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1:2), varID=incr_Id(v)),&
                     'nf90_def_var for incr failed in lnetf_mod') 
             endif
          
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  incr_Id(v),&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate for incr failed in lnetf_mod')             
#endif
             call LIS_verify(nf90_put_att(ftn,incr_Id(v),&
                  "standard_name",standard_name),&
                  'nf90_put_att for incr failed in lnetf_mod')
          end do
          deallocate(state_objs)

!--------------------------------------------------------------------------
!  Localization diagnostics - variable definitions
!--------------------------------------------------------------------------
          write(unit=finst, fmt='(i2.2)') k

!--------------------------------------------------------------------------
!  n_local_obs: Number of local observations within localization radius
!--------------------------------------------------------------------------
          varname = "n_local_obs_"//trim(finst)
          standard_name = "Number_of_local_observations_for_DA_instance_"//&
               trim(finst)

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=n_local_obs_Id),&
                  'nf90_def_var for n_local_obs failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=n_local_obs_Id),&
                  'nf90_def_var for n_local_obs failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               n_local_obs_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for n_local_obs failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,n_local_obs_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for n_local_obs failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,n_local_obs_Id,&
               "long_name","Number of observations within localization radius"),&
               'nf90_put_att for n_local_obs long_name failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,n_local_obs_Id,&
               "units","count"),&
               'nf90_put_att for n_local_obs units failed in lnetf_mod')

!--------------------------------------------------------------------------
!  mean_loc_weight: Mean Gaspari-Cohn localization weight
!--------------------------------------------------------------------------
          varname = "mean_loc_weight_"//trim(finst)
          standard_name = "Mean_localization_weight_for_DA_instance_"//&
               trim(finst)

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=mean_loc_weight_Id),&
                  'nf90_def_var for mean_loc_weight failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=mean_loc_weight_Id),&
                  'nf90_def_var for mean_loc_weight failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               mean_loc_weight_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for mean_loc_weight failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,mean_loc_weight_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for mean_loc_weight failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,mean_loc_weight_Id,&
               "long_name","Mean Gaspari-Cohn localization weight (0=edge, 1=center)"),&
               'nf90_put_att for mean_loc_weight long_name failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,mean_loc_weight_Id,&
               "units","dimensionless"),&
               'nf90_put_att for mean_loc_weight units failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,mean_loc_weight_Id,&
               "valid_range", (/0.0, 1.0/)),&
               'nf90_put_att for mean_loc_weight valid_range failed in lnetf_mod')

!--------------------------------------------------------------------------
!  max_loc_weight: Maximum Gaspari-Cohn localization weight (closest obs)
!--------------------------------------------------------------------------
          varname = "max_loc_weight_"//trim(finst)
          standard_name = "Max_localization_weight_for_DA_instance_"//&
               trim(finst)

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=max_loc_weight_Id),&
                  'nf90_def_var for max_loc_weight failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=max_loc_weight_Id),&
                  'nf90_def_var for max_loc_weight failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               max_loc_weight_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for max_loc_weight failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,max_loc_weight_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for max_loc_weight failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,max_loc_weight_Id,&
               "long_name","Maximum Gaspari-Cohn weight (closest observation)"),&
               'nf90_put_att for max_loc_weight long_name failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,max_loc_weight_Id,&
               "units","dimensionless"),&
               'nf90_put_att for max_loc_weight units failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,max_loc_weight_Id,&
               "valid_range", (/0.0, 1.0/)),&
               'nf90_put_att for max_loc_weight valid_range failed in lnetf_mod')

!--------------------------------------------------------------------------
!  min_loc_weight: Minimum Gaspari-Cohn localization weight (farthest obs)
!--------------------------------------------------------------------------
          varname = "min_loc_weight_"//trim(finst)
          standard_name = "Min_localization_weight_for_DA_instance_"//&
               trim(finst)

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=min_loc_weight_Id),&
                  'nf90_def_var for min_loc_weight failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=min_loc_weight_Id),&
                  'nf90_def_var for min_loc_weight failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               min_loc_weight_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for min_loc_weight failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,min_loc_weight_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for min_loc_weight failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,min_loc_weight_Id,&
               "long_name","Minimum Gaspari-Cohn weight (farthest observation in radius)"),&
               'nf90_put_att for min_loc_weight long_name failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,min_loc_weight_Id,&
               "units","dimensionless"),&
               'nf90_put_att for min_loc_weight units failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,min_loc_weight_Id,&
               "valid_range", (/0.0, 1.0/)),&
               'nf90_put_att for min_loc_weight valid_range failed in lnetf_mod')

!--------------------------------------------------------------------------
!  eff_sample_size: Effective sample size (N_eff) for filter degeneracy
!--------------------------------------------------------------------------
          varname = "eff_sample_size_"//trim(finst)
          standard_name = "Effective_sample_size_for_DA_instance_"//&
               trim(finst)

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=eff_sample_size_Id),&
                  'nf90_def_var for eff_sample_size failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=eff_sample_size_Id),&
                  'nf90_def_var for eff_sample_size failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               eff_sample_size_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for eff_sample_size failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for eff_sample_size failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
               "long_name","Effective sample size N_eff = 1/sum(w_i^2)"),&
               'nf90_put_att for eff_sample_size long_name failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
               "units","count"),&
               'nf90_put_att for eff_sample_size units failed in lnetf_mod')
          call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
               "comment","Higher values (close to N_ens) indicate healthy filter; Low values indicate weight collapse"),&
               'nf90_put_att for eff_sample_size comment failed in lnetf_mod')

          call LIS_verify(nf90_enddef(ftn),&
               'nf90_enddef failed in lnetf_mod')
       endif
       
       do v=1,LIS_rc%nstvars(k)
          call LIS_writevar_incr(ftn,n,k,incr_id(v), &
               lnetf_struc(n,k)%anlys_incr(v,:),v)
       enddo

       !----------------------------------------------------------------------
       ! Write localization diagnostic variables
       ! These are tile-based values that must be expanded to patch-level
       ! (replicate each tile value to all ensemble members of that tile)
       ! Layout: [tile1_ens1, tile1_ens2, ..., tile1_ensN, tile2_ens1, ...]
       !----------------------------------------------------------------------
       npatch = LIS_surfaceModel_DAgetStateSpaceSize(n,k)
       nensem = LIS_rc%nensem(n)
       ntiles = npatch / nensem

       allocate(diag_expanded(npatch))

       ! Expand n_local_obs: tile values -> all ensemble members
       do i = 1, ntiles
          do j = 1, nensem
             diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%n_local_obs(i)
          enddo
       enddo
       call LIS_writevar_incr(ftn,n,k,n_local_obs_id, diag_expanded, 1)

       ! Expand mean_loc_weight: tile values -> all ensemble members
       do i = 1, ntiles
          do j = 1, nensem
             diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%mean_loc_weight(i)
          enddo
       enddo
       call LIS_writevar_incr(ftn,n,k,mean_loc_weight_id, diag_expanded, 1)

       ! Expand max_loc_weight: tile values -> all ensemble members
       do i = 1, ntiles
          do j = 1, nensem
             diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%max_loc_weight(i)
          enddo
       enddo
       call LIS_writevar_incr(ftn,n,k,max_loc_weight_id, diag_expanded, 1)

       ! Expand min_loc_weight: tile values -> all ensemble members
       do i = 1, ntiles
          do j = 1, nensem
             diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%min_loc_weight(i)
          enddo
       enddo
       call LIS_writevar_incr(ftn,n,k,min_loc_weight_id, diag_expanded, 1)

       ! Expand eff_sample_size_tile: tile values -> all ensemble members
       do i = 1, ntiles
          do j = 1, nensem
             diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%eff_sample_size_tile(i)
          enddo
       enddo
       call LIS_writevar_incr(ftn,n,k,eff_sample_size_id, diag_expanded, 1)

       deallocate(diag_expanded)

       if(LIS_masterproc) then
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in lnetf_mod')
       endif
    endif

  end subroutine writeAnalysisIncr

!BOP
! 
! !ROUTINE: pruneVarname
! \label{pruneVarname_lnetf}
!
! !INTERFACE:
  subroutine pruneVarname(varname)
! !ARGUMENTS: 
    character(len=*), intent(INOUT) :: varname
! 
! !DESCRIPTION:
! 
!  This routine generates a filename based on the names of the state
!  prognostic variables
! 
!  \begin{description}
!  \item[varname]  name of the variable filename
!  \end{description}
!EOP
    character*100                   :: temp
    character*1                     :: ftemp(100)
    character*1                     :: fdir(9)
    integer                         :: i,c
    
    write(unit=temp,fmt='(A100)') varname
    read(unit=temp,fmt='(100a1)') (ftemp(i),i=1,100)

    write(unit=temp,fmt='(A9)') 'EnkF/'
    read(unit=temp,fmt='(9a1)') (fdir(i),i=1,9)
    
    do i=100,1,-1
       if(ftemp(i).ne.(' ')) then 
          c = i
          exit
       endif
    enddo

    do i=1,c
       if(ftemp(i).eq.(' ')) ftemp(i) = '_'
    enddo
    write(unit=temp,fmt='(100a1)')(fdir(i),i=1,9),(ftemp(i),i=1,c)
    read(unit=temp,fmt='(a100)') varname

  end subroutine pruneVarname

!BOP
! 
! !ROUTINE: assemble_obs_cov
! \label{assemble_obs_cov_lnetf}
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
!
! assemble measurements error covariance (reichle, 27 Jul 2005)
! 
! \begin{description}
!  \item[nob] number of observation types
!  \item[N\_obs] number of observations
!  \item[obs\_param] observation settings
!  \item[Observations] Observations data type
!  \item[Obs\_cov]     observation error covariance
! \end{description}
!EOP
    
    integer :: i, j, i_species, j_species   !! inum, jnum
    
    real :: fac, xcorr_tmp, ycorr_tmp
    
    ! -------------------------------------------------------------
    
    ! assemble measurement error covariance 
    
    ! initialize

    Obs_cov = 0.

    ! diagonal elements

    do i=1,N_obs

       Obs_cov(i,i) = Observations(i)%std**2
       
    end do
    
    ! off-diagonal elements
    
    do i=1,N_obs
       do j=(i+1),N_obs
          
          i_species = Observations(i)%species
          j_species = Observations(j)%species
          
          ! have non-zero correlation only between observations of same type
          
          if (i_species == j_species) then
             
             xcorr_tmp = obs_param(i_species)%xcorr
             ycorr_tmp = obs_param(i_species)%ycorr
             
             ! check for zero correlation distance 
             
             if (xcorr_tmp>0. .and. ycorr_tmp>0.) then
                
                ! compute correlation between observation locations
                
                !!inum = Observations(i)%catnum
                !!jnum = Observations(j)%catnum 
                
                ! compute Gaussian correlation
                
                !!fac =  & 
                !!  ((tile_coord(inum)%com_lon-tile_coord(jnum)%com_lon)**2 &
                !!  /xcorr_tmp**2 )                                         &
                !!  +                                                       &
                !!  ((tile_coord(inum)%com_lat-tile_coord(jnum)%com_lat)**2 & 
                !!  /ycorr_tmp**2 )               
                
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
! !ROUTINE: lnetf_final
! \label{lnetf_final}
! 
! !INTERFACE: 
  subroutine lnetf_final
! 
! !DESCRIPTION: 
!  This method performs the finalization for all LNETF
!  related structures and subroutines. 
!
!EOP
    deallocate(lnetf_struc)
  end subroutine lnetf_final


!BOP
! 
! !ROUTINE: getObsPert
! \label{getObsPert_lnetf}
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
! 
!  This routine retrieves the perturbations to be applied to the
!  observations from the Observation perturbation state
! 
! The arguments are: 
! \begin{description}
!  \item[OBS\_Pert\_State]  Observation perturbation state 
!  \item[dim1,dim2]         dimensions of the perturbations array
!  \item[pert]              perturbation array
! \end{description}
!EOP
    real, pointer         :: obs_temp(:,:)
    integer               :: i 
    integer               :: obs_state_count
    integer               :: status
    character*100,allocatable     :: obs_state_objs(:)
    type(ESMF_Field), allocatable :: obs_field(:)

    
    call ESMF_StateGet(OBS_Pert_State,itemCount=obs_state_count,rc=status)
    call LIS_verify(status, &
         'ESMF_StateGet failed in getObsPert')
    
    allocate(obs_state_objs(obs_state_count))
    allocate(obs_field(obs_state_count))
    
    call ESMF_StateGet(OBS_Pert_State,itemNameList=obs_state_objs,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed in getObsPert')        

    do i=1,obs_state_count
       call ESMF_StateGet(OBS_Pert_State,obs_state_objs(i),&
            obs_field(i),rc=status)
       call LIS_verify(status,&
            'ESMF_StateGet failed in getObsPert')
       call ESMF_FieldGet(obs_field(i),localDE=0,farrayPtr=obs_temp,rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in getObsPert')
       pert((i-1)*gsize + 1:i*gsize,:) = obs_temp(:,:)
    enddo
    deallocate(obs_state_objs)
    deallocate(obs_field)

  end subroutine getObsPert

!BOP
! 
! !ROUTINE: generateObservations
! \label{generateObservations}
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
!   'obs\_type' datastructure. 
!  
!   \begin{description}
!   \item[n]  index of the nest
!   \item[Nobjs] Number of observation types
!   \item[OBS\_State] ESMF State containing observations
!   \item[Observations] Observations datastructure being returned. 
!   \end{description}
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
    integer, allocatable               :: gid1(:)
    real,    allocatable               :: obsstd(:)
    integer                        :: pert_type
    integer, allocatable               :: obsassimflag(:)
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
       call LIS_verify(status, &
            'ESMF_StateGet failed in generateObservations')

       call ESMF_StateGet(LIS_OBS_Pert_State,&
            "Observation"//vid(1)//vid(2),&
            pertField,rc=status)
       call LIS_verify(status,&
            'ESMF_StateGet failed in generateObservations')
       
       call ESMF_AttributeGet(pertField, "Perturbation Type",&
            pert_type, rc=status)
       call LIS_verify(status, &
            "Perturbation Type attribute not found in lnetf_Mod")

       call ESMF_AttributeGet(LIS_OBS_State,name="Number Of Observations",&
            value = counts,rc=status)
       call LIS_verify(status, &
            'Number of Observations attribute not found: lnetf_mod')

       call ESMF_FieldGet(valfield,localDE=0,farrayPtr=value1,rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in generateObservations')
       
       allocate(gid1(counts))
       allocate(obsassimflag(counts))
       allocate(obsstd(counts))

       if(counts.gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",obsstd,&
               rc=status)
          call LIS_verify(status,&
               'Standard Deviation attribute not found: lnetf_mod')

          call ESMF_AttributeGet(valfield,"Grid Number",gid1,&
               rc=status)
          call LIS_verify(status, &
               'ESMF_AttributeGet for Grid Number in generateObservations')
          
          call ESMF_AttributeGet(valfield,"Assimilation Flag",&
               obsassimflag,rc=status)
          
          call LIS_verify(status,&
               'ESMF_AttributeGet for Assimilation Flag in generateObservations')
       endif

       value(count1: count1+(counts-1))  = value1(:)          
       gid(count1:count1+(counts-1))     = gid1(:)       
       pertyp(count1:count1+(counts-1))  = pert_type
       if(pert_type.eq.0) then 
          std(count1:count1+(counts-1))     = obsstd(:)
       elseif(pert_type.eq.1) then 
          do g = count1,count1+(counts-1)
             if(value1(g).ne.-9999.0) then 
                std(g) = obsstd(g)*value1(g)
                if(std(g).eq.0) std(g) = 0.1
             else
                std(g) = -9999.0
             endif
          enddo
       endif
       species(count1:count1+(counts-1)) = typ
       assimflag(count1:count1+(counts-1)) = obsassimflag(:)
       count1 = count1 + counts
       typ    = typ+1
       
       deallocate(gid1)
       deallocate(obsassimflag)
       deallocate(obsstd)
    enddo

!----------------------------------------------------------------------------
! Map obs_state to Observations object
!----------------------------------------------------------------------------
    gsize = N_obs_size/Nobjs
    do t=1,N_obs_size
       gval = t - nint(real((t-1)/gsize))*gsize
       Observations(t)%species = species(t)
       Observations(t)%catnum = gid(t)

       col = LIS_obs_domain(n,k)%col(gval)
       row = LIS_obs_domain(n,k)%row(gval)
       Observations(t)%lon = LIS_obs_domain(n,k)%lon(&
            col+(row-1)*LIS_rc%obs_lnc(k))
       Observations(t)%lat = LIS_obs_domain(n,k)%lat(&
            col+(row-1)*LIS_rc%obs_lnc(k))
       Observations(t)%value = value(t)
       Observations(t)%std = std(t)
       Observations(t)%pert_type = pertyp(t)
       if(assimflag(t).eq.1) then 
          Observations(t)%assim = .true.
       else
          Observations(t)%assim = .false.
       endif
    enddo
    
  end subroutine generateObservations

!BOP
! 
! !ROUTINE: generateObsparam
! \label{generateObsparam_lnetf}
!
! !INTERFACE: 
  subroutine generateObsparam(Nobjs, OBS_Pert_State, obs_param)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer                       :: Nobjs
    type(ESMF_State)              :: OBS_Pert_State
    type(obs_param_type)          :: obs_param(Nobjs)
!
! !DESCRIPTION: 
!   This routine obtains the metadata information associated with 
!   the observations and puts them into the obs\_param datastructure
!  
!   The arguments are: 
!   \begin{description}
!    \item[Nobjs]            Number of observation types \newline
!    \item[OBS\_Pert\_State] ESMF state containing observation 
!                            perturbations
!    \item[obs\_param]       obs\_param datastructure being 
!                            returned
!   \end{description}
!
!EOP    
    integer                       :: i 
    character*100                 :: temp
    character*1                   :: vid(2)
    type(ESMF_Field)              :: pertField
    integer                       :: status
    real                          :: std_normal_max(Nobjs)
    real                          :: xcorr(Nobjs)
    real                          :: ycorr(Nobjs)

    do i=1,Nobjs
       write(unit=temp,fmt='(i2.2)') i
       read(unit=temp,fmt='(2a1)') vid
       
       call ESMF_StateGet(OBS_Pert_State,"Observation"//vid(1)//vid(2),&
            pertField,rc=status)
       call LIS_verify(status, &
            'ESMF_StateGet failed in generateObsParam')
       
       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max(i),&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet: Std Normal Max failed in generateObsparam')
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr(i),&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet: X Correlation Scale failed in generateObsparam')
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr(i),&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet: Y Correlation Scale failed in generateObsparam')        
       
       obs_param(i)%species = i
       obs_param(i)%assim = .TRUE.
       obs_param(i)%std_normal_max = std_normal_max(i)
       obs_param(i)%xcorr = xcorr(i)
       obs_param(i)%ycorr = ycorr(i)
    enddo
  end subroutine generateObsparam


!BOP
!
! !ROUTINE: print_increment_statistics
! \label{print_increment_statistics_lnetf}
!
! !INTERFACE:
  subroutine print_increment_statistics(n, k)
! !USES:
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n  ! nest index
    integer, intent(in) :: k  ! DA instance index
!
! !DESCRIPTION:
!   Print statistics of analysis increments for quality control.
!   Outputs min, max, mean, and std for each state variable.
!EOP
    integer :: v, i, count_valid
    real :: incr_min, incr_max, incr_mean, incr_std
    real :: incr_sum, incr_sumsq
    real :: val
    integer :: total_size

    write(LIS_logunit,*) '[INFO] =========================================='
    write(LIS_logunit,*) '[INFO] Analysis Increment Statistics:'
    write(LIS_logunit,*) '[INFO] =========================================='

    total_size = size(lnetf_struc(n,k)%anlys_incr, 2)

    do v = 1, LIS_rc%nstvars(k)
       incr_min = 1.0e30
       incr_max = -1.0e30
       incr_sum = 0.0
       incr_sumsq = 0.0
       count_valid = 0

       do i = 1, total_size
          val = lnetf_struc(n,k)%anlys_incr(v, i)
          ! Skip undefined values
          if (abs(val - LIS_rc%udef) > 1.0e-5 .and. abs(val) < 1.0e10) then
             count_valid = count_valid + 1
             if (val < incr_min) incr_min = val
             if (val > incr_max) incr_max = val
             incr_sum = incr_sum + val
             incr_sumsq = incr_sumsq + val*val
          endif
       end do

       if (count_valid > 0) then
          incr_mean = incr_sum / real(count_valid)
          if (count_valid > 1) then
             incr_std = sqrt((incr_sumsq - count_valid*incr_mean*incr_mean) / &
                            real(count_valid - 1))
          else
             incr_std = 0.0
          endif

          write(LIS_logunit,'(A,I2,A)') ' [INFO] State Variable ', v, ':'
          write(LIS_logunit,'(A,I10)')  '        Valid points: ', count_valid
          write(LIS_logunit,'(A,E12.4)')  '        Min:  ', incr_min
          write(LIS_logunit,'(A,E12.4)')  '        Max:  ', incr_max
          write(LIS_logunit,'(A,E12.4)')  '        Mean: ', incr_mean
          write(LIS_logunit,'(A,E12.4)')  '        Std:  ', incr_std

          ! Check for reasonable increment values (soil moisture: typically < 0.1 m3/m3)
          if (abs(incr_max) > 0.5 .or. abs(incr_min) > 0.5) then
             write(LIS_logunit,*) '        [WARN] Large increments detected!'
          endif
       else
          write(LIS_logunit,'(A,I2,A)') ' [INFO] State Variable ', v, ': No valid increments'
       endif
    end do

    write(LIS_logunit,*) '[INFO] =========================================='

  end subroutine print_increment_statistics

!BOP
!
! !ROUTINE: lnetf_multisource_increments
! \label{lnetf_multisource_increments}
!
! !INTERFACE:
  subroutine lnetf_multisource_increments(n)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
!
! !DESCRIPTION:
!  This routine computes analysis increments for multi-source simultaneous
!  data assimilation. Observations from all DA instances are combined into
!  a single observation vector with species tags, then a single LNETF
!  analysis is performed.
!
!  Key difference from sequential mode:
!  - Sequential: DA1 update -> DA2 update (order matters)
!  - Simultaneous: [DA1 + DA2 obs] -> single analysis (optimal)
!
!  The observation error covariance R is block-diagonal:
!    R = [R_DA1    0   ]
!        [  0    R_DA2 ]
!
!EOP
    integer :: k, kk, status
    integer :: N_obs_total, N_obs_actual
    integer :: Nobjs_k, N_obs_size_k, Nobs_k
    integer :: state_size, N_state, N_ens
    integer :: i, j, m, mm, t, jj
    real    :: dx, dy, xcompact, ycompact
    real    :: sigma_km, xcompact_km, ycompact_km
    real    :: km_per_deg_lat, km_per_deg_lon, mean_lat
    real, parameter :: DEG_TO_RAD = 3.14159265358979323846 / 180.0
    logical :: data_status_any
    logical :: data_status_k

    type(obs_type), allocatable :: Observations_all(:)
    type(obs_type), allocatable :: Observations_filtered(:)
    type(obs_type), allocatable :: Observations_k(:)
    real, allocatable :: Obs_pred_all(:,:), Obs_pert_all(:,:)
    real, allocatable :: Obs_pred_filtered(:,:), Obs_pert_filtered(:,:)
    real, allocatable :: Obs_pred_k(:,:), Obs_pert_k(:,:)

    integer :: obs_count, species_offset
    real, allocatable :: state_incr(:,:), state_tmp(:,:), stvar(:,:)
    real, allocatable :: obs_lon_arr(:), obs_lat_arr(:)
    real, allocatable :: lats(:), lons(:)
    real              :: max_dist, dist, tile_lon, tile_lat
    integer           :: tileid, N_local_obs, local_count, N_selected_obs
    integer, allocatable :: local_obs_idx(:)
    logical           :: use_km_loc, obspred_flag, assim

    ! Local observation arrays for lnetf_analysis
    type(obs_type), allocatable     :: obs_da(:)
    type(obs_param_type), allocatable :: obs_param(:)
    real, allocatable               :: obspred_da(:,:), obspert_da(:,:)
    real, allocatable               :: Obs_cov(:,:)

    ! Diagnostic variables
    integer :: tiles_with_obs, tiles_no_obs, total_local_obs_count
    integer :: min_local_obs, max_local_obs
    real    :: tile_n_eff, tile_mean_loc_weight, tile_max_loc_weight, tile_min_loc_weight

    ! Innovation computation variables
    real    :: Obs_pred_i_mean, innov_val
    real    :: std_innov_arr(1)  ! Array for row_variance output

    write(LIS_logunit,*) '[INFO] =============================================='
    write(LIS_logunit,*) '[INFO] LNETF Multi-Source Simultaneous DA'
    write(LIS_logunit,*) '[INFO]   Number of DA instances: ', LIS_rc%ndas
    write(LIS_logunit,*) '[INFO] =============================================='

    ! Check if any DA instance has data
    data_status_any = .false.
    do k = 1, LIS_rc%ndas
       call ESMF_AttributeGet(LIS_OBS_State(n,k), name="Data Update Status",&
            value=data_status_k, rc=status)
       if(data_status_k) data_status_any = .true.
       write(LIS_logunit,*) '[DEBUG] DA instance ', k, ' data_status = ', data_status_k
    enddo

    if(.not. data_status_any) then
       write(LIS_logunit,*) '[INFO] No observations from any DA instance, skipping multi-source analysis'
       return
    endif

    ! Get common parameters from first DA instance
    k = 1
    call LIS_getDomainResolutions(n, dx, dy)
    state_size = LIS_surfaceModel_DAgetStateSpaceSize(n, k)
    N_state = LIS_rc%nstvars(k)
    N_ens = LIS_rc%nensem(n)

    ! Initialize combined increments to zero for ALL DA instances
    ! This ensures all diagnostics outputs have valid data
    do kk = 1, LIS_rc%ndas
       lnetf_struc(n,kk)%anlys_incr = 0.0
       ! Initialize localization diagnostic arrays to undefined
       ! CRITICAL: Must initialize every cycle to prevent stale values from previous time steps
       lnetf_struc(n,kk)%n_local_obs = LIS_rc%udef
       lnetf_struc(n,kk)%mean_loc_weight = LIS_rc%udef
       lnetf_struc(n,kk)%max_loc_weight = LIS_rc%udef
       lnetf_struc(n,kk)%min_loc_weight = LIS_rc%udef
       lnetf_struc(n,kk)%eff_sample_size_tile = LIS_rc%udef
       ! Initialize innovation arrays if allocated (winnov flag controlled)
       if(allocated(lnetf_struc(n,kk)%innov)) then
          lnetf_struc(n,kk)%innov = LIS_rc%udef
          lnetf_struc(n,kk)%norm_innov = LIS_rc%udef
          lnetf_struc(n,kk)%forecast_var = LIS_rc%udef
          lnetf_struc(n,kk)%anlys_res = LIS_rc%udef
       endif
    enddo

    !--------------------------------------------------------------------------
    ! Phase 1: Count total observations from all DA instances
    !--------------------------------------------------------------------------
    N_obs_total = 0
    do kk = 1, LIS_rc%ndas
       call ESMF_StateGet(LIS_OBS_State(n,kk), itemCount=Nobjs_k, rc=status)
       call ESMF_AttributeGet(LIS_OBS_State(n,kk), name="Number Of Observations", &
            value=N_obs_size_k, rc=status)
       N_obs_total = N_obs_total + Nobjs_k * N_obs_size_k
    enddo

    write(LIS_logunit,*) '[INFO] Total observation grid points (all sources): ', N_obs_total

    if(N_obs_total == 0) then
       write(LIS_logunit,*) '[INFO] No observations in combined array, skipping analysis'
       return
    endif

    !--------------------------------------------------------------------------
    ! Phase 2: Allocate combined observation arrays
    !--------------------------------------------------------------------------
    allocate(Observations_all(N_obs_total))
    allocate(Obs_pred_all(N_obs_total, N_ens))
    allocate(Obs_pert_all(N_obs_total, N_ens))

    !--------------------------------------------------------------------------
    ! Phase 3: Collect observations from all DA instances
    !--------------------------------------------------------------------------
    obs_count = 0
    species_offset = 0

    do kk = 1, LIS_rc%ndas
       call ESMF_AttributeGet(LIS_OBS_State(n,kk), name="Data Update Status", &
            value=data_status_k, rc=status)

       if(.not. data_status_k) then
          write(LIS_logunit,*) '[INFO] DA instance ', kk, ' has no data, skipping'
          cycle
       endif

       call ESMF_StateGet(LIS_OBS_State(n,kk), itemCount=Nobjs_k, rc=status)
       call ESMF_AttributeGet(LIS_OBS_State(n,kk), name="Number Of Observations", &
            value=N_obs_size_k, rc=status)
       Nobs_k = Nobjs_k * N_obs_size_k

       if(Nobs_k == 0) cycle

       ! Get observations for this DA instance
       allocate(Observations_k(Nobs_k))
       allocate(Obs_pred_k(Nobs_k, N_ens))
       allocate(Obs_pert_k(Nobs_k, N_ens))

       call generateObservations(n, kk, Nobjs_k, Nobs_k, &
            LIS_OBS_State(n,kk), LIS_OBS_Pert_State(n,kk), Observations_k)
       call LIS_surfaceModel_DAGetObsPred(n, kk, Obs_pred_k)
       call getObsPert(LIS_OBS_Pert_State(n,kk), N_obs_size_k, Nobs_k, N_ens, Obs_pert_k)

       ! Copy to combined arrays with species offset
       do i = 1, Nobs_k
          obs_count = obs_count + 1
          Observations_all(obs_count) = Observations_k(i)
          ! Assign species based on DA instance (for block-diagonal R)
          Observations_all(obs_count)%species = kk
          Obs_pred_all(obs_count, :) = Obs_pred_k(i, :)
          Obs_pert_all(obs_count, :) = Obs_pert_k(i, :)
       enddo

       write(LIS_logunit,*) '[INFO] DA instance ', kk, ': collected ', Nobs_k, ' obs (species=', kk, ')'

       !-----------------------------------------------------------------------
       ! Compute and store innovations for this DA instance (for diagnostics output)
       ! This must be done BEFORE filtering to preserve spatial order
       !-----------------------------------------------------------------------
       if(LIS_rc%winnov(kk).eq.1) then
          do i = 1, Nobs_k
             if(Observations_k(i)%assim) then
                ! Compute innovation = observation - mean(model prediction)
                Obs_pred_i_mean = sum(Obs_pred_k(i,:)) / real(N_ens)
                innov_val = Observations_k(i)%value - Obs_pred_i_mean

                ! Compute forecast variance for normalization
                call row_variance(1, N_ens, Obs_pred_k(i,:), std_innov_arr(1))
                lnetf_struc(n,kk)%forecast_var(i) = std_innov_arr(1)
                std_innov_arr = std_innov_arr + (Observations_k(i)%std)**2
                std_innov_arr = sqrt(std_innov_arr)

                ! Store innovation and normalized innovation
                lnetf_struc(n,kk)%innov(i) = innov_val
                if(std_innov_arr(1) > 1.0e-10) then
                   lnetf_struc(n,kk)%norm_innov(i) = innov_val / std_innov_arr(1)
                else
                   lnetf_struc(n,kk)%norm_innov(i) = LIS_rc%udef
                endif
             else
                lnetf_struc(n,kk)%innov(i) = LIS_rc%udef
                lnetf_struc(n,kk)%norm_innov(i) = LIS_rc%udef
                lnetf_struc(n,kk)%forecast_var(i) = LIS_rc%udef
             endif
          enddo
       endif

       deallocate(Observations_k)
       deallocate(Obs_pred_k)
       deallocate(Obs_pert_k)
    enddo

    write(LIS_logunit,*) '[INFO] Total combined observations: ', obs_count

    !--------------------------------------------------------------------------
    ! Phase 4: Filter to keep only assim=.true. observations
    !--------------------------------------------------------------------------
    N_obs_actual = 0
    do i = 1, obs_count
       if(Observations_all(i)%assim) N_obs_actual = N_obs_actual + 1
    enddo

    write(LIS_logunit,*) '[INFO] Observations with assimflag=1: ', N_obs_actual

    if(N_obs_actual == 0) then
       write(LIS_logunit,*) '[INFO] No valid observations after QC filtering'
       deallocate(Observations_all)
       deallocate(Obs_pred_all)
       deallocate(Obs_pert_all)
       return
    endif

    !--------------------------------------------------------------------------
    ! Phase 5: TRUE SIMULTANEOUS LNETF analysis with combined observations
    !--------------------------------------------------------------------------
    ! All observations from multiple sources are combined with species tags:
    !   - DA instance 1 (AMSR2): species = 1
    !   - DA instance 2 (SMAP):  species = 2
    !
    ! The R matrix is block-diagonal: R = [R_DA1  0; 0  R_DA2]
    ! This is automatically achieved because assemble_obs_cov() sets
    ! correlation = 0 when species differ.
    !--------------------------------------------------------------------------

    write(LIS_logunit,*) '[INFO] =============================================='
    write(LIS_logunit,*) '[INFO] TRUE SIMULTANEOUS LNETF Analysis'
    write(LIS_logunit,*) '[INFO]   Total combined observations: ', obs_count
    write(LIS_logunit,*) '[INFO]   Valid obs (assimflag=1): ', N_obs_actual
    write(LIS_logunit,*) '[INFO]   Species mapping:'
    do kk = 1, LIS_rc%ndas
       write(LIS_logunit,*) '[INFO]     DA instance ', kk, ' -> species ', kk
    enddo
    write(LIS_logunit,*) '[INFO] =============================================='

    ! Filter observations to keep only assim=.true.
    allocate(Observations_filtered(N_obs_actual))
    allocate(Obs_pred_filtered(N_obs_actual, N_ens))
    allocate(Obs_pert_filtered(N_obs_actual, N_ens))

    j = 1
    do i = 1, obs_count
       if(Observations_all(i)%assim) then
          Observations_filtered(j) = Observations_all(i)
          Obs_pred_filtered(j, :) = Obs_pred_all(i, :)
          Obs_pert_filtered(j, :) = Obs_pert_all(i, :)
          j = j + 1
       endif
    enddo

    ! Get observation lat/lon arrays
    allocate(obs_lon_arr(N_obs_actual))
    allocate(obs_lat_arr(N_obs_actual))
    do jj = 1, N_obs_actual
       obs_lon_arr(jj) = Observations_filtered(jj)%lon
       obs_lat_arr(jj) = Observations_filtered(jj)%lat
    enddo

    ! Print observation range
    write(LIS_logunit,*) '[INFO] Combined observation location range:'
    write(LIS_logunit,*) '[INFO]   lon: min=', minval(obs_lon_arr), ' max=', maxval(obs_lon_arr)
    write(LIS_logunit,*) '[INFO]   lat: min=', minval(obs_lat_arr), ' max=', maxval(obs_lat_arr)

    !--------------------------------------------------------------------------
    ! Set up localization parameters (using k=1 settings)
    !--------------------------------------------------------------------------
    k = 1
    use_km_loc = lnetf_struc(n,k)%use_km_localization

    if(use_km_loc) then
       sigma_km = lnetf_struc(n,k)%localization_scale_km
       max_dist = 6.0 * sigma_km  ! Search radius = 6*sigma (compact support = 3*sigma)
       xcompact_km = 3.0 * sigma_km
       ycompact_km = 3.0 * sigma_km
       ! Convert to degrees for lnetf_analysis (approximate at mean latitude)
       mean_lat = (minval(obs_lat_arr) + maxval(obs_lat_arr)) / 2.0
       km_per_deg_lat = 111.0
       km_per_deg_lon = 111.0 * cos(mean_lat * DEG_TO_RAD)
       xcompact = xcompact_km / km_per_deg_lon
       ycompact = ycompact_km / km_per_deg_lat
       write(LIS_logunit,*) '[INFO] km-based localization:'
       write(LIS_logunit,*) '[INFO]   sigma = ', sigma_km, ' km'
       write(LIS_logunit,*) '[INFO]   max search distance = ', max_dist, ' km'
    else
       xcompact = lnetf_struc(n,k)%localization_factor * dx
       ycompact = lnetf_struc(n,k)%localization_factor * dy
       max_dist = 2.0 * sqrt(xcompact**2 + ycompact**2)
    endif

    !--------------------------------------------------------------------------
    ! Retrieve state variables (using k=1 since model state is shared)
    !--------------------------------------------------------------------------
    allocate(stvar(N_state, state_size))
    allocate(state_incr(N_state, state_size))
    allocate(state_tmp(N_state, state_size))
    allocate(lats(state_size))
    allocate(lons(state_size))

    call LIS_surfaceModel_DASetFreshIncrementsStatus(n, k, .false.)
    call LIS_surfaceModel_DAGetStateVar(n, k)
    call LIS_surfaceModel_DAScaleStateVar(n, k)
    call LIS_surfaceModel_DAextractStateVector(n, k, N_state, state_size, stvar)
    call LIS_surfaceModel_getlatlons(n, k, state_size, lats, lons)

    state_incr = stvar
    state_tmp = stvar
    lnetf_struc(n,k)%anlys_incr = 0.0

    ! Initialize localization diagnostics
    lnetf_struc(n,k)%n_local_obs = LIS_rc%udef
    lnetf_struc(n,k)%mean_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%max_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%min_loc_weight = LIS_rc%udef
    lnetf_struc(n,k)%eff_sample_size_tile = LIS_rc%udef

    ! Diagnostic counters
    tiles_with_obs = 0
    tiles_no_obs = 0
    total_local_obs_count = 0
    min_local_obs = 999999
    max_local_obs = 0

    ! Get observation parameters (use combined from all DA instances)
    allocate(obs_param(LIS_rc%ndas))
    do kk = 1, LIS_rc%ndas
       call ESMF_StateGet(LIS_OBS_Pert_State(n,kk), itemCount=Nobjs_k, rc=status)
       if(Nobjs_k > 0) then
          call generateObsparam(1, LIS_OBS_Pert_State(n,kk), obs_param(kk:kk))
       endif
    enddo

    !--------------------------------------------------------------------------
    ! LOCAL LNETF: Main tile loop (TRUE SIMULTANEOUS)
    ! For each tile, find observations from ALL sources within radius
    !--------------------------------------------------------------------------
    write(LIS_logunit,*) '[INFO] Starting local LNETF analysis over ', &
         state_size/N_ens, ' tiles...'

    do i = 1, state_size/N_ens
       obspred_flag = .true.
       tileid = (i-1)*N_ens + 1
       tile_lon = lons(tileid)
       tile_lat = lats(tileid)

       ! Count observations within localization radius (from ALL sources)
       N_local_obs = 0
       do jj = 1, N_obs_actual
          if(use_km_loc) then
             dist = haversine_km(tile_lon, tile_lat, obs_lon_arr(jj), obs_lat_arr(jj))
          else
             dist = sqrt((obs_lon_arr(jj)-tile_lon)**2 + (obs_lat_arr(jj)-tile_lat)**2)
          endif
          if(dist < max_dist .and. Observations_filtered(jj)%assim) then
             N_local_obs = N_local_obs + 1
          endif
       enddo

       ! Skip tile if no nearby observations
       if(N_local_obs == 0) then
          tiles_no_obs = tiles_no_obs + 1
          state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0
          lnetf_struc(n,k)%anlys_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0
          lnetf_struc(n,k)%n_local_obs(i) = 0.0
          lnetf_struc(n,k)%mean_loc_weight(i) = 0.0
          lnetf_struc(n,k)%max_loc_weight(i) = 0.0
          lnetf_struc(n,k)%min_loc_weight(i) = 0.0
          lnetf_struc(n,k)%eff_sample_size_tile(i) = 0.0
          cycle
       endif

       ! Update statistics
       tiles_with_obs = tiles_with_obs + 1
       total_local_obs_count = total_local_obs_count + N_local_obs
       if(N_local_obs < min_local_obs) min_local_obs = N_local_obs
       if(N_local_obs > max_local_obs) max_local_obs = N_local_obs

       ! Build index of nearby observations
       allocate(local_obs_idx(N_local_obs))
       local_count = 0
       do jj = 1, N_obs_actual
          if(use_km_loc) then
             dist = haversine_km(tile_lon, tile_lat, obs_lon_arr(jj), obs_lat_arr(jj))
          else
             dist = sqrt((obs_lon_arr(jj)-tile_lon)**2 + (obs_lat_arr(jj)-tile_lat)**2)
          endif
          if(dist < max_dist .and. Observations_filtered(jj)%assim) then
             local_count = local_count + 1
             local_obs_idx(local_count) = jj
          endif
       enddo

       ! Prepare local observation arrays
       N_selected_obs = N_local_obs
       assim = .true.

       allocate(obs_da(N_selected_obs))
       allocate(obspred_da(N_selected_obs, N_ens))
       allocate(obspert_da(N_selected_obs, N_ens))
       allocate(Obs_cov(N_selected_obs, N_selected_obs))

       do kk = 1, N_local_obs
          jj = local_obs_idx(kk)
          obs_da(kk) = Observations_filtered(jj)
          obspred_da(kk,:) = Obs_pred_filtered(jj,:)
          obspert_da(kk,:) = Obs_pert_filtered(jj,:)

          do mm = 1, N_ens
             if(obspred_da(kk,mm) == LIS_rc%udef) obspred_flag = .false.
          enddo
       enddo

       ! Assemble observation covariance (block-diagonal by species)
       call assemble_obs_cov_multisource(LIS_rc%ndas, N_selected_obs, obs_da, Obs_cov)

       if(assim .and. obspred_flag) then
          ! Debug output for first few tiles
          if(tiles_with_obs <= 5) then
             write(LIS_logunit,*) '[DEBUG-SIMUL] Tile ', i, ': loc=', tile_lon, tile_lat
             write(LIS_logunit,*) '[DEBUG-SIMUL]   Local obs from ALL sources: ', N_local_obs
          endif

          ! Call LNETF analysis with combined observations
          call lnetf_analysis(i, N_state, N_selected_obs, N_ens, &
               obs_da, obspred_da, obspert_da, Obs_cov, &
               state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens), &
               (/tile_lon/), (/tile_lat/), xcompact, ycompact, &
               tile_n_eff, tile_mean_loc_weight, tile_max_loc_weight, tile_min_loc_weight)

          ! Store diagnostics
          lnetf_struc(n,k)%n_local_obs(i) = real(N_local_obs)
          lnetf_struc(n,k)%mean_loc_weight(i) = tile_mean_loc_weight
          lnetf_struc(n,k)%max_loc_weight(i) = tile_max_loc_weight
          lnetf_struc(n,k)%min_loc_weight(i) = tile_min_loc_weight
          lnetf_struc(n,k)%eff_sample_size_tile(i) = tile_n_eff

          ! Compute increment = X^a - X^f
          state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
               state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) - &
               state_tmp(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens)
       else
          state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0
          lnetf_struc(n,k)%n_local_obs(i) = real(N_local_obs)
          lnetf_struc(n,k)%mean_loc_weight(i) = 0.0
          lnetf_struc(n,k)%max_loc_weight(i) = 0.0
          lnetf_struc(n,k)%min_loc_weight(i) = 0.0
          lnetf_struc(n,k)%eff_sample_size_tile(i) = 0.0
       endif

       ! Store increment in structure
       lnetf_struc(n,k)%anlys_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens) = &
            state_incr(:, (i-1)*N_ens+1:(i-1)*N_ens+N_ens)

       deallocate(obs_da)
       deallocate(obspred_da)
       deallocate(obspert_da)
       deallocate(Obs_cov)
       deallocate(local_obs_idx)
    enddo

    !--------------------------------------------------------------------------
    ! Print analysis summary
    !--------------------------------------------------------------------------
    write(LIS_logunit,*) '[INFO] =============================================='
    write(LIS_logunit,*) '[INFO] TRUE SIMULTANEOUS LNETF Analysis Complete'
    write(LIS_logunit,*) '[INFO]   Tiles with observations: ', tiles_with_obs
    write(LIS_logunit,*) '[INFO]   Tiles without observations: ', tiles_no_obs
    if(tiles_with_obs > 0) then
       write(LIS_logunit,*) '[INFO]   Avg local obs per tile: ', &
            real(total_local_obs_count)/real(tiles_with_obs)
       write(LIS_logunit,*) '[INFO]   Min local obs: ', min_local_obs
       write(LIS_logunit,*) '[INFO]   Max local obs: ', max_local_obs
    endif
    write(LIS_logunit,*) '[INFO] =============================================='

    !--------------------------------------------------------------------------
    ! Cleanup
    !--------------------------------------------------------------------------
    deallocate(Observations_all)
    deallocate(Obs_pred_all)
    deallocate(Obs_pert_all)
    deallocate(Observations_filtered)
    deallocate(Obs_pred_filtered)
    deallocate(Obs_pert_filtered)
    deallocate(obs_lon_arr)
    deallocate(obs_lat_arr)
    deallocate(stvar)
    deallocate(state_incr)
    deallocate(state_tmp)
    deallocate(lats)
    deallocate(lons)
    deallocate(obs_param)

    write(LIS_logunit,*) '[INFO] Multi-source SIMULTANEOUS LNETF completed'
    write(LIS_logunit,*) '[INFO] =============================================='

  end subroutine lnetf_multisource_increments


!BOP
!
! !ROUTINE: assemble_obs_cov_multisource
! \label{assemble_obs_cov_multisource}
!
! !INTERFACE:
  subroutine assemble_obs_cov_multisource(ndas, N_obs, Observations, Obs_cov)
! !USES:
    use lnetf_types, only: obs_type

! !ARGUMENTS:
    integer, intent(in)                              :: ndas, N_obs
    type(obs_type), dimension(N_obs), intent(in)     :: Observations
    real, intent(out), dimension(N_obs,N_obs)        :: Obs_cov
!
! !DESCRIPTION:
!
! Assemble observation error covariance for multi-source simultaneous DA.
! Creates a block-diagonal R matrix where:
!   - Diagonal elements: variance from each observation's std
!   - Off-diagonal within same species: spatial correlation (optional)
!   - Off-diagonal between different species: ZERO (block-diagonal)
!
! For multi-source DA, different satellites (AMSR2, SMAP) have different
! species tags, ensuring no cross-correlation between observation sources.
!
! This simplified version uses diagonal R (no within-species spatial correlation)
! for computational efficiency in local analysis where observations are
! already localized spatially.
!
! \begin{description}
!  \item[ndas] number of DA instances (sources)
!  \item[N\_obs] number of observations
!  \item[Observations] Observations data type with species tag
!  \item[Obs\_cov] observation error covariance (output)
! \end{description}
!EOP

    integer :: i, j, i_species, j_species

    ! -------------------------------------------------------------
    ! Initialize R matrix to zero
    ! -------------------------------------------------------------
    Obs_cov = 0.0

    ! -------------------------------------------------------------
    ! Diagonal elements: variance = std^2
    ! This is the primary error term for each observation
    ! -------------------------------------------------------------
    do i = 1, N_obs
       Obs_cov(i,i) = Observations(i)%std**2
    end do

    ! -------------------------------------------------------------
    ! Off-diagonal elements for multi-source:
    ! For simplicity in local LNETF analysis, we use diagonal R
    ! (no spatial correlation within species). This is appropriate
    ! because:
    !   1. Observations are already localized to nearby tiles
    !   2. Spatial correlation is implicitly handled by localization
    !   3. Different species (satellites) are always uncorrelated
    !
    ! If within-species correlation is needed, uncomment the loop below:
    ! -------------------------------------------------------------
    !
    ! do i = 1, N_obs
    !    do j = (i+1), N_obs
    !       i_species = Observations(i)%species
    !       j_species = Observations(j)%species
    !
    !       ! Only correlate within same species (same satellite)
    !       if (i_species == j_species) then
    !          ! Compute spatial correlation factor
    !          ! fac = exp(-0.5 * distance^2 / correlation_length^2)
    !          ! Obs_cov(i,j) = Observations(i)%std * Observations(j)%std * fac
    !          ! Obs_cov(j,i) = Obs_cov(i,j)
    !       endif
    !    end do
    ! end do

  end subroutine assemble_obs_cov_multisource


!BOP
!
! !ROUTINE: lnetf_multisource_update
! \label{lnetf_multisource_update}
!
! !INTERFACE:
  subroutine lnetf_multisource_update(n)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
!
! !DESCRIPTION:
!  This routine applies the analysis increments computed by
!  lnetf_multisource_increments to the model state. Since all DA instances
!  operate on the same model state, we only apply updates from k=1.
!
!EOP
    integer :: k

    integer :: kk, status
    logical :: data_status_k

    write(LIS_logunit,*) '[INFO] Applying multi-source LNETF increments...'

    ! Apply increments from combined analysis (stored in k=1)
    k = 1
    call lnetf_update(n, k)

    ! Set Data Assimilate Status for ALL DA instances that had data
    ! This ensures diagnostics output works for all instances
    do kk = 1, LIS_rc%ndas
       call ESMF_AttributeGet(LIS_OBS_State(n,kk), name="Data Update Status", &
            value=data_status_k, rc=status)
       if(data_status_k) then
          call ESMF_AttributeSet(LIS_OBS_State(n,kk), "Data Assimilate Status", &
               .true., rc=status)
          write(LIS_logunit,*) '[INFO] Set Data Assimilate Status = true for DA instance ', kk
       endif
    enddo

    write(LIS_logunit,*) '[INFO] Multi-source LNETF update completed'

  end subroutine lnetf_multisource_update

!BOP
!
! !ROUTINE: lnetf_multisource_diagnostics
! \label{lnetf_multisource_diagnostics}
!
! !INTERFACE:
  subroutine lnetf_multisource_diagnostics(n)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n
!
! !DESCRIPTION:
!  This subroutine generates the LNETF related diagnostics for multi-source
!  simultaneous data assimilation and outputs to consolidated files.
!  For multi-source mode, this produces:
!    - Single increment file with variables for each observation source
!      (e.g., anlys_incr_SM_Layer1_obs01, anlys_incr_SM_Layer1_obs02)
!    - Single innovation file with variables for each observation source
!      (e.g., innov_obs01, innov_obs02, ninnov_obs01, ninnov_obs02)
!    - Single spread file with combined ensemble spread
!
!EOP
    integer                :: kk, k, v
    integer                :: status
    logical                :: assim_status, any_assim
    logical                :: alarmCheck
    character*3            :: fda
    character*2            :: finst

    ! Check if any DA instance had assimilation
    any_assim = .false.
    do kk = 1, LIS_rc%ndas
       call ESMF_AttributeGet(LIS_OBS_State(n,kk),name="Data Assimilate Status",&
            value=assim_status,rc=status)
       if(assim_status) any_assim = .true.
    enddo

    ! Check alarm for first DA instance (they should all be synchronized)
    write(fda,'(i3.3)') 1
    alarmCheck = LIS_isAlarmRinging(LIS_rc,"LIS DA output "//trim(fda))

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    if(any_assim) then
!--------------------------------------------------------------------------
! Write consolidated innovation file with variables for each observation source
! * Normalized innovations: ninnov_obs01, ninnov_obs02, ...
! * Raw innovations: innov_obs01, innov_obs02, ...
! * Analysis residuals: analysis_residual_obs01, ...
! * Forecast variance: forecast_sigma_obs01, ...
!--------------------------------------------------------------------------
       if(.not.lnetf_struc(n,1)%fileopen.and.LIS_masterproc) then
          call LIS_create_output_directory("LNETF")
       endif

       call writeMultisourceInnovation(n)
       call writeMultisourceAnalysisIncr(n)
    endif

    if(alarmCheck) then
       if(.not.lnetf_struc(n,1)%fileopen.and.LIS_masterproc) then
          call LIS_create_output_directory("LNETF")
       endif

       call writeMultisourceEnsembleSpread(n)
    endif

#endif
  end subroutine lnetf_multisource_diagnostics


!BOP
!
! !ROUTINE: writeMultisourceInnovation
! \label{writeMultisourceInnovation}
!
! !INTERFACE:
  subroutine writeMultisourceInnovation(n)
! !ARGUMENTS:
    integer,  intent(in)    :: n
!
! !DESCRIPTION:
!  This routine writes the innovation values for all observation sources
!  to a single consolidated NetCDF file. Each observation source gets
!  its own set of variables (e.g., innov_obs01, innov_obs02).
!
!EOP
    integer                :: ftn
    integer                :: kk
    character(len=LIS_CONST_PATH_LEN) :: innovfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3), ninnov_Id, innov_id
    integer                :: ares_Id, forecast_sigma_id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    character(len=10)       :: cdate
    character(len=14)       :: cdate1
    character(len=LIS_CONST_PATH_LEN)  :: dname
    integer                :: status
    logical                :: assim_status

    shuffle = 1
    deflate = 1
    deflate_level = 9

    ! Create multi-source innovation filename (without .a01, use .ms for multisource)
    if(LIS_masterproc) then
       write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn

       if(LIS_rc%wstyle.eq."4 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
          dname = trim(dname)//trim(cdate)//'/'
          write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
               LIS_rc%yr, LIS_rc%mo, LIS_rc%da
          dname = trim(dname)//trim(cdate)
          innovfile = trim(dname)//'/LIS_DA_LNETF_'//trim(cdate1)//&
               '_innov.ms.d01.nc'
       elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
          dname = trim(dname)//trim(cdate)//'/'
          innovfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_innov.ms.d01.nc'
       else
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          innovfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_innov.ms.d01.nc'
       endif

#if (defined USE_NETCDF4)
       status = nf90_create(path=innovfile,cmode=nf90_hdf5,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(innovfile)//&
            ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
       status = nf90_create(path=innovfile,cmode=nf90_clobber,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(innovfile)//&
            ' failed in lnetf_Mod')
#endif

       ! Use dimensions from first DA instance (assume same obs grid)
       call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%obs_gnc(1),&
            dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
       call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%obs_gnr(1),&
            dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')

       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att failed for missing_value in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"title", "Multi-source LNETF Innovations"),&
            'nf90_put_att failed for title in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"num_obs_sources", LIS_rc%ndas),&
            'nf90_put_att failed for num_obs_sources in lnetf_mod')

       ! Define variables for each observation source
       do kk = 1, LIS_rc%ndas
          write(unit=finst, fmt='(i2.2)') kk

          ! Check if this DA instance had data
          call ESMF_AttributeGet(LIS_OBS_State(n,kk),name="Data Assimilate Status",&
               value=assim_status,rc=status)

          if(LIS_rc%winnov(kk).eq.1 .and. assim_status) then

             ! Normalized innovations
             varname = "ninnov_obs"//trim(finst)
             standard_name = "Normalized_innovations_for_obs_source_"//trim(finst)

             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=ninnov_Id),&
                  'nf90_def_var failed for '//trim(varname)//' in lnetf_mod')
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  ninnov_Id,&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate failed for ninnov in lnetf_mod')
#endif
             call LIS_verify(nf90_put_att(ftn,ninnov_Id,&
                  "standard_name",standard_name),&
                  'nf90_put_att failed for ninnov in lnetf_mod')

             ! Raw innovations
             varname = "innov_obs"//trim(finst)
             standard_name = "Innovations_for_obs_source_"//trim(finst)

             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=innov_Id),&
                  'nf90_def_var failed for innov')
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  innov_Id,&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate for innov failed in lnetf_mod')
#endif
             call LIS_verify(nf90_put_att(ftn,innov_Id,&
                  "standard_name",standard_name),&
                  'nf90_put_att failed for innov in lnetf_mod')

             ! Analysis residuals
             varname = "analysis_residual_obs"//trim(finst)
             standard_name = "Analysis_residuals_for_obs_source_"//trim(finst)

             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=ares_Id),&
                  'nf90_def_var failed for '//trim(varname)//' in lnetf_mod')
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  ares_Id,&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate failed for analysis_residual in lnetf_mod')
#endif
             call LIS_verify(nf90_put_att(ftn,ares_Id,&
                  "standard_name",standard_name),&
                  'nf90_put_att failed for analysis_residual in lnetf_mod')

             ! Forecast variance
             varname = "forecast_sigma_obs"//trim(finst)
             standard_name = "Forecast_variance_for_obs_source_"//trim(finst)

             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=forecast_sigma_Id),&
                  'nf90_def_var for forecast_sigma failed in lnetf_mod')
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  forecast_sigma_Id,&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var for forecast_sigma failed in lnetf_mod')
#endif
             call LIS_verify(nf90_put_att(ftn,forecast_sigma_Id,&
                  "standard_name",standard_name),&
                  'nf90_put_att for forecast_sigma failed in lnetf_mod')
          endif
       enddo

       call LIS_verify(nf90_enddef(ftn),&
            'nf90_enddef failed in lnetf_mod')
    endif

    ! Write data for each observation source
    do kk = 1, LIS_rc%ndas
       call ESMF_AttributeGet(LIS_OBS_State(n,kk),name="Data Assimilate Status",&
            value=assim_status,rc=status)

       if(LIS_rc%winnov(kk).eq.1 .and. assim_status) then
          ! Get variable IDs and write data
          write(unit=finst, fmt='(i2.2)') kk

          if(LIS_masterproc) then
             status = nf90_inq_varid(ftn, "ninnov_obs"//trim(finst), ninnov_Id)
             status = nf90_inq_varid(ftn, "innov_obs"//trim(finst), innov_Id)
             status = nf90_inq_varid(ftn, "analysis_residual_obs"//trim(finst), ares_Id)
             status = nf90_inq_varid(ftn, "forecast_sigma_obs"//trim(finst), forecast_sigma_Id)
          endif

          call LIS_writevar_innov(ftn,n, kk, ninnov_id, &
               lnetf_struc(n,kk)%norm_innov)
          call LIS_writevar_innov(ftn,n, kk, innov_id, &
               lnetf_struc(n,kk)%innov)
          call LIS_writevar_innov(ftn,n, kk, ares_id, &
               lnetf_struc(n,kk)%anlys_res)
          call LIS_writevar_innov(ftn,n, kk, forecast_sigma_id, &
               lnetf_struc(n,kk)%forecast_var)
       endif
    enddo

    if(LIS_masterproc) then
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in lnetf_mod')
       write(LIS_logunit,*) '[INFO] Wrote multi-source innovation file: '//trim(innovfile)
    endif

  end subroutine writeMultisourceInnovation


!BOP
!
! !ROUTINE: writeMultisourceAnalysisIncr
! \label{writeMultisourceAnalysisIncr}
!
! !INTERFACE:
  subroutine writeMultisourceAnalysisIncr(n)
! !ARGUMENTS:
    integer,  intent(in)    :: n
!
! !DESCRIPTION:
!  This routine writes the analysis increments for multi-source DA
!  to a single consolidated NetCDF file. The increments are actually
!  the same for all sources (since it's simultaneous assimilation),
!  but we also write localization diagnostics with source-specific naming.
!
!EOP
    integer                :: ftn
    integer                :: kk, k, v
    character(len=LIS_CONST_PATH_LEN) :: incrfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: incr_id(20)  ! Max state variables
    integer                :: n_local_obs_id, mean_loc_weight_id
    integer                :: max_loc_weight_id, min_loc_weight_id
    integer                :: eff_sample_size_id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    character(len=10)       :: cdate
    character(len=14)       :: cdate1
    character(len=LIS_CONST_PATH_LEN)  :: dname
    integer                :: status
    real, allocatable      :: stvar(:,:)
    character*100, allocatable :: state_objs(:)

    ! Temporary arrays for expanding tile-level diagnostics to patch-level
    real, allocatable      :: diag_expanded(:)
    integer                :: ntiles, npatch, i, j, nensem
    logical                :: assim_status

    ! Use first DA instance for parameters
    k = 1

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Assimilate Status",&
         value=assim_status,rc=status)
    if(.not.assim_status) return

    if(LIS_rc%wensems(k).ne.1) return

    shuffle = 1
    deflate = 1
    deflate_level = 9

    ! Create multi-source increment filename
    if(LIS_masterproc) then
       write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn

       if(LIS_rc%wstyle.eq."4 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
          dname = trim(dname)//trim(cdate)//'/'
          write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
               LIS_rc%yr, LIS_rc%mo, LIS_rc%da
          dname = trim(dname)//trim(cdate)
          incrfile = trim(dname)//'/LIS_DA_LNETF_'//trim(cdate1)//&
               '_incr.ms.d01.nc'
       elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
          dname = trim(dname)//trim(cdate)//'/'
          incrfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_incr.ms.d01.nc'
       else
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          incrfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_incr.ms.d01.nc'
       endif

#if (defined USE_NETCDF4)
       status = nf90_create(path=incrfile,cmode=nf90_hdf5,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(incrfile)//&
            ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
       status = nf90_create(path=incrfile,cmode=nf90_clobber,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(incrfile)//&
            ' failed in lnetf_Mod')
#endif

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_dim(ftn,'ngrid',&
               LIS_rc%glbngrid_red(n),&
               dimID(1)),'nf90_def_dim for ngrid failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')
       endif

       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att for missing_value failed in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"title", "Multi-source LNETF Analysis Increments"),&
            'nf90_put_att failed for title in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"num_obs_sources", LIS_rc%ndas),&
            'nf90_put_att failed for num_obs_sources in lnetf_mod')

       !----------------------------------------------------------------------
       ! Analysis increments - use combined increments from DA instance 1
       ! (In simultaneous mode, all increments are computed together)
       !----------------------------------------------------------------------
       allocate(state_objs(LIS_rc%nstvars(k)))
       call LIS_surfaceModel_DAgetStateVarNames(n,k,state_objs)

       do v = 1, LIS_rc%nstvars(k)
          varname = "anlys_incr_"//trim(state_objs(v))//"_ms"
          standard_name = "Analysis_incr_multisource_"//trim(state_objs(v))

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=incr_Id(v)),&
                  'nf90_def_var for incr failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=incr_Id(v)),&
                  'nf90_def_var for incr failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               incr_Id(v),&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for incr failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,incr_Id(v),&
               "standard_name",standard_name),&
               'nf90_put_att for incr failed in lnetf_mod')
       end do
       deallocate(state_objs)

       !----------------------------------------------------------------------
       ! Localization diagnostics - n_local_obs (combined from all sources)
       !----------------------------------------------------------------------
       varname = "n_local_obs_ms"
       standard_name = "Number_of_local_observations_multisource"

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1), varID=n_local_obs_Id),&
               'nf90_def_var for n_local_obs failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=n_local_obs_Id),&
               'nf90_def_var for n_local_obs failed in lnetf_mod')
       endif
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,&
            n_local_obs_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate for n_local_obs failed in lnetf_mod')
#endif
       call LIS_verify(nf90_put_att(ftn,n_local_obs_Id,&
            "standard_name",standard_name),&
            'nf90_put_att for n_local_obs failed in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,n_local_obs_Id,&
            "long_name","Combined observations from all sources within localization radius"),&
            'nf90_put_att for n_local_obs long_name failed in lnetf_mod')

       !----------------------------------------------------------------------
       ! mean_loc_weight
       !----------------------------------------------------------------------
       varname = "mean_loc_weight_ms"
       standard_name = "Mean_localization_weight_multisource"

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1), varID=mean_loc_weight_Id),&
               'nf90_def_var for mean_loc_weight failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=mean_loc_weight_Id),&
               'nf90_def_var for mean_loc_weight failed in lnetf_mod')
       endif
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,&
            mean_loc_weight_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate for mean_loc_weight failed in lnetf_mod')
#endif
       call LIS_verify(nf90_put_att(ftn,mean_loc_weight_Id,&
            "standard_name",standard_name),&
            'nf90_put_att for mean_loc_weight failed in lnetf_mod')

       !----------------------------------------------------------------------
       ! max_loc_weight
       !----------------------------------------------------------------------
       varname = "max_loc_weight_ms"
       standard_name = "Max_localization_weight_multisource"

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1), varID=max_loc_weight_Id),&
               'nf90_def_var for max_loc_weight failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=max_loc_weight_Id),&
               'nf90_def_var for max_loc_weight failed in lnetf_mod')
       endif
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,&
            max_loc_weight_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate for max_loc_weight failed in lnetf_mod')
#endif
       call LIS_verify(nf90_put_att(ftn,max_loc_weight_Id,&
            "standard_name",standard_name),&
            'nf90_put_att for max_loc_weight failed in lnetf_mod')

       !----------------------------------------------------------------------
       ! min_loc_weight
       !----------------------------------------------------------------------
       varname = "min_loc_weight_ms"
       standard_name = "Min_localization_weight_multisource"

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1), varID=min_loc_weight_Id),&
               'nf90_def_var for min_loc_weight failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=min_loc_weight_Id),&
               'nf90_def_var for min_loc_weight failed in lnetf_mod')
       endif
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,&
            min_loc_weight_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate for min_loc_weight failed in lnetf_mod')
#endif
       call LIS_verify(nf90_put_att(ftn,min_loc_weight_Id,&
            "standard_name",standard_name),&
            'nf90_put_att for min_loc_weight failed in lnetf_mod')

       !----------------------------------------------------------------------
       ! eff_sample_size
       !----------------------------------------------------------------------
       varname = "eff_sample_size_ms"
       standard_name = "Effective_sample_size_multisource"

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1), varID=eff_sample_size_Id),&
               'nf90_def_var for eff_sample_size failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=eff_sample_size_Id),&
               'nf90_def_var for eff_sample_size failed in lnetf_mod')
       endif
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_deflate(ftn,&
            eff_sample_size_Id,&
            shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate for eff_sample_size failed in lnetf_mod')
#endif
       call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
            "standard_name",standard_name),&
            'nf90_put_att for eff_sample_size failed in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,eff_sample_size_Id,&
            "long_name","Effective sample size N_eff = 1/sum(w_i^2)"),&
            'nf90_put_att for eff_sample_size long_name failed in lnetf_mod')

       call LIS_verify(nf90_enddef(ftn),&
            'nf90_enddef failed in lnetf_mod')
    endif

    ! Write analysis increments (use DA instance 1 - same for all in simultaneous mode)
    do v=1,LIS_rc%nstvars(k)
       call LIS_writevar_incr(ftn,n,k,incr_id(v), &
            lnetf_struc(n,k)%anlys_incr(v,:),v)
    enddo

    !----------------------------------------------------------------------
    ! Write localization diagnostic variables
    !----------------------------------------------------------------------
    npatch = LIS_surfaceModel_DAgetStateSpaceSize(n,k)
    nensem = LIS_rc%nensem(n)
    ntiles = npatch / nensem

    allocate(diag_expanded(npatch))

    ! Expand n_local_obs
    do i = 1, ntiles
       do j = 1, nensem
          diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%n_local_obs(i)
       enddo
    enddo
    call LIS_writevar_incr(ftn,n,k,n_local_obs_id, diag_expanded, 1)

    ! Expand mean_loc_weight
    do i = 1, ntiles
       do j = 1, nensem
          diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%mean_loc_weight(i)
       enddo
    enddo
    call LIS_writevar_incr(ftn,n,k,mean_loc_weight_id, diag_expanded, 1)

    ! Expand max_loc_weight
    do i = 1, ntiles
       do j = 1, nensem
          diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%max_loc_weight(i)
       enddo
    enddo
    call LIS_writevar_incr(ftn,n,k,max_loc_weight_id, diag_expanded, 1)

    ! Expand min_loc_weight
    do i = 1, ntiles
       do j = 1, nensem
          diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%min_loc_weight(i)
       enddo
    enddo
    call LIS_writevar_incr(ftn,n,k,min_loc_weight_id, diag_expanded, 1)

    ! Expand eff_sample_size_tile
    do i = 1, ntiles
       do j = 1, nensem
          diag_expanded((i-1)*nensem + j) = lnetf_struc(n,k)%eff_sample_size_tile(i)
       enddo
    enddo
    call LIS_writevar_incr(ftn,n,k,eff_sample_size_id, diag_expanded, 1)

    deallocate(diag_expanded)

    if(LIS_masterproc) then
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in lnetf_mod')
       write(LIS_logunit,*) '[INFO] Wrote multi-source increment file: '//trim(incrfile)
    endif

  end subroutine writeMultisourceAnalysisIncr


!BOP
!
! !ROUTINE: writeMultisourceEnsembleSpread
! \label{writeMultisourceEnsembleSpread}
!
! !INTERFACE:
  subroutine writeMultisourceEnsembleSpread(n)
! !ARGUMENTS:
    integer,  intent(in)    :: n
!
! !DESCRIPTION:
!  This routine writes the ensemble spread (standard deviation)
!  for multi-source simultaneous DA to a single consolidated file.
!
!EOP
    integer                :: ftn
    integer                :: kk, k, v
    character(len=LIS_CONST_PATH_LEN) :: spreadfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: ensspread_id(20)  ! Max state variables
    character*100          :: varname, vardimname, standard_name
    character(len=10)       :: cdate
    character(len=14)       :: cdate1
    character(len=LIS_CONST_PATH_LEN)  :: dname
    integer                :: status
    integer                :: state_size
    real, allocatable      :: stvar(:,:)
    character*100, allocatable :: state_objs(:)
    logical                :: any_alarm

    ! Use first DA instance for parameters
    k = 1

    if(LIS_rc%wensems(k).ne.1) return

    shuffle = 1
    deflate = 1
    deflate_level = 9

    state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)

    ! Create multi-source spread filename
    if(LIS_masterproc) then
       write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn

       if(LIS_rc%wstyle.eq."4 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
          dname = trim(dname)//trim(cdate)//'/'
          write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
               LIS_rc%yr, LIS_rc%mo, LIS_rc%da
          dname = trim(dname)//trim(cdate)
          spreadfile = trim(dname)//'/LIS_DA_LNETF_'//trim(cdate1)//&
               '_spread.ms.d01.nc'
       elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
          dname = trim(dname)//trim(cdate)//'/'
          spreadfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_spread.ms.d01.nc'
       else
          dname = trim(LIS_rc%odir)//'/'
          dname = trim(dname)//'LNETF/'
          spreadfile = trim(dname)//'LIS_DA_LNETF_'//trim(cdate1)//&
               '_spread.ms.d01.nc'
       endif

#if (defined USE_NETCDF4)
       status = nf90_create(path=spreadfile,cmode=nf90_hdf5,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(spreadfile)//&
            ' failed in lnetf_Mod')
#endif
#if (defined USE_NETCDF3)
       status = nf90_create(path=spreadfile,cmode=nf90_clobber,&
            ncid = ftn)
       call LIS_verify(status,&
            'creating netcdf file '//trim(spreadfile)//&
            ' failed in lnetf_Mod')
#endif

       if(LIS_rc%wopt.eq."1d gridspace") then
          call LIS_verify(nf90_def_dim(ftn,'ngrid',&
               LIS_rc%glbngrid_red(n),&
               dimID(1)),'nf90_def_dim for ngrid failed in lnetf_mod')
       elseif(LIS_rc%wopt.eq."2d gridspace") then
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),'nf90_def_dim for east_west failed in lnetf_mod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),'nf90_def_dim for north_south failed in lnetf_mod')
       endif

       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att for missing_value failed in lnetf_mod')
       call LIS_verify(nf90_put_att(ftn,&
            NF90_GLOBAL,"title", "Multi-source LNETF Ensemble Spread"),&
            'nf90_put_att failed for title in lnetf_mod')

       !----------------------------------------------------------------------
       ! Ensemble spread - meta data
       !----------------------------------------------------------------------
       allocate(state_objs(LIS_rc%nstvars(k)))
       call LIS_surfaceModel_DAgetStateVarNames(n,k,state_objs)

       do v = 1, LIS_rc%nstvars(k)
          varname = "ensspread_"//trim(state_objs(v))//"_ms"
          standard_name = "Ensemble_spread_multisource_"//trim(state_objs(v))

          if(LIS_rc%wopt.eq."1d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1), varID=ensspread_Id(v)),&
                  'nf90_def_var for ensspread failed in lnetf_mod')

          elseif(LIS_rc%wopt.eq."2d gridspace") then
             call LIS_verify(nf90_def_var(ftn,varname,&
                  nf90_float,&
                  dimids = dimID(1:2), varID=ensspread_Id(v)),&
                  'nf90_def_var for ensspread failed in lnetf_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               ensspread_Id(v),&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for ensspread failed in lnetf_mod')
#endif
          call LIS_verify(nf90_put_att(ftn,ensspread_Id(v),&
               "standard_name",standard_name),&
               'nf90_put_att for ensspread failed in lnetf_mod')
          call LIS_verify(nf90_enddef(ftn),&
               'nf90_enddef failed in lnetf_mod')
       end do
       deallocate(state_objs)
    endif

    allocate(stvar(LIS_rc%nstvars(k), state_size))

    call LIS_surfaceModel_DAGetStateVar(n,k)

    call LIS_surfaceModel_DAextractStateVector(n,k,&
         LIS_rc%nstvars(k),state_size, stvar)

    do v=1,LIS_rc%nstvars(k)
       call LIS_writevar_spread(ftn,n,k,ensspread_id(v), &
            stvar(v,:),v)
    enddo

    deallocate(stvar)

    if(LIS_masterproc) then
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in lnetf_mod')
       write(LIS_logunit,*) '[INFO] Wrote multi-source spread file: '//trim(spreadfile)
    endif

  end subroutine writeMultisourceEnsembleSpread

end module lnetf_Mod
