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
! !MODULE: letkf_types
!
! !DESCRIPTION:
!   This module contains type definitions for the LETKF algorithm.
!   The types are identical to those used in EnKF and LNETF to maintain
!   consistency across filter implementations.
!
! !REVISION HISTORY:
!   28 Jan 2026: Initial implementation (copied from enkf_types.F90)
!EOP
module letkf_types

  implicit none

  private

  public :: obs_type
  public :: obs_param_type
  public :: update_region_type

  ! --------------------------------------------------------------------------
  ! Observation type: holds individual observation data
  ! --------------------------------------------------------------------------
  type :: obs_type
     integer :: species        ! Type/source of observation (e.g., SMAP=1, ASCAT=2)
     real    :: lon            ! Longitude of observation [degrees]
     real    :: lat            ! Latitude of observation [degrees]
     real    :: value          ! Observation value
     real    :: std            ! Observation error standard deviation
     logical :: assim          ! Flag: .true. if observation should be assimilated
     integer :: pert_type      ! Perturbation type: 0=additive, 1=multiplicative
  end type obs_type

  ! --------------------------------------------------------------------------
  ! Observation parameter type: holds configuration for each obs type
  ! --------------------------------------------------------------------------
  type :: obs_param_type
     integer :: species        ! Observation type identifier
     character(LEN=100) :: path ! Path to observation data
     real    :: std            ! Default error standard deviation
     real    :: cross_corr     ! Cross-correlation parameter
     real    :: corr_len       ! Correlation length [degrees or km]
  end type obs_param_type

  ! --------------------------------------------------------------------------
  ! Update region type: geographic bounds for DA update
  ! --------------------------------------------------------------------------
  type :: update_region_type
     real :: lat_min           ! Minimum latitude [degrees]
     real :: lat_max           ! Maximum latitude [degrees]
     real :: lon_min           ! Minimum longitude [degrees]
     real :: lon_max           ! Maximum longitude [degrees]
  end type update_region_type

end module letkf_types
