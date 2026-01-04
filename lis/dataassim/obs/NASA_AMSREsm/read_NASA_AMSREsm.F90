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
!BOP
! !ROUTINE: read_NASA_AMSREsm
! \label{read_NASA_AMSREsm}
!
! !REVISION HISTORY:
!  06Nov2007: Bailing Li; Initial Specification
!  17Jun2010: Sujay Kumar; Updated for use with NASA AMSRE Version 6,
!                          added LSM-based QC, generic handling of obs scaling
!  22Dec2025: Gyuyeon Choi; Modified to read AMSR2 L2 HDF5 swath data
!
! !INTERFACE:
subroutine read_NASA_AMSREsm(n, k, OBS_State, OBS_Pert_State)
! !USES:
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod,    only : LIS_rc, LIS_domain, &
       LIS_masterproc, LIS_npes, LIS_masterproc, LIS_localPet
  use LIS_DAobservationsMod, only : LIS_obs_domain, LIS_checkForValidObs
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock, LIS_tick
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use NASA_AMSREsm_Mod, only : NASA_AMSREsm_struc
  use LIS_pluginIndices
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!
!  reads the AMSR2 L2 soil moisture observations
!  from HDF5 files and packages it
!  into an ESMF State with certain predefined
!  attributes
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter     :: MAX_SM_VALUE=0.55
  type(ESMF_Field)    :: smField

  integer             :: iret

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  real                :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                :: sm_current(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))
  real                :: obs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  character(len=LIS_CONST_PATH_LEN) :: name_hourly
  character(len=LIS_CONST_PATH_LEN) :: name_prev
  character(len=LIS_CONST_PATH_LEN) :: amsr2_filename(100)
  character(len=2)  :: fhr
  character(len=2)  :: fhr_prev
  character(len=8)  :: prev_date
  integer           :: prev_hr, prev_yr, prev_mo, prev_da

  logical             :: alarmCheck

  logical             :: readflag
  integer             :: status
  integer             :: valid_gindex_count, smobs_idx, first_data_found, stored_count, skipped_ocean
  integer             :: slash_pos, pattern_start, validObsFlag_temp
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  integer             :: t,c,r,p,jj,igd

  ! Time variables for SMAP-like time window filtering
  real*8              :: timenow, time1, time3
  integer             :: yr, mo, da, hr, mn, ss, doy
  real                :: gmt
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer             :: fname_idx

  integer             :: col, row, i, gridid
  integer             :: binval
  real                :: cdf_obsval
  real                :: smvalue
  real                :: model_delta(LIS_rc%obs_ngrid(k))
  real                :: obs_delta(LIS_rc%obs_ngrid(k))

  integer             :: fnd
  integer             :: ierr
  integer             :: ftn, nfiles

  ! write(LIS_logunit,*) '[DEBUG] ========================================='
  ! write(LIS_logunit,*) '[DEBUG] read_NASA_AMSREsm CALLED!'
  ! write(LIS_logunit,*) '[DEBUG] Current LIS time: ', LIS_rc%yr, '/', LIS_rc%mo, '/', LIS_rc%da, ' ', LIS_rc%hr, ':', LIS_rc%mn
  ! write(LIS_logunit,*) '[DEBUG] ========================================='
  ! call flush(LIS_logunit)

  smobs    = LIS_rc%udef
  obs_unsc = LIS_rc%udef

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)

  ! write(LIS_logunit,*) '[DEBUG] Data directory: ', trim(smobsdir)
  ! call flush(LIS_logunit)

  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)
  data_upd = .false.
!-------------------------------------------------------------------------
!   Read the data at 0z daily (following SMAP/AMSR-E pattern)
!-------------------------------------------------------------------------

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "NASA AMSR-E read alarm")

  ! write(LIS_logunit,*) '[DEBUG] NASA_AMSREsm: alarmCheck = ', alarmCheck
  ! write(LIS_logunit,*) '[DEBUG] NASA_AMSREsm: startMode = ', NASA_AMSREsm_struc(n)%startMode
  ! write(LIS_logunit,*) '[DEBUG] Alarm name check complete'
  ! call flush(LIS_logunit)

  ! Read data if alarm is ringing OR on first call (startMode)
  ! Following SMAP L2 pattern
  if(alarmCheck .or. NASA_AMSREsm_struc(n)%startMode) then
     NASA_AMSREsm_struc(n)%startMode = .false.

     ! write(LIS_logunit,*) '[DEBUG] ========================================='
     ! write(LIS_logunit,*) '[DEBUG] Reading AMSR2 L2 data (alarm or startMode)'
     ! write(LIS_logunit,*) '[DEBUG] obs_lnc(k) = ', LIS_rc%obs_lnc(k)
     ! write(LIS_logunit,*) '[DEBUG] obs_lnr(k) = ', LIS_rc%obs_lnr(k)
     ! write(LIS_logunit,*) '[DEBUG] obs_ngrid(k) = ', LIS_rc%obs_ngrid(k)
     ! write(LIS_logunit,*) '[DEBUG] ========================================='
     ! call flush(LIS_logunit)

     ! Initialize observation arrays to udef (critical for swath data)
     NASA_AMSREsm_struc(n)%smobs = LIS_rc%udef
     NASA_AMSREsm_struc(n)%smtime = 0.0
     NASA_AMSREsm_struc(n)%smqc = 0

     ! AMSR2 L2: Read files for CURRENT and PREVIOUS hour
     ! Time window [T-30min, T) may need observations from previous hour
     ! Example: 01:00 DA needs [00:30, 01:00) -> 00:52 file is from hour 00
     call NASA_AMSREsm_filename(name,smobsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr)

     ! Find "GW1AM2_" position to build hourly pattern
     pattern_start = index(trim(name), 'GW1AM2_')

     ! Calculate previous hour (handle midnight crossing)
     if(LIS_rc%hr .eq. 0) then
        prev_hr = 23
        ! Need previous day - calculate manually
        if(LIS_rc%da .eq. 1) then
           ! First day of month - need previous month
           if(LIS_rc%mo .eq. 1) then
              prev_yr = LIS_rc%yr - 1
              prev_mo = 12
              prev_da = 31  ! December always has 31 days
           else
              prev_yr = LIS_rc%yr
              prev_mo = LIS_rc%mo - 1
              ! Days in previous month
              select case(prev_mo)
              case(1,3,5,7,8,10,12)
                 prev_da = 31
              case(4,6,9,11)
                 prev_da = 30
              case(2)
                 if(mod(prev_yr,4).eq.0 .and. (mod(prev_yr,100).ne.0 .or. mod(prev_yr,400).eq.0)) then
                    prev_da = 29  ! Leap year
                 else
                    prev_da = 28
                 endif
              end select
           endif
        else
           prev_yr = LIS_rc%yr
           prev_mo = LIS_rc%mo
           prev_da = LIS_rc%da - 1
        endif
     else
        prev_hr = LIS_rc%hr - 1
        prev_yr = LIS_rc%yr
        prev_mo = LIS_rc%mo
        prev_da = LIS_rc%da
     endif

     ! Build pattern for current hour
     write(fhr,'(i2.2)') LIS_rc%hr
     name_hourly = trim(name(1:pattern_start+6)) // &  ! up to "GW1AM2_"
                   trim(name(pattern_start+7:pattern_start+14)) // &  ! YYYYMMDD
                   trim(fhr) // '*_*_L2SGSMCLF*.h5'

     ! Build pattern for previous hour
     write(fhr_prev,'(i2.2)') prev_hr
     write(prev_date,'(i4.4,i2.2,i2.2)') prev_yr, prev_mo, prev_da
     name_prev = trim(smobsdir)//'/'//prev_date(1:4)//'/'//prev_date(5:6)// &
                 '/GW1AM2_'//trim(prev_date)//trim(fhr_prev)//'*_*_L2SGSMCLF*.h5'

     write(LIS_logunit,*) '[INFO] Searching for AMSR2 L2 files (current hour): ', trim(name_hourly)
     write(LIS_logunit,*) '[INFO] Searching for AMSR2 L2 files (prev hour): ', trim(name_prev)
     call flush(LIS_logunit)

     ! Create list of files matching BOTH patterns (current + previous hour)
     call system('ls '//trim(name_prev)//' '//trim(name_hourly)//' > AMSR2_filelist_lis.dat 2>/dev/null')

     ftn = LIS_getNextUnitNumber()
     open(ftn, file="./AMSR2_filelist_lis.dat", status='old', iostat=ierr)

     if(ierr .eq. 0) then
        readflag = .true.
        nfiles = 0
        do i = 1, 100
           read(ftn, '(a)', iostat=ierr) amsr2_filename(i)
           if(ierr .ne. 0) exit
           nfiles = nfiles + 1
           write(LIS_logunit,*) '[INFO] Found AMSR2 L2 file: ', trim(amsr2_filename(i))
           call flush(LIS_logunit)
        enddo
        close(ftn)
        call LIS_releaseUnitNumber(ftn)

        if(nfiles .gt. 0) then
           write(LIS_logunit,*) '[INFO] Found ', nfiles, ' AMSR2 L2 files'
           call flush(LIS_logunit)

           ! Calculate time window BEFORE reading files
           ! For DA at time T, use observations from [T-timestep, T)
           yr = LIS_rc%yr
           mo = LIS_rc%mo
           da = LIS_rc%da
           hr = LIS_rc%hr
           mn = LIS_rc%mn
           ss = LIS_rc%ss
           call LIS_tick(time1, doy, gmt, yr, mo, da, hr, mn, ss, -LIS_rc%ts)

           yr = LIS_rc%yr
           mo = LIS_rc%mo
           da = LIS_rc%da
           hr = LIS_rc%hr
           mn = LIS_rc%mn
           ss = LIS_rc%ss
           call LIS_tick(time3, doy, gmt, yr, mo, da, hr, mn, ss, 0.0)

           write(LIS_logunit,*) '[DEBUG-TIME] Pre-filter time window: [', time1, ',', time3, ')'
           call flush(LIS_logunit)

           ! Read only files within time window
           do i = 1, nfiles
              ! Extract filename from full path
              fname_idx = index(amsr2_filename(i), '/', back=.true.)
              if(fname_idx .gt. 0) then
                 filename = amsr2_filename(i)(fname_idx+1:len_trim(amsr2_filename(i)))
              else
                 filename = amsr2_filename(i)
              endif

              ! Parse timestamp from filename: GW1AM2_YYYYMMDDHHMM_...
              ! Example: GW1AM2_201901312344_171A_L2SGSMCLF3300300.h5
              ! Extract FULL timestamp (year, month, day, hour, minute) from filename
              read(filename(8:11), '(i4)') yr
              read(filename(12:13), '(i2)') mo
              read(filename(14:15), '(i2)') da
              read(filename(16:17), '(i2)') hr
              read(filename(18:19), '(i2)') mn
              ss = 0

              ! Calculate timenow using ACTUAL observation time from filename
              call LIS_tick(timenow, doy, gmt, yr, mo, da, hr, mn, ss, 0.0)

              ! CHECK: Skip files outside time window [time1, time3)
              ! This prevents earlier files from overwriting later valid observations
              if(timenow .lt. time1 .or. timenow .ge. time3) then
                 write(LIS_logunit,*) '[INFO] SKIPPING file outside time window: ', trim(filename)
                 write(LIS_logunit,*) '[INFO]   File time: ', yr,'/',mo,'/',da,' ',hr,':',mn, &
                      ' (timenow=', timenow, ')'
                 call flush(LIS_logunit)
                 cycle  ! Skip this file
              endif

              write(LIS_logunit,*) '[DEBUG-TIME] Reading file: ', trim(filename)
              write(LIS_logunit,*) '[DEBUG-TIME] Observation time from filename: ', &
                   yr, '/', mo, '/', da, ' ', hr, ':', mn
              write(LIS_logunit,*) '[DEBUG-TIME] timenow = ', timenow, ' doy=', doy, ' gmt=', gmt
              call flush(LIS_logunit)

              ! Pass timenow to read_AMSREsm (like SMAP does)
              call read_AMSREsm(n, k, amsr2_filename(i), timenow)
           enddo

           write(LIS_logunit,*) '[INFO] Total observations after reading: ', &
                count(NASA_AMSREsm_struc(n)%smobs(:,1).ne.LIS_rc%udef)
           write(LIS_logunit,*) '[DEBUG-TIME] Total smtime non-zero: ', &
                count(NASA_AMSREsm_struc(n)%smtime(:,1).ne.0.0)
           call flush(LIS_logunit)

           ! QC masking after reading all files
           call maskObs_basedonQC(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
                NASA_AMSREsm_struc(n)%smobs, &
                NASA_AMSREsm_struc(n)%smqc)
        else
           readflag = .false.
           write(LIS_logunit,*) '[WARN] No AMSR2 L2 files found for: ', trim(name)
           call flush(LIS_logunit)
        endif
     else
        readflag = .false.
        write(LIS_logunit,*) '[WARN] No AMSR2 L2 files found for: ', trim(name)
        call flush(LIS_logunit)
     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  fnd = 0
  sm_current = LIS_rc%udef

  ! Calculate time window for observation filtering
  ! For DA at time T, use observations from [T-timestep, T)
  ! This ensures observations BEFORE the analysis time are used
  ! time1 = current LIS time - timestep (window start)
  ! time3 = current LIS time (window end)
  yr = LIS_rc%yr
  mo = LIS_rc%mo
  da = LIS_rc%da
  hr = LIS_rc%hr
  mn = LIS_rc%mn
  ss = LIS_rc%ss
  call LIS_tick(time1, doy, gmt, yr, mo, da, hr, mn, ss, -LIS_rc%ts)

  yr = LIS_rc%yr
  mo = LIS_rc%mo
  da = LIS_rc%da
  hr = LIS_rc%hr
  mn = LIS_rc%mn
  ss = LIS_rc%ss
  call LIS_tick(time3, doy, gmt, yr, mo, da, hr, mn, ss, 0.0)

  write(LIS_logunit,*) '[DEBUG-TIME] Time window calculation:'
  write(LIS_logunit,*) '[DEBUG-TIME]   LIS_rc time: ', &
       LIS_rc%yr,'/',LIS_rc%mo,'/',LIS_rc%da,' ',LIS_rc%hr,':',LIS_rc%mn,':',LIS_rc%ss
  write(LIS_logunit,*) '[DEBUG-TIME]   time1 = ', time1
  write(LIS_logunit,*) '[DEBUG-TIME]   time3 = ', time3
  write(LIS_logunit,*) '[DEBUG-TIME]   timestep (LIS_rc%ts) = ', LIS_rc%ts
  call flush(LIS_logunit)

  ! Apply time window filtering
  call maskObs_basedonTime(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASA_AMSREsm_struc(n)%smobs, &
       NASA_AMSREsm_struc(n)%smtime, smobs, fnd, time1, time3)

  ! Map smobs to sm_current and obs_unsc
  ! write(LIS_logunit,*) '[DEBUG] Mapping smobs to sm_current...'
  ! write(LIS_logunit,*) '[DEBUG] smobs array size: ', size(smobs)
  ! write(LIS_logunit,*) '[DEBUG] smobs valid count: ', count(smobs.ne.LIS_rc%udef)
  ! write(LIS_logunit,*) '[DEBUG] sm_current shape: ', shape(sm_current)
  ! call flush(LIS_logunit)

  valid_gindex_count = 0
  first_data_found = 0

  ! Check gindex values in the region where we have data (around c=273-281, r=117-121)
  ! write(LIS_logunit,*) '[DEBUG] Checking gindex in data region (c=270-285, r=115-125):'
  ! do r = 115, 125
  !    do c = 270, 285
  !       if(c .ge. 1 .and. c .le. LIS_rc%obs_lnc(k) .and. &
  !          r .ge. 1 .and. r .le. LIS_rc%obs_lnr(k)) then
  !          smobs_idx = c+LIS_rc%obs_lnc(k)*(r-1)
  !          if(smobs(smobs_idx).ne.LIS_rc%udef) then
  !             write(LIS_logunit,'(A,I4,A,I4,A,I8,A,I8,A,F12.6)') &
  !                  '  c=', c, ' r=', r, ' idx=', smobs_idx, &
  !                  ' gindex=', LIS_obs_domain(n,k)%gindex(c,r), ' SM=', smobs(smobs_idx)
  !          endif
  !       endif
  !    enddo
  ! enddo
  ! call flush(LIS_logunit)

  do r =1,LIS_rc%obs_lnr(k)
     do c =1,LIS_rc%obs_lnc(k)
        smobs_idx = c+LIS_rc%obs_lnc(k)*(r-1)

        ! Debug: check first few points with data
        ! if(smobs(smobs_idx).ne.LIS_rc%udef .and. first_data_found .lt. 5) then
        !    write(LIS_logunit,*) '[DEBUG] Data at (', c, ',', r, ') idx=', smobs_idx, &
        !         ' gindex=', LIS_obs_domain(n,k)%gindex(c,r), ' SM=', smobs(smobs_idx)
        !    first_data_found = first_data_found + 1
        ! endif

        if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1) then
           valid_gindex_count = valid_gindex_count + 1
           sm_current(c,r) = smobs(smobs_idx)
           if(sm_current(c,r).ne.LIS_rc%udef) then
              obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = sm_current(c,r)
           endif
        else
           ! If data exists but gindex is -1, still copy to sm_current
           ! This is the FIX for the issue
           if(smobs(smobs_idx).ne.LIS_rc%udef) then
              sm_current(c,r) = smobs(smobs_idx)
           endif
        endif
     enddo
  enddo

  ! write(LIS_logunit,*) '[DEBUG] Valid gindex count: ', valid_gindex_count
  ! call flush(LIS_logunit)

  ! write(LIS_logunit,*) '[DEBUG] sm_current valid count: ', count(sm_current.ne.LIS_rc%udef)
  ! write(LIS_logunit,*) '[DEBUG] obs_unsc valid count: ', count(obs_unsc.ne.LIS_rc%udef)
  ! call flush(LIS_logunit)

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!  TEMPORARILY DISABLED: CDF files are not compatible with AMSR2 L2 data
!-------------------------------------------------------------------------
  ! write(LIS_logunit,*) '[DEBUG] Before CDF scaling...'
  ! write(LIS_logunit,*) '[DEBUG] scal = ', NASA_AMSREsm_struc(n)%scal
  ! write(LIS_logunit,*) '[DEBUG] fnd = ', fnd
  ! call flush(LIS_logunit)

  if(NASA_AMSREsm_struc(n)%scal.ne.0.and.fnd.ne.0) then  ! CDF scaling enabled
     ! write(LIS_logunit,*) '[DEBUG] Entering CDF scaling...'
     ! call flush(LIS_logunit)

     do t=1,LIS_rc%obs_ngrid(k)
        model_delta(t) = NASA_AMSREsm_struc(n)%model_xrange(t,2,1)-&
             NASA_AMSREsm_struc(n)%model_xrange(t,1,1)
        obs_delta(t) = NASA_AMSREsm_struc(n)%obs_xrange(t,2,1)-&
             NASA_AMSREsm_struc(n)%obs_xrange(t,1,1)
        if(model_delta(t).eq.0) then
           model_delta(t) = 1.0
        endif
        if(obs_delta(t).eq.0) then
           obs_delta(t) = 1.0
        endif
     enddo

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              ! Add bounds check for CDF array access
              igd = LIS_obs_domain(n,k)%gindex(c,r)
              if(igd .lt. 1 .or. igd .gt. LIS_rc%obs_ngrid(k)) then
                 ! Index out of bounds, skip this point
                 sm_current(c,r) = LIS_rc%udef
                 cycle
              endif

              smvalue = sm_current(c,r)
              binval = 0

              if((smvalue.ne.LIS_rc%udef).and.&
                   (smvalue.ne.9999.0)) then
                 ! Add safety checks for CDF values to avoid floating point exceptions
                 if(ieee_is_nan(NASA_AMSREsm_struc(n)%obs_xrange(igd,1,1)) .or. &
                    ieee_is_nan(NASA_AMSREsm_struc(n)%obs_xrange(igd,2,1))) then
                    sm_current(c,r) = LIS_rc%udef
                    cycle
                 endif

                 smvalue = min(smvalue, NASA_AMSREsm_struc(n)%obs_xrange(igd,2,1))
                 smvalue = max(smvalue, NASA_AMSREsm_struc(n)%obs_xrange(igd,1,1))

                 do i=1,NASA_AMSREsm_struc(n)%nbins
                    if(ieee_is_nan(NASA_AMSREsm_struc(n)%obs_cdf(igd,i,1)) .or. &
                       ieee_is_nan(NASA_AMSREsm_struc(n)%obs_cdf(igd,i+1,1))) then
                       exit  ! Skip if CDF has NaN values
                    endif

                    if(smvalue.ge.NASA_AMSREsm_struc(n)%obs_cdf(igd,i,1).and.&
                         (smvalue.lt.NASA_AMSREsm_struc(n)%obs_cdf(igd,i+1,1))) then
                       binval = i
                       exit
                    endif
                 enddo

                 if(binval.ne.0) then
                    cdf_obsval = NASA_AMSREsm_struc(n)%model_cdf(igd,binval,1)
                    if(.not. ieee_is_nan(cdf_obsval)) then
                       sm_current(c,r) = cdf_obsval
                    else
                       sm_current(c,r) = LIS_rc%udef
                    endif
                 else
                    sm_current(c,r) = LIS_rc%udef
                 endif

              endif
           endif
        enddo
     enddo

  endif

  ! write(LIS_logunit,*) '[DEBUG] After CDF scaling (if applied)'
  ! write(LIS_logunit,*) '[DEBUG] sm_current valid count after CDF: ', count(sm_current.ne.LIS_rc%udef)
  ! call flush(LIS_logunit)

!-------------------------------------------------------------------------
!  Update obsl with CDF-scaled data
!-------------------------------------------------------------------------
  ! write(LIS_logunit,*) '[DEBUG] Updating obsl with CDF-scaled data...'
  ! call flush(LIS_logunit)

  obsl = LIS_rc%udef
  stored_count = 0
  skipped_ocean = 0

  do r = 1, LIS_rc%obs_lnr(k)
     do c = 1, LIS_rc%obs_lnc(k)
        if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
           ! Valid land point - store in obsl using gindex
           if(sm_current(c, r) .ne. LIS_rc%udef) then
              obsl(LIS_obs_domain(n, k)%gindex(c, r)) = sm_current(c, r)
              stored_count = stored_count + 1
              ! if(stored_count .le. 5) then
              !    write(LIS_logunit,*) '[DEBUG] Stored at gindex', LIS_obs_domain(n, k)%gindex(c, r), &
              !         ' (c=', c, ' r=', r, ') SM=', sm_current(c, r)
              ! endif
           endif
        else if (sm_current(c, r) .ne. LIS_rc%udef) then
           ! Data exists but gindex=-1 (ocean or outside domain)
           ! Cannot use for assimilation - skip
           skipped_ocean = skipped_ocean + 1
        endif
     enddo
  enddo

  write(LIS_logunit,*) '[INFO] Stored in obsl (land points):', stored_count
  write(LIS_logunit,*) '[INFO] Skipped (ocean/outside domain):', skipped_ocean

  ! Update fnd to reflect actual stored count
  fnd = stored_count

  write(LIS_logunit,*) '[INFO] Observations before LSM QC: ', count(obsl.ne.LIS_rc%udef)
  call flush(LIS_logunit)

  ! ! Print statistics of obsl values before LSM QC
  ! if(count(obsl.ne.LIS_rc%udef) .gt. 0) then
  !    write(LIS_logunit,*) '[DEBUG] obsl min/max before LSM QC: ', &
  !         minval(obsl, mask=obsl.ne.LIS_rc%udef), '/', &
  !         maxval(obsl, mask=obsl.ne.LIS_rc%udef)
  !    write(LIS_logunit,*) '[DEBUG] obsl sample values (first 10 valid):'
  !    t = 0
  !    do i = 1, min(LIS_rc%obs_ngrid(k), 100)
  !       if(obsl(i) .ne. LIS_rc%udef) then
  !          write(LIS_logunit,*) '  obsl(', i, ') = ', obsl(i)
  !          t = t + 1
  !          if(t .ge. 10) exit
  !       endif
  !    enddo
  ! endif
  ! call flush(LIS_logunit)

!-------------------------------------------------------------------------
!  Apply LSM based QC and screening of observations
!  RE-ENABLED: LSM QC threshold adjusted (shdfac: 0.9 -> 0.95)
!-------------------------------------------------------------------------
  if(fnd .gt. 0) then
     write(LIS_logunit,*) '[INFO] Applying LSM QC (relaxed threshold)...'
     call flush(LIS_logunit)

     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_NASA_AMSREsmobsId)//char(0), n, k, OBS_state)

     write(LIS_logunit,*) '[INFO] Observations after LSM QC: ', count(obsl.ne.LIS_rc%udef)

     if(count(obsl.ne.LIS_rc%udef) .eq. 0) then
        write(LIS_logunit,*) '[WARN] LSM QC removed ALL observations!'
        write(LIS_logunit,*) '[WARN] May need further LSM QC threshold adjustment'
     endif
     call flush(LIS_logunit)
  endif

  ! write(LIS_logunit,*) '[DEBUG] Calling LIS_checkForValidObs...'
  ! call flush(LIS_logunit)

  ! LIS_checkForValidObs returns binary flag (0 or 1), not count
  validObsFlag_temp = 0

  call LIS_checkForValidObs(n, k, obsl, validObsFlag_temp, sm_current)

  ! If no valid observations, set fnd to 0
  if(validObsFlag_temp .eq. 0) then
     write(LIS_logunit,*) '[WARN] No valid observations after all QC steps'
     fnd = 0
  endif

  write(LIS_logunit,*) '[INFO] Final observation count: ', fnd
  call flush(LIS_logunit)

  if (fnd .eq. 0) then
     data_upd_flag_local = .false.
  else
     data_upd_flag_local = .true.
  endif

#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local, 1, &
       MPI_LOGICAL, data_upd_flag(:), &
       1, MPI_LOGICAL, LIS_mpi_comm, status)
  data_upd = any(data_upd_flag)
#else
  data_upd = data_upd_flag_local
#endif

  ! write(LIS_logunit,*) '[DEBUG] data_upd = ', data_upd
  ! call flush(LIS_logunit)

  if (data_upd) then
     ! write(LIS_logunit,*) '[DEBUG] Setting attributes...'
     ! call flush(LIS_logunit)

     ! ! Debug: print LIS_rc%udef value
     ! write(LIS_logunit,*) '[DEBUG] LIS_rc%udef = ', LIS_rc%udef

     ! ! Debug: count non-udef values in obsl
     ! write(LIS_logunit,*) '[DEBUG] obsl values != LIS_rc%udef: ', &
     !      count(obsl(1:LIS_rc%obs_ngrid(k)) .ne. LIS_rc%udef)
     ! write(LIS_logunit,*) '[DEBUG] obsl values != -9999.0: ', &
     !      count(obsl(1:LIS_rc%obs_ngrid(k)) .ne. -9999.0)
     ! call flush(LIS_logunit)

     ! ! Print first few obsl values
     ! write(LIS_logunit,*) '[DEBUG] First 10 obsl values:'
     ! do t = 1, min(10, LIS_rc%obs_ngrid(k))
     !    write(LIS_logunit,*) '  obsl(', t, ') = ', obsl(t)
     ! enddo
     ! call flush(LIS_logunit)

     do t = 1, LIS_rc%obs_ngrid(k)
        gid(t) = t
        ! FIXED: Use LIS_rc%udef instead of hardcoded -9999.0
        if (obsl(t) .ne. LIS_rc%udef) then
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo

     write(LIS_logunit,*) '[INFO] Observations with assimflag=1: ', count(assimflag.eq.1)
     call flush(LIS_logunit)

     call ESMF_AttributeSet(OBS_State, "Data Update Status", &
          .true., rc=status)
     call LIS_verify(status)

     if (LIS_rc%obs_ngrid(k) .gt. 0) then
        call ESMF_AttributeSet(smField, "Grid Number", &
             gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error setting Grid Number attribute')

        call ESMF_AttributeSet(smField, "Assimilation Flag", &
             assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error setting Assimilation Flag attribute')

        call ESMF_AttributeSet(smfield, "Unscaled Obs", &
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')
     endif

  else
     call ESMF_AttributeSet(OBS_State, "Data Update Status", &
          .false., rc=status)
     call LIS_verify(status)
  endif

end subroutine read_NASA_AMSREsm


!BOP
! !ROUTINE: read_AMSREsm
! \label{read_AMSREsm}
!
! !INTERFACE:
subroutine read_AMSREsm(n, k, name, time)
! !USES:
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logmod,       only : LIS_logunit
  use LIS_timeMgrMod,   only : LIS_tick
  use NASA_AMSREsm_Mod, only : NASA_AMSREsm_struc
  use LIS_DAobservationsMod, only : LIS_obs_domain
#if (defined USE_HDF5)
  use hdf5
#endif
  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)  :: k
  character(len=*)     :: name
  real*8, intent(in)   :: time
!
! !DESCRIPTION:
!  Reads AMSR2 L2 soil moisture swath data from HDF5 files
!  and stores observation time for temporal filtering
!
!EOP
  real                 :: sb_rqc(NASA_AMSREsm_struc(n)%mo)

#if (defined USE_HDF5)
  ! AMSR2 L2 HDF5 variables
  integer(hid_t)        :: file_id
  integer(hid_t)        :: dset_id_sm, dset_id_lat, dset_id_lon, dset_id_qc
  integer(hsize_t), dimension(3) :: dims_sm
  integer(hsize_t), dimension(2) :: dims_latlon, dims_qc
  integer(hsize_t), dimension(3) :: maxdims_sm
  integer(hsize_t), dimension(2) :: maxdims_latlon, maxdims_qc
  integer(hid_t)        :: dspace_id

  ! AMSR2 L2 swath dimensions (will be read from file)
  integer               :: swath_nr, swath_nc

  ! Data arrays
  real,allocatable      :: sm_swath(:,:,:)      ! (1, swath_nc, swath_nr)
  integer*1,allocatable :: qc_swath(:,:,:)      ! (1, swath_nc, swath_nr) - QC flag (3D)
  real,allocatable      :: lat_swath(:,:)       ! (swath_nc, swath_nr)
  real,allocatable      :: lon_swath(:,:)       ! (swath_nc, swath_nr)

  ! 1D arrays for interpolation
  real,allocatable      :: sm_1d(:)
  real,allocatable      :: lat_1d(:)
  real,allocatable      :: lon_1d(:)
  logical*1,allocatable :: li(:)

  integer               :: mi, status, iret
  integer               :: t, t_best, i, j, idx, igd, valid_count, r, c
  integer               :: obs_c, obs_r, obs_idx, map_count
  integer               :: obs_lnc_global, obs_lnr_global  ! Global grid dimensions
  integer               :: qc_filtered, range_filtered
  integer, allocatable  :: count_interp(:)
  real                  :: udef, scale_factor
  real                  :: sm_value_fraction
  real                  :: lat_target, lon_target, min_dist, dist, dlat, dlon
  real                  :: obs_dlon, obs_dlat, obs_lon_start, obs_lat_start
  real                  :: obs_lat, obs_lon  ! For grid cell calculations
  real, parameter       :: MAX_DIST_THRESHOLD = 16.0  ! 4 degrees squared (obs grid ~1.25 deg)
  real,allocatable      :: sm_interp(:)

  ! For parsing filename timestamp
  integer               :: yr, mo, da, hr, mn, ss, doy
  real                  :: gmt
  real*8                :: timenow

  character*100         :: sm_field_name, lat_field_name, lon_field_name, qc_field_name
  character*100         :: filename
  integer               :: fname_idx
  integer(hid_t), save  :: my_real = -1
  integer(hid_t), save  :: my_int1 = -1
  logical, save         :: hdf5_initialized = .false.

  ! DEBUG: Print received time parameter
  write(LIS_logunit,*) '[DEBUG-TIME] read_AMSREsm called with time = ', time
  call flush(LIS_logunit)

  ! Field names for AMSR2 L2
  sm_field_name  = "Geophysical Data"
  lat_field_name = "Latitude of Observation Point"
  lon_field_name = "Longitude of Observation Point"
  qc_field_name  = "Pixel Data Quality"

  write(LIS_logunit,*) '[INFO] Reading AMSR2 L2 file: ', trim(name)
  call flush(LIS_logunit)

  udef = -9999.0
  scale_factor = 0.1  ! AMSR2 L2 scale factor
  igd = 1  ! AMSR2 L2 uses single index (not ascending/descending like AMSR-E)

  ! Extract filename from full path
  fname_idx = index(name, '/', back=.true.)
  if(fname_idx .gt. 0) then
     filename = name(fname_idx+1:len_trim(name))
  else
     filename = name
  endif

  ! Parse timestamp from filename: GW1AM2_YYYYMMDDHHMM_...
  ! Example: GW1AM2_201901312344_171A_L2SGSMCLF3300300.h5
  read(filename(8:11), '(i4)') yr
  read(filename(12:13), '(i2)') mo
  read(filename(14:15), '(i2)') da
  read(filename(16:17), '(i2)') hr
  read(filename(18:19), '(i2)') mn
  ss = 0

  ! write(LIS_logunit,*) '[INFO] Observation time: ', yr, '/', mo, '/', da, ' ', hr, ':', mn, ' UTC'

  ! Initialize HDF5 interface (only once per process)
  if(.not. hdf5_initialized) then
     call h5open_f(status)
     if(status .ne. 0) then
        write(LIS_logunit,*) '[ERR] Failed to initialize HDF5 interface, status=', status
        return
     endif
     hdf5_initialized = .true.
  endif

  ! Copy datatype for HDF5 reading (only once)
  if (my_real == -1) then
     call h5tcopy_f(H5T_NATIVE_REAL, my_real, status)
     if(status .ne. 0) then
        write(LIS_logunit,*) '[ERR] Failed to copy H5T_NATIVE_REAL, status=', status
        return
     endif
  endif

  if (my_int1 == -1) then
     ! Use H5T_STD_I8LE for 1-byte integer (QC data is uint8)
     call h5tcopy_f(H5T_STD_I8LE, my_int1, status)
     if(status .ne. 0) then
        write(LIS_logunit,*) '[ERR] Failed to copy H5T_STD_I8LE, status=', status
        return
     endif
  endif

  ! Open HDF5 file first to get dimensions
  call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to open file: ', trim(name)
     return
  endif

  ! Open Geophysical Data dataset to get dimensions
  call h5dopen_f(file_id, trim(sm_field_name), dset_id_sm, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to open Geophysical Data dataset'
     call h5fclose_f(file_id, status)
     return
  endif

  ! Get dataspace and dimensions
  call h5dget_space_f(dset_id_sm, dspace_id, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to get dataspace'
     call h5dclose_f(dset_id_sm, status)
     call h5fclose_f(file_id, status)
     return
  endif

  call h5sget_simple_extent_dims_f(dspace_id, dims_sm, maxdims_sm, status)
  if(status .eq. -1) then
     write(LIS_logunit,*) '[ERR] Failed to get dimensions'
     call h5sclose_f(dspace_id, status)
     call h5dclose_f(dset_id_sm, status)
     call h5fclose_f(file_id, status)
     return
  endif

  ! Extract actual dimensions from file
  ! HDF5 Fortran returns dimensions in reverse order: C (1979,243,1) -> Fortran (1,243,1979)
  swath_nr = int(dims_sm(3))  ! Last dimension = first in C order = 1979
  swath_nc = int(dims_sm(2))  ! Middle dimension = 243
  mi = swath_nr * swath_nc

  write(LIS_logunit,*) '[INFO] Swath size: ', swath_nr, ' x ', swath_nc, ' (', mi, ' points)'

  ! Close dataset and dataspace
  call h5sclose_f(dspace_id, status)
  call h5dclose_f(dset_id_sm, status)

  ! Allocate arrays now that we know the dimensions
  ! HDF5 Fortran API uses C order: (1979,243,1) stays as (1979,243,1) in Fortran
  allocate(sm_swath(swath_nr, swath_nc, 1))
  allocate(qc_swath(swath_nr, swath_nc, 1))
  allocate(lat_swath(swath_nr, swath_nc))
  allocate(lon_swath(swath_nr, swath_nc))
  allocate(sm_1d(mi))
  allocate(lat_1d(mi))
  allocate(lon_1d(mi))
  allocate(li(mi))
  ! sm_interp will be allocated later for observation domain size

  ! Read Geophysical Data (soil moisture)
  call h5dopen_f(file_id, trim(sm_field_name), dset_id_sm, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to open Geophysical Data dataset'
     goto 100
  endif

  ! HDF5 Fortran API uses C order: match file dimensions (1979,243,1)
  dims_sm(1) = swath_nr  ! 1979
  dims_sm(2) = swath_nc  ! 243
  dims_sm(3) = 1
  ! Use my_real to let HDF5 convert int16 to real automatically
  call h5dread_f(dset_id_sm, my_real, sm_swath, dims_sm, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to read Geophysical Data, status=', status
     call h5dclose_f(dset_id_sm, status)
     goto 100
  endif
  ! write(LIS_logunit,*) '[DEBUG] SM data read, min/max: ', minval(sm_swath), maxval(sm_swath)
  ! write(LIS_logunit,*) '[DEBUG] Positive SM count: ', count(sm_swath > 0)
  call h5dclose_f(dset_id_sm, status)

  ! Read Latitude
  call h5dopen_f(file_id, trim(lat_field_name), dset_id_lat, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to open Latitude dataset'
     goto 100
  endif

  ! HDF5 Fortran API uses C order: match file dimensions (1979,243)
  dims_latlon(1) = swath_nr  ! 1979
  dims_latlon(2) = swath_nc  ! 243
  call h5dread_f(dset_id_lat, my_real, lat_swath, dims_latlon, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to read Latitude'
     call h5dclose_f(dset_id_lat, status)
     goto 100
  endif
  call h5dclose_f(dset_id_lat, status)

  ! Read Longitude
  call h5dopen_f(file_id, trim(lon_field_name), dset_id_lon, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to open Longitude dataset'
     goto 100
  endif

  call h5dread_f(dset_id_lon, my_real, lon_swath, dims_latlon, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[ERR] Failed to read Longitude'
     call h5dclose_f(dset_id_lon, status)
     goto 100
  endif
  call h5dclose_f(dset_id_lon, status)

  ! Read Pixel Data Quality (QC flag)
  ! QC data is stored as uint8 (H5T_STD_U8LE) with shape (1, 1978, 243)
  ! QC=0 means normal retrieval (good data), QC!=0 means various issues
  call h5dopen_f(file_id, trim(qc_field_name), dset_id_qc, status)
  if(status .ne. 0) then
     write(LIS_logunit,*) '[WARN] Failed to open Pixel Data Quality dataset, using all data'
     qc_swath = 0  ! If QC not available, assume all good
  else
     ! Read QC using same dimensions as SM (3D: 1 x swath_nc x swath_nr)
     ! Use my_int1 (1-byte type) to match integer*1 array
     call h5dread_f(dset_id_qc, my_int1, qc_swath, dims_sm, status)
     if(status .ne. 0) then
        write(LIS_logunit,*) '[WARN] Failed to read Pixel Data Quality, using all data'
        qc_swath = 0
     ! else
     !    write(LIS_logunit,*) '[INFO] QC data read successfully'
     !    write(LIS_logunit,*) '[INFO] QC=0 count (good data): ', count(qc_swath.eq.0)
     !    write(LIS_logunit,*) '[INFO] QC!=0 count (flagged): ', count(qc_swath.ne.0)
     endif
     call h5dclose_f(dset_id_qc, status)
  endif

  ! Close file
  call h5fclose_f(file_id, status)

  ! Convert 2D swath to 1D arrays for interpolation
  ! Arrays are in C order: sm_swath(i, j, 1), lat_swath(i, j)
  ! Apply QC filtering: only use pixels with QC flag = 0 (normal retrieval)

  idx = 0
  valid_count = 0
  qc_filtered = 0
  range_filtered = 0

  !! DEBUG: Print first few lat/lon values to verify reading
  !write(LIS_logunit,*) '[DEBUG] First 5 lat/lon values from swath:'
  !do i = 1, min(5, swath_nr)
  !   do j = 1, min(5, swath_nc)
  !      write(LIS_logunit,'(A,I4,A,I4,A,F10.4,A,F10.4,A,F10.2)') &
  !           '  [', i, ',', j, '] lat=', lat_swath(i,j), ' lon=', lon_swath(i,j), &
  !           ' SM=', sm_swath(i,j,1)
  !   enddo
  !enddo
  call flush(LIS_logunit)

  do i = 1, swath_nr
     do j = 1, swath_nc
        idx = idx + 1
        ! Apply QC filter and scale factor
        ! QC = 0: normal retrieval, no warnings
        ! QC != 0: possible precipitation, invalid L1, errors, etc.
        if(qc_swath(i,j,1) .eq. 0 .and. sm_swath(i,j,1) .gt. 0) then
           sm_value_fraction = sm_swath(i,j,1) * scale_factor / 100.0  ! Convert % to fraction

           ! Apply relaxed physical range check (0.0 < SM < 1.0 or 100%)
           ! Satellite QC flag already filters bad retrievals, so trust QC flag
           ! and allow full physical range
           if(sm_value_fraction .ge. 0.001 .and. sm_value_fraction .le. 1.0) then
              sm_1d(idx) = sm_value_fraction
              lat_1d(idx) = lat_swath(i,j)
              lon_1d(idx) = lon_swath(i,j)
              li(idx) = .true.
              valid_count = valid_count + 1
           else
              ! Outside physical range (negative or >100%) - reject
              li(idx) = .false.
              range_filtered = range_filtered + 1
           endif
        else
           ! Failed QC check
           qc_filtered = qc_filtered + 1
           sm_1d(idx) = udef
           lat_1d(idx) = udef
           lon_1d(idx) = udef
           li(idx) = .false.
        endif
     enddo
  enddo

  write(LIS_logunit,*) '[INFO] Swath filtering: Total=', mi, ' Valid=', valid_count, &
       ' QC_filtered=', qc_filtered, ' Range_filtered=', range_filtered
  call flush(LIS_logunit)

  ! Perform interpolation to OBSERVATION domain (not model domain!)
  ! Map each swath point to observation grid cell using lat/lon calculation

  ! Get observation grid parameters from LIS obs_gridDesc
  obs_dlon = LIS_rc%obs_gridDesc(k,9)
  obs_dlat = LIS_rc%obs_gridDesc(k,10)

  ! For L2 swath data: Use global grid as intermediate buffer
  ! This ensures all swath points can be mapped regardless of location
  ! Memory usage: ~1440*720*4bytes = 4MB (acceptable)
  obs_lon_start = -180.0 + obs_dlon/2.0  ! West edge (-179.875 for 0.25 deg)
  obs_lat_start = -90.0 + obs_dlat/2.0   ! South pole (-89.875 for 0.25 deg)
  obs_lnc_global = nint(360.0 / obs_dlon)  ! Full longitude coverage: 360/0.25 = 1440
  obs_lnr_global = nint(180.0 / obs_dlat)  ! Full latitude coverage: 180/0.25 = 720

  write(LIS_logunit,*) '[INFO] Mapping swath to global grid: ', obs_lnc_global, ' x ', obs_lnr_global
  write(LIS_logunit,*) '[INFO] Global grid range: lon=[', -180.0, ',', 180.0, &
                       '] lat=[', -90.0, ',', 90.0, ']'
  write(LIS_logunit,*) '[INFO] LIS obs domain: lon=[', LIS_rc%obs_gridDesc(k,5), ',', &
                       LIS_rc%obs_gridDesc(k,8), '] lat=[', LIS_rc%obs_gridDesc(k,4), ',', &
                       LIS_rc%obs_gridDesc(k,7), ']'
  call flush(LIS_logunit)

  ! Reallocate sm_interp for GLOBAL observation domain
  if(allocated(sm_interp)) deallocate(sm_interp)
  allocate(sm_interp(obs_lnc_global*obs_lnr_global))
  sm_interp = udef

  ! Also track number of points mapped to each cell for averaging
  allocate(count_interp(obs_lnc_global*obs_lnr_global))
  count_interp = 0

  ! For each valid swath point, calculate which obs grid cell it belongs to
  map_count = 0
  do idx = 1, mi
     if(li(idx) .and. lat_1d(idx).ne.udef .and. lon_1d(idx).ne.udef) then
        ! Calculate observation grid column and row
        obs_c = nint((lon_1d(idx) - obs_lon_start) / obs_dlon) + 1
        obs_r = nint((lat_1d(idx) - obs_lat_start) / obs_dlat) + 1

        ! Check if within valid GLOBAL grid bounds
        if(obs_c .ge. 1 .and. obs_c .le. obs_lnc_global .and. &
           obs_r .ge. 1 .and. obs_r .le. obs_lnr_global) then

           obs_idx = (obs_r-1)*obs_lnc_global + obs_c

           ! Accumulate values for averaging
           if(count_interp(obs_idx) .eq. 0) then
              sm_interp(obs_idx) = sm_1d(idx)
              count_interp(obs_idx) = 1
           else
              ! Running average
              sm_interp(obs_idx) = (sm_interp(obs_idx)*count_interp(obs_idx) + sm_1d(idx)) / &
                                   real(count_interp(obs_idx)+1)
              count_interp(obs_idx) = count_interp(obs_idx) + 1
           endif

           map_count = map_count + 1
           !! Debug first few mappings
           !if(map_count .le. 10) then
           !   write(LIS_logunit,'(A,I6,A,F8.3,A,F9.3,A,F8.5,A,I5,A,I4,A,I4,A)') &
           !        '[DEBUG] Swath', idx, ': (', lat_1d(idx), ',', lon_1d(idx), &
           !        ') SM=', sm_1d(idx), ' -> grid(', obs_c, ',', obs_r, ') idx=', obs_idx, ')'
           !endif
        else
           !! DEBUG: Print first few out-of-bounds points
           if(map_count .eq. 0 .and. idx .le. 20) then
              write(LIS_logunit,'(A,I6,A,F8.3,A,F9.3,A,I5,A,I4,A,I5,A,I4,A)') &
                   '[DEBUG] OUT-OF-BOUNDS Swath', idx, ': (', lat_1d(idx), ',', lon_1d(idx), &
                   ') -> grid(', obs_c, ',', obs_r, ') limits=(', obs_lnc_global, ',', obs_lnr_global, ')'
           endif
        endif
     endif
  enddo

  write(LIS_logunit,*) '[INFO] Swath-to-grid mapping: ', map_count, ' points -> ', &
       count(sm_interp .ne. udef), ' unique cells'
  deallocate(count_interp)
  !call flush(LIS_logunit)

  ! Now sm_interp is in GLOBAL observation domain coordinates
  ! Extract the subset that matches LIS obs domain and copy to smobs array
  valid_count = 0

  ! DEBUG: Print global grid info
  !write(LIS_logunit,*) '[DEBUG] Global grid info:'
  !write(LIS_logunit,*) '  obs_lat_start=', obs_lat_start
  !write(LIS_logunit,*) '  obs_lon_start=', obs_lon_start
  !write(LIS_logunit,*) '  obs_dlat=', obs_dlat
  !write(LIS_logunit,*) '  obs_dlon=', obs_dlon
  !write(LIS_logunit,*) '  obs_lnr_global=', obs_lnr_global
  !write(LIS_logunit,*) '  obs_lnc_global=', obs_lnc_global
  !write(LIS_logunit,*) '[DEBUG] LIS obs domain:'
  !write(LIS_logunit,*) '  LIS_rc%obs_gridDesc(k,4)=', LIS_rc%obs_gridDesc(k,4), ' (south)'
  !write(LIS_logunit,*) '  LIS_rc%obs_gridDesc(k,5)=', LIS_rc%obs_gridDesc(k,5), ' (west)'
  !write(LIS_logunit,*) '  LIS_rc%obs_gridDesc(k,7)=', LIS_rc%obs_gridDesc(k,7), ' (north)'
  !write(LIS_logunit,*) '  LIS_rc%obs_gridDesc(k,8)=', LIS_rc%obs_gridDesc(k,8), ' (east)'
  !write(LIS_logunit,*) '  LIS_rc%obs_lnr(k)=', LIS_rc%obs_lnr(k)
  !write(LIS_logunit,*) '  LIS_rc%obs_lnc(k)=', LIS_rc%obs_lnc(k)
  call flush(LIS_logunit)

  ! Map from global grid to LIS obs grid
  ! LIS obs grid may be subset (e.g., southern hemisphere only)
  ! Calculate which rows/cols of global grid correspond to LIS obs grid
  do r = 1, LIS_rc%obs_lnr(k)
     do c = 1, LIS_rc%obs_lnc(k)
        ! Calculate lat/lon of this LIS obs grid cell
        obs_lat = LIS_rc%obs_gridDesc(k,4) + (r-1)*obs_dlat
        obs_lon = LIS_rc%obs_gridDesc(k,5) + (c-1)*obs_dlon

        ! Find corresponding cell in global grid
        obs_r = nint((obs_lat - obs_lat_start) / obs_dlat) + 1
        obs_c = nint((obs_lon - obs_lon_start) / obs_dlon) + 1

        ! Check bounds and copy
        if(obs_r .ge. 1 .and. obs_r .le. obs_lnr_global .and. &
           obs_c .ge. 1 .and. obs_c .le. obs_lnc_global) then

           obs_idx = (obs_r-1)*obs_lnc_global + obs_c
           t = (r-1)*LIS_rc%obs_lnc(k) + c

           if(sm_interp(obs_idx) .ne. udef) then
              ! Store observation value AND time (following SMAP pattern)
              NASA_AMSREsm_struc(n)%smobs(t, igd) = sm_interp(obs_idx)
              NASA_AMSREsm_struc(n)%smtime(t, igd) = time
              valid_count = valid_count + 1
              ! DEBUG: Print first few mappings
              if(valid_count .le. 5) then
                 write(LIS_logunit,*) '[DEBUG-TIME] Storing obs #',valid_count,&
                      ': t=',t,' igd=',igd,' smtime=',time
              endif
              if(valid_count .le. 2) then
                 write(LIS_logunit,'(A,I4,A,I4,A,I4,A,I4,A,F8.3,A,F8.3,A,F8.5)') &
                      '[DEBUG] Mapped: LIS(', c, ',', r, ') -> Global(', obs_c, ',', obs_r, &
                      ') lat=', obs_lat, ' lon=', obs_lon, ' SM=', sm_interp(obs_idx)
              endif
           endif
        endif
     enddo
  enddo

  write(LIS_logunit,*) '[INFO] Added ', valid_count, ' valid observations to global array'
  call flush(LIS_logunit)

100 continue

  ! Clean up arrays
  if(allocated(sm_swath)) deallocate(sm_swath)
  if(allocated(qc_swath)) deallocate(qc_swath)
  if(allocated(lat_swath)) deallocate(lat_swath)
  if(allocated(lon_swath)) deallocate(lon_swath)
  if(allocated(sm_1d)) deallocate(sm_1d)
  if(allocated(lat_1d)) deallocate(lat_1d)
  if(allocated(lon_1d)) deallocate(lon_1d)
  if(allocated(li)) deallocate(li)
  if(allocated(sm_interp)) deallocate(sm_interp)

#else
  write(LIS_logunit,*) '[ERR] LIS not compiled with HDF5 support'
  write(LIS_logunit,*) '[ERR] Recompile with -DUSE_HDF5'
  call LIS_endrun()
#endif

end subroutine read_AMSREsm


!BOP
! !ROUTINE: maskObs_basedonQC
! \label{maskObs_basedonQC}
!
! !INTERFACE:
subroutine maskObs_basedonQC(npts, smc, qc)
! !USES:
  use LIS_coreMod
  use LIS_logMod,  only : LIS_logunit

  implicit none
! !ARGUMENTS:
  integer             :: npts
  real                :: smc(npts)
  integer*2           :: qc(npts)
!
! !DESCRIPTION:
!  This subroutine masks the observations based on QC flags
!  For AMSR2 L2, QC masking is minimal since we filter during reading
!
!EOP
  integer             :: i

  ! For AMSR2 L2, most QC is done during file reading
  ! Here we just check for any remaining invalid values
  do i=1,npts
     if(qc(i) .ne. 0) then
        smc(i) = LIS_rc%udef
     endif
  enddo

end subroutine maskObs_basedonQC


!BOP
! !ROUTINE: maskObs_basedonTime
! \label{maskObs_basedonTime}
!
! !INTERFACE:
subroutine maskObs_basedonTime(npts, smc, smtime, rsmc, fnd, time1, time3)
! !USES:
  use LIS_coreMod
  use LIS_logMod

  implicit none
! !ARGUMENTS:
  integer             :: npts
  real                :: smc(npts,2)
  real*8              :: smtime(npts,2)
  real                :: rsmc(npts)
  integer             :: fnd
  real*8              :: time1  ! Current LIS time
  real*8              :: time3  ! Current LIS time + timestep
!
! !DESCRIPTION:
!  This subroutine masks observations based on time window
!  Following SMAP L2 pattern to avoid double counting:
!  - Only use observations where (smtime - time1) >= 0 AND < (time3 - time1)
!  - This ensures each observation is used only once per timestep
!
!EOP
  integer             :: i
  real*8              :: dt
  integer             :: valid_count_igd1, valid_count_igd2
  integer             :: filtered_count
  integer             :: first_valid_debug

  rsmc = LIS_rc%udef
  fnd = 0
  valid_count_igd1 = 0
  valid_count_igd2 = 0
  filtered_count = 0

  write(LIS_logunit,*) '[DEBUG-TIME] maskObs_basedonTime: Checking first few observations:'
  write(LIS_logunit,*) '[DEBUG-TIME]   time1=',time1,' time3=',time3,' window=',time3-time1

  ! Apply time window filtering (following SMAP pattern)
  ! dt is not defined as absolute value to avoid double counting
  first_valid_debug = 0
  do i = 1, npts
     ! Debug: print first few to check what's happening
     if(i .le. 5) then
        write(LIS_logunit,*) '[DEBUG-TIME]   obs #',i,': smc(i,1)=',smc(i,1),&
             ' smtime(i,1)=',smtime(i,1),' udef=',LIS_rc%udef
     endif

     ! Try igd=1 first
     if(smc(i,1) .ne. LIS_rc%udef .and. smtime(i,1) .ne. 0.0) then
        dt = smtime(i,1) - time1
        ! Debug: print first few VALID observations
        if(first_valid_debug .lt. 5) then
           first_valid_debug = first_valid_debug + 1
           write(LIS_logunit,*) '[DEBUG-TIME] Valid obs #',first_valid_debug,&
                ' at index ',i,': smtime=',smtime(i,1),' dt=',dt,&
                ' window=',time3-time1,' passes?',(dt .ge. 0.0 .and. dt .lt. (time3-time1))
        endif
        if(dt .ge. 0.0 .and. dt .lt. (time3 - time1)) then
           rsmc(i) = smc(i,1)
           fnd = fnd + 1
           valid_count_igd1 = valid_count_igd1 + 1
        else
           filtered_count = filtered_count + 1
        endif
     ! If igd=1 doesn't pass, try igd=2
     elseif(smc(i,2) .ne. LIS_rc%udef .and. smtime(i,2) .ne. 0.0) then
        dt = smtime(i,2) - time1
        if(dt .ge. 0.0 .and. dt .lt. (time3 - time1)) then
           rsmc(i) = smc(i,2)
           fnd = fnd + 1
           valid_count_igd2 = valid_count_igd2 + 1
        else
           filtered_count = filtered_count + 1
        endif
     endif
  enddo

  write(LIS_logunit,*) '[INFO] maskObs_basedonTime: Window [', time1, ',', time3, ')'
  write(LIS_logunit,*) '[INFO] maskObs_basedonTime: ', fnd, ' observations within time window'
  write(LIS_logunit,*) '[INFO]   - igd=1: ', valid_count_igd1, ' igd=2: ', valid_count_igd2
  write(LIS_logunit,*) '[INFO]   - Filtered (outside window): ', filtered_count
  call flush(LIS_logunit)

end subroutine maskObs_basedonTime


!BOP
! !ROUTINE: NASA_AMSREsm_filename
! \label{NASA_AMSREsm_filename}
!
! !INTERFACE:
subroutine NASA_AMSREsm_filename(name, ndir, yr, mo, da, hr)
! !USES:
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS:
  character(len=*)  :: name
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
!
! !DESCRIPTION:
!  This subroutine creates the AMSR2 L2 filename pattern based on the time and date
!  Pattern: GW1AM2_YYYYMMDDHH*_*_L2SGSMCLF*.h5 (hour-specific)
!
!  The arguments are:
!  \begin{description}
!  \item[name] name of the AMSR2 L2 soil moisture filename pattern
!  \item[ndir] name of the AMSR2 L2 soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[hr]  current hour
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo, fda, fhr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  ! AMSR2 L2 filename pattern: GW1AM2_YYYYMMDDHHMM_OOOX_L2SGSMCLF3300300.h5
  ! We use wildcard for minutes to match all passes within the hour
  name = trim(ndir)//'/'//trim(fyr)//'/'//trim(fmo)//'/'//'GW1AM2_'&
         //trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'*_*_L2SGSMCLF*.h5'

  write(LIS_logunit,*) '[DEBUG] AMSR2 L2 filename pattern: ', trim(name)

end subroutine NASA_AMSREsm_filename
