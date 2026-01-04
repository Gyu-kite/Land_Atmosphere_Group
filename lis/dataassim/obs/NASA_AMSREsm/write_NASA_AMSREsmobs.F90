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
!
! !ROUTINE: write_NASA_AMSREsmobs
! \label{write_NASA_AMSREsmobs}
!
! !REVISION HISTORY:
! 25Jan2008: Sujay Kumar; Initial Specification
! 22Dec2025: Updated to use obs_ngrid(k) for AMSR2 L2 support
!
! !INTERFACE:
subroutine write_NASA_AMSREsmobs(n, k, OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:

  integer,     intent(in)  :: n
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
! writes the transformed (interpolated/upscaled/reprojected)
! NASA AMSRE observations to a file
!
!EOP
  type(ESMF_Field)         :: smField
  logical                  :: data_update
  real, pointer            :: smobs(:)
  real                     :: smobs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  write(LIS_logunit,*) '[DEBUG-WRITE] Entering write_NASA_AMSREsmobs, n=', n, ' k=', k
  write(LIS_logunit,*) '[DEBUG-WRITE] obs_ngrid(k)=', LIS_rc%obs_ngrid(k)
  call flush(LIS_logunit)

  call ESMF_AttributeGet(OBS_State, "Data Update Status", &
       data_update, rc=status)
  call LIS_verify(status)
  write(LIS_logunit,*) '[DEBUG-WRITE] data_update=', data_update
  call flush(LIS_logunit)

  if(data_update) then

     write(LIS_logunit,*) '[DEBUG-WRITE] Getting Observation01 field...'
     call flush(LIS_logunit)
     call ESMF_StateGet(OBS_State, "Observation01",smField, &
          rc=status)
     call LIS_verify(status, 'Error in StateGet Observation01')
     write(LIS_logunit,*) '[DEBUG-WRITE] Got smField, status=', status
     call flush(LIS_logunit)

     call ESMF_FieldGet(smField, localDE=0, farrayPtr=smobs, rc=status)
     call LIS_verify(status, 'Error in FieldGet smobs')
     write(LIS_logunit,*) '[DEBUG-WRITE] Got smobs pointer, status=', status
     call flush(LIS_logunit)

     write(LIS_logunit,*) '[DEBUG-WRITE] Checking obs_ngrid(k)...'
     call flush(LIS_logunit)
     if(LIS_rc%obs_ngrid(k).gt.0) then
        write(LIS_logunit,*) '[DEBUG-WRITE] Getting Unscaled Obs attribute...'
        call flush(LIS_logunit)
        call ESMF_AttributeGet(smfield,"Unscaled Obs",smobs_unsc,&
             itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error in AttributeGet-Unscaled Obs')
        write(LIS_logunit,*) '[DEBUG-WRITE] Got Unscaled Obs, status=', status
        call flush(LIS_logunit)
     endif

     if(LIS_masterproc) then
        write(LIS_logunit,*) '[DEBUG-WRITE] Getting file unit number...'
        call flush(LIS_logunit)
        ftn = LIS_getNextUnitNumber()
        write(LIS_logunit,*) '[DEBUG-WRITE] ftn=', ftn
        call flush(LIS_logunit)
        call NASA_AMSRE_smobsname(n,k,obsname)
        write(LIS_logunit,*) '[DEBUG-WRITE] obsname=', trim(obsname)
        call flush(LIS_logunit)

        call LIS_create_output_directory('DAOBS')
        write(LIS_logunit,*) '[DEBUG-WRITE] Opening file...'
        call flush(LIS_logunit)
        open(ftn,file=trim(obsname), form='unformatted')
        write(LIS_logunit,*) '[DEBUG-WRITE] File opened'
        call flush(LIS_logunit)
     endif

     write(LIS_logunit,*) '[DEBUG-WRITE] Writing smobs_unsc...'
     call flush(LIS_logunit)
     call LIS_writevar_gridded_obs(ftn,n,k,smobs_unsc)
     write(LIS_logunit,*) '[DEBUG-WRITE] Writing smobs...'
     call flush(LIS_logunit)
     call LIS_writevar_gridded_obs(ftn,n,k,smobs)
     write(LIS_logunit,*) '[DEBUG-WRITE] Write complete'
     call flush(LIS_logunit)

     if(LIS_masterproc) then
        call LIS_releaseUnitNumber(ftn)
        write(LIS_logunit,*) '[DEBUG-WRITE] File unit released'
        call flush(LIS_logunit)
     endif

  endif

  write(LIS_logunit,*) '[DEBUG-WRITE] write_NASA_AMSREsmobs completed'
  call flush(LIS_logunit)

end subroutine write_NASA_AMSREsmobs

!BOP
! !ROUTINE: NASA_AMSRE_smobsname
! \label{NASA_AMSRE_smobsname}
!
! !INTERFACE:
subroutine NASA_AMSRE_smobsname(n,k,obsname)
! !USES:
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS:
  integer                 :: n
  integer                 :: k
  character(len=*)        :: obsname
!
! !DESCRIPTION:
!
!EOP

  character(len=12) :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'

end subroutine NASA_AMSRE_smobsname
