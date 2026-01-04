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
!
! !ROUTINE: read_era5
! \label{read_era5}
! 
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
subroutine read_era5(n, kk,order, year, month, day, hour, findex,          &
     fname, ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use era5_forcingMod, only : era5_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: kk
  integer, intent(in)          :: order
  integer, intent(in)          :: year
  integer, intent(in)          :: month
  integer, intent(in)          :: day
  integer, intent(in)          :: hour
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: fname
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ERA5 data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
!  ZS: surface orography (m)
!  UREF: reference height for the wind (m)
!  ZREF: reference height for the temperature (m)
!
!  Tair: Atmospheric temperature (K) (2m)
!  Qair: Atmospheric humidity (kg/kg)
!  PSurf: Atmospheric pressure (Pa)
!  Rainf: Rain (kg/m2/s)
!  Snowf: Snow (kg/m2/s)
!  Wind: Wind speed (m/s) (10m)
!  Wind_DIR: Wind direction (degrees from N, clockwise)
!  LWdown: Long-wave radiation (W/m2)
!  DIR_SWdown: direct short-wave radiation (W/m2)
!  SCA_SWdown: diffuse short-wave radiation (W/m2)
!  CO2air: near surface CO2 concentration (kg/m3)
!   
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour ERA5 analysis file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP
  
  integer   :: ftn
  integer   :: tmpId, qId, windId, lwdId,psId,rainfId,snowfID,dirSWId,difSWId
  integer   :: tindex
  logical   :: file_exists
  integer   :: c,r,t,k,l,iret
  integer   :: mo,rec_size
  logical   :: read_lnd
  logical   :: read_flag
  
  real, allocatable      :: tair3d(:,:,:)
  real, allocatable      :: qair3d(:,:,:)
  real, allocatable      :: swd3d(:,:,:)
  real, allocatable      :: lwd3d(:,:,:)
  real, allocatable      :: wind3d(:,:,:)
  real, allocatable      :: ps3d(:,:,:)
  real, allocatable      :: rainf3d(:,:,:)
    
  real, allocatable      :: tair(:)
  real, allocatable      :: qair(:)
  real, allocatable      :: swd(:)
  real, allocatable      :: lwd(:)
  real, allocatable      :: wind(:)
  real, allocatable      :: ps(:)
  real, allocatable      :: rainf(:)
  
  integer :: ntime, nlat, nlon
  integer :: time_step_interval  ! 3시간 간격
  integer :: hour_index          ! 3시간별 시간 인덱스

  integer :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
! __________________________________________________________________________

  read_flag = .false.
  time_step_interval = 3  ! 3시간 간격

  ! Check if it is the switch of a month
  if(order.eq.1) then
     if(era5_struc(n)%mo1.ne.month) then
        era5_struc(n)%mo1 = month
        read_flag = .true.
     endif
  else
     if(era5_struc(n)%mo2.ne.month) then
        era5_struc(n)%mo2 = month
        read_flag = .true.
     endif
  endif
  ferror = 1

  if(read_flag) then 
#if (defined USE_NETCDF4) 

     if((mod(year,4) .eq. 0 .and. mod(year, 100).ne.0) &!leap year
          .or.(mod(year,400) .eq.0)) then 
        days(2) = 29
     else 
        days(2) = 28
     endif
     
     mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

     ! 3시간별 데이터: 하루 8개 시간대 (0,3,6,9,12,15,18,21)
     rec_size = days(month) * 24 / time_step_interval  ! 월별 총 시간 개수

! Read single layer file (*slv) fields:
     inquire(file=fname,exist=file_exists) 
     if(file_exists) then 
        write(LIS_logunit,*)'[INFO] Reading ERA5 3-hourly file (bookend,', order,' -',trim(fname), ')'

        call LIS_verify(nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_era5')

        ! 차원 정보 읽기
        call LIS_verify(nf90_inq_dimid(ftn, 'time', tmpId), &
             'nf90_inq_dimid failed for time')
        call LIS_verify(nf90_inquire_dimension(ftn, tmpId, len=ntime), &
             'nf90_inquire_dimension failed for time')
        
        call LIS_verify(nf90_inq_dimid(ftn, 'lat', tmpId), &
             'nf90_inq_dimid failed for lat')
        call LIS_verify(nf90_inquire_dimension(ftn, tmpId, len=nlat), &
             'nf90_inquire_dimension failed for lat')
             
        call LIS_verify(nf90_inq_dimid(ftn, 'lon', tmpId), &
             'nf90_inq_dimid failed for lon')
        call LIS_verify(nf90_inquire_dimension(ftn, tmpId, len=nlon), &
             'nf90_inquire_dimension failed for lon')

        write(LIS_logunit,*) '[INFO] ERA5 dimensions: time=', ntime, ' lat=', nlat, ' lon=', nlon
        ! 예상 차원과 비교
        write(LIS_logunit,*) '[DEBUG] Expected dimensions: time=248, lat=721, lon=1440'
        write(LIS_logunit,*) '[DEBUG] Actual dimensions: time=', ntime, ' lat=', nlat, ' lon=', nlon
        ! 차원이 예상과 다르면 경고
        if(ntime /= 248 .or. nlat /= 721 .or. nlon /= 1440) then
           write(LIS_logunit,*) '[WARN] Dimensions differ from expected values!'
        endif
        ! 3차원 배열 할당
        allocate(ps3d(ntime, nlat, nlon))
        allocate(tair3d(ntime, nlat, nlon))
        allocate(qair3d(ntime, nlat, nlon))
        allocate(wind3d(ntime, nlat, nlon))
        allocate(rainf3d(ntime, nlat, nlon))
        allocate(swd3d(ntime, nlat, nlon))
        allocate(lwd3d(ntime, nlat, nlon))
        write(LIS_logunit,*) '[DEBUG] Allocated 3D arrays successfully'

        ! 1차원 배열 할당
        allocate(ps(nlat*nlon))
        allocate(tair(nlat*nlon))
        allocate(qair(nlat*nlon))
        allocate(wind(nlat*nlon))
        allocate(rainf(nlat*nlon))
        allocate(swd(nlat*nlon))
        allocate(lwd(nlat*nlon))
        write(LIS_logunit,*) '[DEBUG] Allocated 1D arrays successfully'

        ! PSurf 변수 읽기 (범위 지정으로 안전하게)
        write(LIS_logunit,*) '[DEBUG] Looking for PSurf variable...'
        iret = nf90_inq_varid(ftn,'PSurf',psId)
        if(iret /= nf90_noerr) then
           write(LIS_logunit,*) '[ERR] Cannot find PSurf variable'
           call LIS_endrun()
        endif

        write(LIS_logunit,*) '[DEBUG] Attempting to read PSurf data with explicit ranges...'
        
        ! 시간별로 나누어 읽기
        do l = 1, ntime
           write(LIS_logunit,*) '[DEBUG] Reading time slice', l, 'of', ntime
           
           ! 한 번에 하나의 시간 슬라이스만 읽기
           iret = nf90_get_var(ftn, psId, ps3d(l:l,:,:), &
                start=[l,1,1], count=[1,nlat,nlon])
           
           if (iret /= nf90_noerr) then
              write(LIS_logunit,*) '[ERR] Failed to read time slice', l
              write(LIS_logunit,*) '[ERR] NetCDF error:', trim(nf90_strerror(iret))
              call LIS_endrun()
           endif
           
           if (mod(l, 50) == 0) then  ! 50개마다 진행상황 출력
              write(LIS_logunit,*) '[DEBUG] Read', l, 'time slices successfully'
           endif
        enddo
        
        write(LIS_logunit,*) '[DEBUG] Successfully read all PSurf data'
        
        ! NaN 처리
        do l=1,ntime
           do r=1,nlat
              do c=1,nlon
                 if (ps3d(l,r,c) /= ps3d(l,r,c)) then  ! NaN 체크
                    ps3d(l,r,c) = LIS_rc%udef
                 endif
              enddo
           enddo
        enddo


        !iret = nf90_inq_varid(ftn,'PSurf',psId)
        !if(iret /= nf90_noerr) then
        !   write(LIS_logunit,*) '[ERR] Cannot find PSurf variable'
        !   call LIS_endrun()
        !endif
        !! 변수 차원 정보 확인
        !call LIS_verify(nf90_inquire_variable(ftn, psId, ndims=iret), &
        !     'nf90_inquire_variable failed for pressure')
        !write(LIS_logunit,*) '[DEBUG] PSurf variable has', iret, 'dimensions'
        !
        !write(LIS_logunit,*) '[DEBUG] Attempting to read PSurf data...'
        !
        !iret = nf90_get_var(ftn, psId, ps3d)
        !if (iret /= nf90_noerr) then
        !   write(LIS_logunit,*) '[ERR] nf90_get_var failed for PSurf'  
        !   write(LIS_logunit,*) '[ERR] NetCDF error code:', iret
        !   write(LIS_logunit,*) '[ERR] NetCDF error message:', trim(nf90_strerror(iret))
        !   write(LIS_logunit,*) '[ERR] Array dimensions: ', ntime, nlat, nlon
        !   call LIS_endrun()
        !endif
        !
        !write(LIS_logunit,*) '[DEBUG] Successfully read PSurf data'
        !! NaN 값을 LIS_rc%udef로 변환
        !do l=1,ntime
        !   do r=1,nlat
        !      do c=1,nlon
        !         ! NaN 체크 및 변환 (isnan 함수 사용)
        !         if (ps3d(l,r,c) /= ps3d(l,r,c)) then  ! NaN 체크 방법
        !            ps3d(l,r,c) = LIS_rc%udef
        !         endif
        !      enddo
        !   enddo
        !enddo

        write(LIS_logunit,*) '[DEBUG] Looking for Tair variable...'
        call LIS_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
             'nf90_inq_varid failed for Tair in read_era5')
        call LIS_verify(nf90_get_var(ftn,tmpId, tair3d),&
             'nf90_get_var failed for tair in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read Tair data'

        write(LIS_logunit,*) '[DEBUG] Looking for Qair variable...'
        call LIS_verify(nf90_inq_varid(ftn,'Qair',qId), &
             'nf90_inq_varid failed for Qair in read_era5')
        call LIS_verify(nf90_get_var(ftn,qId, qair3d),&
             'nf90_get_var failed for qair in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read Qair data'

        write(LIS_logunit,*) '[DEBUG] Looking for Wind variable...'
        call LIS_verify(nf90_inq_varid(ftn,'Wind',windId), &
             'nf90_inq_varid failed for Wind in read_era5')
        call LIS_verify(nf90_get_var(ftn,windId, wind3d),&
             'nf90_get_var failed for wind in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read Wind data'

        write(LIS_logunit,*) '[DEBUG] Looking for Rainf variable...'
        call LIS_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
             'nf90_inq_varid failed for Rainf in read_era5')
        call LIS_verify(nf90_get_var(ftn,rainfId, rainf3d),&
             'nf90_get_var failed for rainf in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read Rainf data'

        ! SWdown 읽기 (DIR_SWdown 대신)
        write(LIS_logunit,*) '[DEBUG] Looking for SWdown variable...'
        call LIS_verify(nf90_inq_varid(ftn,'SWdown',dirSWId), &
             'nf90_inq_varid failed for SWdown in read_era5')
        call LIS_verify(nf90_get_var(ftn,dirSWId, swd3d),&
             'nf90_get_var failed for swd in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read SWdown data'

        write(LIS_logunit,*) '[DEBUG] Looking for LWdown variable...'
        call LIS_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
             'nf90_inq_varid failed for LWdown in read_era5')
        call LIS_verify(nf90_get_var(ftn,lwdId, lwd3d),&
             'nf90_get_var failed for lwd in read_era5')
        write(LIS_logunit,*) '[DEBUG] Successfully read LWdown data'

        write(LIS_logunit,*) '[DEBUG] All variables read successfully'
        ! 배열 초기화
        if(order.eq.1) then
           era5_struc(n)%ps1     = LIS_rc%udef
           era5_struc(n)%tair1   = LIS_rc%udef
           era5_struc(n)%qair1   = LIS_rc%udef
           era5_struc(n)%wind1   = LIS_rc%udef
           era5_struc(n)%rainf1  = LIS_rc%udef
           era5_struc(n)%swd1    = LIS_rc%udef
           era5_struc(n)%lwd1    = LIS_rc%udef
        else
           era5_struc(n)%ps2     = LIS_rc%udef
           era5_struc(n)%tair2   = LIS_rc%udef
           era5_struc(n)%qair2   = LIS_rc%udef
           era5_struc(n)%wind2   = LIS_rc%udef
           era5_struc(n)%rainf2  = LIS_rc%udef
           era5_struc(n)%swd2    = LIS_rc%udef
           era5_struc(n)%lwd2    = LIS_rc%udef
        endif

        ! 전체 시간대 처리 (timeinterp_era5에서 시간 보간용)
        do l=1,ntime
           ! 3차원에서 1차원으로 변환
           do r=1,nlat
              do c=1,nlon
                 k = c+(r-1)*nlon
                 ps(k) = ps3d(l,r,c)
                 tair(k) = tair3d(l,r,c)
                 qair(k) = qair3d(l,r,c)
                 wind(k) = wind3d(l,r,c)
                 rainf(k) = rainf3d(l,r,c)
                 swd(k) = swd3d(l,r,c)
                 lwd(k) = lwd3d(l,r,c)
              enddo
           enddo

           ! 시간 인덱스 조정 (3시간 → 1시간 확장)
           ! l=1 → 1,2,3시간, l=2 → 4,5,6시간, ...
           do t=1,3
              hour_index = (l-1)*3 + t
              if(hour_index <= 745) then  ! 배열 크기 내에서만
                 if(order.eq.1) then 
                    call interp_era5_var(n,findex,month,ps,     .false., era5_struc(n)%ps1(:,hour_index))
                    call interp_era5_var(n,findex,month,tair,   .false., era5_struc(n)%tair1(:,hour_index))
                    call interp_era5_var(n,findex,month,qair,   .false., era5_struc(n)%qair1(:,hour_index))
                    call interp_era5_var(n,findex,month,wind,   .false., era5_struc(n)%wind1(:,hour_index))
                    call interp_era5_var(n,findex,month,rainf,  .true.,  era5_struc(n)%rainf1(:,hour_index))
                    call interp_era5_var(n,findex,month,swd,    .false., era5_struc(n)%swd1(:,hour_index))
                    call interp_era5_var(n,findex,month,lwd,    .false., era5_struc(n)%lwd1(:,hour_index))
                 else
                    call interp_era5_var(n,findex,month,ps,     .false., era5_struc(n)%ps2(:,hour_index))
                    call interp_era5_var(n,findex,month,tair,   .false., era5_struc(n)%tair2(:,hour_index))
                    call interp_era5_var(n,findex,month,qair,   .false., era5_struc(n)%qair2(:,hour_index))
                    call interp_era5_var(n,findex,month,wind,   .false., era5_struc(n)%wind2(:,hour_index))
                    call interp_era5_var(n,findex,month,rainf,  .true.,  era5_struc(n)%rainf2(:,hour_index))
                    call interp_era5_var(n,findex,month,swd,    .false., era5_struc(n)%swd2(:,hour_index))
                    call interp_era5_var(n,findex,month,lwd,    .false., era5_struc(n)%lwd2(:,hour_index))
                 endif
              endif
           enddo
        enddo

        ! 메모리 해제
        deallocate(ps3d, tair3d, qair3d, wind3d, rainf3d, swd3d, lwd3d)
        deallocate(ps, tair, qair, wind, rainf, swd, lwd)

     else
        write(LIS_logunit,*) '[ERR] ',trim(fname)//' does not exist'
        call LIS_endrun()
     endif
#endif
  endif

  ! 기존 시간 인덱스 계산 유지 (timeinterp_era5에서 사용)
  tindex = (day - 1)*24 + hour + 1
  
  ! 기존 assign_processed_era5forc 호출 유지
  if(order.eq.1) then 
     call assign_processed_era5forc(n,kk,order,1,era5_struc(n)%tair1(:,tindex))
     call assign_processed_era5forc(n,kk,order,2,era5_struc(n)%qair1(:,tindex))
     call assign_processed_era5forc(n,kk,order,3,era5_struc(n)%swd1(:,tindex))
     call assign_processed_era5forc(n,kk,order,4,era5_struc(n)%lwd1(:,tindex))
     call assign_processed_era5forc(n,kk,order,5,era5_struc(n)%wind1(:,tindex))
     call assign_processed_era5forc(n,kk,order,6,era5_struc(n)%ps1(:,tindex))
     call assign_processed_era5forc(n,kk,order,7,era5_struc(n)%rainf1(:,tindex))
  else
     call assign_processed_era5forc(n,kk,order,1,era5_struc(n)%tair2(:,tindex))
     call assign_processed_era5forc(n,kk,order,2,era5_struc(n)%qair2(:,tindex))
     call assign_processed_era5forc(n,kk,order,3,era5_struc(n)%swd2(:,tindex))
     call assign_processed_era5forc(n,kk,order,4,era5_struc(n)%lwd2(:,tindex))
     call assign_processed_era5forc(n,kk,order,5,era5_struc(n)%wind2(:,tindex))
     call assign_processed_era5forc(n,kk,order,6,era5_struc(n)%ps2(:,tindex))
     call assign_processed_era5forc(n,kk,order,7,era5_struc(n)%rainf2(:,tindex))
  endif


end subroutine read_era5


!BOP
! 
! !ROUTINE: interp_era5_var
! \label{interp_era5_var}
! 
! !INTERFACE: 
subroutine interp_era5_var(n,findex, month, input_var, &
     pcp_flag, output_var)

! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use era5_forcingMod, only : era5_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(era5_struc(n)%npts)
  real,    intent(out)   :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical, intent(in)    :: pcp_flag

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a ERA5 field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lb(era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
! _____________________________________________________________

  input_size = era5_struc(n)%ncold*era5_struc(n)%nrold
  output_var = LIS_rc%udef

!-----------------------------------------------------------------------    
! Apply corrections
!-----------------------------------------------------------------------  
     
  lb = .false.
  do r=1,era5_struc(n)%nrold
     do c=1,era5_struc(n)%ncold           
        k= c+(r-1)*era5_struc(n)%ncold
        !if(era5_struc(n)%g2p(c,r).gt.0) then 
        !   f(k) = input_var(era5_struc(n)%g2p(c,r))
        !   lb(k) = .true.
        !else
        !   f(k) = LIS_rc%udef
        !   lb(k) = .false.
        if(k <= size(input_var) .and. input_var(k) /= LIS_rc%udef) then
           f(k) = input_var(k)
           lb(k) = .true.
        else
           f(k) = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo
  
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     
  if(pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0) then 
     call LIS_generatePcpClimoRatioField(n,findex,"ERA5",&
          month, & 
          input_size, &
          f, &
          lb)     
  endif
          
  if(pcp_flag.and.&
       trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%w112,era5_struc(n)%w122,&
          era5_struc(n)%w212,era5_struc(n)%w222,&
          era5_struc(n)%n112,era5_struc(n)%n122,&
          era5_struc(n)%n212,era5_struc(n)%n222,&
          LIS_rc%udef, iret)
     
  elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
       trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%w111,era5_struc(n)%w121,&
          era5_struc(n)%w211,era5_struc(n)%w221,&
          era5_struc(n)%n111,era5_struc(n)%n121,&
          era5_struc(n)%n211,era5_struc(n)%n221,&
          LIS_rc%udef, iret)
     
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),era5_struc(n)%mi,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%n113,LIS_rc%udef,iret)
  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LIS_rc%met_interp(findex))//&
          ' not supported for ERA5'
     call LIS_endrun()
  endif
  
  if( pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0 ) then 
     
     call LIS_pcpClimoDownscaling(n, findex, month,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), output_var(:), lo)
     
  endif

end subroutine interp_era5_var

!BOP
! 
! !ROUTINE: assign_processed_era5forc
! \label{assign_processed_era5forc}
! 
! !INTERFACE: 
subroutine assign_processed_era5forc(n,kk,order,var_index,era5forc)
! !USES: 
  use LIS_coreMod
  use era5_forcingMod, only : era5_struc
!
! !DESCRIPTION: 
!  This routine assigns the interpolated ERA5 forcing data
!  to the module data structures to be used later for 
!  time interpolation 
!
!EOP
  implicit none

  integer :: n
  integer :: kk
  integer :: order
  integer :: var_index
  real    :: era5forc(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  

  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              era5_struc(n)%metdata1(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           elseif(order.eq.2) then 
              era5_struc(n)%metdata2(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_era5forc



