!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_gldas
! \label{read_gldas}
!
! !REVISION HISTORY:
!  19 Sept 2008: Sujay Kumar: Initial Implementation
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  24 Jan 2019: Sanghee Jun; Switched to the use netcdf
! !INTERFACE:
subroutine read_gldas( order, n, findex, name, ferror, try )
! !USES:  
  use LIS_coreMod,         only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_timeMgrMod,      only : LIS_get_nstep
  use LIS_metforcingMod,   only : LIS_forc
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_warning, &
                                  LIS_endrun
  use gldas_forcingMod,    only : gldas_struc

!#if (defined USE_GRIBAPI)
!  use grib_api
!#endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: order    
  integer, intent(in)           :: n
  integer, intent(in)           :: findex
  character(len=*),  intent(in) :: name
  integer, intent(out)          :: ferror 
  integer, intent(inout)        :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GLDAS data and spatially interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 3 hour GLDAS forecast file
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \item[try]
!    index of the tries (in case of missing data)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gldas](\ref{interp_gldas}) \newline
!    spatially interpolates a GLDAS variable
!  \item[fillgaps\_gldas](\ref{fillgaps_gldas}) \newline
!    fills the data gaps due to mismatches in landmask of LIS
!    domain and GLDAS mask.
!  \end{description}
!EOP
!==== Local Variables=======================
!arambarius
!  integer, parameter         :: iv_total=8
  integer, parameter         :: iv_total=7
  integer                    :: ftn, varid
  integer                    :: i, j, num, c,r,t
  integer                    :: nforce, ngldas
  logical                    :: file_exists
  integer                    :: var_index
  integer                    :: iv
  logical                    :: var_found,pcp_flag
  logical                    :: var_status(10)
  logical*1, allocatable         :: lb(:)
  real,  allocatable             :: f(:)
  real,  allocatable             :: datain(:,:)
  real                       :: missingValue
  real                       :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                    :: kk,nvars
  integer                    :: igrib
  integer                    :: pds7_1(10), pds6_1(10)
  integer                    :: pds5(10), pds7(10), pds6(10)
  integer                    :: pds5_val, pds7_val, pds6_val
  integer                    :: jpds5, jpds7, jpds6
  integer                    :: rc,status,iret
  character(20), dimension(iv_total), parameter :: svars=(/  &
       'Tair       ',    &
       'Qair      ',    & 
       'SWdown     ',    &
       'LWdown     ',    &
       'Wind      ',    &
       'PSurf      ',    &
!       For GLDAS
       'Rainf      '    &    
!       For ERA5
!       'Rainf_tavg      '    &    
! modified arambarius
!       'Snowf_tavg      '    &   
         /)

!=== End Variable Definition =======================

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
!#if(defined USE_GRIBAPI) 
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
!  pds5 = (/ 011,051,204,205,033,034,001,132,131,63 /) !parameter
!  pds6 = (/ 001,001,001,001,001,001,001,001,001,001 /)
!  pds7 = (/ 000,000,000,000,000,000,000,000,000,000 /) !htlev2

!  pds6_1 = (/ 001,001,001,001,105,105,001,001,001,001 /)
!  pds7_1 = (/ 000,000,000,000,010,010,000,000,000,000 /) !htlev2

  ngldas = (gldas_struc(n)%ncold*gldas_struc(n)%nrold)
  nforce = LIS_rc%met_nf(findex)
  ferror = 1  



!  iv_total = 9
  inquire (file=name, exist=file_exists)
  if (file_exists) then      
     ! Open netCDF file.
     iret = nf90_open(trim(name), nf90_NoWrite, ftn)
     write(LIS_logunit,*) "nf90_open",iret 

     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             'Could not open file: ',trim(name)
        ferror = 0
        return
     endif


     allocate(lb(gldas_struc(n)%ncold*gldas_struc(n)%nrold))
     allocate(f(gldas_struc(n)%ncold*gldas_struc(n)%nrold))
     allocate(datain(gldas_struc(n)%ncold,gldas_struc(n)%nrold))


     do kk=1,iv_total
        var_index=kk
        !if (kk .ge. 6) var_index=kk+1
        
        !initialize data
        f = LIS_rc%udef 
        lb = .false.
        
        !open and read file 
        iret = nf90_inq_varid(ftn, trim(svars(kk)), varid)
        write(LIS_logunit,*) "nf90_inq_varid", ftn, trim(svars(kk)), varid, iret
        iret = nf90_get_att(ftn, varid, "_FillValue", missingValue)
        iret = nf90_get_var(ftn, varid, datain, &
                 start=(/1,1,1,1/), &
                 count=(/gldas_struc(n)%ncold,gldas_struc(n)%nrold,1,1/))

       if(iret/=0) then
         if(LIS_masterproc) then 
            write(LIS_logunit,*)'[ERR] Problem get variable: ',svars(kk),iret
            write(LIS_logunit,*)'[ERR]  Stopping...'
            call LIS_endrun
         endif
         deallocate(lb)
         deallocate(f)
         deallocate(datain)
         call LIS_endrun
       end if
       !write(LIS_logunit,*) datain
       !=== Transferring current data to 1-D array for interpolation
       num=0
       do i=1,gldas_struc(n)%nrold
         do j=1,gldas_struc(n)%ncold
            num = num + 1
            f(num) = datain(j,i)
         enddo
       enddo

       !write(LIS_logunit,*) f 
       do t=1,ngldas
           if(f(t).ne.missingValue) then
                lb(t) = .true.
         
                ! 변수별 범위 검사 및 수정
                if(kk.eq.1 .and. (f(t) < 200.0 .or. f(t) > 350.0)) then  ! Tair
                    f(t) = 288.0
                endif
                if(kk.eq.2 .and. (f(t) < 0.0 .or. f(t) > 0.1)) then  ! Qair
                    f(t) = 0.01
                endif
                if(kk.eq.3 .and. f(t) < 0.0) then  ! SWdown
                    f(t) = 0.0
                endif
                if(kk.eq.4 .and. f(t) < 0.0) then  ! LWdown
                    f(t) = 200.0
                endif
                if(kk.eq.5 .and. f(t) < 0.0) then  ! Wind
                    f(t) = 0.5
                endif
                if(kk.eq.6 .and. (f(t) < 30000.0 .or. f(t) > 110000.0)) then  ! Psurf
                    f(t) = 101325.0
                endif
                if(kk.eq.7 .and. f(t) < 0.0) then  ! Rainf
                    f(t) = 0.0
                endif
            endif
           enddo 
          
           where ( f == missingValue )
              f = LIS_rc%udef
           endwhere
           pcp_flag = .false. 
           if(var_index.eq.7) pcp_flag = .true. 

           call interp_gldas(n, findex, pcp_flag,ngldas,&
                f,lb,LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

           call fillgaps_gldas(n,1,varfield)

           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gldas_struc(n)%metdata1(var_index,&
                            LIS_domain(n)%gindex(c,r)) =&
                            varfield(c,r)
                    else
                       gldas_struc(n)%metdata2(var_index,&
                            LIS_domain(n)%gindex(c,r)) = &
                            varfield(c,r)
                    endif
                 endif
              enddo
           enddo

     enddo
    ! Close netCDF file.
    iret=nf90_close(ftn)

     deallocate(lb)
     deallocate(f)   
     deallocate(datain) 
         

  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(name)
     ferror = 0
  endif

#endif


end subroutine read_gldas


!BOP
! !ROUTINE: interp_gldas
! \label{interp_gldas}
!
! !INTERFACE:
subroutine interp_gldas(n,findex,pcp_flag,ngldas,f,lb,lis_gds,nc,nr, &
                        varfield)
! !USES:
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use gldas_forcingMod, only : gldas_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngldas
  real, intent(out)   :: f(ngldas)
  logical*1           :: lb(ngldas)
  real                :: lis_gds(50)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine spatially interpolates a given GLDAS field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[findex]
!  index of the forcing source
! \item[pcp\_flag]
!  flag indicating if precip variables are being interpolated
! \item[ngldas]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo

  real, dimension(nc*nr) :: lis1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------  
  if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gldas_struc(n)%w111,gldas_struc(n)%w121,&
          gldas_struc(n)%w211,gldas_struc(n)%w221,&
          gldas_struc(n)%n111,gldas_struc(n)%n121,&
          gldas_struc(n)%n211,gldas_struc(n)%n221,LIS_rc%udef, iret)
  elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
     if (pcp_flag)then     
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gldas_struc(n)%w112,gldas_struc(n)%w122,&
             gldas_struc(n)%w212,gldas_struc(n)%w222,&
             gldas_struc(n)%n112,gldas_struc(n)%n122,&
             gldas_struc(n)%n212,gldas_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gldas_struc(n)%w111,gldas_struc(n)%w121,&
             gldas_struc(n)%w211,gldas_struc(n)%w221,&
             gldas_struc(n)%n111,gldas_struc(n)%n121,&
             gldas_struc(n)%n211,gldas_struc(n)%n221,LIS_rc%udef,iret)
     endif
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gldas_struc(n)%n113,LIS_rc%udef,iret)
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GLDAS & LDAS. For LDAS land 
! points not included in GLDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gldas
