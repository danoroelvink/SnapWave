module snapwave_date
   implicit none
   contains
   ! 
   subroutine time_difference(datespw,datesim,dtsec)
       !
       integer ijul1, ijul2     
       integer yyyy1,mm1,dd1,hh1,mn1,ss1,yyyy2,mm2,dd2,hh2,mn2,ss2
       integer dtsec,sec1,sec2
       !
       character*15  :: datespw
       character*15  :: datesim
       !
       read(datespw,'(I4,2I2,1X,3I2)')yyyy1,mm1,dd1,hh1,mn1,ss1
       read(datesim,'(I4,2I2,1X,3I2)')yyyy2,mm2,dd2,hh2,mn2,ss2
       !
       ! Compute Julian days
       !
       ijul1 = julian_date (yyyy1,mm1,dd1)
       ijul2 = julian_date (yyyy2,mm2,dd2)
       !
       sec1  = hh1*3600 + mn1*60 + ss1
       sec2  = hh2*3600 + mn2*60 + ss2
       dtsec = (ijul2 - ijul1)*86400 + sec2 - sec1
       !
   end subroutine
   
   FUNCTION julian_date (yyyy, mm, dd) RESULT (julian)
       ! converts calendar date to Julian date
       ! cf Fliegel & Van Flandern, CACM 11(10):657, 1968
       ! example: julian_date(1970,1,1)=2440588
       INTEGER,INTENT(IN) :: yyyy,mm,dd
       INTEGER :: julian
       julian = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + &
       367*(mm-2-((mm-14)/12)*12)/12- &
       3*((yyyy + 4900 + (mm - 14)/12)/100)/4
   END FUNCTION julian_date
   !
   function date_after_seconds(tref, t) result(str)
       !
       ! Calendar date/time 't' seconds after reference date tref ('yyyymmdd HHMMSS'),
       ! returned as 'yyyy-mm-dd HH:MM:SS'
       !
       character(len=*), intent(in) :: tref
       real*8, intent(in)           :: t
       character(len=19)            :: str
       integer :: y0, mo0, d0, h0, mi0, s0, jd, l, n, i, j, k
       integer(kind=8) :: secs, days, rem
       !
       read(tref,'(I4,2I2,1X,3I2)') y0, mo0, d0, h0, mi0, s0
       secs = int(h0*3600 + mi0*60 + s0, 8) + nint(t, 8)
       days = floor(dble(secs)/86400.d0, 8)
       rem  = secs - days*86400
       jd   = julian_date(y0, mo0, d0) + int(days)
       ! inverse of julian_date (Fliegel & Van Flandern)
       l = jd + 68569
       n = 4*l/146097
       l = l - (146097*n + 3)/4
       i = 4000*(l + 1)/1461001
       l = l - 1461*i/4 + 31
       j = 80*l/2447
       k = l - 2447*j/80
       l = j/11
       j = j + 2 - 12*l
       i = 100*(n - 49) + i + l
       write(str,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') i, j, k, &
             int(rem/3600), int(mod(rem,3600_8)/60), int(mod(rem,60_8))
       !
   end function date_after_seconds
   !
   function date_to_iso8601(date_string) result (date_iso8601)
       character(len=*), intent(in) :: date_string
       character*256 :: date_iso8601
       integer yyyy,mm,dd,hh,no_nodes,ss
       ! assume an input format of yyyymmdd HHMMSS
       read(date_string,'(I4,2I2,1X,3I2)')yyyy,mm,dd,hh,no_nodes,ss 
       ! return an iso8601 formatted date
       !                       yyyy -  mm   -  dd      HH   :  MM   :  SS
       write(date_iso8601, '(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)')  yyyy, '-', mm, '-', dd, ' ', hh, ':', no_nodes, ':', ss  !in DFM/FEWS without the T
       !write(date_iso8601, '(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)')  yyyy, '-', mm, '-', dd, 'T', hh, ':', no_nodes, ':', ss   
       !
   end function   
   !
   
   function convert_fewsdate(treftimefews, tbnd, ntbnd) result (tbndout)
       !
       !use sfincs_input
       !
       character*41 :: treftimefews, trefstr
       !
       integer ijul1, ijul2, ntbnd     
       integer yyyy1,mm1,dd1,hh1,mn1,ss1,yyyy2,mm2,dd2,hh2,mn2,ss2, tmp
       integer dtsec,sec1,sec2
       !
       real*4, dimension(:),     allocatable :: tbnd             
       real*4, dimension(:),     allocatable :: tbndout       
       !       
       !write(*,*)'trefstr',trefstr
       read(trefstr,'(I4,2I2,1X,3I2)')yyyy1,mm1,dd1,hh1,mn1,ss1
       !
       ! Fews time input: minutes since 1970-01-01 00:00:00.0 +0000       
       ! yyyy -  mm   -  dd      HH   :  MM   :  SS       
       read(treftimefews, '(A14,I4,A1,I2,A1,I2,A1,I2,A1,I2,A1,I2)')tmp,yyyy2,tmp,mm2,tmp,dd2,tmp,hh2,tmp,mn2,tmp,ss2 !assume time in UTC, whole seconds 
       !
       ijul1 = julian_date (yyyy1,mm1,dd1)
       ijul2 = julian_date (yyyy2,mm2,dd2)
       !
       sec1  = hh1*3600 + mn1*60 + ss1
       sec2  = hh2*3600 + mn2*60 + ss2
       dtsec = (ijul2 - ijul1)*86400 + sec2 - sec1
       !       
       allocate(tbndout(ntbnd))
       !
       tbndout = (tbnd * 60) + dtsec  ! time from fews is in minutes, then correct for use wrt sfincs treftime
       !       
   end function   
   
   !
end module