module snapwave_obspoints
   implicit none
contains
   ! 
   subroutine read_obs_points()
   !
   ! Reads obs files
   !
   use snapwave_data
   use interp
   !
   implicit none
   !
   real*4 dummy
   !
   integer m, n, stat, j1, j2, jdq
   !
   character(len=256)        :: line
   character(len=256)        :: line2
   !
   real*4, dimension(2)      :: value
   !
   ! Read observation points
   !
   nobs = 0
   !
   if (obsfile(1:4) /= 'none') then
      ! 
      write(*,*)'Reading observation points ...'
      !
      open(500, file=trim(obsfile))       
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nobs = nobs + 1
      enddo
      rewind(500)
      allocate(xobs(nobs))
      allocate(yobs(nobs))
      allocate(nameobs(nobs))     
      !
      !
      value(1) = 0.0
      value(2) = 0.0     
      !
      do n = 1, nobs
         read(500,'(a)')line
                           
         j1=index(line,"'")
         jdq=index(line,'"')
         if (j1 == 0 .and. jdq==0) then! no name supplied, give standard name
            j2 = 12
            nameobs(n) = ''
            write(nameobs(n)(1:j2), '(A8,I0.4)') 'station_', n       
         elseif (j1>0) then ! name supplied,         
            line2 = adjustl(trim(line(j1+1:256)))
            j2=index(line2,"'")      
            nameobs(n) = adjustl(trim(line2(1:j2-1)))
         else
            line2 = adjustl(trim(line(jdq+1:256)))
            j2=index(line2,'"')      
            nameobs(n) = adjustl(trim(line2(1:j2-1)))            
         endif 
         !
         read(line,*)(value(m), m = 1, 2)         
         xobs(n) = value(1)
         yobs(n) = value(2)
         ! 
      enddo             
      close(500)
      !
      ! Determine indices and weights of observation points
      !
      allocate(irefobs(4,nobs))
      allocate(nrefobs(no_nodes))
      allocate(wobs(4,nobs))
     !
      call make_map_fm (x, y, face_nodes, no_nodes, no_faces, xobs, yobs, nobs, wobs, irefobs, nrefobs)
!
      ! Allocate arrays output variables at observation points
      allocate(hm0obs(nobs)) 
      allocate(tpobs(nobs)) 
      allocate(hm0igobs(nobs))
      allocate(dwobs(nobs))
      allocate(dfobs(nobs))
      allocate(stobs(nobs))
      allocate(swobs(nobs))
      allocate(hm0xobs(nobs))
      allocate(hm0yobs(nobs))
      allocate(wdobs(nobs))
      
      !
      hm0obs = 0.d0
      tpobs = 0.d0
      hm0igobs = 0.d0
      dwobs = 0.d0
      dfobs = 0.d0
      stobs = 0.d0
      swobs = 0.d0
      hm0xobs = 0.d0
      hm0yobs = 0.d0
      wdobs = 0.d0
   !
   endif
   !
   end subroutine
   !
   subroutine update_obs_points ()
   !
   use snapwave_data
   use interp
   !
   if (nobs>0) then
      buf=H*sqrt(2.)
      call grmap(buf, no_nodes, hm0obs, nobs, irefobs, wobs, 4,  0)
      if (ig==1) then
         buf=H_ig*sqrt(2.)
         call grmap(buf, no_nodes, hm0igobs, nobs, irefobs, wobs, 4,  0)  
      endif
      buf=Tp
      call grmap(buf, no_nodes, tpobs, nobs, irefobs, wobs, 4,  0)
      buf=Dw
      call grmap(buf, no_nodes, dwobs, nobs, irefobs, wobs, 4,  0)
      buf=Df
      call grmap(buf, no_nodes, dfobs, nobs, irefobs, wobs, 4,  0)
      if (wind==1) then
         buf=SwE
         call grmap(buf, no_nodes, swobs, nobs, irefobs, wobs, 4,  0)  
         buf=SwA
         call grmap(buf, no_nodes, stobs, nobs, irefobs, wobs, 4,  0)
      endif
      buf=H*cos(thetam)
      call grmap(buf, no_nodes, hm0xobs, nobs, irefobs, wobs, 4,  0)
      buf=H*sin(thetam)
      call grmap(buf, no_nodes, hm0yobs, nobs, irefobs, wobs, 4,  0)
      wdobs=mod(270.-atan2(hm0yobs,hm0xobs)*180./pi+360.,360.)
   endif
   !
   end subroutine
   !
end module
