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
   real*4                    ::dummy
   !
   integer                   ::m, n, stat, j1, j2, jdq
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
      allocate(zsobs(nobs)) 
      allocate(tpobs(nobs)) 
      allocate(hm0igobs(nobs))
      allocate(dwobs(nobs))
      allocate(dfobs(nobs))
      allocate(stobs(nobs))
      allocate(swobs(nobs))
      allocate(hm0xobs(nobs))
      allocate(hm0yobs(nobs))
      allocate(wdobs(nobs))
      allocate(dirsprobs(nobs))
      
      !
      hm0obs = 0.d0
      zsobs = 0.d0
      tpobs = 0.d0
      hm0igobs = 0.d0
      dwobs = 0.d0
      dfobs = 0.d0
      stobs = 0.d0
      swobs = 0.d0
      hm0xobs = 0.d0
      hm0yobs = 0.d0
      wdobs = 0.d0
      dirsprobs = 0.d0
    !
    endif
   !
   end subroutine
   !
   subroutine update_obs_points ()
    !
    use snapwave_data
    !
   implicit none

   integer :: iobs, ip, k
   real*4  :: spread_deg
   real*4  :: weight
   real*4  :: sqrt2

    if (nobs>0) then
      sqrt2 = sqrt(2.0)
      do iobs = 1, nobs
         if (irefobs(1, iobs) > 0) then
            hm0obs(iobs) = 0.0
            zsobs(iobs) = 0.0
            tpobs(iobs) = 0.0
            dwobs(iobs) = 0.0
            dfobs(iobs) = 0.0
            hm0xobs(iobs) = 0.0
            hm0yobs(iobs) = 0.0
            dirsprobs(iobs) = 0.0
            if (ig == 1) hm0igobs(iobs) = 0.0
            if (wind == 1) then
               swobs(iobs) = 0.0
               stobs(iobs) = 0.0
            end if
            !
            do ip = 1, 4
               k = max(irefobs(ip, iobs), 1)
               weight = real(wobs(ip, iobs), kind(1.0))
               hm0obs(iobs) = hm0obs(iobs) + weight*H(k)
               zsobs(iobs) = zsobs(iobs) + weight*(depth(k) + zb(k))
               tpobs(iobs) = tpobs(iobs) + weight*Tp(k)
               dwobs(iobs) = dwobs(iobs) + weight*Dw(k)
               dfobs(iobs) = dfobs(iobs) + weight*Df(k)
               hm0xobs(iobs) = hm0xobs(iobs) + weight*H(k)*cos(thetam(k))
               hm0yobs(iobs) = hm0yobs(iobs) + weight*H(k)*sin(thetam(k))
               call obs_directional_spreading(ee(:, k), spread_deg)
               dirsprobs(iobs) = dirsprobs(iobs) + weight*spread_deg
               if (ig == 1) hm0igobs(iobs) = hm0igobs(iobs) + weight*H_ig(k)
               if (wind == 1) then
                  swobs(iobs) = swobs(iobs) + weight*SwE(k)
                  stobs(iobs) = stobs(iobs) + weight*SwA(k)
               end if
            end do
            !
            hm0obs(iobs) = hm0obs(iobs)*sqrt2
            if (ig == 1) hm0igobs(iobs) = hm0igobs(iobs)*sqrt2
            wdobs(iobs)=mod(270.-atan2(hm0yobs(iobs),hm0xobs(iobs))*rad2deg+360.,360.)
         end if
      end do
    endif
    !
    end subroutine
    !
   subroutine obs_directional_spreading(ee_point, spread_deg)
      !
      use snapwave_data, only: sector, ntheta, dtheta, deg2rad, rad2deg, thetamean, FILL_VALUE
      !
      implicit none
      !
      integer, parameter :: sp = kind(1.0)
      integer, parameter :: dp = kind(1.0d0)
      !
      real(sp), intent(in)  :: ee_point(:)
      real(sp), intent(out) :: spread_deg
      !
      integer  :: j
      real(dp) :: theta_deg, theta_rad
      real(dp) :: energy, weight
      real(dp) :: m0, a1, b1, r1
      real(dp) :: dtheta_deg
      real(dp) :: offset
      !
      if (ntheta <= 0 .or. dtheta <= 0.0_dp) then
         spread_deg = FILL_VALUE
         return
      end if
      !
      dtheta_deg = dtheta*rad2deg
      offset = 0.5_dp
      m0 = 0.0_dp
      a1 = 0.0_dp
      b1 = 0.0_dp
      !
      do j = 1, ntheta
         theta_deg = thetamean*rad2deg - sector/2.0 + (real(j - 1, dp) + offset)*dtheta_deg
         theta_rad = theta_deg*deg2rad
         energy = ee_point(j)
         if (energy < 0.0_dp) energy = 0.0_dp
         weight = energy*dtheta
         m0 = m0 + weight
         a1 = a1 + weight*cos(theta_rad)
         b1 = b1 + weight*sin(theta_rad)
      end do
      !
      if (m0 > 0.0_dp) then
         a1 = a1/m0
         b1 = b1/m0
         r1 = sqrt(a1*a1 + b1*b1)
         if (r1 > 1.0_dp) r1 = 1.0_dp
         if (r1 < 0.0_dp) r1 = 0.0_dp
         spread_deg = sqrt(2.0_dp*(1.0_dp - r1))*rad2deg
      else
         spread_deg = FILL_VALUE
      end if
      !
   end subroutine obs_directional_spreading
   !
end module
