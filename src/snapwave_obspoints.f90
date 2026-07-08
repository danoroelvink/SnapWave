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
   integer, parameter        :: dp = kind(1.0d0)
   real(dp)                  :: dummy
   !
   integer                   ::m, n, stat, j1, j2, jdq
   !
   character(len=256)        :: line
   character(len=256)        :: line2
   !
   real(dp), dimension(2)    :: value
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
      value(1) = 0.0_dp
      value(2) = 0.0_dp
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
      allocate(wdobs(nobs))
      allocate(dirsprobs(nobs))
      
      !
      hm0obs = FILL_VALUE
      zsobs = FILL_VALUE
      tpobs = FILL_VALUE
      hm0igobs = FILL_VALUE
      dwobs = FILL_VALUE
      dfobs = FILL_VALUE
      stobs = FILL_VALUE
      swobs = FILL_VALUE
      wdobs = FILL_VALUE
      dirsprobs = FILL_VALUE
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

   integer, parameter :: sp = kind(1.0), dp = kind(1.0d0)
   integer :: iobs, ip, k, itheta
   real(sp) :: weight_sp
   real(sp) :: sqrt2
   real(dp) :: weight
   real(dp) :: hm0x_sum, hm0y_sum
   real(dp) :: m0_obs, a1_obs, b1_obs
   real(dp) :: energy_weight, r1_obs
   real(dp) :: dtheta_dp, rad2deg_dp
   real(dp) :: energy_bin, theta_bin

    if (nobs>0) then
      sqrt2 = sqrt(2.0_sp)
      dtheta_dp = real(dtheta, dp)
      rad2deg_dp = real(rad2deg, dp)
      do iobs = 1, nobs
         hm0obs(iobs) = FILL_VALUE
         zsobs(iobs) = FILL_VALUE
         tpobs(iobs) = FILL_VALUE
         hm0igobs(iobs) = FILL_VALUE
         dwobs(iobs) = FILL_VALUE
         dfobs(iobs) = FILL_VALUE
         stobs(iobs) = FILL_VALUE
         swobs(iobs) = FILL_VALUE
         wdobs(iobs) = FILL_VALUE
         dirsprobs(iobs) = FILL_VALUE
         if (irefobs(1, iobs) > 0) then
            hm0obs(iobs) = 0.0
            zsobs(iobs) = 0.0
            tpobs(iobs) = 0.0
            dwobs(iobs) = 0.0
            dfobs(iobs) = 0.0
            dirsprobs(iobs) = 0.0
            hm0x_sum = 0.0_dp
            hm0y_sum = 0.0_dp
            m0_obs = 0.0_dp
            a1_obs = 0.0_dp
            b1_obs = 0.0_dp
            if (ig == 1) hm0igobs(iobs) = 0.0
            if (wind == 1) then
               swobs(iobs) = 0.0
               stobs(iobs) = 0.0
            end if
            !
            do ip = 1, 4
               k = max(irefobs(ip, iobs), 1)
               weight = wobs(ip, iobs)
               weight_sp = real(weight, sp)
               hm0obs(iobs) = hm0obs(iobs) + weight_sp*H(k)
               zsobs(iobs) = zsobs(iobs) + weight_sp*(depth(k) + zb(k))
               tpobs(iobs) = tpobs(iobs) + weight_sp*Tp(k)
               dwobs(iobs) = dwobs(iobs) + weight_sp*Dw(k)
               dfobs(iobs) = dfobs(iobs) + weight_sp*Df(k)
               hm0x_sum = hm0x_sum + weight*real(H(k), dp)*cos(real(thetam(k), dp))
               hm0y_sum = hm0y_sum + weight*real(H(k), dp)*sin(real(thetam(k), dp))
               do itheta = 1, ntheta
                  energy_bin = max(real(ee(itheta, k), dp), 0.0_dp)
                  theta_bin = real(theta(itheta), dp)
                  energy_weight = weight*energy_bin*dtheta_dp
                  m0_obs = m0_obs + energy_weight
                  a1_obs = a1_obs + energy_weight*cos(theta_bin)
                  b1_obs = b1_obs + energy_weight*sin(theta_bin)
               end do
               if (ig == 1) hm0igobs(iobs) = hm0igobs(iobs) + weight_sp*H_ig(k)
               if (wind == 1) then
                  swobs(iobs) = swobs(iobs) + weight_sp*SwE(k)
                  stobs(iobs) = stobs(iobs) + weight_sp*SwA(k)
               end if
            end do
            !
            hm0obs(iobs) = hm0obs(iobs)*sqrt2
            if (ig == 1) hm0igobs(iobs) = hm0igobs(iobs)*sqrt2
            wdobs(iobs)=real(mod(270.0_dp - atan2(hm0y_sum,hm0x_sum)*rad2deg_dp + &
                                 360.0_dp, 360.0_dp), sp)
            if (m0_obs > 0.0_dp) then
               a1_obs = a1_obs/m0_obs
               b1_obs = b1_obs/m0_obs
               r1_obs = sqrt(a1_obs*a1_obs + b1_obs*b1_obs)
               if (r1_obs > 1.0_dp) r1_obs = 1.0_dp
               if (r1_obs < 0.0_dp) r1_obs = 0.0_dp
               dirsprobs(iobs) = real(sqrt(2.0_dp*(1.0_dp - r1_obs))*rad2deg_dp, sp)
            else
               dirsprobs(iobs) = FILL_VALUE
            end if
         end if
      end do
    endif
    !
    end subroutine
end module
