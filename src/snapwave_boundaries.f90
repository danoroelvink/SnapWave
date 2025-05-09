module snapwave_boundaries
   !
contains   
   !
subroutine read_boundary_data()
   !
   ! Reads bnd, bhs, spw etc. files
   !Boundary conditions > lever tp, Hs en wave dir aand (ook dirspr?) en zet om naar spectrum op de randcellen (ntheta, nsigma) > JONSWAP met gamma als user waarde
   !Of tijdseries als bij FM > block per punt? > voorkeur van Maarten. tim-file
   
   !Voor b.c., sp2 files van XBeach code kopieren?    > een binarire file aanhouden > alvast een version header toevoegen oid
   !Maken matlab converter van sp2 naar binair
   !
   use snapwave_data
   !
   implicit none
   !
   nwbnd = 0
   ntwbnd = 0
   itwbndlast = 2
   itwindbndlast = 2
   !
   if (jonswapfile /= '') then
      !
      ! Read data from timeseries in single point
      !
      call read_boundary_data_singlepoint()
      !
   else
      !
      ! Read space- and time-varying data
      !
      call read_boundary_data_timeseries()
      !
   endif   
   !
   ! Compute reference table between wave nwbd boundary support points and grid boundary points 
   !
   call find_boundary_indices()
   !
end subroutine read_boundary_data
!
subroutine read_boundary_data_singlepoint()
   !
   ! Reads jonswap and water level file for single point
   ! t(s) Hm0 Tp dir ms zst
   !
   use snapwave_data   
   !
   implicit none
   !
   integer :: irec
   integer :: nrec
   integer :: ier
   real*4  :: dum
   !   
   write(*,*)'Reading boundary file ', trim(jonswapfile), ' ...'
   !
   ! Read jonswap wave time series
   !
   open(11,file=jonswapfile)
   irec=0
   ier=0
   do while (ier==0)
      read(11,*,iostat=ier)dum
      irec = irec + 1
   enddo
   rewind(11)
   nrec = irec - 1
   !
   allocate(t_bwv(nrec))   
   allocate(hs_bwv(1, nrec))
   allocate(tp_bwv(1, nrec))
   allocate(wd_bwv(1, nrec))
   allocate(ds_bwv(1, nrec))
   allocate(zs_bwv(1, nrec))
   !   
   do irec=1,nrec
      read(11,*)t_bwv(irec),hs_bwv(1, irec),tp_bwv(1, irec),wd_bwv(1, irec),ds_bwv(1, irec),zs_bwv(1, irec)      
   enddo
   !
   wd_bwv = (270.0 - wd_bwv)*pi/180
   ds_bwv = ds_bwv*pi/180.

   !
   close(11)
   !
   nwbnd = 1
   !
   allocate(hst_bwv(nwbnd))
   allocate(tpt_bwv(nwbnd))
   allocate(wdt_bwv(nwbnd))
   allocate(dst_bwv(nwbnd))
   allocate(zst_bwv(nwbnd)) 
   allocate(eet_bwv(ntheta, nwbnd)) 
   !
   ntwbnd = nrec
   !
end subroutine read_boundary_data_singlepoint
!
subroutine read_boundary_data_timeseries()
   !
   ! Reads bnd, bhs, spw etc. files
   !Boundary conditions > lever tp, Hs en wave dir aand (ook dirspr?) en zet om naar spectrum op de randcellen (ntheta, nsigma) > JONSWAP met gamma als user waarde
   !
   use snapwave_data
   !
   implicit none
   !
   integer n, itb, ib, stat
   !
   real*4 dummy
   !
   ! Read wave boundaries
   !
   if (bndfile(1:4) /= 'none') then    ! Normal ascii input files
      ! temporarily use this input of hs/tp/wavdir/dirspr in separate files, later change to DFM type tim-files
      !
      write(*,*)'Reading wave boundary locations ...'
      !
      open(500, file=trim(bndfile)) !as in bwvfile of SFINCS
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nwbnd = nwbnd + 1
      enddo
      rewind(500)
      allocate(x_bwv(nwbnd))
      allocate(y_bwv(nwbnd))
      do n = 1, nwbnd
         read(500,*)x_bwv(n),y_bwv(n)
      enddo
      close(500)
      !
      ! Read wave boundaries
      !
      write(*,*)'Reading wave boundaries ...'
      !
      ! Wave time series
      !
      ! First find times in bhs file
      !
      open(500, file=trim(bhsfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntwbnd = ntwbnd + 1
      enddo
      close(500)
      !
      allocate(t_bwv(ntwbnd))
      !
      allocate(hst_bwv(nwbnd))
      allocate(tpt_bwv(nwbnd))
      allocate(wdt_bwv(nwbnd))
      allocate(dst_bwv(nwbnd))
      allocate(zst_bwv(nwbnd)) 
      allocate(eet_bwv(ntheta,nwbnd)) 
      !
      ! Hs (significant wave height)
      ! Times in btp and bwd files must be the same as in bhs file!
      !
      open(500, file=trim(bhsfile))
      allocate(hs_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(hs_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Tp (peak period)
      !
      open(500, file=trim(btpfile))
      allocate(tp_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(tp_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Wd (wave direction)
      !
      open(500, file=trim(bwdfile))
      allocate(wd_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(wd_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Convert to cartesian, going-to, radians
      wd_bwv = (270.0 - wd_bwv)*pi/180.
      !
      ! Ds (directional spreading)
      !
      open(500, file=trim(bdsfile))
      allocate(ds_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(ds_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      ds_bwv = ds_bwv*pi/180.
      !
      ! zs (water level)
      !
      open(500, file=trim(bzsfile))
      allocate(zs_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(zs_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      write(*,*)'   Input boundary points found: ',nwbnd
      !
   endif
   !
end subroutine read_boundary_data_timeseries
!
subroutine find_boundary_indices()
   !
   use snapwave_data
   !
   implicit none
   !
   ! For each grid boundary point (kcs=2) :
   !
   ! Determine indices and weights of boundary points
   ! For tide and surge, these are the indices and weights of the points in the bnd file
   ! For waves, these are the indices and weights of the points in the cst file
   !
   integer k, ib1, ib2, ic
   !
   real xgb, ygb, dst1, dst2, dst
   !
   if (nwbnd>0 .and. nb>0) then  !nbnd
      !
      ! Count number of boundary points
      !
      ! Allocate boundary arrays
      !
      ! Water levels at boundary
      !
      ! Wave arrays
      !
      allocate(ind1_bwv_cst(nb))
      allocate(ind2_bwv_cst(nb))
      allocate(fac_bwv_cst(nb))
      !
      ! Find two closest boundary condition points for each boundary point
      ! And the two closest coastline points
      !
      nb = 0
      !
      ! Loop through all grid points
      !
      do k = 1, no_nodes
         !
         ! Check if this point is a boundary point
         !
         if (msk(k) == 2) then
            !
            nb = nb + 1
            nmindbnd(nb) = k
            !
            xgb = x(k)
            ygb = y(k)
            !
            ! Indices and weights for wave boundaries
            !
            if (nwbnd>1) then
               !
               dst1 = 1.0e10
               dst2 = 1.0e10
               ib1 = 0
               ib2 = 0
               !
               ! Loop through all water level boundary points
               !
               do ic = 1, nwbnd
                  !
                  ! Compute distance of this point to grid boundary point
                  !
                  dst = sqrt((x_bwv(ic) - xgb)**2 + (y_bwv(ic) - ygb)**2)
                  !
                  if (dst<dst1) then
                     !
                     ! Nearest point found
                     !
                     dst2 = dst1
                     ib2  = ib1
                     dst1 = dst
                     ib1  = ic
                     !
                  elseif (dst<dst2) then
                     !
                     ! Second nearest point found
                     !
                     dst2 = dst
                     ib2  = ic
                     !
                  endif
               enddo
               !
               ind1_bwv_cst(nb)  = ib1
               ind2_bwv_cst(nb)  = ib2
               fac_bwv_cst(nb) = dst2/(dst1 + dst2)
               !
            else
               !
               ind1_bwv_cst(nb)  = 1
               ind2_bwv_cst(nb)  = 1
               fac_bwv_cst(nb)   = 1.0
               !
            endif
            !
         endif
      enddo
   endif
   !
end subroutine find_boundary_indices  
!   
subroutine update_boundary_conditions(t)
   !
   ! Update all wave boundary conditions
   !
   use snapwave_data
   !
   implicit none
   !
   real*8, intent(in)        :: t
   !
   write(*,*)'Update boundary conditions ...'
   !   
   ! Update boundary conditions at boundary points
   !
   call update_boundary_points(t)
   !
   ! Update wind forcing 
   !
   if (wind) then
      call update_wind_field(t)
   endif
   !
   ! Update the upwind neighbours on 360 dir-grid given the updated wind direction
   !
   !call update_upwind_neigbours_wind()   
   !
   ! Make directional grid around boundary mean wave/wind direction
   !
   thetamean = wdmean_bwv
   if (ntwbnd > 0) then
      call make_theta_grid(wdmean_bwv)
   else
      thetamean=u10dmean
      call make_theta_grid(u10dmean)
   endif 
   !
   ! Build spectra on the boundary support points
   !
   call build_boundary_support_points_spectra()
   !
   ! Update boundary conditions at grid points
   !
   call update_boundaries()
   !
end subroutine update_boundary_conditions  
!
!subroutine update_upwind_neigbours_wind()
!   !
!   use snapwave_data
!   use snapwave_domain, only: find_upwind_neighbours_1dir
!   !
!   implicit none
!   !
!   real*4, dimension(ntheta360) :: windspread360k
!   integer                      :: itheta, k, count
!   !
!   write(*,*)'Update upwind neighbours wind'
!   !
!   ! Find the upwind neigbours in each cell for wind direction
!   !
!   call find_upwind_neighbours_1dir(x,y,no_nodes,sferic,kp,np, u10dir, wu10, prevu10, dsu10)
!   !
!   ! approximate grid resolution to closed upwind boundary with half of the ds360 average
!   !
!   do k = 1, no_nodes
!       if (dsu10(k) == 0.0) then
!          !
!          count = 0
!          do itheta = 1, ntheta360
!              !
!              if (ds360(itheta, k) > 0.0) then
!                  count = count + 1
!              endif
!              !
!          enddo
!          !
!          dsu10(k) = sum(ds360(:, k)) / count / 2.0
!          !
!       endif 
!   enddo 
!   !
!   ! Update the distribution array of wind input
!   ! 
!   do k = 1, no_nodes
!     windspread360k = (cos(theta360-u10dir(k)))**2.0
!     where(cos(theta360-u10dir(k))<0.0) windspread360k = 0.0
!     windspread360k = (windspread360k/sum(windspread360k))/dtheta ! normalized and converted to input per rad
!     windspread360(:,k) = windspread360k
!   enddo
!   !
!end subroutine update_upwind_neigbours_wind
!
subroutine update_boundary_points(t)
   !
   ! Update boundary conditions at boundary points
   !
   use snapwave_data
   !
   implicit none
   !
   real*8, intent(in)  :: t
   !   
   integer ib, itb
   !
   real*4  :: tbfac
   real*4  :: hs, tps, wd, dsp, zst
   !
   ! Interpolate boundary conditions in timeseries to boundary points
   !
   do itb = itwbndlast, ntwbnd ! Loop in time
      !
      if (t_bwv(itb)>t .or. itb==ntwbnd) then
         !
         tbfac  = (t - t_bwv(itb - 1))/(t_bwv(itb) - t_bwv(itb - 1))
         !
         do ib = 1, nwbnd ! Loop along boundary points
            !
            hs    = hs_bwv(ib, itb - 1) + (hs_bwv(ib, itb) - hs_bwv(ib, itb - 1))*tbfac
            tps   = tp_bwv(ib, itb - 1) + (tp_bwv(ib, itb) - tp_bwv(ib, itb - 1))*tbfac
            dsp    = ds_bwv(ib, itb - 1) + (ds_bwv(ib, itb) - ds_bwv(ib, itb - 1))*tbfac    !dirspr
            zst    = zs_bwv(ib, itb - 1) + (zs_bwv(ib, itb) - zs_bwv(ib, itb - 1))*tbfac  
            !
            call weighted_average(wd_bwv(ib, itb - 1), wd_bwv(ib, itb), 1.0 - tbfac, 2, wd)  !wavdir
            !
            hst_bwv(ib) = hs
            tpt_bwv(ib) = tps                  
            wdt_bwv(ib) = wd
            dst_bwv(ib) = dsp
            zst_bwv(ib) = zst
            !
         enddo
         !
         itwbndlast = itb
         exit
         !
      endif
   enddo
   !
   ! Now generate wave spectra at the boundary points
   !
   ! Average wave period and direction to determine theta grid
   !
   if (ntwbnd > 0) then
      tpmean_bwv = sum(tpt_bwv)/size(tpt_bwv)
      zsmean_bwv = sum(zst_bwv)/size(zst_bwv)
      depth = max(zsmean_bwv-zb,hmin)
      wdmean_bwv = atan2(sum(sin(wdt_bwv)*hst_bwv)/sum(hst_bwv),sum(cos(wdt_bwv)*hst_bwv)/sum(hst_bwv))  
   else
      depth = max(zsini-zb,hmin)
   endif
   !         
end subroutine update_boundary_points
!
subroutine update_wind_field(t)
   !
   ! Update wind field at all grid cells
   !
   use snapwave_data
   !
   implicit none
   !
   real*8, intent(in)  :: t
   !   
   integer itb, k
   real*4, dimension(:), allocatable :: windspread360k
   !
   real*4  :: tbfac
   real*4  :: u10k, u10dirk
   !
   allocate(windspread360k(ntheta360))
   windspread360k=0.0
   write(*,*)'Update wind field ...'
   !
   ! Interpolate boundary conditions in timeseries to boundary points
   !
   if (ntu10bnd==1) then
       !
       ! Just one 
       !
       do k=1, no_nodes
           u10(k) = u10_bwv(1, k)
           u10dir(k) = u10dir_bwv(1, k)
       enddo
       !
   else
       
      do itb = itwindbndlast, ntu10bnd ! Loop in time
         !
         if (t_u10_bwv(itb)>=t .or. itb==ntu10bnd) then
            !
            tbfac  = (t - t_u10_bwv(itb - 1))/(t_u10_bwv(itb) - t_u10_bwv(itb - 1))
            !
            do k = 1, no_nodes ! Loop along boundary points
               !
               u10k    = u10_bwv(itb - 1, k) + (u10_bwv(itb, k) - u10_bwv(itb - 1, k))*tbfac
               !
               call weighted_average(u10dir_bwv(itb - 1, k), u10dir_bwv(itb, k), 1.0 - tbfac, 2, u10dirk)  
               !
               u10(k) = u10k
               u10dir(k) = u10dirk                  
               !
            enddo
            !
            itwindbndlast = itb
            exit
            !
         endif
         ! 
      enddo
      !
   endif
   ! average wind direction
   ! 
   u10dmean = atan2(sum(sin(u10dir)*u10),sum(cos(u10dir)*u10))  
   !
   ! Initialize the distribution array of wind input
   ! 
   do k = 1,no_nodes
      windspread360k = (cos(theta360-u10dir(k)))**2.0
      where(cos(theta360-u10dir(k))<0.0) windspread360k = 0.0
      windspread360k = (windspread360k/sum(windspread360k))/dtheta ! normalized and converted to input per rad
      windspread360(:,k) = windspread360k
   enddo
   !   
end subroutine update_wind_field
!    
subroutine make_theta_grid(central_theta)
   !
   ! make theta grid based on boundary mean wave direction
   !
   use snapwave_data
   !
   implicit none
   !
   real, intent(in) :: central_theta
   integer k, itheta, ind
   !
   ! Determine theta grid and adjust w, prev and ds tables
   !
   ! Definition of directional grid
   !   
   ind=nint(central_theta/dtheta)-ntheta/2;
   do itheta = 1, ntheta
      i360(itheta)=mod2(itheta+ind,ntheta360)
   enddo
   !
   do itheta = 1, ntheta
      !
      theta(itheta) = theta360(i360(itheta))
      !
      do k = 1, no_nodes
         w(1, itheta, k)    = w360(1, i360(itheta), k)
         w(2, itheta, k)    = w360(2, i360(itheta), k)
         prev(1, itheta, k) = prev360(1, i360(itheta), k)
         prev(2, itheta, k) = prev360(2, i360(itheta), k)
         ds(itheta, k)      = ds360(i360(itheta), k)
         !
         windspreadfac(itheta, k) = windspread360(i360(itheta), k)
      enddo
      !
   enddo  
   !
   if (wind==1) then
       !
       ! initialization of distribution array of wind input
       !
       do k = 1, no_nodes
          windspreadfac(:,k)=(cos(theta-u10dir(k)))**mwind
          where(cos(theta-u10dir(k))<0.0) windspreadfac(:,k)=0.0
          if (sum(windspreadfac(:,k))>0.0) then      
             windspreadfac(:,k) = (windspreadfac(:,k)/sum(windspreadfac(:,k)))/dtheta
          else
             windspreadfac(:,k)=0.0
          endif  
       enddo
       !
   endif
   !
end subroutine make_theta_grid
!
subroutine build_boundary_support_points_spectra()
   !
   ! Update directional spectra on boundary points from time series
   !
   use snapwave_data
   !
   implicit none
   !
   integer :: ib
   !
   real*4  :: E0, ms
   !  
   do ib = 1, nwbnd ! Loop along boundary points
      E0   = 0.0625*rho*g*hst_bwv(ib)**2
      ms = 1.0/dst_bwv(ib)**2-1.0
      dist = sign(1.0,cos(theta - thetamean))*abs(cos(theta - thetamean))**ms
      where (abs(mod(pi+theta - thetamean,2.0*pi)-pi)>0.999*pi/2.0) dist = 0.0
      eet_bwv(:,ib) = dist/sum(dist)*E0/dtheta
   enddo
   !
end subroutine build_boundary_support_points_spectra
!
subroutine update_boundaries()
   !
   ! Update values at boundary points
   !
   use snapwave_data
   !
   implicit none
   !
   integer ib, i, k
   !      
   ! Set wave parameters in all boundary points on grid
   !
   ! Loop through grid boundary points
   ! Now for all grid boundary points do spatial interpolation of the wave spectra from boundary points at polygon
   !
   if (writetestfiles) open(113,file='testbnd.txt')
   do ib = 1, nb
      !
      k = nmindbnd(ib)
      !   
      do i = 1, ntheta
         ! 
         ee(i,k) = eet_bwv(i,ind1_bwv_cst(ib))*fac_bwv_cst(ib)  + eet_bwv(i,ind2_bwv_cst(ib))*(1.0 - fac_bwv_cst(ib))
         !                 
      enddo
      !
      if (wind/=1) then
          Tp(k) = tpmean_bwv
      else
          Tp(k) = Tpt_bwv(ind1_bwv_cst(ib))*fac_bwv_cst(ib)  + Tpt_bwv(ind2_bwv_cst(ib))*(1.0 - fac_bwv_cst(ib))
      endif
      !
      if (writetestfiles) write(113,*)ib,k,(ee(i,k),i=1,ntheta)
   enddo
   if (writetestfiles) close(113)
   !
   ! Update initial wave period to average bwv if not specified in .inp
   !
   if (wind/=1) then
       Tpini = tpmean_bwv
   endif
   !
end subroutine update_boundaries   
!
subroutine read_wind_data()
   !
   use snapwave_data
   use snapwave_domain, only: read_interpolate_map_input
   !
   implicit none
   !
   real*4  :: u10_0, u10dir_0
   integer :: it, k
   logical                                     :: fileu10_exists, fileu10dir_exists
   !
   if (windlistfile /= '') then
      !
      ! Read data from windlistfile: timeseries of spatially uniform / varying winds
      !
      call read_wind_data_from_list()
      !
   else
      !
      ! only one wind boundary condition given, read from .inp
      ! 
      ntu10bnd = 1
      !
      allocate(t_u10_bwv(ntu10bnd))   
      allocate(u10_bwv(ntu10bnd, no_nodes))
      allocate(u10dir_bwv(ntu10bnd, no_nodes))
      
      inquire(file=trim(u10str), exist=fileu10_exists)
      inquire(file=trim(u10dirstr), exist=fileu10dir_exists)
      !
      if (fileu10_exists) then
         !
         call read_interpolate_map_input(u10str,no_nodes, x, y, sferic, u10_bwv(1, :))
         !
      else
         !
         write(*,*)'   u10 has uniform value of ', trim(u10str)
         read( u10str, * )  u10_0
         u10_bwv(1, :) = u10_0
         !
      endif
      !
      if (fileu10dir_exists) then
          !
          call read_interpolate_map_input(u10dirstr,no_nodes, x, y, sferic, u10dir_bwv(1,:))
          !
      else
          !
          write(*,*)'   u10 has uniform value of ', trim(u10dirstr)
          read( u10dirstr, * )  u10dir_0 
          u10dir_bwv(1, :) = u10dir_0 
          !
      endif
      !
      do it = 1,ntu10bnd
          do k = 1,no_nodes
             u10dir_bwv(it, k) = mod(270.0-u10dir_bwv(it, k), 360.0)*(4.0*atan(1.0))/180.0 ! from nautical coming from in degrees to cartesian going to in radians
          enddo
      enddo
      !
   endif  
   
end subroutine read_wind_data
!
subroutine read_wind_data_from_list()
   !
   use snapwave_data
   use snapwave_domain, only: read_interpolate_map_input
   !
   implicit none
   !
   integer :: irec
   integer :: ier
   integer :: it
   real*4  :: dum
   real*4                              :: u10dir_0, u10_0
   character*256                       :: u10str_it, u10dirstr_it
   logical                             :: fileu10_exists, fileu10dir_exists
   !   
   write(*,*)'Reading windlist file ', trim(windlistfile), ' ...'
   !
   ! Find out how many consecutive wind fields are given
   !
   open(11,file=windlistfile)
   irec=0
   ier=0
   do while (ier==0)
      read(11,*,iostat=ier)dum
      irec = irec + 1
   enddo
   !
   ntu10bnd = irec - 1
   !
   allocate(t_u10_bwv(ntu10bnd))   
   allocate(u10_bwv(ntu10bnd, no_nodes))
   allocate(u10dir_bwv(ntu10bnd, no_nodes))
   !   
   ! Read all timestamps line by line 
   !
   rewind(11)
   do it = 1, ntu10bnd
      ! 
      !read(11,*)t_u10_bwv(it),u10str_it,u10dirstr_it   
      read(11,*)t_u10_bwv(it),u10str_it,u10dirstr_it   
      !
      inquire(file=trim(u10str_it), exist=fileu10_exists)
      inquire(file=trim(u10dirstr_it), exist=fileu10dir_exists)
      !
      ! read wind magnitude
      !
      if (fileu10_exists) then
         ! 
         call read_interpolate_map_input(u10str_it,no_nodes, x, y, sferic, u10_bwv(it,:))
         !
      else   
         !
         ! convert read value to double, and assign uniformly to grid
         !
         read(u10str_it, * )  u10_0
         u10_bwv(it, :) = u10_0
         !
      endif
      !
      ! read wind direction
      !
      if (fileu10dir_exists) then
         ! 
         call read_interpolate_map_input(u10dirstr_it,no_nodes, x, y, sferic, u10dir_bwv(it,:))
         !
      else   
         !
         ! convert read value to double, and assign uniformly to grid
         !
         read(u10dirstr_it, * )  u10dir_0
         u10dir_bwv(it, :) = u10dir_0
         !
      endif
      !   
   enddo
   !
   close(11)
   !
   u10dir_bwv=mod(270.0-u10dir_bwv, 360.0)*pi/180.0 ! from nautical coming from in degrees to cartesian going to in radians
   !
end subroutine read_wind_data_from_list
!
!subroutine read_wind_field()
!   !
!   use snapwave_data
!   use snapwave_domain, only: read_interpolate_map_input, find_upwind_neighbours_1dir
!   !
!   implicit none
!   !
!   real*4, dimension(ntheta360) :: windspread360k
!   real*4                       :: u10_0, u10dir_0
!   integer                      :: file_exists, k, itheta, count
!
!   !
!   ! Spatially-uniform/varying wind field magnitude
!   !
!   inquire(file=trim(u10str), exist=file_exists)
!   if (file_exists) then
!       call read_interpolate_map_input(u10str,no_nodes, x, y, sferic, u10)
!   else   
!       ! convert read value to double, and assign uniformly to grid
!       write(*,*)'u10 has uniform value of ', trim(u10str)
!       read( u10str, '(f10.4)' )  u10_0
!       u10=u10_0
!   endif
!   !
!   ! Spatially-uniform/varying wind field direction
!   !
!   inquire(file=trim(u10dirstr), exist=file_exists)
!   if (file_exists) then
!       call read_interpolate_map_input(u10dirstr,no_nodes, x, y, sferic, u10dir)
!   else   
!       ! convert read value to double, and assign uniformly to grid
!       write(*,*)'u10dir has uniform value of ', trim(u10dirstr)
!       read(u10dirstr, '(f10.4)' )  u10dir_0
!       u10dir=u10dir_0
!   endif
!   u10dir=mod(270.0-u10dir, 360.0)*(4.0*atan(1.0))/180.0 ! from nautical coming from in degrees to cartesian going to in radians
!   !
!   ! average wind direction
!   ! 
!   u10dmean = atan2(sum(sin(u10dir)*u10)/sum(u10),sum(cos(u10dir)*u10)/sum(u10))  
!   !
!   ! Initialize the distribution array of wind input
!   ! 
!   do k = 1,no_nodes
!      windspread360k = (cos(theta360-u10dir(k)))**2.0
!      where(cos(theta360-u10dir(k))<0.0) windspread360k = 0.0
!      windspread360k = (windspread360k/sum(windspread360k))/dtheta ! normalized and converted to input per rad
!      windspread360(:,k) = windspread360k
!   enddo
!   !
!   ! Find the upwind neigbours in each cell for wind direction
!   !
!   call find_upwind_neighbours_1dir(x, y, no_nodes, sferic, kp, np, u10dir, wu10,prevu10, dsu10)
!   !
!   ! approximate grid resolution to closed upwind boundary with half of the ds360 average
!   !
!   do k = 1, no_nodes
!       if (dsu10(k) == 0.0) then
!          !
!          count = 0
!          do itheta = 1, ntheta360
!              !
!              if (ds360(itheta, k) > 0.0) then
!                  count = count + 1
!              endif
!              !
!          enddo
!          !
!          dsu10(k) = sum(ds360(:, k)) / count / 2
!          !
!       endif 
!   enddo 
!   !
!end subroutine read_wind_field
!   
subroutine weighted_average(val1,val2,fac,iopt,val3) ! as in sfincs_boundaries.f90
   !
   implicit none
   !
   integer, intent(in)    :: iopt
   real*4,  intent(in)    :: val1
   real*4,  intent(in)    :: val2
   real*4,  intent(in)    :: fac
   real*4,  intent(out)   :: val3
   !
   real*4                 :: u1
   real*4                 :: v1
   real*4                 :: u2
   real*4                 :: v2
   real*4                 :: u
   real*4                 :: v
   !
   if (iopt==1) then
      !
      ! Regular
      !
      val3 = val1*fac  + val2*(1.0 - fac)
      !
   else
      !
      ! Angles (input must be in radians!)
      !
      u1 = cos(val1)
      v1 = sin(val1)
      u2 = cos(val2)
      v2 = sin(val2)
      !
      u = u1*fac  + u2*(1.0 - fac)
      v = v1*fac  + v2*(1.0 - fac)
      !
      val3 = atan2(v, u)
      !
   endif
   !
end subroutine weighted_average
!   
function mod2 (a,b) result (c)
   integer a,b,c
   !
   c=mod(a,b)
   if (c==0) c = b
   if (c<0)  c = c + b
   
end function
!   
!
!
end module
