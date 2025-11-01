module snapwave_domain
   implicit none

contains

subroutine initialize_snapwave_domain()
   !
   use snapwave_data
   use snapwave_results
   use snapwave_ncoutput
   use interp
   use m_alloc
   !
   implicit none
   !
   ! Local input variables
   !
   integer                                     :: k, j, ino, ic, m, n, itheta 
   integer*4, dimension(:),     allocatable    :: indices
   real*8,    dimension(:),     allocatable    :: theta360d0 ! wave angles,sine and cosine of wave angles
   real*8,    dimension(:,:),   allocatable    :: ds360d0
   real*8,    dimension(:,:,:), allocatable    :: w360d0
   character*2                                 :: ext
   logical                                     :: file_exists
   real*8, parameter                           :: dmiss=-99999.0
   !
   ! First set some constants
   !
   np   = 22 ! why?
   cosrot = cos(rotation*pi/180)
   sinrot = sin(rotation*pi/180)
   !
   write(*,*)'Initializing SnapWave domain ...'
   !
   j=index(gridfile,'.')      
   ext = gridfile(j+1:j+2)
   if (ext=='nc') then
      call nc_read_net()
      allocate(kp(np,no_nodes))
   !
   elseif (mskfile/='') then
      !
      ! Read structured index and mask files
      !
      ! Read index file (first time, only read number of active points)
      !
      open(unit = 500, file = trim(indfile), form = 'unformatted', access = 'stream')
      read(500)no_nodes
      close(500)
      !
      ! Allocate arrays
      write(*,*)'   ', no_nodes,' active SnapWave points'
      !
      allocate(x(no_nodes))
      allocate(y(no_nodes))
      allocate(xs(no_nodes))
      allocate(ys(no_nodes))
      allocate(zb(no_nodes))
      allocate(msk(no_nodes))
      allocate(index_v_n(no_nodes))
      allocate(index_v_m(no_nodes))
      allocate(index_v_nm(no_nodes))
      allocate(index_g_nm(nmax, mmax))
      allocate(indices(no_nodes))
      !
      ! Reading of grid-based variables
      !
      allocate(kp(np,no_nodes))
      !
      ! Read index file
      !
      write(*,*)'   Reading ',trim(indfile)
      open(unit = 500, file = trim(indfile), form = 'unformatted', access = 'stream')
      read(500)no_nodes
      read(500)indices
      close(500)
      !
      ! Read binary depth file
      !
      write(*,*)'   Reading ',trim(depfile)
      open(unit = 500, file = trim(depfile), form = 'unformatted', access = 'stream')
      read(500)zb
      close(500)
      !
      ! Read binary mask file
      !
      write(*,*)'   Reading ',trim(mskfile)
      open(unit = 500, file = trim(mskfile), form = 'unformatted', access = 'stream')
      read(500)msk
      close(500)
      !
      ! Grid indices
      !
      index_g_nm = 0
      !
      do k = 1, no_nodes
         !
         m = int((indices(k) - 1)/(nmax - 2)) + 1
         n = indices(k) - (m - 1)*(nmax - 2)
         m = m + 1
         n = n + 1
         index_v_m(k)    = m
         index_v_n(k)    = n
         index_v_nm(k)   = k
         index_g_nm(n, m) = k
         !
         ! Grid coordinates
         !
         x(k) = x0 + cosrot*(1.0*(m - 1) - 0.5)*dx - sinrot*(1.0*(n - 1) - 0.5)*dy
         y(k) = y0 + sinrot*(1.0*(m - 1) - 0.5)*dx + cosrot*(1.0*(n - 1) - 0.5)*dy
         !
      enddo
      !
      ! Connected nodes
      !
      ! Count cells
      !
      no_faces = 0
      !
      do k = 1, no_nodes
         !
         m = index_v_m(k)
         n = index_v_n(k)
         !
         if (index_g_nm(n    , m + 1)>0 .and. index_g_nm(n + 1, m + 1)>0 .and. index_g_nm(n + 1, m   )>0) then
            no_faces = no_faces + 1
         endif            
         !
      enddo
      !
      write(*,*)no_faces,' active SnapWave cells'
      !
      ! Allocate cells
      !
      allocate(face_nodes(4, no_faces))
      face_nodes = 0
      !
      ic = 0
      do k = 1, no_nodes
         !
         m = index_v_m(k)
         n = index_v_n(k)
         !
         if (index_g_nm(n    , m + 1)>0 .and. index_g_nm(n + 1, m + 1)>0 .and. index_g_nm(n + 1, m   )>0) then
            !
            ic = ic + 1
            !            
            face_nodes(1, ic) = k
            face_nodes(2, ic) = index_g_nm(n    , m + 1)
            face_nodes(3, ic) = index_g_nm(n + 1, m + 1)
            face_nodes(4, ic) = index_g_nm(n + 1, m    )       
            !
         endif            
         !
      enddo      
      !
      !
   else
      !
      ! Read snapwave unstructured mesh from ASCII file 
      !
      ! Read grid input
      !
      write(*,*)'   Reading mesh ...'
      !
      open(11,file=gridfile)
      read(11,*)no_nodes, no_faces !, no_edges
   !   allocate(xloc(no_nodes))
   !   allocate(yloc(no_nodes))
      allocate(x(no_nodes))
      allocate(y(no_nodes))
      allocate(xs(no_nodes))
      allocate(ys(no_nodes))
   !   allocate(x00(no_nodes))
   !   allocate(y00(no_nodes))
      allocate(zb(no_nodes))
      allocate(msk(no_nodes))
      allocate(face_nodes(4,no_faces))
     ! allocate(edge_nodes(2,no_edges))
      !
      ! Reading of grid-based variables
      !
      allocate(kp(np,no_nodes))
      !
      xmn = 1.0e9
      ymn = 1.0e9
      do k = 1, no_nodes
         read(11,*)x(k),y(k),zb(k),msk(k)
         zb(k) = max(zb(k), -200.0)
   !      xmn=min(xmn, x(k))
   !      ymn=min(ymn, y(k))
   !      xs(k) = 1.0*x(k)
   !      ys(k) = 1.0*y(k)
      enddo
   !   x = x - xmn
   !   y = y - ymn
      !
      do k = 1, no_faces
         read(11,*)(face_nodes(j, k), j = 1, 4)
      enddo
      !
      !do k = 1, no_edges
      !   read(11,*)(edge_nodes(j, k), j = 1, 2)
      !enddo
      close(11)
      !
   endif
   !
   zb = posdwn * zb
   !
   do k=1,no_faces
       if (face_nodes(4,k)==0) face_nodes(4,k) = -999
   enddo 
   !
   ! keep on also if ja_vegetation==0, so array Dveg is initialized with zeroes
   if (ja_vegetation==1) then
      call veggie_init()   
   else
      allocate(veg_Cd(no_nodes, no_secveg))
      allocate(veg_ah(no_nodes, no_secveg))
      allocate(veg_bstems(no_nodes,  no_secveg))
      allocate(veg_Nstems(no_nodes,  no_secveg))
   endif
   !
   ntheta360 = nint(360./dtheta)
   ntheta    = nint(sector/dtheta)
   !
   ! Allocation of spatial arrays
   !
   !allocate(inner(no_nodes))
   allocate(depth(no_nodes))
   allocate(dhdx(no_nodes))
   allocate(dhdy(no_nodes))
   allocate(kwav(no_nodes))
   allocate(nwav(no_nodes))
   allocate(C(no_nodes))
   allocate(Cg(no_nodes))
   allocate(sinhkh(no_nodes))
   allocate(Hmx(no_nodes))
   allocate(fw(no_nodes))
   allocate(H(no_nodes))
   allocate(Tp(no_nodes))
   allocate(kwav_ig(no_nodes))
   allocate(nwav_ig(no_nodes))
   allocate(C_ig(no_nodes))
   allocate(Cg_ig(no_nodes))
   allocate(sinhkh_ig(no_nodes))
   allocate(Hmx_ig(no_nodes))
   allocate(fw_ig(no_nodes))
   allocate(H_ig(no_nodes))
   allocate(Dveg(no_nodes))
   allocate(Dw(no_nodes))
   allocate(F(no_nodes))
   allocate(Fx(no_nodes))
   allocate(Fy(no_nodes))
   allocate(u10(no_nodes))
   allocate(u10dir(no_nodes))
   allocate(Df(no_nodes))
   allocate(thetam(no_nodes))
   allocate(ctheta(ntheta,no_nodes))
   allocate(ctheta_ig(ntheta,no_nodes))
   allocate(ctheta360(ntheta360,no_nodes))
   allocate(w(2,ntheta,no_nodes))
   allocate(w360(2,ntheta360,no_nodes))
   allocate(w360d0(2,ntheta360,no_nodes))
   allocate(prev(2,ntheta,no_nodes))
   allocate(prev360(2,ntheta360,no_nodes))
   allocate(ds(ntheta,no_nodes))
   allocate(ds360(ntheta360,no_nodes))
   allocate(ds360d0(ntheta360,no_nodes))
   allocate(i360(ntheta))
   allocate(neumannconnected(no_nodes))
   allocate(theta(ntheta))
   allocate(dist(ntheta))
   allocate(theta360d0(ntheta360))
   allocate(theta360(ntheta360))
   allocate(windspread360(ntheta360, no_nodes))
   allocate(windspreadfac(ntheta, no_nodes))
   allocate(ee(ntheta,no_nodes))
   allocate(ee_ig(ntheta,no_nodes))
   allocate(aa    (ntheta,no_nodes))
   allocate(sig   (no_nodes))
   allocate(WsorE (ntheta,no_nodes))
   allocate(WsorA (ntheta,no_nodes))
   allocate(SwE   (no_nodes))
   allocate(SwA   (no_nodes))
   !
   ! Spatially-uniform/varying bottom friction coefficients
   !
   ! (1) Fw
   !
   inquire(file=trim(fwstr), exist=file_exists)
   if (file_exists) then
      call read_interpolate_map_input(fwstr,no_nodes, x, y, sferic, fw)
   else   
      ! convert read value to double, and assign uniformly to grid
      write(*,*)'   fw has uniform value of ', trim(fwstr)
      read( fwstr, * )  fw0
      fw=fw0
   endif
   !
   where(zsini-zb>fwcutoff) fw = 0d0
   !
   ! (2) fw_ig
   !
   inquire(file=trim(fw_igstr), exist=file_exists)
   if (file_exists) then
      call read_interpolate_map_input(fw_igstr,no_nodes, x, y, sferic, fw_ig)
   else   
      ! convert read value to double, and assign uniformly to grid
      write(*,*)'   fw_ig has uniform value of ', trim(fw_igstr)
      read( fw_igstr, * )  fw0_ig
      fw_ig=fw0_ig
   endif
   !
   where(zsini-zb>fwcutoff) fw_ig = 0d0
   !
   ! Initialization of reference tables
   !
   ! Definition of directional grid (dtheta is user input)
   !
   do itheta = 1, ntheta360
      theta360d0(itheta) = 0.0 + 0.5*dtheta + (itheta - 1)*dtheta
      theta360(itheta)   = 1.0*theta360d0(itheta)
   enddo
   !
   theta360d0 = theta360d0*pi/180.0
   !
   kp   = 0
   w    = 0.0
   prev   = 0
   ds     = 0.0
   ee     = 0.0
   ee_ig  = 0.0
   ds360d0=0.d0
   w360d0 =0.d0
   prev360=0
   aa     = 0.0
   sig    = 0.0
   WsorE  = 0.0
   WsorA  = 0.0
   SwE    = 0.0
   SwA    = 0.0
   windspreadfac = 0.0
   dist = 0.0
   !
   if (upwfile=='') then   
      !
      write(*,*)'   Getting surrounding points ...'
      !
      call fm_surrounding_points(x, y, zb, no_nodes, sferic, face_nodes, no_faces, kp, np, dhdx, dhdy)
      !   
      write(*,*)'   Finding upwind neighbors ...'
      call find_upwind_neighbours(x, y, no_nodes, sferic, theta360d0, ntheta360, kp, np, w360d0, prev360, ds360d0)    
      !
      ! Distances and weights over 360 degrees (single precision)   
      !
      ds360 = ds360d0*1.0
      w360  = w360d0*1.0
      !
      !
      ! Read polygon outlining valid boundary points
      !
      call read_boundary_enclosure ()

      if (n_bndenc > 0) then
         !
         do k=1,no_nodes
             do itheta=1,ntheta360
                 if (ds360d0(itheta,k)==0.d0) then
                     call ipon(x_bndenc,y_bndenc,n_bndenc,x(k),y(k),ino)
                     if (ino>0) msk(k)=2
                 endif
             enddo
         enddo
         !
      endif
      !
      ! If present read neumann polyline
      !
      call read_neumann_boundary()
      if (n_neu>0) then
         !
         call neuboundaries(x,y,no_nodes,x_neu,y_neu,n_neu,tol,neumannconnected)
         !
         do k=1,no_nodes
           if (neumannconnected(k)>0) then
              if (msk(k)==1) then
                ! k is inner and can be neumannconnected
                ! inner(neumannconnected(k))= .false.
                msk(neumannconnected(k)) = 3
              else
                ! we don't allow neumannconnected links if the node is an open boundary
                neumannconnected(k) = 0  
              endif
           endif
         enddo
         !
      else
          neumannconnected = 0
      endif
      !
      if (writetestfiles) then
         open(112,file='masktest.txt')
         do k=1,no_nodes
             write(112,*)x(k),y(k),dble(msk(k))
         enddo
         close(112)
      endif
      !
      nb = 0
      nnmb = 0
      do k = 1, no_nodes
         !if (msk(k)==1) then
         !   inner(k) = .true.
         !else
         !   inner(k) = .false.
         !endif
         if (msk(k)==2) nb = nb + 1
         if (msk(k)==3) nnmb = nnmb + 1
      enddo
      !
      allocate(nmindbnd(nb))
      !
      if (nnmb>0) then
         allocate(neubnd(nnmb))
      endif
      !
      nb = 0
      nnmb = 0
      do k = 1, no_nodes
         if (msk(k)==2) then  ! ML: this used to be >1 but since neumann connected points are now msk(k)=3 changed this to ==2
            nb = nb + 1
            nmindbnd(nb) = k
         elseif (msk(k)==3) then  
            nnmb = nnmb + 1
            neubnd(nnmb) = k   
         endif   
      enddo   
      !
      ! Write upwind neighbors to file
      !
      upwfile = 'snapwave.upw'
      if (upwfile/='') then
         open(145,file=upwfile,access='stream', form='unformatted')
         write(145)prev360
         write(145)w360
         write(145)ds360
         write(145)dhdx
         write(145)dhdy
         close(145)
      endif
      !
   else
      !
      write(*,*)'   Reading upwind neighbors file ...'
      !
      open(unit=145, file=trim(upwfile), form='unformatted', access='stream')
      read(145)prev360
      read(145)w360
      read(145)ds360
      read(145)dhdx
      read(145)dhdy
      close(145)
      !
   endif   
   !
   theta360 = theta360*pi/180.0
   dtheta   = dtheta*pi/180.0
   !
   write(*,*)'Finished initializing SnapWave.'
   !
end subroutine initialize_snapwave_domain

subroutine read_interpolate_map_input(tmpstr,no_nodes, x, y, sferic, proj)
   !
   ! reads samples from file and interpolates them onto the mesh x,y
   !
   use interp, only: triintfast
   use m_alloc
   !
   implicit none
   !
   character*256, intent(in)                   :: tmpstr
   integer, intent(in)                         :: no_nodes
   real*8, dimension(no_nodes), intent(in)     :: x, y
   integer, intent(in)                         :: sferic
   real*4, dimension(no_nodes), intent(out)    :: proj
   real*8,    dimension(:),     allocatable    :: xsam, ysam, zsam
   real*8,    dimension(:),     allocatable    :: work
   real*4                                      :: t1, t2
   integer                                     :: nsam, k,jsferic, ino
   !
   allocate(work(no_nodes))
   !
   ! determine no of sample points
   !
   write(*,*)'   Reading samples from ',trim(tmpstr)
   open(unit = 500, file = trim(tmpstr))
   nsam = 0
   !
   do
     read(500,*,iostat=ino)
     if (ino/=0) exit
     nsam = nsam + 1
   enddo
   !
   rewind(500)
   write(*,*) nsam, 'samples read from ',trim(tmpstr)
   !
   ! allocate and read samples
   !
   call realloc(xsam,nsam,keepExisting=.false.,fill=0d0)
   call realloc(ysam,nsam,keepExisting=.false.,fill=0d0)
   call realloc(zsam,nsam,keepExisting=.false.,fill=0d0)
   do k=1,nsam
      read(500,*) xsam(k), ysam(k), zsam(k)  
   enddo
   close(500)
   !
   ! Interpolate samples onto mesh
   !
   write(*,*) '   Interpolation of samples on mesh started...'
   call timer(t1)
   jsferic=sferic    ! not superfluous, jsferic gets changed inside triintfast routines
   work=0d0
   call triintfast(xsam, ysam, zsam, nsam, 1 ,x ,y, work, no_nodes, jsferic, 0d0, 1.1d0)
   proj=real(work,kind(1.0))
   call timer(t2)
   write(*,*) '   ',nsam, ' samples interpolated on ', no_nodes, ' mesh nodes in ', t2-t1, ' s.'
   !
end subroutine read_interpolate_map_input
!   
subroutine find_upwind_neighbours_1dir(x,y,no_nodes,sferic,kp,np,dir,w,prev,ds)
   !
   ! finds the upwind neigbours for a specific direction, for example the mean wave direction or the wind direction (mapfields)
   !
   implicit none
   !
   real*8,  dimension(no_nodes),          intent(in)        :: x,y                 ! x, y coordinates of grid
   integer,                               intent(in)        :: no_nodes            ! length of x,y;
   integer,                               intent(in)        :: sferic              ! sferic (1) or cartesian (0) grid
   integer,                               intent(in)        :: np                  ! number of boundary points
   integer, dimension(np,no_nodes),       intent(in)        :: kp                  ! grid indices of surrounfing points per grid point
   real*4, dimension(no_nodes),           intent(in)        :: dir                 ! mapfield of direction (e.g. mean wave direction or wind direction)
   !
   real*4,  dimension(2, no_nodes) ,      intent(out)       :: w                   ! per grid point, weight of upwind points
   integer, dimension(2, no_nodes),       intent(out)       :: prev                ! per grid point, indices of upwind points
   real*4,  dimension(no_nodes),          intent(out)       :: ds                  ! upwind distance to intersection point 
   !
   real*8,  dimension(2)                              :: xsect,ysect,ww
   real*8                                             :: pi,dss,xi,yi
   integer                                            :: ind1,ind2
   integer                                            :: ip,nploc
   integer                                            :: k 
   real*8                                             :: circumf_eq=40075017.,circumf_pole=40007863.
   !
   ! Find upwind neighbours for the mean wave direction for each cell in an unstructured grid x,y (1d
   ! vectors)
   !
   pi = 4*atan(1.0)
   !
   do k = 1, no_nodes
      !
      call findloc(kp(:,k), np, 0, nploc)
      nploc = nploc - 1
      !
      if (kp(1,k)/=0) then
         !
         do ip = 1, nploc - 1
            !
            ind1 = kp(ip,k)
            ind2 = kp(ip+1,k)
            if (sferic==0) then
               xsect=[x(ind1), x(ind2)]
               ysect=[y(ind1), y(ind2)]
               call intersect_angle(x(k), y(k), dir(k) + pi, xsect, ysect, ww, dss, xi, yi)
            else
               xsect=[x(ind1)-x(k), x(ind2)-x(k)]*circumf_eq/360.0*cosd(y(k))
               ysect=[y(ind1)-y(k), y(ind2)-y(k)]*circumf_pole/360.0
               call intersect_angle(0.d0, 0.d0, dir(k) + pi, xsect, ysect, ww, dss, xi, yi)
            endif  
            !
            if (dss/=0) then
               w(1,k) = ww(1)
               w(2,k) = ww(2)
               ds(k)  = dss
               prev(1, k) = kp(ip  ,k  )
               prev(2, k) = kp(ip+1,k )
               exit
            endif
            !
         enddo
         !
         if (dss==0) then
            prev(1, k) = 1
            prev(2, k) = 1
            w(1, k)    = 0.0
            w(2, k)    = 0.0 
            ds(k)      = 0.0 ! very big number
         endif
         !
      endif
      !
   enddo
   !
end subroutine find_upwind_neighbours_1dir
!
subroutine find_upwind_neighbours(x,y,no_nodes,sferic,theta,ntheta,kp,np,w,prev,ds)
   !
   implicit none
   !
   integer,                               intent(in)        :: no_nodes,ntheta,np  ! length of x,y; length of theta; length of neighbouring points
   integer,                               intent(in)        :: sferic              ! sferic (1) or cartesian (0) grid
   real*8,  dimension(no_nodes),          intent(in)        :: x,y                 ! x, y coordinates of grid
   integer, dimension(np,no_nodes),       intent(in)        :: kp                  ! grid indices of surrounfing points per grid point
   real*8,  dimension(ntheta),            intent(in)        :: theta               ! array of wave angles
   real*8,  dimension(2,ntheta,no_nodes), intent(out)       :: w                   ! per grid point and direction, weight of upwind points
   integer, dimension(2,ntheta,no_nodes), intent(out)       :: prev                ! per grid point and direction, indices of upwind points
   real*8,  dimension(ntheta,no_nodes),   intent(out)       :: ds                  ! upwind distance to intersection point for each direction
   !
   real*8,  dimension(2)                              :: xsect,ysect,ww
   real*8                                             :: pi,dss,xi,yi
   integer                                            :: ind1,ind2
   integer                                            :: ip,nploc
   integer                                            :: k, itheta
   real*8                                             :: circumf_eq=40075017.,circumf_pole=40007863.
   real*8, parameter                                  :: epsdist=0.001
   !
   ! Find upwind neighbours for each cell in an unstructured grid x,y (1d
   ! vectors) given vector of directions theta
   !
   pi = 4.0*atan(1.0)
   !
   do k = 1, no_nodes
      call findloc(kp(:,k), np, 0, nploc)
      nploc = nploc - 1
      if (kp(1,k)/=0) then
         do itheta = 1, ntheta
            do ip = 1, nploc - 1
               ind1 = kp(ip,k)
               ind2 = kp(ip+1,k)
               if (sferic==0) then
                  xsect=[x(ind1), x(ind2)]
                  ysect=[y(ind1), y(ind2)]
                  call intersect_angle(x(k), y(k), theta(itheta) + pi, xsect, ysect, ww, dss, xi, yi)
               else
                  xsect=[x(ind1)-x(k), x(ind2)-x(k)]*circumf_eq/360.0*cosd(y(k))
                  ysect=[y(ind1)-y(k), y(ind2)-y(k)]*circumf_pole/360.0
                  call intersect_angle(0.d0, 0.d0, theta(itheta) + pi, xsect, ysect, ww, dss, xi, yi)
               endif               
               if (dss>epsdist) then
                  w(1, itheta,k) = ww(1)
                  w(2, itheta,k) = ww(2)
                  ds(itheta, k) = dss
                  prev(1, itheta, k) = kp(ip  ,k  )
                  prev(2, itheta, k) = kp(ip+1,k )
                  exit
               endif
            enddo
            if (dss<=epsdist) then
!               write(*,*)k,x(k),y(k)
               prev(1, itheta, k) = 1
               prev(2, itheta, k) = 1
               w(1, itheta, k)    = 0.0
               w(2, itheta, k)    = 0.0
            endif
         enddo
      endif
   enddo

end subroutine find_upwind_neighbours
   
subroutine intersect_angle(x0, y0, phi, x, y, W, ds, xi, yi)
   !
   implicit none
   !
   real*8, intent(in)               :: x0, y0, phi
   real*8, dimension(2),intent(in)  :: x, y
   real*8, dimension(2),intent(out) :: W
   real*8, intent(out)              :: ds, xi, yi
   real*8                           :: eps, m, a, b, n, L, d1, d2, err
   !
   eps = 1.0e-2
   !
   if (abs(x(2) - x(1))>eps) then
      m  = (y(2) - y(1))/(x(2) - x(1))
      a  = y(1) - m*x(1)
      n  = tan(phi)
      b  = y0 - n*x0
      xi = (b - a)/(m - n)
      yi = a + m*xi
   else
      yi = (x(1) - x0)*tan(phi) + y0
      xi = x(1)
   endif
   !
   L   = hypot(x(2) - x(1), y(2) - y(1))
   d1  = hypot(xi - x(1), yi - y(1))
   d2  = hypot(xi - x(2), yi - y(2))
   ds  = hypot(xi - x0, yi - y0)
   err = hypot((x0 - xi) + ds*cos(phi), (y0 - yi) + ds*sin(phi))
   if (abs(L - d1 - d2)<eps .and. err<eps) then
      W(1) = d2/L
      W(2) = d1/L
   else
      W(1) = 0.0
      W(2) = 0.0
      ds   = 0.0
   endif
   !
end subroutine intersect_angle
   
   

subroutine fm_surrounding_points(xn,yn,zn,no_nodes,sferic,face_nodes,no_faces,kp,np,dhdx,dhdy)
   !
   implicit none
   !
   real*8,  dimension(no_nodes),                 intent(in)  :: xn,yn              ! coordinates of network nodes
   real*4,  dimension(no_nodes),                 intent(in)  :: zn                 ! coordinates of network nodes
   integer,                                      intent(in)  :: no_nodes           ! number of network nodes
   integer,                                      intent(in)  :: sferic
   integer, dimension(4,no_faces),               intent(in)  :: face_nodes         ! node numbers connected to each cell
   integer,                                      intent(in)  :: no_faces           ! number of cells
   integer, dimension(np,no_nodes),              intent(out) :: kp                 ! sorted surrounding node numbers for each node
   integer,                                      intent(in)  :: np                 ! max. number of surrounding nodes
   real*4,  dimension(no_nodes),                 intent(out) :: dhdx,dhdy          ! bed slopes in each node
   !
   ! Local variables
   !
   integer, dimension(:),   allocatable         :: no_connected_cells
   integer, dimension(:,:), allocatable         :: connected_cells
   integer, dimension(4)                        :: kpts, edge
   integer, dimension(np)                       :: surr_points
   integer, dimension(np)                       :: surr_pts
   real*4,  dimension(np)                       :: xp,yp,zp              ! x,y,z of sorted surrounding nodes for each node
   real*8                                       :: circumf_eq=40075017.,circumf_pole=40007863.
   real*4                                       :: dzdx,dzdy
   real*8                                       :: sy2, syz, sx2, sxz
   !
   integer                                      :: k, inode, knode, kcell, kn, ip, isp, j, jj, next, isp2
   !
   allocate(no_connected_cells(no_nodes))
   allocate(connected_cells(12,no_nodes))
   no_connected_cells=0
   connected_cells=0
   kp=0
   !
   do k=1,no_faces
      do inode=1,4
         knode=face_nodes(inode,k)
         if (knode /= -999) then
            no_connected_cells(knode)=no_connected_cells(knode)+1
            connected_cells(no_connected_cells(knode),knode)=k
         endif
      enddo
   enddo
   !
   do kn=1,no_nodes                      ! for each node
      surr_points=0
      surr_pts=0
      isp=1
      do kcell=1,no_connected_cells(kn)       ! for each cell connected to node
         k=connected_cells(kcell,kn)          ! get cell number
         kpts=face_nodes(:,k)                 ! all nodes in that cell
         if (kpts(4)==-999) then
            jj=0
            do j=1,3
               if (kpts(j)/=kn) then
                  jj=jj+1
                  edge(jj)=kpts(j)
               endif
            enddo
            surr_points(isp:isp+1)=edge(1:2)
            isp=isp+2
         else
            ip=minloc(abs(kpts-kn),1)
            if (ip==1) then
               edge=[kpts(2),kpts(3),kpts(3),kpts(4)]
            elseif (ip==2) then
               edge=[kpts(3),kpts(4),kpts(4),kpts(1)]
            elseif (ip==3) then
               edge=[kpts(4),kpts(1),kpts(1),kpts(2)]
            else
               edge=[kpts(1),kpts(2),kpts(2),kpts(3)]
            endif
            surr_points(isp:isp+3)=edge(1:4)
            isp=isp+4
         endif
      enddo
      isp=isp-1  ! number of surrounding points
      !   now connect the edges
      if (isp>=2) then
         surr_pts(1:2)=surr_points(1:2)
         surr_points=[surr_points(3:),0,0]
         isp2=2
         isp=isp-2
         do while (isp>=2)
            call findloc(surr_points,isp,surr_pts(isp2),next)
            if (next/=-1) then
               if (mod(next,2)==1) then
                  surr_pts(1:isp2+1)=[surr_pts(1:isp2),surr_points(next+1)]
                  surr_points=[surr_points(1:next-1),surr_points(next+2:),0,0]
               else
                  surr_pts(1:isp2+1)=[surr_pts(1:isp2),surr_points(next-1)]
                  surr_points=[surr_points(1:next-2),surr_points(next+1:),0,0]
               endif
            else
               call findloc(surr_points,isp,surr_pts(1),next)
               if (next/=-1) then
                  if (mod(next,2)==1) then
                     surr_pts(1:isp2+1)=[surr_points(next+1),surr_pts(1:isp2)]
                     surr_points=[surr_points(1:next-1),surr_points(next+2:),0,0]
                  else
                     surr_pts(1:isp2+1)=[surr_points(next-1),surr_pts(1:isp2)]
                     surr_points=[surr_points(1:next-2),surr_points(next+1:),0,0]
                  endif
               else
                  surr_pts=0
               endif
            endif
            isp2=isp2+1
            isp=isp-2
         end do
      else
         surr_pts=0
      endif
      kp(:,kn)=surr_pts
      sxz=0.0
      syz=0.0
      sx2=0.0
      sy2=0.0
      if (sum(surr_pts)>0) then
         do j=1,isp2-1
            xp(j)=xn(surr_pts(j))-xn(kn)
            yp(j)=yn(surr_pts(j))-yn(kn)
            zp(j)=zn(surr_pts(j))      
            if (sferic==1) then
               xp(j)=xp(j)*circumf_eq/360.0*cosd(yn(kn))
               yp(j)=yp(j)*circumf_pole/360.0
            endif
         enddo
         call plane_fit(xp,yp,zp,isp2-1, dzdx,dzdy)
         dhdx(kn)=-dzdx
         dhdy(kn)=-dzdy
      else
          dhdx(kn)=0.
          dhdy(kn)=0.
      endif
      !write(*,*)'fmgsp 2',kn, no_nodes,isp,isp2

   enddo
   !
end subroutine

subroutine slope_plane_fit (xn,yn,fn,no_nodes,sferic,kp,np,dfdx,dfdy)
   !
   implicit none
   !
   real*8,  dimension(no_nodes),                 intent(in)  :: xn,yn              ! coordinates of network nodes
   real*4,  dimension(no_nodes),                 intent(in)  :: fn                 ! variable f for which to determine gradients
   integer,                                      intent(in)  :: no_nodes           ! number of network nodes
   integer,                                      intent(in)  :: sferic
   integer, dimension(np,no_nodes),              intent(in)  :: kp                 ! sorted surrounding node numbers for each node
   integer,                                      intent(in)  :: np                 ! max. number of surrounding nodes
   real*4,  dimension(no_nodes),                 intent(out) :: dfdx,dfdy          ! gradients of f in each node

   real*8                                       :: circumf_eq=40075017.,circumf_pole=40007863.
   integer, dimension(np)                       :: surr_pts
   real*4,  dimension(np)                       :: xp,yp,fp
   integer                                      :: kn,indx,j

   do kn=1,no_nodes
      surr_pts=kp(:,kn)
      if (sum(surr_pts)>0) then
         call findloc(surr_pts,np,0,indx)
         do j=1,indx-1
            xp(j)=xn(surr_pts(j))-xn(kn)
            yp(j)=yn(surr_pts(j))-yn(kn)
            fp(j)=fn(surr_pts(j))      
            if (sferic==1) then
               xp(j)=xp(j)*circumf_eq/360.0*cosd(yn(kn))
               yp(j)=yp(j)*circumf_pole/360.0
            endif
         enddo
         call plane_fit(xp,yp,fp,indx-1, dfdx(kn),dfdy(kn))
      endif
   enddo
end subroutine slope_plane_fit

subroutine findloc(a,n,b,indx)
   !
   implicit none
   !
   integer, dimension(n), intent(in) :: a
   integer, intent(in)               :: n,b
   integer, intent(out)              :: indx
   integer                           :: i
   !
   indx=-1
   do i=1,n
      if (a(i)==b) then
         indx=i
         return
      endif
   enddo
end subroutine findloc
    

subroutine linear_interp(x, y, n, xx, yy)
   !
   ! NAME
   !    linear_interp
   ! SYNOPSIS
   !    Interpolate linearly into an array Y given the value of X.
   !    Return both the interpolated value and the position in the
   !    array.
   !
   ! ARGUMENTS
   !    - X: independent array (sorted in ascending order)
   !    - Y: array whose values are a function of X
   !    - XX: specified value, interpolating in first array
   !    - YY: interpolation result
   !    - INDINT: (optional) index in array (x(indint) <= xx <= x(indint+1))
   !
   ! SOURCE
   !
      integer,              intent(in) :: n
      real*4, dimension(n), intent(in) :: x
      real*4, dimension(n), intent(in) :: y
      real*4, intent(in)               :: xx
      real*4, intent(out)              :: yy
      !****
      !
      ! CODE: linear interpolation
      !
      !
      !
      real*4           :: a,  b, dyy
      integer          :: j, indint

      yy = 0.0
      indint=0
      if (n.le.0) return
      !
      ! *** N GREATER THAN 0
      !
      if (n.eq.1) then
         yy = y(1)
         return
      endif

      call binary_search( x, n, xx, j )

      if ( j .le. 0 ) then
         yy = y(1)
      elseif ( j .ge. n ) then
         yy = y(n)
      else
         a = x (j+1)
         b = x (j)
         if ( a .eq. b ) then
            dyy = 0.0
         else
            dyy = (y(j+1) - y(j)) / (a - b)
         endif
         yy = y(j) + (xx - x(j)) * dyy
      endif

!      if ( present(indint) ) then
!         indint = j
!      endif

      return

end subroutine linear_interp

subroutine binary_search(xx, n, x, j)
   !****f* Interpolation/binary_search
   !
   ! NAME
   !    binary_search
   ! SYNOPSIS
   !    Perform a binary search in an ordered real array
   !    to find the largest entry equal or lower than a given value:
   !    Given an array XX of length N, given value X, return a value J
   !    such that X is between XX(J) en XX (J+1)
   !
   !    XX must be monotonic, either decreasing or increasing
   !    J=0 or J=N indicates X is out of range.
   !
   ! ARGUMENTS
   !    - XX: ordered array of values
   !    - X: value to be found
   !    - J: index such that X is between XX(J) and XX(J+1)
   !
   ! SOURCE
   !
   integer,              intent(in) :: N
   real*4, dimension(N), intent(in) :: xx
   real*4, intent(in)               :: x
   integer, intent(out)             :: j
   !****
   !
   ! CODE: binary search in (real) arrays
   !
   ! Requirement:
   !    Parameter wp set to the proper kind
   !
   ! Subroutine from 'Numerical recipes' Fortran  edition.
   ! Given an array XX of length N, given value X, return a value J
   ! such that X is between XX(J) en XX (J+1)
   ! XX must be monotonic, either decreasing or increasin
   ! J=0 or J=N indicates X is out of range.


   !
   ! Local variables
   !
   integer   jl, ju, jm
   logical   l1, l2

   jl = 0
   ju = n+1
10 if (ju-jl .gt. 1) then
      jm = (ju+jl)/2
      l1 = xx(n) .gt. xx(1)
      l2 = x .gt. xx(jm)
      if ( (l1.and.l2) .or. (.not. (l1 .or. l2)) ) then
         jl = jm
      else
         ju = jm
      endif
      goto 10
   endif

   j = jl

   return

end subroutine binary_search

!subroutine boundaries(x,y,no_nodes,xb,yb,nb,tol,bndpts,nobndpts,bndindx,bndweight)
!   !
!   implicit none
!   !
!   real*4,  dimension(no_nodes),        intent(in)    :: x,y
!   integer,                             intent(in)    :: no_nodes
!   real*4,  dimension(nb),              intent(in)    :: xb,yb
!   integer,                             intent(in)    :: nb
!   real*4,                              intent(in)    :: tol
!   integer                                            :: nobndpts
!   integer, dimension(nobndpts)                       :: bndpts
!   integer, dimension(2,nobndpts)                     :: bndindx
!   real*4,  dimension(2,nobndpts)                     :: bndweight
!   !
!   integer                               :: ib,k,ibnd
!   real*4                                :: alpha, cosa,sina, x1, y1, xend
!   !
!   ibnd=0
!   do ib=1,nb-1
!    if (xb(ib).ne.-999.and. xb(ib+1).ne.-999) then
!      alpha=atan2(yb(ib+1)-yb(ib),xb(ib+1)-xb(ib))
!      cosa=cos(alpha)
!      sina=sin(alpha)
!      xend=(xb(ib+1)-xb(ib))*cosa+(yb(ib+1)-yb(ib))*sina
!      do k=1,no_nodes
!         x1= (x(k)-xb(ib))*cosa+(y(k)-yb(ib))*sina
!         y1=-(x(k)-xb(ib))*sina+(y(k)-yb(ib))*cosa
!         if (x1>=0. .and. x1<=xend) then
!            if (abs(y1)<tol) then
!               ibnd=ibnd+1
!               bndpts(ibnd)=k
!               bndindx(1,ibnd)=ib
!               bndindx(2,ibnd)=ib+1
!               bndweight(2,ibnd)=x1/xend
!               bndweight(1,ibnd)=1.-bndweight(2,ibnd)
!            endif
!         endif
!      enddo
!    endif
!   enddo
!   nobndpts=ibnd
!   
!end subroutine boundaries            
  
subroutine neuboundaries(x,y,no_nodes,xneu,yneu,n_neu,tol,neumannconnected)
   !
   implicit none
   !
   integer, intent(in)                        :: no_nodes
   integer, intent(in)                        :: n_neu
   real*8, dimension(no_nodes), intent(in)    :: x,y
   real*8, dimension(n_neu), intent(in)       :: xneu,yneu
   real*4, intent(in)                         :: tol
   integer, dimension(no_nodes), intent(out)  :: neumannconnected
   !
   integer                                    :: ib,k,kmin, k2
   real*8                                     :: alpha, cosa,sina, distmin, x1,y1,x2,y2, xend
   !
   neumannconnected=0
   do ib=1,n_neu-1
      if (xneu(ib).ne.-999.and.xneu(ib+1).ne.-999) then 
         alpha=atan2(yneu(ib+1)-yneu(ib),xneu(ib+1)-xneu(ib))
         cosa=cos(alpha)
         sina=sin(alpha)
         xend=(xneu(ib+1)-xneu(ib))*cosa+(yneu(ib+1)-yneu(ib))*sina
         do k=1,no_nodes
            x1= (x(k)-xneu(ib))*cosa+(y(k)-yneu(ib))*sina
            y1=-(x(k)-xneu(ib))*sina+(y(k)-yneu(ib))*cosa
            if (x1>=0.d0 .and. x1<=xend) then
               if (abs(y1)<tol) then
                  ! point k is on the neumann boundary
                  distmin=1d10
                  kmin=0
                  do k2=1,no_nodes
                     x2= (x(k2)-xneu(ib))*cosa+(y(k2)-yneu(ib))*sina
                     y2=-(x(k2)-xneu(ib))*sina+(y(k2)-yneu(ib))*cosa
                     if (abs(x2-x1)<tol .and. (k2.ne.k)) then
                        if (abs(y2-y1)<distmin) then
                           kmin=k2
                           distmin=abs(y2-y1)
                        endif
                     endif
                  enddo
                  if (kmin>0) then
                     neumannconnected(kmin)=k
                     write(*,*)kmin,k
                  endif
               endif
            endif
         enddo
      endif
   enddo
   !
end subroutine neuboundaries
  !
subroutine read_neumann_boundary()
   !
   ! Reads neu file
   !
   use snapwave_data
   !
   implicit none
   !
   integer n, stat
   !
   real*4 dummy
   !
   ! Read polylines
   if (neumannfile(1:4) /= 'none') then    ! Normal ascii input files
      write(*,*)neumannfile(1:4)

      write(*,*)'   Reading neumann boundary locations ...'
      !
      open(500, file=trim(neumannfile)) !as in bwvfile of SFINCS
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         n_neu = n_neu + 1
      enddo
      rewind(500)
      if (n_neu==0) then
            write(*,*)'  neumann boundary enclosure file ',trim(neumannfile),' is empty or does not exist'
      else           
         allocate(x_neu(n_neu))
         allocate(y_neu(n_neu))
         do n = 1, n_neu
            read(500,*)x_neu(n),y_neu(n)
         enddo
         close(500) 
      endif
   else
      n_neu=0 ! There is no neumann boundary specified
   endif
      !
end subroutine read_neumann_boundary
   
subroutine read_boundary_enclosure()
   !
   ! Reads enc file
   !
   use snapwave_data
   !
   implicit none
   !
   integer n, stat
   !
   real*4 dummy
   !
   ! Read wave boundaries
   !
   if (encfile(1:4) /= 'none') then    
      write(*,*)'   Reading wave boundary enclosure ...'
      open(500, file=trim(encfile)) !as in bwvfile of SFINCS
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         n_bndenc = n_bndenc+1
      enddo
      rewind(500)
      if (n_bndenc==0) then
          write(*,*)'  boundary enclosure file ',trim(encfile),' is empty or does not exist'
      endif
      allocate(x_bndenc(n_bndenc))
      allocate(y_bndenc(n_bndenc))
      do n = 1, n_bndenc
         read(500,*)x_bndenc(n),y_bndenc(n)
      enddo
      close(500)
      !
   else
      n_bndenc = 0 
   endif
   !
end subroutine read_boundary_enclosure
   
subroutine timer(t)
   !
   implicit none
   !
   real*4,intent(out)               :: t
   !
   integer*4                        :: count,count_rate,count_max
   !
   call system_clock (count,count_rate,count_max)
   t = dble(count)/count_rate
   !
end subroutine timer


subroutine veggie_init()
	use snapwave_data
    use snapwave_input, only: read_int_input, read_real_input, read_char_input
    use snapwave_ncoutput, only: nc_check_err
    use netcdf  
    !
    implicit none    
    !
    real*4,  dimension(:),       allocatable    :: xveg, yveg                  ! x,y coordinates of vegetation grid (single precision)
    !
    integer                                     :: i,j,k,fid,ier,is,m,ind
    integer                                     :: ierror
    integer                                     :: idfile, iddim_nsec
    integer                                     :: idvar_nsec, idvar_veg_ah, idvar_veg_bstems, idvar_veg_Nstems, idvar_veg_Cd
    character*2                                 :: ext 
    character(NF90_MAX_NAME)                    :: string
    
    write(*,*)'Initializing vegetation input settings '
	!
    no_secvegmax = 3
    ! Get maximum number of vegetation sections within model domain - needed to set size of Cd, ah, bv and Nv matrix
    no_secveg = MIN(no_secvegmax, no_secveg)
    !
    j=index(vegmapfile,'.')      
	ext = vegmapfile(j+1:j+2)
    if (ext == 'nc') then
       !
       ! open vegmapfile
       !
       ierror = nf90_open(vegmapfile, NF90_NOWRITE, idfile)       
       call nc_check_err(ierror, "opening file", vegmapfile)
       
       ierror = nf90_inq_dimid(idfile, 'nsec', iddim_nsec) 
       call nc_check_err(ierror, "inq_dimid nsec", vegmapfile)
       
       ierror = nf90_inquire_dimension(idfile, iddim_nsec, string, no_secveg)
       call nc_check_err(ierror, "inq_dim nsec", vegmapfile)

       ! check if variables are available on the netcdf
       ierror = nf90_inq_varid(idfile, 'veg_ah'    , idvar_veg_ah      )                 
       call nc_check_err(ierror, "inq_varid veg_ah", vegmapfile)
       ierror = nf90_inq_varid(idfile, 'veg_bstems', idvar_veg_bstems  )                 
       call nc_check_err(ierror, "inq_varid veg_bstems", vegmapfile)
       ierror = nf90_inq_varid(idfile, 'veg_Nstems', idvar_veg_Nstems  )                 
       call nc_check_err(ierror, "inq_varid veg_Nstems", vegmapfile)
       ierror = nf90_inq_varid(idfile, 'veg_Cd'    , idvar_veg_Cd      )                 
       call nc_check_err(ierror, "inq_varid veg_Cd", vegmapfile)
       
	   ! Allocate vegetation specific variables
	   allocate(veg_Cd(no_nodes, no_secveg))
	   allocate(veg_ah(no_nodes, no_secveg))
	   allocate(veg_bstems(no_nodes,  no_secveg))
	   allocate(veg_Nstems(no_nodes,  no_secveg))
       
       ! read variables from file
       ierror = nf90_get_var(idfile, idvar_veg_ah         , veg_ah             , start=(/ 1 /)   , count=(/ no_nodes, no_secveg /) )
       call nc_check_err(ierror, "get_var veg_ah", vegmapfile)
       ierror = nf90_get_var(idfile, idvar_veg_bstems     , veg_bstems         , start=(/ 1 /)   , count=(/ no_nodes, no_secveg /) )
       call nc_check_err(ierror, "get_var veg_bstems", vegmapfile)
       ierror = nf90_get_var(idfile, idvar_veg_Nstems     , veg_Nstems         , start=(/ 1 /)   , count=(/ no_nodes,no_secveg /) )
       call nc_check_err(ierror, "get_var veg_Nstems", vegmapfile) ! should it be read as an integer?
       ierror = nf90_get_var(idfile, idvar_veg_Cd         , veg_Cd         , start=(/ 1 /)   , count=(/ no_nodes, no_secveg /) )
       call nc_check_err(ierror, "get_var veg_Cd", vegmapfile)
        
    else 
        write(*,*)'COULD NOT READ VEG INPUT FILE, PLEASE SUPPLY NETCDF FILE'
    endif
        
    end subroutine veggie_init

  integer function count_lines(u,f)
      !
      ! returns number of lines in file f
      ! to be called by all, xmaster will count
      ! result is broadcasted to all
      !
      !use filefunctions
      character(len=*), intent(in) :: f
      integer lines,u,ierr
      !fid = create_new_fid()
      open(unit=u,file=f)
      lines = 0
      ierr = 0
      do while(ierr == 0)
          read(u,*,iostat=ierr)
          if (ierr == 0) lines = lines+1
      enddo
      close(u)
      count_lines = lines

   end function count_lines

   subroutine plane_fit(x,y,z,n, dzdx,dzdy)
    implicit none
    integer :: n                             ! Number of points
    real(4), dimension(n) :: x, y, z         ! Coordinates
    real(4), dimension(n, 3) :: XX           ! Designmatrix
    real(4), dimension(3, 3) :: XT_X         ! X-transpose * X
    real(4), dimension(3) :: XT_z, coeffs    ! X-transpose * z and solution [a, b, c]
    real(4) :: mean_slope                    ! Mean slope
    real(4) :: dzdx,dzdy                     ! Slope in x and y direction
    integer :: i, j, k

    ! Build the designmatrix XX
    do i = 1, n
        XX(i, 1) = x(i)
        XX(i, 2) = y(i)
        XX(i, 3) = 1.0d0
    end do

    ! Compute XT_X = XX^T * XX
    XT_X = 0.0d0
    do i = 1, 3
        do j = 1, 3
            do k = 1, n
                XT_X(i, j) = XT_X(i, j) + XX(k, i) * XX(k, j)
            end do
        end do
    end do

    ! Compute XT_z = XX^T * z
    XT_z = 0.0d0
    do i = 1, 3
        do j = 1, n
            XT_z(i) = XT_z(i) + XX(j, i) * z(j)
        end do
    end do

    ! Solve XT_X * coeffs = XT_z for coeffs (a, b, c)
    call solve_linear_system(XT_X, XT_z, coeffs)

    ! Coefficients of the plane
    dzdx = coeffs(1)
    dzdy = coeffs(2)

    ! Mean slope
    mean_slope = sqrt(coeffs(1)**2 + coeffs(2)**2)

   end subroutine plane_fit
   
       subroutine solve_linear_system(A, b, x)
        real(4), dimension(3, 3), intent(in) :: A
        real(4), dimension(3), intent(in) :: b
        real(4), dimension(3), intent(out) :: x
        real(4), dimension(3, 3) :: A_inv
        integer :: i, j, k
        real(4) :: factor

        ! Kopieer matrix A
        A_inv = A

        ! Initialiseer x met b
        x = b

        ! Gaussiaanse eliminatie
        do i = 1, 3
            factor = A_inv(i, i)
            if (abs(factor) >= 1.d-10) then
               A_inv(i, :) = A_inv(i, :) / factor
               x(i) = x(i) / factor
            endif

            do j = 1, 3
                if (i /= j) then
                    factor = A_inv(j, i)
                    A_inv(j, :) = A_inv(j, :) - factor * A_inv(i, :)
                    x(j) = x(j) - factor * x(i)
                end if
            end do
        end do
       end subroutine solve_linear_system

end module