module interp
   use mathconsts
   use physicalconsts
   use m_ec_triangle
   use kdtree2Factory
   use m_alloc
   
   implicit none
   
   private
   public linear_interp, interp_in_cyclic_function, interp_using_trapez_rule,ipon
   public grmap2, grmap_sg, mkmap, grmap, linear_interp_2d, make_map, make_map_fm, mkmap_step, trapezoidal
   public triintfast
   save
contains
   pure subroutine linear_interp_2d(X,nx,Y,ny,Z,xx,yy,zz,method,exception)

   implicit none
   ! input/output
   integer,intent(in)                    :: nx,ny
   real*8,dimension(nx),intent(in)       :: X
   real*8,dimension(ny),intent(in)       :: Y
   real*8,dimension(nx,ny),intent(in)    :: Z
   real*8,intent(in)                     :: xx,yy
   real*8,intent(out)                    :: zz
   character(len=*),intent(in)           :: method
   real*8,intent(in)                     :: exception
   ! internal
   integer,dimension(4)                  :: ind
   real*8,dimension(2)                   :: yint
   real*8                                :: modx,mody,disx,disy
   logical                               :: interpX,interpY

   ! does the interpolation point fall within the data?
   if (xx>=minval(X) .and. xx<=maxval(X)) then
      interpX = .true.
   else
      interpX = .false.
   endif
   if (yy>=minval(Y) .and. yy<=maxval(Y)) then
      interpY = .true.
   else
      interpY = .false.
   endif

   if (interpX .and. interpY) then
      ! find rank position xx in X direction
      ind(1) = minloc(X,1,X>=xx)
      ind(2) = maxloc(X,1,X<=xx)
      ! find rank position yy in Y direction
      ind(3) = minloc(Y,1,Y>=yy)
      ind(4) = maxloc(Y,1,Y<=yy)
      ! distance between X points and Y points
      disx = X(ind(1))-X(ind(2))
      disy = Y(ind(3))-Y(ind(4))
      ! relative position of (xx,yy) on disx,disy
      if (disx>0.d0) then
         modx = (xx-X(ind(2)))/disx
      else
         modx = 0.d0   ! xx corresponds exactly to X point
      endif
      if (disy>0.d0) then
         mody = (yy-Y(ind(4)))/disy
      else
         mody = 0.d0   ! yy corresponds exactly to Y point
      endif
      ! interpolate the correct Y value, based on two nearest X intersects
      ! if disy==0 then only single interpolation needed. This could also be done for X,
      ! but since this is used by waveparams and the interp angles (Y) are more often equal
      ! this is probably faster.
      if (disy>0.d0) then
         yint(1) = (1.d0-modx)*Z(ind(2),ind(3))+modx*Z(ind(1),ind(3))
         yint(2) = (1.d0-modx)*Z(ind(2),ind(4))+modx*Z(ind(1),ind(4))
         zz = (1.d0-mody)*yint(2)+mody*yint(1)
      else
         zz = (1.d0-modx)*Z(ind(2),ind(3))+modx*Z(ind(1),ind(3))
      endif
   else
      select case (method)
      case ('interp')
         zz = exception
      case ('extendclosest')
         ! find closest X point
         ind(1) = minloc(abs(X-xx),1)
         ! find closest Y point
         ind(3) = minloc(abs(Y-yy),1)
         ! external value given that of closest point
         zz = Z(ind(1),ind(3))
         case default
         zz = exception
      endselect
   endif

   end subroutine linear_interp_2d

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
   pure subroutine linear_interp(x, y, n, xx, yy, indint)
   
   implicit none
   
   integer,              intent(in) :: n
   real*8, dimension(n), intent(in) :: x
   real*8, dimension(n), intent(in) :: y
   real*8, intent(in)               :: xx
   real*8, intent(out)              :: yy
   integer, intent(out), optional   :: indint
   !****
   !
   ! CODE: linear interpolation
   !
   !
   !
   real*8           :: a,  b, dyy
   integer          :: j

   yy = 0.0d0
   if ( present(indint) ) then
      indint = 0
   endif

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
         dyy = 0.0d0
      else
         dyy = (y(j+1) - y(j)) / (a - b)
      endif
      yy = y(j) + (xx - x(j)) * dyy
   endif

   if ( present(indint) ) then
      indint = j
   endif

   return

   end subroutine linear_interp

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
   pure subroutine binary_search(xx, n, x, j)
   
   implicit none
   
   integer,              intent(in) :: N
   real*8, dimension(N), intent(in) :: xx
   real*8, intent(in)               :: x
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

   subroutine mkmap(code      ,x1        ,y1        ,m1        ,n1        , &
      & x2        ,y2        ,n2        ,xs        ,ys        , &
      & nrx       ,nry       ,iflag     ,nrin      ,w         , &
      & iref      ,iprint    ,covered   ,xymiss)
   !!--description-----------------------------------------------------------------
   !
   !
   !     SUBROUTINE MKMAP
   !     Interpolation of curvilinear, numerically ordered grid (grid1)
   !     to random points, with weighting points (grid2).
   !
   !     J.A. Roelvink
   !     Deltares
   !     24-2-1992 (MKMAP
   !
   !     Given: numerically ordered grid M1*N1
   !     with coordinates X1 (1:M1,1:N1)
   !                  and Y1 (1:M1,1:N1)
   !
   !     Also given: random points X2(1:N2)
   !                           and Y2(1:N2)
   !
   !     To be determined:  weighting factors and pointers for bilinear interpolation
   !     Weighting factors and locations of the requested (random) points in the
   !     ordered grid are saved in resp.
   !     W(1:4,1:N2) and Iref(1:4,1:N2)
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                    , intent(in)  :: iprint
   integer                    , intent(in)  :: m1
   integer                    , intent(in)  :: n1
   integer                                  :: n2
   integer , dimension(4 , n2), intent(out) :: iref
   integer , dimension(m1, n1), intent(in)  :: code
   integer , dimension(n2)                  :: covered !  0: target point is not covered by source grid (default)
   !  1: target point is covered by   valid points of source grid
   ! -1: target point is covered by invalid points of source grid
   integer , dimension(n2)                  :: iflag
   integer , dimension(n2)                  :: nrin
   integer , dimension(n2)                  :: nrx
   integer , dimension(n2)                  :: nry
   real                       , intent(in)  :: xymiss  ! missing value in grid 1
   real*8  , dimension(4 , n2)              :: w       ! Remains single precision!
   real*8  , dimension(m1, n1), intent(in)  :: x1
   real*8  , dimension(m1, n1), intent(in)  :: y1
   real*8  , dimension(n2)                  :: x2
   real*8  , dimension(n2)                  :: xs
   real*8  , dimension(n2)                  :: y2
   real*8  , dimension(n2)                  :: ys
   !
   ! Local variables
   !
   integer                           :: i
   integer                           :: i1
   integer                           :: i2
   integer                           :: ier
   integer                           :: iin
   integer                           :: inout
   integer                           :: ip
   integer                           :: j1
   integer                           :: lomaxx
   integer                           :: lominx
   integer                           :: lomnx
   integer                           :: lomny
   integer                           :: m1max
   integer                           :: m1min
   integer                           :: n1max
   integer                           :: n1min
   integer                           :: nin
   real*8                           :: eps
   real*8                            :: xpmax
   real*8                            :: xpmean
   real*8                            :: xpmin
   real*8                            :: ypmax
   real*8                            :: ypmean
   real*8                            :: ypmin
   real*8  , dimension(5)            :: xp
   real*8  , dimension(5)            :: yp
   real*8  , dimension(4)            :: hbuff
   !
   !! executable statements -------------------------------------------------------
   !
   eps = 0.00001
   !
   !     Initialise tables
   !
   if (iprint==1) write (*, *) 'in mkmap', m1, n1, n2
   lomnx = 1
   lomny = 1
   m1min = 1
   m1max = m1
   n1min = 1
   n1max = n1
   nrx=0
   nry=0
   iflag=0
   nrin=0
   xs=0.d0
   ys=0.d0
   iref=0
   w=0.d0
   !
   ! Sort X2 en Y2
   !
   call sort(n2        ,x2        ,xs        ,nrx       )
   call sort(n2        ,y2        ,ys        ,nry       )
   !
   ! Loop over all cels of grid1
   !
   !
   do j1 = n1min, n1max - 1
      do i1 = m1min, m1max - 1
         !
         ! Cell definition
         !
         xp(1) = x1(i1, j1)
         xp(2) = x1(i1 + 1, j1)
         xp(3) = x1(i1 + 1, j1 + 1)
         xp(4) = x1(i1, j1 + 1)
         yp(1) = y1(i1, j1)
         yp(2) = y1(i1 + 1, j1)
         yp(3) = y1(i1 + 1, j1 + 1)
         yp(4) = y1(i1, j1 + 1)
         if (     (xp(1) > xymiss-eps .and. xp(1) < xymiss+eps) &
            & .or. (xp(2) > xymiss-eps .and. xp(2) < xymiss+eps) &
            & .or. (xp(3) > xymiss-eps .and. xp(3) < xymiss+eps) &
            & .or. (xp(4) > xymiss-eps .and. xp(4) < xymiss+eps) &
            & .or. (yp(1) > xymiss-eps .and. yp(1) < xymiss+eps) &
            & .or. (yp(2) > xymiss-eps .and. yp(2) < xymiss+eps) &
            & .or. (yp(3) > xymiss-eps .and. yp(3) < xymiss+eps) &
            & .or. (yp(4) > xymiss-eps .and. yp(4) < xymiss+eps) ) cycle
         !
         ! Determine minimum and maximum X and Y of the cell
         !
         xpmin =  1.e10
         xpmax = -1.e10
         ypmin =  1.e10
         ypmax = -1.e10
         do ip = 1, 4
            xpmin = min(xp(ip), xpmin)
            xpmax = max(xp(ip), xpmax)
            ypmin = min(yp(ip), ypmin)
            ypmax = max(yp(ip), ypmax)
         enddo
         xpmean = .5*(xpmin + xpmax)
         ypmean = .5*(ypmin + ypmax)
         !
         ! First selection of points of grid2 that can be located in the cell
         !
         ! Find centre of the cell in tables Xs and Ys
         !
         call hunt(xs        ,n2        ,xpmean    ,lomnx     )
         call hunt(ys        ,n2        ,ypmean    ,lomny     )
         !
         ! For points with X-values between Xpmin and Xpmax set: iFlag(i)=1
         !
         lominx = lomnx
         lomaxx = lomnx
         do i = lomnx, 1, -1
            if (xs(i)>=xpmin) then
               lominx        = i
               iflag(nrx(i)) = 1
            else
               exit
            endif
         enddo
         do i = lomnx + 1, n2
            if (xs(i)<=xpmax) then
               lomaxx        = i
               iflag(nrx(i)) = 1
            else
               exit
            endif
         enddo
         !
         ! For the points with Y-values between Ypmin and Ypmax,
         ! that also lie between Xpmin and Xpmax: Save them in NrIn
         !
         iin = 1
         do i = lomny, 1, -1
            if (ys(i)>=ypmin) then
               nrin(iin) = nry(i)*iflag(nry(i))
               iin       = iin + iflag(nry(i))
            else
               exit
            endif
         enddo
         do i = lomny + 1, n2
            if (ys(i)<=ypmax) then
               nrin(iin) = nry(i)*iflag(nry(i))
               iin       = iin + iflag(nry(i))
            else
               exit
            endif
         enddo
         nin = iin - 1
         !
         ! Put iFlag back to 0
         !
         do i = lominx, lomaxx
            if (i/=0) iflag(nrx(i)) = 0
         enddo
         !
         ! Check whether selected points of grid2 lie within the cell
         ! using subroutine IPON; if so, determine weights W of the surrounding
         ! values in grid1 using subroutine INTRP4. Save the weights in Wtab
         ! The reference to grid1 is saved in arrays Iref and Jref.
         !
         do iin = 1, nin
            i2 = nrin(iin)
            inout = -1
            call ipon(xp        ,yp        ,4         ,x2(i2)    , &
               &          y2(i2)    ,inout     )
            if (inout>=0) then
               !
               ! Check on point with nonsense information on source grid 1;
               ! This depends on agreements
               !
               if (      (code(i1    , j1    )==1 .or. code(i1    , j1    )==3) &
                  & .and. (code(i1 + 1, j1    )==1 .or. code(i1 + 1, j1    )==3) &
                  & .and. (code(i1 + 1, j1 + 1)==1 .or. code(i1 + 1, j1 + 1)==3) &
                  & .and. (code(i1    , j1 + 1)==1 .or. code(i1    , j1 + 1)==3) ) then
                  !
                  ! grid2 point is covered by 4 valid grid1 points
                  !
                  covered(i2) = 1
                  call bilin5(xp, yp, x2(i2), y2(i2), hbuff, ier)
                  w(:, i2) = real(hbuff(:),4)
                  !
                  iref(1, i2) = i1 + (j1 - 1)*m1
                  iref(2, i2) = i1 + 1 + (j1 - 1)*m1
                  iref(3, i2) = i1 + 1 + j1*m1
                  iref(4, i2) = i1 + j1*m1
               elseif (code(i1, j1)>0 .and. code(i1 + 1, j1)>0 .and. &
                  & code(i1 + 1, j1 + 1)>0 .and. code(i1, j1 + 1)>0) then
                  !
                  ! Grid2 point is covered by 4 valid grid1 points, containing
                  ! one or more boundary points:
                  ! - do not use grid1 values
                  ! - use extrapolation values
                  !
                  covered(i2) = 0
               else
                  !
                  ! Grid2 point is not covered by 4 valid grid1 points:
                  ! - do not use grid1 values
                  ! - do not use extrapolation values
                  !
                  covered(i2) = -1
               endif
            endif
         enddo
      enddo
   enddo
   end subroutine mkmap

   subroutine make_map(code, x1, y1, m1, n1, x2, y2, n2, w, iref, &
      &                   iprint    ,covered   ,xymiss)
   !!--description-----------------------------------------------------------------
   !
   !
   !     SUBROUTINE MKMAP
   !     Interpolation of curvilinear, numerically ordered grid (grid1)
   !     to random points, with weighting points (grid2).
   !
   !     J.A. Roelvink
   !     Deltares
   !     24-2-1992 (MKMAP
   !
   !     Given: numerically ordered grid M1*N1
   !     with coordinates X1 (1:M1,1:N1)
   !                  and Y1 (1:M1,1:N1)
   !
   !     Also given: random points X2(1:N2)
   !                           and Y2(1:N2)
   !
   !     To be determined:  weighting factors and pointers for bilinear interpolation
   !     Weighting factors and locations of the requested (random) points in the
   !     ordered grid are saved in resp.
   !     W(1:4,1:N2) and Iref(1:4,1:N2)
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                    , intent(in)  :: iprint
   integer                    , intent(in)  :: m1
   integer                    , intent(in)  :: n1
   integer                    , intent(in)  :: n2
   integer , dimension(4 , n2), intent(out) :: iref
   real*8  , dimension(4 , n2), intent(out) :: w       ! Remains single precision!
   integer , dimension(m1, n1), intent(in)  :: code
   integer , dimension(n2)    , intent(out) :: covered !  0: target point is not covered by source grid (default)
   !  1: target point is covered by   valid points of source grid
   ! -1: target point is covered by invalid points of source grid
   real                       , intent(in)  :: xymiss  ! missing value in grid 1
   real*8  , dimension(m1, n1), intent(in)  :: x1
   real*8  , dimension(m1, n1), intent(in)  :: y1
   real*8  , dimension(n2), intent(in)      :: x2
   real*8  , dimension(n2), intent(in)      :: y2
   !
   ! Local variables
   !
   integer, dimension(:),allocatable             :: iflag
   integer, dimension(:),allocatable             :: nrin
   integer, dimension(:),allocatable             :: nrx
   integer, dimension(:),allocatable             :: nry
   real*8,  dimension(:),allocatable             :: xs
   real*8,  dimension(:),allocatable             :: ys

   integer                           :: i
   integer                           :: i1
   integer                           :: i2
   integer                           :: ier
   integer                           :: iin
   integer                           :: inout
   integer                           :: ip
   integer                           :: j1
   integer                           :: lomaxx
   integer                           :: lominx
   integer                           :: lomnx
   integer                           :: lomny
   integer                           :: m1max
   integer                           :: m1min
   integer                           :: n1max
   integer                           :: n1min
   integer                           :: nin
   real*8                           :: eps
   real*8                            :: xpmax
   real*8                            :: xpmean
   real*8                            :: xpmin
   real*8                            :: ypmax
   real*8                            :: ypmean
   real*8                            :: ypmin
   real*8  , dimension(5)            :: xp
   real*8  , dimension(5)            :: yp
   real*8  , dimension(4)            :: hbuff
   !
   allocate(iflag(n2))
   allocate(nrin(n2))
   allocate(nrx(n2))
   allocate(nry(n2))
   allocate(xs(n2))
   allocate(ys(n2))
   !! executable statements -------------------------------------------------------
   !
   eps = 0.00001
   !
   !     Initialise tables
   !
   if (iprint==1) write (*, *) 'in mkmap', m1, n1, n2
   lomnx = 1
   lomny = 1
   m1min = 1
   m1max = m1
   n1min = 1
   n1max = n1
   nrx=0
   nry=0
   iflag=0
   nrin=0
   xs=0.d0
   ys=0.d0
   iref=0
   w=0.d0
   !
   ! Sort X2 en Y2
   !
   call sort(n2        ,x2        ,xs        ,nrx       )
   call sort(n2        ,y2        ,ys        ,nry       )
   !
   ! Loop over all cels of grid1
   !
   !
   do j1 = n1min, n1max - 1
      do i1 = m1min, m1max - 1
         !
         ! Cell definition
         !
         xp(1) = x1(i1, j1)
         xp(2) = x1(i1 + 1, j1)
         xp(3) = x1(i1 + 1, j1 + 1)
         xp(4) = x1(i1, j1 + 1)
         yp(1) = y1(i1, j1)
         yp(2) = y1(i1 + 1, j1)
         yp(3) = y1(i1 + 1, j1 + 1)
         yp(4) = y1(i1, j1 + 1)
         if (     (xp(1) > xymiss-eps .and. xp(1) < xymiss+eps) &
            & .or. (xp(2) > xymiss-eps .and. xp(2) < xymiss+eps) &
            & .or. (xp(3) > xymiss-eps .and. xp(3) < xymiss+eps) &
            & .or. (xp(4) > xymiss-eps .and. xp(4) < xymiss+eps) &
            & .or. (yp(1) > xymiss-eps .and. yp(1) < xymiss+eps) &
            & .or. (yp(2) > xymiss-eps .and. yp(2) < xymiss+eps) &
            & .or. (yp(3) > xymiss-eps .and. yp(3) < xymiss+eps) &
            & .or. (yp(4) > xymiss-eps .and. yp(4) < xymiss+eps) ) cycle
         !
         ! Determine minimum and maximum X and Y of the cell
         !
         xpmin =  1.e10
         xpmax = -1.e10
         ypmin =  1.e10
         ypmax = -1.e10
         do ip = 1, 4
            xpmin = min(xp(ip), xpmin)
            xpmax = max(xp(ip), xpmax)
            ypmin = min(yp(ip), ypmin)
            ypmax = max(yp(ip), ypmax)
         enddo
         xpmean = .5*(xpmin + xpmax)
         ypmean = .5*(ypmin + ypmax)
         !
         ! First selection of points of grid2 that can be located in the cell
         !
         ! Find centre of the cell in tables Xs and Ys
         !
         call hunt(xs        ,n2        ,xpmean    ,lomnx     )
         call hunt(ys        ,n2        ,ypmean    ,lomny     )
         !
         ! For points with X-values between Xpmin and Xpmax set: iFlag(i)=1
         !
         lominx = lomnx
         lomaxx = lomnx
         do i = lomnx, 1, -1
            if (xs(i)>=xpmin) then
               lominx        = i
               iflag(nrx(i)) = 1
            else
               exit
            endif
         enddo
         do i = lomnx + 1, n2
            if (xs(i)<=xpmax) then
               lomaxx        = i
               iflag(nrx(i)) = 1
            else
               exit
            endif
         enddo
         !
         ! For the points with Y-values between Ypmin and Ypmax,
         ! that also lie between Xpmin and Xpmax: Save them in NrIn
         !
         iin = 1
         do i = lomny, 1, -1
            if (ys(i)>=ypmin) then
               nrin(iin) = nry(i)*iflag(nry(i))
               iin       = iin + iflag(nry(i))
            else
               exit
            endif
         enddo
         do i = lomny + 1, n2
            if (ys(i)<=ypmax) then
               nrin(iin) = nry(i)*iflag(nry(i))
               iin       = iin + iflag(nry(i))
            else
               exit
            endif
         enddo
         nin = iin - 1
         !
         ! Put iFlag back to 0
         !
         do i = lominx, lomaxx
            if (i/=0) iflag(nrx(i)) = 0
         enddo
         !
         ! Check whether selected points of grid2 lie within the cell
         ! using subroutine IPON; if so, determine weights W of the surrounding
         ! values in grid1 using subroutine INTRP4. Save the weights in Wtab
         ! The reference to grid1 is saved in arrays Iref and Jref.
         !
         do iin = 1, nin
            i2 = nrin(iin)
            inout = -1
            call ipon(xp        ,yp        ,4         ,x2(i2)    , &
               &          y2(i2)    ,inout     )
            if (inout>=0) then
               !
               ! Check on point with nonsense information on source grid 1;
               ! This depends on agreements
               !
               if (      (code(i1    , j1    )==1 .or. code(i1    , j1    )==3) &
                  & .and. (code(i1 + 1, j1    )==1 .or. code(i1 + 1, j1    )==3) &
                  & .and. (code(i1 + 1, j1 + 1)==1 .or. code(i1 + 1, j1 + 1)==3) &
                  & .and. (code(i1    , j1 + 1)==1 .or. code(i1    , j1 + 1)==3) ) then
                  !
                  ! grid2 point is covered by 4 valid grid1 points
                  !
                  covered(i2) = 1
                  call bilin5(xp, yp, x2(i2), y2(i2), hbuff, ier)
                  w(:, i2) = real(hbuff(:),4)
                  !
                  iref(1, i2) = i1 + (j1 - 1)*m1
                  iref(2, i2) = i1 + 1 + (j1 - 1)*m1
                  iref(3, i2) = i1 + 1 + j1*m1
                  iref(4, i2) = i1 + j1*m1
               elseif (code(i1, j1)>0 .and. code(i1 + 1, j1)>0 .and. &
                  & code(i1 + 1, j1 + 1)>0 .and. code(i1, j1 + 1)>0) then
                  !
                  ! Grid2 point is covered by 4 valid grid1 points, containing
                  ! one or more boundary points:
                  ! - do not use grid1 values
                  ! - use extrapolation values
                  !
                  covered(i2) = 0
               else
                  !
                  ! Grid2 point is not covered by 4 valid grid1 points:
                  ! - do not use grid1 values
                  ! - do not use extrapolation values
                  !
                  covered(i2) = -1
               endif
            endif
         enddo
      enddo
   enddo

   end subroutine make_map

   subroutine make_map_fm (x1, y1, face_nodes, no_nodes1, no_faces1, x2, y2, n2, w, iref, nref)
   !!--description-----------------------------------------------------------------
   !
   !
   !     SUBROUTINE MKMAP
   !     Interpolation of unstructured grid (grid1)
   !     to random points, with weighting points (grid2).
   !
   !     J.A. Roelvink
   !     Deltares
   !     24-2-1992 (MKMAP
   !
   !     Given: unstructured grid
   !     with coordinates X1 (1:no_nodes1)
   !                  and Y1 (1:no_nodes1)
   !     and cells connecting nodes face_nodes (4,no_faces1)
   !
   !     Also given: random points X2(1:N2)
   !                           and Y2(1:N2)
   !
   !     To be determined:  weighting factors and pointers for bilinear interpolation
   !     Weighting factors and node numbers of the requested (random) points in the
   !     unstructured grid 1 are saved in resp.
   !     W(1:4,1:N2) and Iref(1:4,1:N2)
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                       , intent(in)  :: no_nodes1
   integer                       , intent(in)  :: no_faces1
   integer                       , intent(in)  :: n2
   integer , dimension(4 , n2)   , intent(out) :: iref
   integer , dimension(no_nodes1), intent(out) :: nref
   real*8  , dimension(4 , n2)   , intent(out) :: w       ! Remains single precision!
   real*8  , dimension(no_nodes1), intent(in)  :: x1
   real*8  , dimension(no_nodes1), intent(in)  :: y1
   integer , dimension(4,no_faces1), intent(in) :: face_nodes
   real*8  , dimension(n2), intent(in)      :: x2
   real*8  , dimension(n2), intent(in)      :: y2
   !
   ! Local variables
   !
   integer, dimension(:),allocatable             :: iflag
   integer, dimension(:),allocatable             :: nrin
   integer, dimension(:),allocatable             :: nrx
   integer, dimension(:),allocatable             :: nry
   real*8,  dimension(:),allocatable             :: xs
   real*8,  dimension(:),allocatable             :: ys

   integer                           :: i
   integer                           :: i1
   integer                           :: i2
   integer                           :: k1
   integer                           :: np
   integer                           :: ier
   integer                           :: iin
   integer                           :: inout
   integer                           :: ip
   integer                           :: j1
   integer                           :: lomaxx
   integer                           :: lominx
   integer                           :: lomnx
   integer                           :: lomny
   integer                           :: nin
   real*8                            :: eps
   real*8                            :: xpmax
   real*8                            :: xpmean
   real*8                            :: xpmin
   real*8                            :: ypmax
   real*8                            :: ypmean
   real*8                            :: ypmin
   real*8  , dimension(4)            :: xp
   real*8  , dimension(4)            :: yp
   real*8  , dimension(4)            :: hbuff
   !
   allocate(iflag(n2))
   allocate(nrin(n2))
   allocate(nrx(n2))
   allocate(nry(n2))
   allocate(xs(n2))
   allocate(ys(n2))
   !! executable statements -------------------------------------------------------
   !
   eps = 0.00001
   !
   !     Initialise tables
   !
   lomnx = 1
   lomny = 1
   nrx=0
   nry=0
   iflag=0
   nrin=0
   xs=0.d0
   ys=0.d0
   iref=0
   nref=0
   w=0.d0
   !
   ! Sort X2 en Y2
   !
   call sort(n2        ,x2        ,xs        ,nrx       )
   call sort(n2        ,y2        ,ys        ,nry       )
   !
   ! Loop over all cells of grid1
   !
   !
   do i1 = 1,no_faces1
      !
      ! Cell definition
      !
      if (face_nodes(4,i1)==-999) then
         np=3
      else
         np=4
      endif

      do ip=1,np
         k1=face_nodes(ip,i1)
         xp(ip)=x1(k1)
         yp(ip)=y1(k1)
      enddo
      !
      ! Determine minimum and maximum X and Y of the cell
      !
      xpmin = minval(xp(1:np))
      xpmax = maxval(xp(1:np))
      ypmin = minval(yp(1:np))
      ypmax = maxval(yp(1:np))
      xpmean = .5*(xpmin + xpmax)
      ypmean = .5*(ypmin + ypmax)
      !
      ! First selection of points of grid2 that can be located in the cell
      !
      ! Find centre of the cell in tables Xs and Ys
      !
      call hunt(xs        ,n2        ,xpmean    ,lomnx     )
      call hunt(ys        ,n2        ,ypmean    ,lomny     )
      !
      ! For points with X-values between Xpmin and Xpmax set: iFlag(i)=1
      !
      lominx = lomnx
      lomaxx = lomnx
      do i = lomnx, 1, -1
         if (xs(i)>=xpmin) then
            lominx        = i
            iflag(nrx(i)) = 1
         else
            exit
         endif
      enddo
      do i = lomnx + 1, n2
         if (xs(i)<=xpmax) then
            lomaxx        = i
            iflag(nrx(i)) = 1
         else
            exit
         endif
      enddo
      !
      ! For the points with Y-values between Ypmin and Ypmax,
      ! that also lie between Xpmin and Xpmax: Save them in NrIn
      !
      iin = 1
      do i = lomny, 1, -1
         if (ys(i)>=ypmin) then
            nrin(iin) = nry(i)*iflag(nry(i))
            iin       = iin + iflag(nry(i))
         else
            exit
         endif
      enddo
      do i = lomny + 1, n2
         if (ys(i)<=ypmax) then
            nrin(iin) = nry(i)*iflag(nry(i))
            iin       = iin + iflag(nry(i))
         else
            exit
         endif
      enddo
      nin = iin - 1
      !
      ! Put iFlag back to 0
      !
      do i = lominx, lomaxx
         if (i/=0) iflag(nrx(i)) = 0
      enddo
      !
      ! Check whether selected points of grid2 lie within the cell
      ! using subroutine IPON; if so, determine weights W of the surrounding
      ! values in grid1 using subroutine INTRP4. Save the weights in Wtab
      ! The reference to grid1 is saved in arrays Iref and Jref.
      !
      do iin = 1, nin
         i2 = nrin(iin)
         inout = -1
         call ipon(xp        ,yp        ,np         ,x2(i2)    , &
            &          y2(i2)    ,inout     )
         if (inout>=0) then
            !
            ! grid2 point is covered by 4 valid grid1 points
            !
            if (np==4) then
               call bilin5(xp, yp, x2(i2), y2(i2), w(:,i2), ier)
               do ip=1,4
                  iref(ip,i2)=face_nodes(ip,i1)
                  nref(iref(ip,i2))=nref(iref(ip,i2))+1
               enddo
            else
               call triangle_intp(xp,yp,x2(i2),y2(i2),w(1:3,i2))
               do ip=1,3
                  iref(ip,i2)=face_nodes(ip,i1)
                  nref(iref(ip,i2))=nref(iref(ip,i2))+1
               enddo
               w(4,i2)=0.d0
               iref(4,i2)=1
            endif
         endif
      enddo
   enddo

   end subroutine make_map_fm

   subroutine mkmap_step (x1        ,y1        ,m1        ,n1        ,    &
      &       alfaz     ,dsu       ,dnv       ,               &
      &       x2        ,y2        ,n2        ,w         ,  iref      )
   !!--description-----------------------------------------------------------------
   !
   !
   !     SUBROUTINE MKMAP_STEP
   !     Interpolation of curvilinear, numerically ordered grid (grid1)
   !     to random points, with weighting points (grid2).
   !
   !     J.A. Roelvink
   !     UNESCO-IHE / Deltares
   !     8 May 2013
   !
   !     Given: numerically ordered grid M1*N1
   !     with (new)coordinates X1 (1:M1,1:N1)
   !                  and Y1 (1:M1,1:N1)
   !
   !     Also given: random points X2(1:N2)
   !                           and Y2(1:N2)
   !                 previous iref
   !
   !     To be determined:  weighting factors and updated pointers for bilinear interpolation
   !     Weighting factors and locations of the requested (random) points in the
   !     ordered grid are saved in resp.
   !     W(1:4,1:N2) and Iref(1:4,1:N2)
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                    , intent(in)  :: m1
   integer                    , intent(in)  :: n1
   integer                                  :: n2
   integer , dimension(4 , n2)              :: iref
   real*8  , dimension(4 , n2)              :: w
   real*8  , dimension(m1, n1), intent(in)  :: x1
   real*8  , dimension(m1, n1), intent(in)  :: y1
   real*8  , dimension(m1, n1), intent(in)  :: alfaz
   real*8  , dimension(m1, n1), intent(in)  :: dsu
   real*8  , dimension(m1, n1), intent(in)  :: dnv
   real*8  , dimension(n2)                  :: x2
   real*8  , dimension(n2)                  :: y2
   real*8                                   :: dx
   real*8                                   :: dy
   real*8                                   :: ds
   real*8                                   :: dn
   real*8                                   :: dsrel
   real*8                                   :: dnrel
   integer                                  :: i1,j1,i2

   do i2=1,n2
      if (iref(1,i2)>0) then
         i1=mod(iref(1,i2)-1,m1)+1
         j1=(iref(1,i2)-i1)/m1+1
         dx=x2(i2)-x1(i1,j1)
         dy=y2(i2)-y1(i1,j1)
         ds= dx*cos(alfaz(i1,j1))+dy*sin(alfaz(i1,j1))
         dn=-dx*sin(alfaz(i1,j1))+dy*cos(alfaz(i1,j1))
         dsrel=ds/dsu(i1,j1)
         dnrel=dn/dnv(i1,j1)
         i1=i1+floor(dsrel)
         j1=j1+floor(dnrel)
         dsrel=dsrel-floor(dsrel)
         dnrel=dnrel-floor(dnrel)
         iref(1, i2) = i1 + (j1 - 1)*m1
         iref(2, i2) = i1 + 1 + (j1 - 1)*m1
         iref(3, i2) = i1 + 1 + j1*m1
         iref(4, i2) = i1 + j1*m1
         w(1, i2) = (1.d0-dsrel)*(1.d0-dnrel)
         w(2, i2) =       dsrel *(1.d0-dnrel)
         w(3, i2) =       dsrel *      dnrel
         w(4, i2) = (1.d0-dsrel)*      dnrel
      endif
   enddo

   end subroutine mkmap_step

   subroutine grmap(f1        ,n1        ,f2        ,n2        ,iref      , &
      & w         ,np        ,iprint    )
   !!--description-----------------------------------------------------------------
   !
   ! compute interpolated values for all points on grid 2
   !
   ! special treatment of points on grid 2 that are outside
   ! grid 1; in that case iref(1,i2)=0 AND w(ip,i2)=0 for all ip
   !
   ! Iref(1,i2)   i1    ifac   F2(i2)*ifac     Result
   !
   !      0        1      1      F2(i2)        Old value is kept
   !    j,j>0      j      0       0.           F2 is initialized
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   implicit none
   !
   ! Global variables
   !
   integer                   , intent(in)  :: iprint
   integer                   , intent(in)  :: n1
   integer                   , intent(in)  :: n2
   integer                   , intent(in)  :: np
   integer, dimension(np, n2), intent(in)  :: iref
   real*4 , dimension(n1)    , intent(in)  :: f1
   real*4 , dimension(n2)                  :: f2
   real*8 , dimension(np, n2), intent(in)  :: w
   !
   ! Local variables
   !
   integer :: i
   integer :: i1
   integer :: i2
   integer :: ip
   !
   !! executable statements -------------------------------------------------------
   !
   if (iprint==1) write (*, *) 'in grmap n1 n2', n1, n2
   do i2 = 1, n2
      i = iref(1, i2)
      if (i>0) then
         f2(i2)=0.d0
         !        i1 = max(i, 1)
         !        ifac = 1 - i/i1
         !        f2(i2) = f2(i2)*ifac
         !
         ! Function values at grid 2 are expressed as weighted average
         ! of function values in Np surrounding points of grid 1
         !
         if (iprint==1 .and. i2<=n2) &
            & write (*, '(1X,A,I6,4(1X,E10.4))') ' i2 w ', i2, (w(ip, i2) , ip = 1,  &
            & np)
         do ip = 1, np
            i = iref(ip, i2)
            i1 = max(i, 1)
            if (iprint==1 .and. i2<=n2) write (*, *) ' i1,f1(i1) ', i1, f1(i1)
            f2(i2) = f2(i2) + w(ip, i2)*f1(i1)
         enddo
      endif
   enddo
   end subroutine grmap

   subroutine grmap2(f1, cellsz1i     , n1        ,f2  ,cellsz2,  n2        ,iref      , &
      & w         ,np     )
   !!--description-----------------------------------------------------------------
   !
   ! compute interpolated values for all points on grid 1 given reference table
   ! for grid 2; this works the other way round from GRMAP. Assumption is that
   ! grid 2 is much finer than grid 1. For each point in grid 2 we know the
   ! surrounding points in grid 1 and the related weights. Instead of using this to
   ! interpolate from 1 to 2 we now integrate from 2 to 1 using the same weights.
   !
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   implicit none
   !
   ! Global variables
   !
   integer                   , intent(in)  :: n1
   integer                   , intent(in)  :: n2
   integer                   , intent(in)  :: np
   integer, dimension(np, n2), intent(in)  :: iref
   real*8 , dimension(n1)                  :: f1
   real*8 , dimension(n1)    , intent(in)  :: cellsz1i   !array with 1/cell size
   real*8                    , intent(in)  :: cellsz2
   real*8 , dimension(n2)    , intent(in)  :: f2
   real*8 , dimension(np, n2), intent(in)  :: w
   !
   ! Local variables
   !
   integer :: i1
   integer :: i2
   integer :: ip
   !
   !! executable statements -------------------------------------------------------
   !
   do ip = 1, np
      do i2=1,n2
         i1 = iref(ip, i2)
         if (i1>0) then
            !f2(i2) = f2(i2) + w(ip, i2)*f1(i1)
            f1(i1) = f1(i1) + w(ip, i2)*f2(i2)*cellsz2*cellsz1i(i1)
         endif
      enddo
   enddo
   end subroutine grmap2

   subroutine grmap_sg(f1, nref, n1, f2, n2, iref, w, np)
   !!--description-----------------------------------------------------------------
   !
   ! compute interpolated values for all points on grid 1 given reference table
   ! for grid 2; this works the other way round from GRMAP. Assumption is that
   ! grid 2 is much finer than grid 1. For each point in grid 2 we know the
   ! surrounding points in grid 1 and the related weights. Instead of using this to
   ! interpolate from 1 to 2 we now integrate from 2 to 1 using the same weights.
   !
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   implicit none
   !
   ! Global variables
   !
   integer                   , intent(in)  :: n1
   integer                   , intent(in)  :: n2
   integer                   , intent(in)  :: np
   integer, dimension(np, n2), intent(in)  :: iref
   integer, dimension(n1)    , intent(in)  :: nref
   real*8 , dimension(n1)                  :: f1
   real*8 , dimension(n2)    , intent(in)  :: f2
   real*8 , dimension(np, n2), intent(in)  :: w
   !
   ! Local variables
   !
   integer :: i1
   integer :: i2
   integer :: ip
   !
   !! executable statements -------------------------------------------------------
   !
   f1=0.d0
   do i2=1,n2
      do ip = 1, np
         i1 = iref(ip, i2)
         if (i1>0) then
            !f1(i1) = f1(i1) + w(ip, i2)*f2(i2)/nref(i1)
            f1(i1) = f1(i1) +            f2(i2)/nref(i1)
         endif
      enddo
   enddo
   end subroutine grmap_sg

   subroutine ipon(xq     ,yq     ,n      ,xp     ,yp     ,inout     )
   !--description----------------------------------------------------------------
   !
   ! Deltares                                                               *
   ! AUTHOR : J.A.ROELVINK                                                  *
   ! DATE   : 22-12-1988                                                    *
   ! DETERMINE WHETHER POINT (xp,yp) LIES IN POLYGON (x,y) OF n POINTS      *
   ! POINT n+1 IS SET EQUAL TO POINT 1                                      *
   ! (ARRAY MUST HAVE DIMENSION n+1 IN MAIN PROGRAMME                       *
   ! inpout = -1 :  OUTSIDE POLYGON                                         *
   ! inpout =  0 :  ON EDGE OF POLYGON                                      *
   ! inpout =  1 :  INSIDE POLYGON                                          *
   ! USED METHOD :         - DRAW A VERTICAL LINE THROUGH (xp,yp)           *
   !                       - DETERMINE NUMBER OF INTERSECTIONS WITH POLYGON *
   !                         UNDER yp : nunder                              *
   !                       - IF nunder IS EVEN, THEN THE POINT LIES OUTSIDE *
   !                         THE POLYGON, OTHERWISE IT LIES INSIDE          *
   !                       - THE EDGE IS TREATED SEPARATELY                 *
   !
   !--pseudo code and references-------------------------------------------------
   ! NONE
   !--declarations---------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer               , intent(out) :: inout
   integer               , intent(in)  :: n
   real*8                , intent(in)  :: xp
   real*8                , intent(in)  :: yp
   real*8  , dimension(*)              :: xq
   real*8  , dimension(*)              :: yq
   !
   ! Local variables
   !
   integer                             :: i
   integer                             :: ierr
   integer                             :: nunder
   real*4                              :: ysn
   real*4  , dimension(:), allocatable :: x
   real*4  , dimension(:), allocatable :: y
   !
   ! executable statements ------------------------------------------------------
   !
   allocate(x(n+1))
   allocate(y(n+1))
   do i = 1, n
      x(i) = real( xq(i)-xp , 4)
      y(i) = real( yq(i)-yp , 4)
   enddo
   x(n + 1) = x(1)
   y(n + 1) = y(1)
   nunder   = 0
   do i = 1, n
      if ((x(i)<0. .and. x(i + 1)>=0.).or.(x(i + 1)<0. .and. x(i)>=0.)) then
         if (y(i)<0. .and. y(i + 1)<0.) then
            nunder = nunder + 1
         elseif ((y(i)<=0. .and. y(i + 1)>=0.) .or.                       &
            & (y(i + 1)<=0. .and. y(i)>=0.)) then
            ysn = (y(i)*x(i + 1) - x(i)*y(i + 1))/(x(i + 1) - x(i))
            if (ysn<0.) then
               nunder = nunder + 1
            elseif (ysn<=0.) then
               !
               ! Edge
               !
               inout = 0
               goto 100
            else
            endif
         else
         endif
      elseif (abs(x(i))<1.0E-8 .and. abs(x(i + 1))<1.0E-8) then
         if ((y(i)<=0. .and. y(i + 1)>=0.).or.(y(i + 1)<=0..and.y(i)>=0.)) &
            & then
         !
         ! Edge
         !
         inout = 0
         goto 100
      endif
   else
   endif
   enddo
   if (mod(nunder, 2)==0) then
      !
      ! Outside
      !
      inout = -1
   else
      !
      ! Inside
      !
      inout = 1
   endif
100 continue
   deallocate(x, stat=ierr)
   deallocate(y, stat=ierr)
   end subroutine ipon

   subroutine hunt(xx        ,n         ,x         ,jlo       )
   !!--description-----------------------------------------------------------------
   ! NONE
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                                  :: jlo
   integer                    , intent(in)  :: n
   real*8              , intent(in)  :: x
   real*8, dimension(n), intent(in)  :: xx
   !
   ! Local variables
   !
   integer :: inc
   integer :: jhi
   integer :: jm
   logical :: ascnd
   !
   !! executable statements -------------------------------------------------------
   !
   ascnd = xx(n)>=xx(1)
   if (jlo<=0 .or. jlo>n) then
      jlo = 0
      jhi = n + 1
      goto 3
   endif
   inc = 1
   if (x>=xx(jlo) .eqv. ascnd) then
1     continue
      jhi = jlo + inc
      if (jhi>n) then
         jhi = n + 1
      elseif (x>=xx(jhi) .eqv. ascnd) then
         jlo = jhi
         inc = inc + inc
         goto 1
      else
      endif
   else
      jhi = jlo
2     continue
      jlo = jhi - inc
      if (jlo<1) then
         jlo = 0
      elseif (x<xx(jlo) .eqv. ascnd) then
         jhi = jlo
         inc = inc + inc
         goto 2
      else
      endif
   endif
3  continue
   if (jhi - jlo==1) then
      return
   endif
   jm = (jhi + jlo)/2
   if (x>xx(jm) .eqv. ascnd) then
      jlo = jm
   else
      jhi = jm
   endif
   goto 3
   end subroutine hunt

   subroutine indexx(n         ,arrin     ,indx      )
   !----- GPL ---------------------------------------------------------------------
   !
   !  Copyright (C)  Stichting Deltares, 2011.
   !
   !  This program is free software: you can redistribute it and/or modify
   !  it under the terms of the GNU General Public License as published by
   !  the Free Software Foundation version 3.
   !
   !  This program is distributed in the hope that it will be useful,
   !  but WITHOUT ANY WARRANTY; without even the implied warranty of
   !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   !  GNU General Public License for more details.
   !
   !  You should have received a copy of the GNU General Public License
   !  along with this program.  If not, see <http://www.gnu.org/licenses/>.
   !
   !  contact: delft3d.support@deltares.nl
   !  Stichting Deltares
   !  P.O. Box 177
   !  2600 MH Delft, The Netherlands
   !
   !  All indications and logos of, and references to, "Delft3D" and "Deltares"
   !  are registered trademarks of Stichting Deltares, and remain the property of
   !  Stichting Deltares. All rights reserved.
   !
   !-------------------------------------------------------------------------------
   !  $Id: interp.F90 5598 2019-12-09 10:21:29Z ridde_mo $
   !  $HeadURL: https://svn.oss.deltares.nl/repos/xbeach/trunk/src/xbeachlibrary/interp.F90 $
   !!--description-----------------------------------------------------------------
   ! NONE
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                       , intent(in)  :: n
   integer, dimension(n)                       :: indx
   real*8   , dimension(n), intent(in)  :: arrin
   !
   ! Local variables
   !
   integer :: i
   integer :: indxt
   integer :: ir
   integer :: j
   integer :: l
   real*8 :: q
   !
   !! executable statements -------------------------------------------------------
   !
   do j = 1, n
      indx(j) = j
   enddo
   l = n/2 + 1
   ir = n
10 continue
   if (l>1) then
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
   else
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
      if (ir==1) then
         indx(1) = indxt
         return
      endif
   endif
   i = l
   j = l + l
20 continue
   if (j<=ir) then
      if (j<ir) then
         if (arrin(indx(j))<arrin(indx(j + 1))) j = j + 1
      endif
      if (q<arrin(indx(j))) then
         indx(i) = indx(j)
         i = j
         j = j + j
      else
         j = ir + 1
      endif
      goto 20
   endif
   indx(i) = indxt
   goto 10
   end subroutine indexx

   subroutine bilin5(xa, ya, x0, y0, w, ier)
   !!--description-----------------------------------------------------------------
   ! NONE
   !!--pseudo code and references--------------------------------------------------
   !
   ! Author: H. Petit
   !
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer               , intent(out) :: ier
   real*8                , intent(in)  :: x0
   real*8                , intent(in)  :: y0
   real*8  , dimension(4), intent(out) :: w
   real*8  , dimension(4), intent(in)  :: xa
   real*8  , dimension(4), intent(in)  :: ya
   !
   ! Local variables
   !
   real*8   :: a
   real*8   :: a21
   real*8   :: a22
   real*8   :: a31
   real*8   :: a32
   real*8   :: a41
   real*8   :: a42
   real*8   :: b
   real*8   :: c
   real*8   :: det
   real*8   :: discr
   real*8   :: eta
   real*8   :: x
   real*8   :: x1
   real*8   :: x2
   real*8   :: x3
   real*8   :: x3t
   real*8   :: x4
   real*8   :: xi
   real*8   :: xt
   real*8   :: y
   real*8   :: y1
   real*8   :: y2
   real*8   :: y3
   real*8   :: y3t
   real*8   :: y4
   real*8   :: yt
   !
   !! executable statements -------------------------------------------------------
   !
   ! read(12,*)x1,y1,f1
   x1 = xa(1)
   y1 = ya(1)
   ! read(12,*)x2,y2,f2
   x2 = xa(2)
   y2 = ya(2)
   ! read(12,*)x3,y3,f3
   x3 = xa(3)
   y3 = ya(3)
   ! read(12,*)x4,y4,f4
   x4 = xa(4)
   y4 = ya(4)
   x  = x0
   y  = y0
   !
   ! The bilinear interpolation problem is first transformed
   ! to the quadrangle with nodes (0,0),(1,0),(x3t,y3t),(0,1)
   ! and required location (xt,yt)
   !
   a21 = x2 - x1
   a22 = y2 - y1
   a31 = x3 - x1
   a32 = y3 - y1
   a41 = x4 - x1
   a42 = y4 - y1
   det = a21*a42 - a22*a41
   if (abs(det) < 1.0e-20) then
      ier = 1
      goto 99999
   endif
   x3t = (  a42*a31      - a41*a32     ) / det
   y3t = ( -a22*a31      + a21*a32     ) / det
   xt  = (  a42*(x - x1) - a41*(y - y1)) / det
   yt  = ( -a22*(x - x1) + a21*(y - y1)) / det
   if ((x3t < 0.0) .or. (y3t < 0.0)) then
      ! write (*, *) 'distorted quadrangle'
      ier = 1
      goto 99999
   endif
   if (abs(x3t - 1.0d0) < 1.0e-7) then
      xi = xt
      if (abs(y3t - 1.0d0) < 1.0e-7) then
         eta = yt
      elseif (abs(1.0d0 + (y3t - 1.0d0)*xt) < 1.0e-6) then
         ! write (*, *) 'extrapolation over too large a distance'
         ier = 1
         goto 99999
      else
         eta = yt / (1.0d0 + (y3t - 1.0d0)*xt)
      endif
   elseif (abs(y3t - 1.0d0) < 1.0e-6) then
      eta = yt
      if (abs(1.0d0 + (x3t - 1.0d0)*yt) < 1.0e-6) then
         ! write (*, *) 'extrapolation over too large a distance'
         ier = 1
         goto 99999
      else
         xi = xt / (1.0d0 + (x3t - 1.0d0)*yt)
      endif
   else
      a     = y3t - 1.0d0
      b     = 1.0d0 + (x3t - 1.0d0)*yt - (y3t - 1.0d0)*xt
      c     = -xt
      discr = b*b - 4.0d0*a*c
      if (discr < 1.0e-6) then
         ! write (*, *) 'extrapolation over too large a distance'
         ier = 1
         goto 99999
      endif
      xi  = (-b + sqrt(discr)) / (2.0d0*a)
      eta = ((y3t - 1.0d0)*(xi - xt) + (x3t - 1.0d0)*yt) / (x3t - 1.0d0)
   endif
   w(1) = (1.0d0-xi) * (1.0d0-eta)
   w(2) =         xi  * (1.0d0-eta)
   w(3) =         xi  *         eta
   w(4) =        eta  * (1.0d0-xi )
   return
99999 continue
   end subroutine bilin5

   subroutine sort(n         ,ra        ,wksp      ,iwksp     )
   !!--description-----------------------------------------------------------------
   ! Sorts an array, routine from Numerical Recipes
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   !
   implicit none
   !
   ! Global variables
   !
   integer                                     :: n
   integer, dimension(n)                       :: iwksp
   real*8   , dimension(n)              :: ra
   real*8   , dimension(n), intent(out) :: wksp
   !
   ! Local variables
   !
   integer :: j
   !
   !! executable statements -------------------------------------------------------
   !
   call indexx(n         ,ra        ,iwksp     )
   do j = 1, n
      wksp(j) = ra(iwksp(j))
   enddo
   end subroutine sort

   subroutine trapezoidal(x,y,n,x1_in,x2_in,integ)
   ! Compute integral over function y(x) from x1 to x2
   ! x is equidistant
   !(c)2014 Dano Roelvink
   implicit none
   integer,              intent(in)    :: n
   real*8, dimension(n), intent(in)    :: x,y          ! arrays x and y(x)
   real*8,               intent(in)    :: x1_in,x2_in  ! integration limits
   real*8,               intent(out)   :: integ	  ! integral
   real*8                              :: x1,y1,x2,y2,Ifirst,Imid,Ilast,dx
   integer                             :: i1,i2,i,i1p1,i2p1
   if (x1_in>x(n).or.x2_in<x(1)) then
      integ=0
   else

      x1=max(x1_in,x(1)+1d-60)
      x2=min(x2_in,x(n)-1d-60)
      dx=(x(n)-x(1))/(n-1)
      i1=floor((x1-x(1))/dx)+1
      i2=floor((x2-x(1))/dx)+1
      i1p1=min(i1+1,n)
      i2p1=min(i2+1,n)
      ! first partial trapezoid
      y1=y(i1)+(x1-x(i1))/dx*(y(i1p1)-y(i1))
      Ifirst=.5*(x(i1p1)-x1)*(y(i1p1)+y1)
      ! middle part
      Imid=.5*y(i1p1)
      do i=i1+2,i2-1
         Imid=Imid+y(i)
      end do
      Imid=Imid+.5*y(i2)
      Imid=Imid*dx
      ! last partial trapezoid
      y2=y(i2)+(x2-x(i2))/dx*(y(i2p1)-y(i2))
      Ilast=.5*(x2-x(i2))*(y2+y(i2))
      integ=Ifirst+Imid+Ilast
   endif
   end subroutine trapezoidal

   subroutine trapezoidal_cyclic(x,y,n,xcycle,x1,x2,integ)
   ! Compute integral over function y(x) from x1 to x2
   ! x is equidistant, function is cyclic so y(x+k*xcycle)=y(x)
   !(c)2014 Dano Roelvink
   implicit none
   integer,              intent(in)    :: n
   real*8, dimension(n), intent(in)    :: x,y          ! arrays x and y(x)
   real*8,               intent(in)    :: x1,x2,xcycle ! integration limits,cycle length
   real*8,               intent(out)   :: integ	    ! integral
   real*8, dimension(:),allocatable    :: xp,yp
   real*8                              :: dx
   integer                             :: ip,indt,np

   dx=x(2)-x(1)
   np=floor(x2/dx)-(floor(x1/dx)+1)+2
   allocate(xp(np))
   allocate(yp(np))
   xp(1)=x1
   do ip=2,np-1
      xp(ip)=(floor(x1/dx)+ip-1)*dx
   end do
   xp(np)=x2
   if (xcycle>0) then
      call interp_in_cyclic_function(x,y,n,xcycle,xp,np,yp)
   else
      do ip=1,np
         call linear_interp (x,y,n,xp(ip),yp(ip),indt)
      end do
   endif
   integ=0.d0
   do ip=1,np-1
      integ=integ+.5*(xp(ip+1)-xp(ip))*(yp(ip+1)+yp(ip))
   end do

   deallocate(xp)
   deallocate(yp)

   end subroutine trapezoidal_cyclic

   subroutine interp_using_trapez_rule(x,y,n,xcycle,xtarg,ytarg,ntarg)
   implicit none
   ! Compute integral over function y(x) from x1 to x2
   ! x is equidistant, function is cyclic so y(x+k*xcycle)=y(x)
   !(c)2014 Dano Roelvink
   integer,              intent(in)        :: n,ntarg
   real*8, dimension(n), intent(in)        :: x,y      ! arrays x and y(x)
   real*8,               intent(in)        :: xcycle   ! cycle length
   real*8, dimension(ntarg), intent(in)    :: xtarg    ! x values to interpolate to
   real*8, dimension(ntarg), intent(out)   :: ytarg    ! y values to interpolate to
   real*8                                  :: dx,x1,x2,integ
   integer                                 :: itarg

   dx=xtarg(2)-xtarg(1)
   do itarg=1,ntarg
      x1=xtarg(itarg)-.5*dx
      x2=xtarg(itarg)+.5*dx
      call trapezoidal_cyclic(x,y,n,xcycle,x1,x2,integ)
      ytarg(itarg)=integ/dx
   end do

   end subroutine interp_using_trapez_rule

   subroutine interp_in_cyclic_function(x,y,n,xcycle,xp,np,yp)
   implicit none
   integer,               intent(in)    :: n
   real*8, dimension(n),  intent(in)    :: x,y          ! arrays x and y(x)
   real*8, dimension(:), allocatable    :: xc,yc        ! complemented cyclic arrays
   real*8,                intent(in)    :: xcycle       ! cycle length
   integer,               intent(in)    :: np
   real*8, dimension(np), intent(in)    :: xp           ! points to interpolate to
   real*8, dimension(np), intent(out)   :: yp	         ! interpolated yp values
   integer                              :: icycle,ip,ileft,iright,nc,i
   real*8                               :: dx,yleft,yright,facleft,facright

   dx=x(2)-x(1)
   icycle=nint(xcycle/dx)
   allocate(xc(icycle+1))
   allocate(yc(icycle+1))
   if (n>icycle+1) then
      ! nonsense; error
   elseif (n==icycle+1) then
      ! e.g. 0, 30, 60,...360
      xc(1:n)=x
      yc(1:n)=y
   elseif (n==icycle) then
      ! e.g. 30, 60 ...360 or 15, 45 ...345 with n=12 and icycle=12
      ! interpolate between n-th and 1st value
      xc(1:n)=x
      yc(1:n)=y
      xc(n+1)=xc(n)+dx
      yc(n+1)=yc(1)
   elseif (n<icycle) then
      ! do not interpolate between n-th and 1st value; set values in gap to 0
      xc(1:n)=x
      yc(1:n)=y
      do i=n+1,icycle
         xc(i)=xc(i-1)+dx
         yc(i)=0.d0
      end do
      xc(icycle+1)=xc(icycle)+dx
      yc(icycle+1)=yc(1)
   endif

   nc=icycle+1

   do ip=1,np
      ileft=floor((xp(ip)-xc(1))/dx)+1
      do while (ileft<1)
         ileft=ileft+icycle
      end do
      do while (ileft>nc)
         ileft=ileft-icycle
      end do
      if (ileft>nc .or. ileft<1) then
         yleft=0
      else
         yleft=yc(ileft)
      endif
      iright=ileft+1
      do while (iright<1)
         iright=iright+icycle
      end do
      do while (iright>nc)
         iright=iright-icycle
      end do
      if (iright>nc .or. iright<1) then
         yright=0
      else
         yright=yc(iright)
      endif
      facright=mod(xp(ip),dx)/dx
      facleft=1.d0-facright
      yp(ip)=facleft*yleft+facright*yright

   end do
   deallocate(xc)
   deallocate(yc)

   end subroutine interp_in_cyclic_function

   subroutine triangle_intp(x,y,xt,yt,w)
   
   implicit none
   
   real*8, dimension(3), intent(in)      :: x,y
   real*8              , intent(in)      :: xt,yt
   real*8, dimension(3), intent(out)     :: w

   real*8                                :: den

   den= (y(2)-y(3))*(x(1)-x(3)) + (x(3)-x(2))*(y(1)-y(3))
   w(1)=((y(2)-y(3))*(xt  -x(3)) + (x(3)-x(2))*(yt  -y(3)))/den
   w(2)=((y(3)-y(1))*(xt  -x(3)) + (x(1)-x(3))*(yt  -y(3)))/den
   w(3)=1.d0-w(1)-w(2)

   end subroutine triangle_intp

! Do triangular interpolation using kdtree search   
subroutine TRIINTfast(XS_IN, YS_IN, ZS_IN, NS,NDIM,X,Y,Z,NXY, jsferic, dmiss, &
                      searchtol)
   implicit none
   !
   ! parameters
   integer,       intent(in)           :: ns, NDIM               !< sample vector dimension
   integer,       intent(in)           :: nxy
   integer,       intent(inout)        :: jsferic                   ! always 0 on return
   real(kind=hp), intent(in)           :: dmiss
   real(kind=hp), intent(in), target   :: XS_IN(ns), YS_IN(ns), ZS_IN(NDIM,ns) !< sample coordinates and values
   real(kind=hp), intent(in)           :: X(:), Y(:)                           !< node coordinates
   real(kind=hp), intent(inout)        :: Z(NDIM,nxy)                          !< node values
   real(kind=hp), intent(in)           :: searchtol
   !
   ! locals
   real(kind=hp)                       :: af
   integer                             :: i, ierr, k, intri
   integer                             :: n, nxx
   integer                             :: naf
   integer                             :: nbf
   integer                             :: ncf
   integer                             :: ncol
   integer                             :: nrfind
   integer                             :: n2
   integer                             :: idim
   real(kind=hp)                       :: rd
   real(kind=hp)                       :: slo(NDIM)
   real(kind=hp)                       :: xp, yp, zp(NDIM), xpmin, xpmax,ypmin, ypmax
   real(kind=hp)                       :: XL(3),YL(3)
   real(kind=hp), allocatable, target  :: xs1(:), ys1(:), zs1(:,:)      ! for store/restore of xs, ys, zs
   real(kind=hp), pointer              :: xs(:), ys(:), zs(:,:)
   real(kind=hp)                       :: xsmin, ysmin, xsmax, ysmax  ! bounding box corner coordinates
   integer                             :: indf(3)
   real(kind=hp)                       :: wf(3)
   integer                             :: NS1           ! for store/restore of NS
   integer                             :: jadum, ierror
   logical                             :: Ldeleteddata  ! for store/restore of xs, ys and zs
   integer                             :: KMOD
   type(kdtree_instance)               :: sampletree
   integer                             :: numsam, ii, k1, jakdtree, jakdtree2, jdla
   real(kind=hp)                       :: R2Search, dist2, cof1, cof2

   jakdtree=1
   jakdtree2 = jakdtree

   if ( jakdtree == 1 ) then
      !       enforce generation of kdtree
      treeglob%itreestat = ITREE_EMPTY
   end if

   !     JSLO=1, OOK SLOPES RD4
   IF (NS  <=  2) RETURN
   jdla=1
   IF (JDLA == 1) THEN

      NXX = NXY
      N2  = NS
      xs => xs_in
      ys => ys_in
      zs => zs_in

      !        determine sample bounding box
      xsmin = minval(xs,mask=xs /= DMISS)
      xsmax = maxval(xs,mask=xs /= DMISS)
      ysmin = minval(ys,mask=ys /= DMISS)
      ysmax = maxval(ys,mask=ys /= DMISS)

      ierr = 1

      call dlaun(xs,ys,N2,3,ierr)   ! generate edgeindex and triedge => moet leveren: numtri,

      if ( ierr /= 0 ) then
         goto 1234
      end if
      JDLA = 0
   ENDIF

   IF (NUMTRI  <  1) THEN
      goto 1234
   ENDIF

   R2search = 0d0
   if( searchtol /= dmiss ) R2search = searchtol**2
   if ( jakdtree2 == 1 .and. R2search > 0d0 ) then             ! build kd-tree for sample points
      sampletree%itreestat = ITREE_EMPTY
      call build_kdtree(sampletree, n2,xs,ys,ierror, jsferic, dmiss)
      if ( ierror /= 0 ) then
         if ( sampletree%itreestat /= ITREE_EMPTY ) call delete_kdtree2(sampletree)
         jakdtree2 = 0
      end if
   end if

   DO N = 1,NXx

      idim = 1 ! DMISS check on first sample dimension only

      XP = x(N)
      YP = y(N)
      RD = z(idim,N)

      INTRI = 0
      ! For triangulation, only consider points that are inside the sample bounding box
      if ( xp >= xsmin .and. xp <= xsmax .and. yp >= ysmin .and. yp <= ysmax ) then
         IF (RD == dmiss) then
            if ( jakdtree == 1 ) then
               call findtri_kdtree(XP,YP,ZP,XS,YS,ZS,N2,NDIM,NRFIND,INTRI,&
                                   ierror,indf,wf, dmiss, numtri, numedge,jsferic)
               if ( ierror /= 0 ) then
                  !                    deallocate
                  if ( treeglob%itreestat /= ITREE_EMPTY ) call delete_kdtree2(treeglob)
                  if ( allocated(imask) ) deallocate(imask)
                  if ( allocated(LNtri) ) deallocate(LNtri)
                  !                    disable kdtree
                  jakdtree = 0
               end if
            end if

            IF (INTRI == 1) THEN

               z(:,n)  = zp

               if (jagetwf == 1) then
                  indxx(1:3,n) = indf(1:3)
                  wfxx (1:3,n) = wf (1:3)
               endif

            ENDIF ! (INTRI == 1)
         ENDIF ! (XP  /=  XYMIS .AND. ... .AND. YP  <=  YMAXS )
      endif


      !!!!!!!!!! give it another try with nearest neighbour opr inverse distance.
      if (intri == 0 .and. R2search > 0d0) then

         if (RD == dmiss) then
            if( jakdtree2 == 1 ) then
               call make_queryvector_kdtree(sampletree, xp, yp, jsferic)
               numsam = kdtree2_r_count( sampletree%tree, sampletree%qv, R2search )
               if ( numsam > 0 ) then
                  call realloc_results_kdtree(sampletree, numsam)
                  call kdtree2_n_nearest(sampletree%tree,sampletree%qv,numsam,sampletree%results)
               end if

               do i = 1,idim
                  ii = 0
                  cof1 = 0d0
                  cof2 = 0d0
                  do k = 1,numsam
                     k1 = sampletree%results(k)%idx
                     if( abs( xp - xs(k1) ) < 1d-6 .and. abs( yp - ys(k1) ) < 1d-6 ) then
                        ii = 1
                        z(i,n) = zs(i,k1)
                        exit
                     endif

                     dist2 = ( xp - xs(k1) )**2 +  ( yp - ys(k1) )**2
                     cof1 = cof1 + zs(i,k1) / dist2
                     cof2 = cof2 + 1d0 / dist2
                  enddo
                  if( ii == 0 .and. numsam > 0 ) then
                     z(i,n) = cof1 / cof2
                  end if
               enddo
            endif
         endif ! RDMISS

      end if ! intri ==0

   enddo ! DO N = 1,NXx

1234 continue

   !     deallocate
   if ( jakdtree == 1 ) then
      if ( treeglob%itreestat    /=  ITREE_EMPTY ) call delete_kdtree2(treeglob  )
      if ( sampletree%itreestat  /=  ITREE_EMPTY ) call delete_kdtree2(sampletree)
      if ( allocated(imask) ) deallocate(imask)
      if ( allocated(LNtri) ) deallocate(LNtri)
      if ( allocated(edgeindx) ) deallocate(edgeindx)
      if ( allocated(triedge) ) deallocate(triedge)
      if ( allocated(indx) ) deallocate(indx)
   end if

   end subroutine triintfast   
   
 subroutine DLAUN(XS,YS,NS,jatri,ierr)
   implicit none
   real(kind=hp), intent(in) :: XS(:), YS(:)
   integer, intent(out)      :: ierr
   integer , intent(in)      :: ns
   integer, intent(in)       :: jatri !< Type of DLaun triangulation: 1: just triangulate,
   !! 3: also produce node-edge-triangle mapping tables
   !! for use in Triangulatesamplestonetwork.

   integer              :: maxtri
   integer, parameter   :: nh= 1   ! SPvdP: too small if jatri == 0
   integer              :: nsm
   real(kind=hp)        :: trisize
   real(kind=hp)        :: XH(NH), YH(NH)
   integer, allocatable :: idum(:)

   !     check memory
   call realloc(idum,50*ns,keepExisting=.false.,fill=0,stat=ierr)

   if ( ierr /= 0 ) then
      write(*,*) 'dlaun: out of memory'
      return
   end if

   deallocate(idum)


   NUMTRI = 0
   NUMEDGE = 0
   IF (NS  <  3) RETURN

   if (jatri /= 1 .and. jatri /= 3) then
      return
   end if

   numtri = -1
   do while ( numtri < 0 )
      NSM    = 6*NS + 10
      call realloc(indx,(/ 3,nsm /), keepExisting=.false.,fill=0)

      if (jatri==3) then
         call realloc(EDGEINDX, (/ 2,2*NSM /), keepExisting=.false., fill=0, stat=ierr)
         call realloc(TRIEDGE, (/ 3,NSM /), keepExisting=.false., fill=0, stat=ierr)
      else
         call realloc(EDGEINDX, (/ 2,1 /), keepExisting=.false., fill=0, stat=ierr)
         call realloc(TRIEDGE, (/ 3,1 /), keepExisting=.false., fill=0, stat=ierr)
      end if
      MAXTRI = NSM !?

      numtri = NSM ! Input value should specify max nr of triangles in indx.
      CALL TRICALL(jatri,XS,YS,NS,INDX,NUMTRI,EDGEINDX,NUMEDGE,TRIEDGE,XH,YH,NH,trisize)
      if ( numtri < 0 ) nsm = -numtri
   end do

   end subroutine DLAUN  
   
 subroutine findtri_kdtree(XP,YP,ZP,XS,YS,ZS,NS,NDIM,NRFIND,INTRI,ierror,ind, wf, dmiss, numtri, numedge, jsferic)
   implicit none

   real(kind=hp)   , intent(in)                :: xp, yp    !< node coordinates
   real(kind=hp)   , intent(out)               :: zp(NDIM)  !< node values
   integer,          intent(in)                :: NS        !< number of samples
   integer,          intent(in)                :: NDIM      !< sample vector dimension
   real(kind=hp)   , intent(in)                :: xs(NS), ys(NS), zs(NDIM,NS) !< sample coordinates and values
   integer,          intent(out)               :: NRFIND    !< triangle index
   integer,          intent(out)               :: intri     !< in triangle (1) or not (0)
   integer,          intent(out)               :: ierror    !< error (1) or not (0)
   integer,          intent(out)               :: ind(3)    !<
   real(kind=hp)   , intent(out)               :: wf(3)     !<
   real(kind=hp)   , intent(in)                :: dmiss
   integer,          intent(in)                :: numtri
   integer,          intent(in)                :: numedge
   integer,          intent(inout)             :: jsferic


   real(kind=hp)   , dimension(:), allocatable :: xx, yy !< triangle circumcenter coordinates, dim(numtri)

   real(kind=hp)   , dimension(3)      :: xv, yv    ! triangle node coordinates
   real(kind=hp)   , dimension(NDIM,3) :: zv        ! triangle node vectors


   real(kind=hp)                       :: xz, yz    ! triangle circumcenter
   real(kind=hp)                       :: SL,SM,XCR,YCR

   integer                             :: i, inod, ii, inext, k
   integer                             :: jacros, iothertriangle
   integer                             :: jsferic_store, ierr
   integer                             :: iedge, k1, k2, numsearched

   real(kind=hp)   , parameter         :: dtol = 1d-8
   integer,          parameter         :: MAXTRICON=100
   integer,          parameter         :: INIT_TRIS_PER_NODE=6

   real(kind=hp)   , external          :: dcosphi
   real(kind=hp)   , parameter         :: dfac = 1.000001d0  ! enlargement factor for pinpok3D

   ierror = 1

   NRFIND = 0

   if ( numtri <= 1 ) then
      !            call qnerror('findtri_kdtree: numtri<=1', ' ', ' ')
      goto 1234
   end if

   ! store jsferic
   jsferic_store = jsferic

   if ( treeglob%itreestat /= ITREE_OK ) then

      !           INITIALIZE kdetree2

      call realloc(xx,numtri, keepExisting=.false.,fill=0d0)
      call realloc(yy,numtri, keepExisting=.false.,fill=0d0)
      call realloc(imask,numtri, keepExisting=.false.,fill=0)
      IDENT = 0

      call realloc(LNtri,(/ 2,numedge /),  keepExisting=.false.,fill=0)
      LNtri = 0

      !           dlaun is not spherical proof: set global jsferic for circumcenter
      jsferic = 0

      do i=1,numtri
         do ii=1,3
            !                 generate edge-triangle connectivity
            iedge = triedge(ii,i)
            if ( LNtri(1,iedge) == 0 ) then
               LNtri(1,iedge) = i
            else
               LNtri(2,iedge) = i
            end if

            inod = indx(ii,i)

            !                 get triangle polygon
            xv(ii) = xs(inod)
            yv(ii) = ys(inod)
            zv(1:NDIM,ii) = zs(1:NDIM,inod)
         end do

         !              compute triangle circumcenter
         xx(i) = sum(xv(1:3)) / 3d0
         yy(i) = sum(yv(1:3)) / 3d0

      end do

      !        restore jsferic
      jsferic = jsferic_store

      call build_kdtree(treeglob, numtri, xx, yy, ierror,jsferic, dmiss)

      !           deallocate
      deallocate(xx, yy)

      if ( ierror /= 0 ) then
         goto 1234
      end if
   end if

   !        update identifier
   IDENT=IDENT+1

   !        fill query vector
   call make_queryvector_kdtree(treeglob,xp,yp,jsferic)

   !        get first triangle
   call kdtree2_n_nearest(treeglob%tree,treeglob%qv,1,treeglob%results)
   inext = treeglob%results(1)%idx

   !        perform a line search
   numsearched = 1
   i = 0

   mainloop:do while ( NRFIND == 0 .and. numsearched <= 2*NUMTRI .and. numsearched > 0 )

      inext = treeglob%results(1)%idx

      numsearched = 0

      do while ( inext > 0 .and. inext <= numtri .and. numsearched <= 2*NUMTRI )   ! numsearched: safety
         i = inext

         !           check current triangle
         do ii=1,3
            inod = indx(ii,i)
            xv(ii) = xs(inod)
            yv(ii) = ys(inod)
            zv(1:NDIM,ii) = zs(1:NDIM,inod)
         end do

!        get a point in the cell
            xz = sum(xv(1:3)) / 3d0
            yz = sum(yv(1:3)) / 3d0

         if ( imask(i) == IDENT ) then
            intri = 0
         else
            numsearched = numsearched+1
            !
            call ipon(xv,yv,3,xp,yp,intri)
            imask(i) = IDENT
         end if

         if ( intri == 1 ) then
            NRFIND = i
            exit mainloop
         end if

         !           dlaun is not spherical proof: set global jsferic
         jsferic = 0

         !           proceed to next triangle, which is adjacent to the edge that is cut by the line from the current triangle to the query point
         inext = 0

         do ii=1,3
            iedge=triedge(ii,i)
            if ( LNtri(2,iedge) == 0 ) cycle

            iothertriangle = LNtri(1,iedge) + LNtri(2,iedge) - i
            if ( imask(iothertriangle) == IDENT ) then
               cycle
            end if

            k1 = edgeindx(1,iedge)
            k2 = edgeindx(2,iedge)
            call CROSS(xz, yz, xp, yp, xs(k1), ys(k1), xs(k2), ys(k2), JACROS,SL,SM,XCR,YCR, jsferic, dmiss)

            !              use tolerance
            if ( jacros == 0 ) then
               if ( sm == dmiss .or. sl == dmiss ) then
                  !                    triangle with small area: sm and sl have not been computed
                  jacros = 1
               else IF (SM  >=  0d0-dtol .AND. SM  <=  1d0+dtol .AND. &
                  SL  >=  0d0-dtol .AND. SL  <=  1d0+dtol) THEN
               JACROS = 1
               ENDIF
            end if

            if ( jacros == 1 ) then  ! only select this edge if it has a second adjacent triangle
               inext = iothertriangle
               exit
            end if
         end do

         !        restore jsferic
         jsferic = jsferic_store

      end do

   end do mainloop

   if ( intri == 0 .and. inext > 0 ) then
      if ( numsearched > 10 ) write(6,"('numsearched= ', I0)") numsearched
   end if

   if (intri == 1) then
      call linear(xv, yv, zv, NDIM, xp, yp, zp, wf, dmiss, jsferic)
      do k = 1,3
         ind(k) = indx(k,nrfind)
      enddo

   endif

   ierror = 0
1234 continue

   return
   end subroutine findtri_kdtree   
 
! basically same subroutine as triangle_intp, but with zp calculated, and safe for tiny determinant 
! consider replacing above with below   
subroutine LINEAR ( X, Y, Z, NDIM, XP, YP, ZP, wf, dmiss, jsferic)
   implicit none
   integer      , intent(in)    :: NDIM   !< sample vector dimension
   real(kind=hp), intent(in)    :: X(:),Y(:),Z(:,:)
   real(kind=hp), intent(out)   :: wf(:)
   real(kind=hp), intent(out)   :: zp(:)
   real(kind=hp), intent(in)    :: dmiss
   integer,       intent(in)    :: jsferic
   real(kind=hp), intent(in)    :: xp
   real(kind=hp), intent(in)    :: yp

   real(kind=hp)    :: a11
   real(kind=hp)    :: a12
   real(kind=hp)    :: a21
   real(kind=hp)    :: a22
   real(kind=hp)    :: b1
   real(kind=hp)    :: b2
   real(kind=hp)    :: det
   real(kind=hp)    :: rlam
   real(kind=hp)    :: rmhu

   ZP  = dmiss
   A11 = getdx(x(1),y(1),x(2),y(2),jsferic)   ! X(2) - X(1)
   A21 = getdy(x(1),y(1),x(2),y(2),jsferic)   ! Y(2) - Y(1)
   A12 = getdx(x(1),y(1),x(3),y(3),jsferic)   ! X(3) - X(1)
   A22 = getdy(x(1),y(1),x(3),y(3),jsferic)   ! Y(3) - Y(1)
   B1  = getdx(x(1),y(1),xp  ,yp , jsferic)   ! XP   - X(1)
   B2  = getdy(x(1),y(1),xp  ,yp , jsferic)   ! YP   - Y(1)

   DET  =   A11 * A22 - A12 * A21
   IF (ABS(DET)  <  1E-12) THEN
      RETURN
   ENDIF

   RLAM = ( A22 * B1  - A12 * B2) / DET
   RMHU = (-A21 * B1  + A11 * B2) / DET
   wf(3) = rmhu
   wf(2) = rlam
   wf(1) = 1d0 - rlam - rmhu

   ZP   = Z(:,1) + RLAM * (Z(:,2) - Z(:,1)) + RMHU * (Z(:,3) - Z(:,1))
  
   return
end subroutine LINEAR   
   

double precision function getdx(x1,y1,x2,y2, jsferic)

   implicit none
   double precision                    :: x1, y1, x2, y2
   double precision                    :: xx1, yy1, xx2, yy2
   double precision                    :: diff1, diff2, dy, r, dx2
   integer, intent(in)                 :: jsferic
   double precision                    :: csphi

   if (jsferic == 1) then

      ! fix for poles
      diff1 = abs(abs(y1)-90d0)
      diff2 = abs(abs(y2)-90d0)
      if ( (diff1.le.dtol_pole .and. diff2.gt.dtol_pole) .or. &
         (diff1.gt.dtol_pole .and. diff2.le.dtol_pole) ) then
         getdx = 0d0
         return
      end if

      xx1   = x1
      xx2   = x2
      if      (xx1 - xx2 >  180d0) then
         xx1 = xx1 - 360d0
      else if (xx1 - xx2 < -180d0) then
         xx1 = xx1 + 360d0
      endif
      xx1   = xx1*degrad_hp
      xx2   = xx2*degrad_hp
      yy1   = y1 *degrad_hp
      yy2   = y2 *degrad_hp
      csphi = dcos(0.5d0*(yy1+yy2))
      getdx = earth_radius*csphi*(xx2-xx1)
   else
      getdx = x2-x1
   endif
end function getdx

double precision function getdy(x1,y1,x2,y2, jsferic)

   implicit none

   double precision :: x1, y1, x2, y2
   double precision :: xx1, yy1,yy2
   integer, intent(in)  :: jsferic

   if (jsferic == 1) then
      yy1   = y1*degrad_hp
      yy2   = y2*degrad_hp
      getdy = earth_radius*(yy2-yy1)
   else
      getdy = y2-y1
   endif
end function getdy


subroutine getdxdy(x1,y1,x2,y2,dx,dy, jsferic)

   implicit none

   double precision, intent(in) :: x1, y1, x2, y2
   double precision, intent(out):: dx, dy
   integer, intent(in)          :: jsferic
   
   double precision             :: dx2, dy2, dum

   if (Jsferic == 1) then
      dx = getdx(x1,y1,x2,y2,jsferic)
      dy = getdy(x1,y1,x2,y2,jsferic)
   else
      dx = x2-x1
      dy = y2-y1
   endif

end subroutine getdxdy


subroutine CROSS(x1, y1, x2, y2, x3, y3, x4, y4, JACROS,SL,SM,XCR,YCR, jsferic, dmiss)

   implicit none

   double precision, intent(in)    :: x1, y1, x2, y2, x3, y3, x4, y4   
   integer, intent(inout)          :: jacros
   double precision, intent(out)   :: sl
   double precision, intent(out)   :: sm
   integer, intent(in)             :: jsferic
   double precision, intent(in)    :: dmiss
   double precision, intent(out)   :: xcr, ycr
   
   double precision                :: x21, y21, x31, y31, x43, y43
   double precision                :: det
   double precision                :: eps

   JACROS = 0
   EPS    = 0.00001d0
   SL     = DMISS
   SM     = DMISS
   !     SL     = LABDA TUSSEN 0 EN 1 OP EERSTE PAAR
   !     Sm     = LABDA TUSSEN 0 EN 1 OP TWEEDE PAAR

   call getdxdy(x1,y1,x2,y2,x21,y21,jsferic)
   call getdxdy(x3,y3,x4,y4,x43,y43,jsferic)
   call getdxdy(x1,y1,x3,y3,x31,y31,jsferic)

   DET =  X43*Y21 - Y43*X21

   EPS = max(EPS*MAXVAL((/ X21,Y21,X43,Y43 /)), tiny(0d0))
   IF (ABS(DET) .LT. EPS) THEN
      xcr=dmiss
      ycr=dmiss
      RETURN
   ELSE
      SM = (Y31*X21 - X31*Y21) / DET
      SL = (Y31*X43 - X31*Y43) / DET
      IF (SM .GE. 0d0 .AND. SM .LE. 1d0 .AND. &
         SL .GE. 0d0 .AND. SL .LE. 1d0) THEN
         JACROS = 1
      ENDIF
      XCR    = X1 + SL*(X2-X1)
      YCR    = Y1 + SL*(Y2-Y1)
   ENDIF
   RETURN

   end subroutine CROSS

end module interp
