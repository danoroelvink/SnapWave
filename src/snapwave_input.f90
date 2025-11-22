module snapwave_input
   implicit none
contains

   subroutine read_snapwave_input()
      !
      ! Reads snapwave.inp
      !
      use snapwave_data
      use snapwave_date
      !
      implicit none
      !
      integer :: dtsec
      integer :: irestart
      integer :: iwritetestfiles
      !
      character(len=256) :: filename
      integer :: ios, ii
      logical :: exists

! List of possible reasonable filenames
      character(len=*), parameter :: possible_names(*) = [ &
                                     'snapwave.inp', 'SnapWave.inp', 'SNAPWAVE.INP', &
                                     'snapwave.INP', 'Snapwave.INP', 'SNAPWAVE.inp']

      write (*, *) 'Reading input file ...'

      do ii = 1, size(possible_names)
         filename = trim(possible_names(ii))
         inquire (file=filename, exist=exists, iostat=ios)
         if (exists .and. ios == 0) exit
      end do

      if (.not. exists) then
         write (*, *) 'ERROR: none of the expected input files were found:'
         do ii = 1, size(possible_names)
            write (*, *) '   - ', trim(possible_names(ii))
         end do
         stop 1
      end if

      open (unit=500, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: found file but failed to open: ', trim(filename)
         stop 1
      end if

      write (*, *) 'Successfully opened: ', trim(filename)
      !
      ! Input section
      !
      call read_int_input(500, 'nmax', nmax, 0)
      call read_int_input(500, 'mmax', mmax, 0)
      call read_real_input(500, 'dx', dx, 0.0)
      call read_real_input(500, 'dy', dy, 0.0)
      call read_real_input(500, 'x0', x0, 0.0)
      call read_real_input(500, 'y0', y0, 0.0)
      call read_real_input(500, 'rotation', rotation, 0.0)
      call read_real_input(500, 'posdwn', posdwn, -1.0)
      call read_char_input(500, 'tref', trefstr, '20000101 000000')
      call read_char_input(500, 'tstart', tstartstr, '20000101 000000')
      call read_char_input(500, 'tstop', tstopstr, '20000101 000000')
      call read_real_input(500, 'timestep', timestep, 3600.0)
      call read_int_input(500, 'niter', niter, 10)
      call read_real_input(500, 'crit', crit, 0.00001)
      call read_real_input(500, 'dt', dt, 36000.0)
      call read_real_input(500, 'gamma', gamma, 0.7)
      call read_real_input(500, 'alpha', alpha, 1.0)
      call read_real_input(500, 'hmin', hmin, 0.1)
      call read_real_input(500, 'gammax', gammax, 0.6)
      call read_char_input(500, 'gridfile', gridfile, '.txt')
      call read_int_input(500, 'sferic', sferic, 0)
      call read_char_input(500, 'fw', fwstr, '0.01')
      call read_char_input(500, 'fwig', fw_igstr, '0.015')
      call read_real_input(500, 'fwcutoff', fwcutoff, 200.0)
      call read_real_input(500, 'tol', tol, 10.0)
      call read_real_input(500, 'dtheta', dtheta, 10.0)
      call read_real_input(500, 'sector', sector, 180.0)
      call read_char_input(500, 'jonswapfile', jonswapfile, '')
      call read_char_input(500, 'windlistfile', windlistfile, '')
      call read_char_input(500, 'bndfile', bndfile, 'none')
      call read_char_input(500, 'encfile', encfile, 'none')
      call read_char_input(500, 'neumannfile', neumannfile, 'none')
      call read_char_input(500, 'bhsfile', bhsfile, '')
      call read_char_input(500, 'btpfile', btpfile, '')
      call read_char_input(500, 'bwdfile', bwdfile, '')
      call read_char_input(500, 'bdsfile', bdsfile, '')
      call read_char_input(500, 'bzsfile', bzsfile, '')
      call read_char_input(500, 'upwfile', upwfile, '')
      call read_char_input(500, 'mskfile', mskfile, '')
      call read_char_input(500, 'indfile', indfile, '')
      call read_char_input(500, 'depfile', depfile, '')
      call read_char_input(500, 'obsfile', obsfile, 'none')
      call read_char_input(500, 'outputformat', outputformat, 'bin')
      call read_char_input(500, 'map_file', map_filename, '')
      call read_char_input(500, 'his_file', his_filename, '')
      call read_int_input(500, 'map_depth', map_dep, 1)
      call read_int_input(500, 'map_Hm0', map_Hm0, 1)
      call read_int_input(500, 'map_Hig', map_Hig, 0)
      call read_int_input(500, 'map_Tp', map_Tp, 1)
      call read_int_input(500, 'map_dir', map_dir, 1)
      call read_int_input(500, 'map_Cg', map_Cg, 0)
      call read_int_input(500, 'map_Dw', map_Dw, 0)
      call read_int_input(500, 'map_Df', map_Df, 0)
      call read_int_input(500, 'map_SwE', map_SwE, 0)
      call read_int_input(500, 'map_SwA', map_SwA, 0)
      call read_int_input(500, 'map_sig', map_sig, 0)
      call read_int_input(500, 'map_u10', map_u10, 0)
      call read_int_input(500, 'map_Dveg', map_Dveg, 0)
      call read_int_input(500, 'writetestfiles', iwritetestfiles, 0)
      call read_int_input(500, 'ja_save_each_iter', ja_save_each_iter, 0)

      call read_int_input(500, 'map_ee', map_ee, 0)
      call read_int_input(500, 'map_ctheta', map_ctheta, 0)
      call read_int_input(500, 'restart', irestart, 0)
      !
      call read_char_input(500, 'u10', u10str, '0.0')
      call read_char_input(500, 'u10dir', u10dirstr, '270.0')
      call read_real_input(500, 'Tpini', Tpini, 1.0)
      call read_int_input(500, 'mwind', mwind, 2)
      call read_real_input(500, 'sigmin', sigmin, 8.0 * atan(1.0) / 25.0)
      call read_real_input(500, 'sigmax', sigmax, 8.0 * atan(1.0) / 1.0)
      call read_int_input(500, 'jadcgdx', jadcgdx, 1)
      call read_real_input(500, 'c_dispT', c_dispT, 1.0)
      call read_real_input(500, 'zsini', zsini, 0.0)
      call read_int_input(500, 'ig', ig, 0)
      call read_int_input(500, 'upwindref', upwindref, 0)
      !
      ! Vegetation input
      !
      call read_int_input(500, 'ja_vegetation', ja_vegetation, 0)
      call read_char_input(500, 'vegmapfile', vegmapfile, '.txt')
      !
      wind = 0
      if ((u10str == '0.0') .and. (windlistfile == '')) then
         !
         wind = 0
         !
         write (*, *) '   Uniform wave period in entire domain.'
         !
      else
         !
         write (*, *) '   Wind growth turned on.'
         wind = 1
      end if
      !
      close (500)
      !
      call time_difference(trefstr, tstartstr, dtsec) ! time difference in seconds between tstart and tref
      tstart = dtsec * 1.0 ! time difference in seconds between tstop and tstart
      call time_difference(trefstr, tstopstr, dtsec)
      tstop = dtsec * 1.0 ! time difference in seconds between tstop and tstart
      !
      mmax = mmax + 2 ! Original mmax and nmax are for number of cells in bathy grid. Add two dummy rows.
      nmax = nmax + 2
      !
      restart = .true.
      writetestfiles = .true.
      if (irestart == 0) restart = .false.
      if (iwritetestfiles == 0) writetestfiles = .false.
      !
   end subroutine

   subroutine read_real_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      real*4, intent(out) :: value
      real*4, intent(in) :: default
      integer j, stat
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) value
            exit
         end if
      end do
      !
   end subroutine

   subroutine read_real_array_input(fileid, keyword, value, default, nr)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      integer, intent(in) :: nr
      real*4, dimension(:), intent(out), allocatable :: value
      real*4, intent(in) :: default
      integer j, stat, m
      !
      allocate (value(nr))
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) (value(m), m=1, nr)
            exit
         end if
      end do
      !
   end subroutine

   subroutine read_int_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      integer, intent(out) :: value
      integer, intent(in) :: default
      integer j, stat
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) value
            exit
         end if
      end do
      !
   end subroutine

   subroutine read_char_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      character(*), intent(in) :: default
      character(*), intent(out) :: value
      integer j, stat
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = adjustl(trim(line(j + 1:256)))
            value = valstr
            exit
         end if
      end do
      !
   end subroutine

end module snapwave_input
