module snapwave_input
   use snapwave_log
   implicit none
   !
   ! Keywords requested by SnapWave; used to report unknown or misspelled keywords
   !
   integer, parameter, private :: max_known = 400
   character(len=64), private  :: known_keys(max_known)
   integer, private            :: nknown = 0
   !
   ! Input parameters as read (name, value, default yes/no), written to the map and his files
   !
   integer, public             :: npar = 0
   character(len=64), public   :: par_name(max_known)
   character(len=256), public  :: par_value(max_known)
   logical, public             :: par_default(max_known)
   !
   private :: find_keyword, log_parameter, register_keyword, check_unused_keywords, suggest_keyword, lower, &
              check_input_consistency
   !
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

      call writelog('ls', 'Reading input file ...')

      do ii = 1, size(possible_names)
         filename = trim(possible_names(ii))
         inquire (file=filename, exist=exists, iostat=ios)
         if (exists .and. ios == 0) exit
      end do

      if (.not. exists) then
         call halt_program('input file snapwave.inp not found in the run directory')
      end if

      open (unit=500, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call halt_program('found input file '//trim(filename)//' but failed to open it')
      end if

      call writelog('ls', '   Opened '//trim(filename))
      call writelog('ls', '')
      call writelog('ls', 'Input parameters (values not given in '//trim(filename)//' are defaults, marked (default)):')
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
      call read_char_input(500, 'map_file', map_filename, '')
      call read_char_input(500, 'his_file', his_filename, '')
      call read_real_input(500, 'map_interval', map_interval, timestep)
      call read_real_input(500, 'his_interval', his_interval, timestep)
      if (map_filename /= '' .and. map_interval <= 0.0) then
         call halt_program('map_interval must be positive')
      end if
      if (his_filename /= '' .and. his_interval <= 0.0) then
         call halt_program('his_interval must be positive')
      end if
      call read_int_input(500, 'map_depth', map_dep, 1)
      call read_int_input(500, 'map_Hm0', map_Hm0, 1)
      call read_int_input(500, 'map_Hig', map_Hig, 0)
      call read_int_input(500, 'map_Tp', map_Tp, 1)
      call read_int_input(500, 'map_dir', map_dir, 1)
      call read_int_input(500, 'map_dirspr', map_dirspr, 0)
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
         call writelog('ls', '   Wind growth off: uniform wave period in entire domain')
         !
      else
         !
         call writelog('ls', '   Wind growth turned on')
         wind = 1
      end if
      !
      ! Report lines in the input file that are not used
      !
      call check_unused_keywords(500, filename)
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
      call check_input_consistency()
      !
   end subroutine

   subroutine find_keyword(fileid, keyword, valstr, found)
      !
      ! Find "keyword = value" in the input file; returns the value string.
      ! Leading blanks before the keyword are allowed; the first occurrence is used.
      !
      integer, intent(in) :: fileid
      character(*), intent(in) :: keyword
      character(len=256), intent(out) :: valstr
      logical, intent(out) :: found
      character(len=256) :: line
      character(len=256) :: keystr
      integer :: j, stat
      !
      call register_keyword(keyword)
      found = .false.
      valstr = ''
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat /= 0) exit
         j = index(line, '=')
         if (j <= 1) cycle
         keystr = adjustl(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = adjustl(line(j + 1:))
            found = .true.
            exit
         end if
      end do
      !
   end subroutine find_keyword

   subroutine log_parameter(keyword, valstr, found)
      !
      ! Echo a parameter to the log: user-specified values also to screen, defaults to log file only
      !
      character(*), intent(in) :: keyword
      character(*), intent(in) :: valstr
      logical, intent(in) :: found
      character(len=20) :: kw
      !
      kw = keyword
      if (npar < max_known) then
         npar = npar + 1
         par_name(npar) = keyword
         par_value(npar) = valstr
         par_default(npar) = .not. found
      end if
      if (found) then
         call writelog('ls', '   '//kw//' = '//trim(valstr))
      else
         call writelog('l', '   '//kw//' = '//trim(valstr)//'   (default)')
      end if
      !
   end subroutine log_parameter

   subroutine read_real_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: valstr
      integer, intent(in) :: fileid
      real*4, intent(out) :: value
      real*4, intent(in) :: default
      logical :: found
      integer :: ios
      !
      value = default
      call find_keyword(fileid, keyword, valstr, found)
      if (found) then
         read (valstr, *, iostat=ios) value
         if (ios /= 0) call halt_program('could not read a real number for keyword '//trim(keyword)// &
                                         ' from "'//trim(valstr)//'"')
      end if
      call log_parameter(keyword, num2str(value), found)
      !
   end subroutine

   subroutine read_real_array_input(fileid, keyword, value, default, nr)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: valstr
      integer, intent(in) :: fileid
      integer, intent(in) :: nr
      real*4, dimension(:), intent(out), allocatable :: value
      real*4, intent(in) :: default
      logical :: found
      integer :: m, ios
      !
      allocate (value(nr))
      !
      value = default
      call find_keyword(fileid, keyword, valstr, found)
      if (found) then
         read (valstr, *, iostat=ios) (value(m), m=1, nr)
         if (ios /= 0) call halt_program('could not read '//num2str(nr)//' real numbers for keyword '// &
                                         trim(keyword)//' from "'//trim(valstr)//'"')
         call log_parameter(keyword, valstr, found)
      else
         call log_parameter(keyword, num2str(default), found)
      end if
      !
   end subroutine

   subroutine read_int_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: valstr
      integer, intent(in) :: fileid
      integer, intent(out) :: value
      integer, intent(in) :: default
      logical :: found
      integer :: ios
      !
      value = default
      call find_keyword(fileid, keyword, valstr, found)
      if (found) then
         read (valstr, *, iostat=ios) value
         if (ios /= 0) call halt_program('could not read an integer for keyword '//trim(keyword)// &
                                         ' from "'//trim(valstr)//'"')
      end if
      call log_parameter(keyword, num2str(value), found)
      !
   end subroutine

   subroutine read_char_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: valstr
      integer, intent(in) :: fileid
      character(*), intent(in) :: default
      character(*), intent(out) :: value
      logical :: found
      !
      value = default
      call find_keyword(fileid, keyword, valstr, found)
      if (found) value = valstr
      call log_parameter(keyword, value, found)
      !
   end subroutine

   subroutine register_keyword(keyword)
      character(*), intent(in) :: keyword
      integer :: i
      do i = 1, nknown
         if (trim(known_keys(i)) == trim(keyword)) return
      end do
      if (nknown < max_known) then
         nknown = nknown + 1
         known_keys(nknown) = keyword
      end if
   end subroutine register_keyword

   subroutine check_unused_keywords(fileid, filename)
      !
      ! Warn about lines in snapwave.inp that SnapWave does not use:
      ! lines without '=', unknown (e.g. misspelled) keywords and keywords given more than once.
      !
      integer, intent(in) :: fileid
      character(*), intent(in) :: filename
      character(len=256) :: line, keystr
      character(len=64)  :: seen(max_known)
      integer :: j, i, stat, iline, nseen
      logical :: known, dup
      !
      nseen = 0
      iline = 0
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat /= 0) exit
         iline = iline + 1
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (scan(line(1:1), '#!%*') > 0) cycle          ! comment lines
         j = index(line, '=')
         if (j == 0) then
            call writewarning(trim(filename)//' line '//num2str(iline)//' has no "=" and is ignored: "'// &
                              trim(line)//'"')
            cycle
         end if
         keystr = adjustl(line(1:j - 1))
         known = .false.
         do i = 1, nknown
            if (trim(known_keys(i)) == trim(keystr)) then
               known = .true.
               exit
            end if
         end do
         if (.not. known) then
            call writewarning(trim(filename)//' line '//num2str(iline)//': unknown keyword "'//trim(keystr)// &
                              '" is ignored'//trim(suggest_keyword(keystr)))
            cycle
         end if
         dup = .false.
         do i = 1, nseen
            if (trim(seen(i)) == trim(keystr)) dup = .true.
         end do
         if (dup) then
            call writewarning(trim(filename)//' line '//num2str(iline)//': keyword "'//trim(keystr)// &
                              '" is given more than once; only the first value is used')
         elseif (nseen < max_known) then
            nseen = nseen + 1
            seen(nseen) = keystr
         end if
      end do
      !
   end subroutine check_unused_keywords

   function suggest_keyword(keystr) result(hint)
      ! suggest a known keyword that differs only in upper/lower case
      character(*), intent(in) :: keystr
      character(len=:), allocatable :: hint
      integer :: i
      hint = ''
      do i = 1, nknown
         if (lower(trim(known_keys(i))) == lower(trim(keystr))) then
            hint = ' (did you mean "'//trim(known_keys(i))//'"?)'
            return
         end if
      end do
   end function suggest_keyword

   function lower(s) result(o)
      character(*), intent(in) :: s
      character(len=len(s)) :: o
      integer :: i, c
      o = s
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
      end do
   end function lower

   subroutine check_input_consistency()
      !
      ! Basic checks on input values and on the existence of input files
      !
      use snapwave_data
      !
      implicit none
      !
      if (timestep <= 0.0) call halt_program('timestep must be positive')
      if (tstop < tstart) call halt_program('tstop ('//trim(tstopstr)//') is before tstart ('//trim(tstartstr)//')')
      if (dtheta <= 0.0) call halt_program('dtheta must be positive')
      if (sector <= 0.0 .or. sector > 360.0) call halt_program('sector must be between 0 and 360 degrees')
      if (hmin <= 0.0) call halt_program('hmin must be positive')
      if (gamma <= 0.0) call halt_program('gamma must be positive')
      if (niter < 1) call halt_program('niter must be at least 1')
      if (sferic /= 0 .and. sferic /= 1) call halt_program('sferic must be 0 or 1')
      !
      ! Grid
      !
      if (mskfile == '') then
         call file_exists_or_halt('gridfile', gridfile)
      else
         call file_exists_or_halt('mskfile', mskfile)
         call file_exists_or_halt('indfile', indfile)
         call file_exists_or_halt('depfile', depfile)
      end if
      !
      ! Wave boundary conditions
      !
      if (jonswapfile /= '') then
         call file_exists_or_halt('jonswapfile', jonswapfile)
      elseif (bndfile(1:4) /= 'none') then
         call file_exists_or_halt('bndfile', bndfile)
         if (bhsfile == '') call halt_program('bndfile is given but bhsfile is missing')
         if (btpfile == '') call halt_program('bndfile is given but btpfile is missing')
         if (bwdfile == '') call halt_program('bndfile is given but bwdfile is missing')
         if (bdsfile == '') call halt_program('bndfile is given but bdsfile is missing')
         if (bzsfile == '') call halt_program('bndfile is given but bzsfile is missing')
         call file_exists_or_halt('bhsfile', bhsfile)
         call file_exists_or_halt('btpfile', btpfile)
         call file_exists_or_halt('bwdfile', bwdfile)
         call file_exists_or_halt('bdsfile', bdsfile)
         call file_exists_or_halt('bzsfile', bzsfile)
      else
         call writewarning('no wave boundary conditions given (neither jonswapfile nor bndfile)')
      end if
      if (encfile(1:4) /= 'none') then
         call file_exists_or_halt('encfile', encfile)
      elseif (jonswapfile /= '' .or. bndfile(1:4) /= 'none') then
         call writewarning('no encfile given: no grid nodes will be marked as wave boundary')
      end if
      if (neumannfile(1:4) /= 'none') call file_exists_or_halt('neumannfile', neumannfile)
      if (obsfile(1:4) /= 'none') call file_exists_or_halt('obsfile', obsfile)
      if (windlistfile /= '') call file_exists_or_halt('windlistfile', windlistfile)
      if (upwfile /= '') call file_exists_or_halt('upwfile', upwfile)
      if (ja_vegetation == 1) call file_exists_or_halt('vegmapfile', vegmapfile)
      !
      call writelog('ls', '')
      call writelog('ls', 'Simulation period: '//trim(tstartstr)//' to '//trim(tstopstr)// &
                          ' (tref '//trim(trefstr)//'), '//num2str(nint((tstop - tstart) / timestep) + 1)// &
                          ' time steps of '//num2str(timestep)//' s')
      call writelog('ls', '')
      !
   end subroutine check_input_consistency

end module snapwave_input
