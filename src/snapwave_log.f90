module snapwave_log
   !
   ! Logging for SnapWave, modelled on the XBeach logging module.
   !
   ! Files (in the run directory):
   !    snapwave_log.txt      everything written with destination 'l'
   !    snapwave_warning.txt  warnings (only created when there is a warning)
   !    snapwave_error.txt    the error that stopped the run (only created on error)
   !
   ! Usage:
   !    call writelog('ls', 'Reading mesh ...')               ! log file + screen
   !    call writelog('l',  '   nodes : '//num2str(no_nodes))  ! log file only
   !    call writewarning('boundary time series ends before tstop')
   !    call halt_program('bndfile '//trim(bndfile)//' contains no points')
   !
   ! Destination characters: l = log file, s = screen, w = warning file, e = error file
   !
   implicit none
   !
   private
   public :: open_logfiles, close_logfiles, writelog, writewarning, halt_program, num2str, &
             file_exists_or_halt, nwarnings
   !
   character(len=*), parameter :: logfilename  = 'snapwave_log.txt'
   character(len=*), parameter :: warnfilename = 'snapwave_warning.txt'
   character(len=*), parameter :: errfilename  = 'snapwave_error.txt'
   !
   integer, save :: lunlog  = 0      ! unit numbers from open(newunit=...) are negative,
   integer, save :: lunwarn = 0      ! so open/closed state is kept in separate flags
   integer, save :: lunerr  = 0
   logical, save :: log_open  = .false.
   logical, save :: warn_open = .false.
   logical, save :: err_open  = .false.
   integer, save :: nwarnings = 0
   integer, parameter :: max_screen_warnings = 50   ! further warnings go to the files only
   integer(kind=8), save :: clock_start = 0
   !
   interface num2str
      module procedure int2str, int8_2str, real4_2str, real8_2str
   end interface num2str
   !
contains
   !
   subroutine open_logfiles()
      !
      integer :: ios
      character(len=8)  :: d
      character(len=10) :: tm
      !
      call system_clock(clock_start)
      !
      open (newunit=lunlog, file=logfilename, status='replace', action='write', iostat=ios)
      log_open = (ios == 0)
      if (.not. log_open) then
         write (*, '(a)') 'WARNING: could not open '//logfilename//', logging to screen only'
      end if
      !
      ! Remove warning and error files from a previous run
      !
      call delete_file(warnfilename)
      call delete_file(errfilename)
      !
      call date_and_time(d, tm)
      call writelog('ls', '**********************************************************************')
      call writelog('ls', '                         SnapWave')
      call writelog('ls', '   Stationary directional wave model for unstructured grids')
      call writelog('ls', '**********************************************************************')
      call writelog('ls', 'Simulation started '//d(1:4)//'-'//d(5:6)//'-'//d(7:8)//' '// &
                           tm(1:2)//':'//tm(3:4)//':'//tm(5:6))
      call writelog('ls', '')
      !
   end subroutine open_logfiles
   !
   subroutine close_logfiles()
      !
      character(len=8)  :: d
      character(len=10) :: tm
      !
      call date_and_time(d, tm)
      call writelog('ls', '')
      call writelog('ls', '**********************************************************************')
      call writelog('ls', 'Simulation finished '//d(1:4)//'-'//d(5:6)//'-'//d(7:8)//' '// &
                           tm(1:2)//':'//tm(3:4)//':'//tm(5:6))
      call writelog('ls', 'Total wall clock time : '//num2str(elapsed_seconds())//' s')
      if (nwarnings > 0) then
         call writelog('ls', 'Number of warnings    : '//num2str(nwarnings)//'  (see '//warnfilename//')')
      else
         call writelog('ls', 'Number of warnings    : 0')
      end if
      call writelog('ls', '**********************************************************************')
      !
      if (log_open) close (lunlog)
      if (warn_open) close (lunwarn)
      log_open = .false.
      warn_open = .false.
      !
   end subroutine close_logfiles
   !
   subroutine writelog(destination, message)
      !
      ! destination: any combination of 'l' (log file), 's' (screen), 'w' (warning file), 'e' (error file)
      !
      character(len=*), intent(in) :: destination
      character(len=*), intent(in) :: message
      integer :: ios
      !
      if (scan(destination, 's') > 0) then
         write (*, '(a)') trim(message)
      end if
      !
      if (scan(destination, 'l') > 0 .and. log_open) then
         write (lunlog, '(a)') trim(message)
         flush (lunlog)                       ! keep the log complete if the program crashes
      end if
      !
      if (scan(destination, 'w') > 0) then
         if (.not. warn_open) then
            open (newunit=lunwarn, file=warnfilename, status='replace', action='write', iostat=ios)
            warn_open = (ios == 0)
         end if
         if (warn_open) then
            write (lunwarn, '(a)') trim(message)
            flush (lunwarn)
         end if
      end if
      !
      if (scan(destination, 'e') > 0) then
         if (.not. err_open) then
            open (newunit=lunerr, file=errfilename, status='replace', action='write', iostat=ios)
            err_open = (ios == 0)
         end if
         if (err_open) then
            write (lunerr, '(a)') trim(message)
            flush (lunerr)
         end if
      end if
      !
   end subroutine writelog
   !
   subroutine writewarning(message)
      !
      character(len=*), intent(in) :: message
      !
      nwarnings = nwarnings + 1
      if (nwarnings <= max_screen_warnings) then
         call writelog('lsw', 'WARNING: '//trim(message))
      else
         if (nwarnings == max_screen_warnings + 1) then
            call writelog('s', 'WARNING: more than '//num2str(max_screen_warnings)// &
                               ' warnings; further warnings only in '//warnfilename)
         end if
         call writelog('lw', 'WARNING: '//trim(message))
      end if
      !
   end subroutine writewarning
   !
   subroutine halt_program(message)
      !
      ! Write an error message to screen, log file and error file, and stop.
      !
      character(len=*), intent(in) :: message
      !
      call writelog('lse', '')
      call writelog('lse', 'ERROR: '//trim(message))
      call writelog('lse', '')
      call writelog('ls', 'SnapWave stopped with an error, see '//errfilename)
      if (log_open) close (lunlog)
      if (warn_open) close (lunwarn)
      if (err_open) close (lunerr)
      stop 1
      !
   end subroutine halt_program
   !
   subroutine file_exists_or_halt(keyword, filename)
      !
      ! Stop with a clear message if an input file given in snapwave.inp does not exist.
      !
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: filename
      logical :: exists
      !
      inquire (file=trim(filename), exist=exists)
      if (.not. exists) then
         call halt_program('file "'//trim(filename)//'" given for keyword '//trim(keyword)//' does not exist')
      end if
      !
   end subroutine file_exists_or_halt
   !
   function elapsed_seconds() result(s)
      real*8 :: s
      integer(kind=8) :: c, rate
      call system_clock(c, rate)
      if (rate > 0) then
         s = dble(c - clock_start) / dble(rate)
      else
         s = 0.d0
      end if
   end function elapsed_seconds
   !
   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      logical :: exists
      integer :: lun, ios
      inquire (file=filename, exist=exists)
      if (exists) then
         open (newunit=lun, file=filename, status='old', iostat=ios)
         if (ios == 0) close (lun, status='delete')
      end if
   end subroutine delete_file
   !
   ! ---- number to string conversion (generic num2str) ----
   !
   function int2str(i) result(str)
      integer, intent(in) :: i
      character(len=:), allocatable :: str
      character(len=24) :: buf
      write (buf, '(i0)') i
      str = trim(buf)
   end function int2str
   !
   function int8_2str(i) result(str)
      integer(kind=8), intent(in) :: i
      character(len=:), allocatable :: str
      character(len=24) :: buf
      write (buf, '(i0)') i
      str = trim(buf)
   end function int8_2str
   !
   function real4_2str(r) result(str)
      real*4, intent(in) :: r
      character(len=:), allocatable :: str
      str = real8_2str(dble(r))
   end function real4_2str
   !
   function real8_2str(r) result(str)
      real*8, intent(in) :: r
      character(len=:), allocatable :: str
      character(len=32) :: buf
      integer :: n
      if (r == 0.d0) then
         buf = '0.0'
      elseif (abs(r) >= 1.d-3 .and. abs(r) < 1.d7) then
         write (buf, '(f32.6)') r
         buf = adjustl(buf)
         ! remove trailing zeros after the decimal point (keep one digit)
         n = len_trim(buf)
         do while (n > 2)
            if (buf(n:n) /= '0' .or. buf(n - 1:n - 1) == '.') exit
            buf(n:n) = ' '
            n = n - 1
         end do
      else
         write (buf, '(es32.5)') r
         buf = adjustl(buf)
      end if
      str = trim(buf)
   end function real8_2str
   !
end module snapwave_log
