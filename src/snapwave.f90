program snapwave
   !
   ! Stand-alone program for implicit stationary wave solver
   ! Based on wave enrgy balance for directionally spread waves, single representative frequency
   ! 4-sweep implicit method
   ! This method can be used for structured 2D grids such as XBeach or unstructured 2D grids such as Delft3D-FM; subroutines find_upwind_neighbours and
   ! solve_energy_balance2Dstat are for unstructured grids. This program reads an unstructured gridf from an ascii mesh grid file.
   !
   ! (c) 2020 Dano Roelvink, IHE Delft
   !
   use snapwave_data
   use snapwave_input
   use snapwave_domain
   use snapwave_boundaries
   use snapwave_solver
   use snapwave_ncoutput
   use snapwave_obspoints
   use snapwave_results
   use snapwave_log
   use snapwave_date, only: date_after_seconds
   use omp_lib
   !
   implicit none
   !
   real*8  :: t
   real*8  :: output_tol
   !
   integer :: it, nt, n_notconv
   integer(kind=8) :: clock0, clock1, clock_rate
   !
   call open_logfiles()                  ! snapwave_log.txt (+ warning/error files when needed)
   !
   call read_snapwave_input()            ! Reads snapwave.inp      
   !
   call initialize_snapwave_domain()     ! Read mesh, finds upwind neighbors, etc.
   !
   call read_obs_points()
   !
   ! Read boundary conditions
   !
   call read_boundary_data()
   !
   ! read wind data if specified
   !
   call read_wind_data()
   !
   ! Initialize NetCDF output
   !
   call ncoutput_init()
   !
   ! Start time loop
   !
   !$omp parallel
   !$omp single
   call writelog('ls', 'Running with '//num2str(omp_get_num_threads())//' OpenMP threads')
   !$omp end single
   !$omp end parallel
   it = 0
   t  = tstart
   map_output_count = 0
   his_output_count = 0
   next_map_output = tstart
   next_his_output = tstart
   output_tol = max(1.0d-6, abs(dble(timestep))*1.0d-6)
   !
   nt = floor((tstop - tstart) / dble(timestep)) + 1
   n_notconv = 0
   call writelog('ls', '')
   call writelog('ls', 'Start time loop: '//num2str(nt)//' time steps')
   call system_clock(clock0, clock_rate)
   do while (t<=tstop)
      !
      ! New time step
      !
      it = it + 1
      !
      call update_boundary_conditions(t) ! includes theta_grid creation
      !
      call compute_wave_field(t)   
      !
      call system_clock(clock1)
      if (.not. solver_converged) n_notconv = n_notconv + 1
      call report_time_step(t, it, nt, dble(clock1 - clock0) / dble(max(clock_rate, 1_8)))
      !
      if (his_filename /= '' .and. nobs > 0 .and. t >= next_his_output - output_tol) then
         his_output_count = his_output_count + 1
         call update_obs_points()
         call ncoutput_update_his(t, his_output_count)
         do while (next_his_output <= t + output_tol)
            next_his_output = next_his_output + dble(his_interval)
         end do
      end if
      !
      if (ja_save_each_iter == 0 .and. map_filename /= '' .and. t >= next_map_output - output_tol) then
         map_output_count = map_output_count + 1
         call ncoutput_update_map(t, map_output_count)
         do while (next_map_output <= t + output_tol)
            next_map_output = next_map_output + dble(map_interval)
         end do
      end if
      !
      t = t + timestep      
      !
   enddo
   !
   call ncoutput_finalize()
   !
   call system_clock(clock1)
   call writelog('ls', '')
   call writelog('ls', 'Time loop finished: '//num2str(it)//' time steps, '// &
                       num2str(dble(clock1 - clock0) / dble(max(clock_rate, 1_8)) / dble(max(it, 1)))//' s per step')
   if (n_notconv > 0) then
      call writelog('ls', 'Wave solver did not converge within niter in '//num2str(n_notconv)//' of '//num2str(it)// &
                          ' time steps (see snapwave_warning.txt); consider increasing niter')
   end if
   !
   call close_logfiles()
   !
contains
   !
   subroutine report_time_step(t, it, nt, elapsed)
      !
      ! One line per time step on screen and in the log: model date and time, progress,
      ! mean boundary wave conditions, solver convergence, wall clock time per step and
      ! estimated time to finish
      !
      real*8,  intent(in) :: t, elapsed
      integer, intent(in) :: it, nt
      !
      character(len=200) :: line
      character(len=60)  :: wavestr, solvstr, etastr
      real*8  :: per_step, eta
      real*4  :: hm0, dirn
      integer :: vals(8), fin_s, fin_d
      !
      ! Boundary conditions (Hm0 averaged over boundary points, Tp/dir/zs as used by the solver)
      !
      if (nwbnd > 0 .and. ntwbnd > 0) then
         hm0  = sum(hst_bwv) / real(nwbnd)
         dirn = modulo(270.0 - wdmean_bwv * rad2deg, 360.0)
         write (wavestr, '(a,f6.2,a,f5.1,a,i4,a,f6.2,a)') 'Hm0', hm0, ' m  Tp', tpmean_bwv, ' s  dir', nint(dirn), &
                                                         '  zs', zsmean_bwv, ' m'
      else
         wavestr = 'no wave boundary conditions'
      end if
      !
      write (solvstr, '(a,i4,a,es8.1,a,f6.1,a)') 'iter', solver_iterations, '  max err', solver_error, &
                                                '  ok', solver_percok, '%'
      !
      ! Wall clock: average time per step so far, estimated time to finish and clock time at finish
      !
      per_step = elapsed / dble(it)
      eta = per_step * dble(max(nt - it, 0))
      call date_and_time(values=vals)
      fin_s = vals(5) * 3600 + vals(6) * 60 + vals(7) + nint(eta)
      fin_d = fin_s / 86400
      fin_s = mod(fin_s, 86400)
      if (fin_d > 0) then
         write (etastr, '(a,a,i2.2,a,i2.2,a,i0,a)') trim(duration_str(eta)), ' (', fin_s / 3600, ':', mod(fin_s, 3600) / 60, &
                                                    ' +', fin_d, 'd)'
      else
         write (etastr, '(a,a,i2.2,a,i2.2,a)') trim(duration_str(eta)), ' (', fin_s / 3600, ':', mod(fin_s, 3600) / 60, ')'
      end if
      !
      write (line, '(a,a,i6,a,i6,f6.1,a,a,a,a,f7.2,a,a)') '   ', date_after_seconds(trefstr, t)//'  step', it, '/', nt, &
             100.0 * real(it) / real(nt), '%  ', trim(wavestr), '  ', trim(solvstr), per_step, ' s/step  ETA ', trim(etastr)
      call writelog('ls', line)
      !
   end subroutine report_time_step
   !
   function duration_str(s) result(str)
      !
      ! e.g. 2d 03h, 3h 25m, 12m 05s, 40s
      !
      real*8, intent(in) :: s
      character(len=16) :: str
      integer :: is
      is = nint(s)
      if (is >= 86400) then
         write (str, '(i0,a,i2.2,a)') is / 86400, 'd ', mod(is, 86400) / 3600, 'h'
      elseif (is >= 3600) then
         write (str, '(i0,a,i2.2,a)') is / 3600, 'h ', mod(is, 3600) / 60, 'm'
      elseif (is >= 60) then
         write (str, '(i0,a,i2.2,a)') is / 60, 'm ', mod(is, 60), 's'
      else
         write (str, '(i0,a)') is, 's'
      end if
   end function duration_str
   !
end program
