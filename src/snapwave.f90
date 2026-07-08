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
   use omp_lib
   !
   implicit none
   !
   real*8  :: t
   real*8  :: output_tol
   !
   integer :: it
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
   write(*,'(A,I2,A)') 'Running with ', omp_get_num_threads(), ' OMP threads.' 
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
   write(*,*)'Start time loop'
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
end program
