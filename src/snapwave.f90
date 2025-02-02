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
   !
   implicit none
   !
   real*8  :: t 
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
   it = 0
   t  = tstart
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
      call update_obs_points()
      !
      if (ja_save_each_iter==0) call ncoutput_update(t, it)
      !
      t = t + timestep      
      !
   enddo
   !
   call ncoutput_finalize()
   !
end program
