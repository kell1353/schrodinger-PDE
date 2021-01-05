! Program: schrodinger
! By: Austin Keller
!-----------------------------------------------------------------------------
! This program solves time dependent Schrodinger equation within a box for different
! zero and harmonic potentials. 
!
! Constructs the initial Gaussina wave function as a function of x and then constructs
! the evolution matrix.
! Using these it contructs the time evolution matrix of the wave function
! which is ultimately written to a .dat file for use in a Jupyter Notebook. 
!
! In addition, this program calulate expectation values and writes them to 
! a .dat file for use in a Jupyter Notebook as well.
!-----------------------------------------------------------------------------
program schrodinger 

use types
use read_write, only : read_input, write_time_evolution, write_expectation_values
use quantum, only : sample_box, construct_initial_wavefunction, construct_time_evolution_matrix, &
    evolve_wave_function, expectation_values

implicit none

real(dp) :: length, delta_t, width, center, k_oscillator
integer :: n_points, n_steps
character(len=1024) :: time_file, density_file, potential_scenario
real(dp), allocatable :: x_vector(:) !will be of size n_points.
real(dp), allocatable :: wave_function(:) !will be of size 2*n_points.
real(dp), allocatable :: evolution_matrix(:,:) !will be of size 2*n_points by 2*n_points
real(dp), allocatable :: time_wave_function(:,:) !will be of size n_points by n_steps + 1 (the +1 is so that you can store the t=0 value)
real(dp), allocatable :: norm(:), position(:), sigma(:) !all of size n_steps + 1 

call read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator, potential_scenario, time_file, density_file)

! allocate arrays
allocate(x_vector(1:n_points))
allocate(wave_function(1:2*n_points))

allocate(norm(1:n_steps+1))
allocate(position(1:n_steps+1))
allocate(sigma(1:n_steps+1))

! allocate matricies
allocate(evolution_matrix(1:2*n_points, 1:2*n_points))
allocate(time_wave_function(1:n_points, 1:n_steps+1))

call sample_box(length, n_points, x_vector)
call construct_initial_wavefunction(x_vector, width, center, wave_function)
call construct_time_evolution_matrix(x_vector, k_oscillator, potential_scenario, delta_t, n_points, length, evolution_matrix)
call evolve_wave_function(wave_function, evolution_matrix, time_wave_function)

call expectation_values(x_vector, time_wave_function, norm, position, sigma)
call write_time_evolution(density_file, x_vector, time_wave_function)
call write_expectation_values(time_file, potential_scenario, width, norm, position, sigma, k_oscillator, center)

end program schrodinger