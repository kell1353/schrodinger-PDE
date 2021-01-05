!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This module contains subroutines intended for reading the provided 
!! namelist, as well as establishing default values, in order to run
!! the program with intended parameters. 
!!
!! In addition, this module also contains two subroutines for writing the 
!! results for given potential scenarios to a file, i.e. the results of 
!! the probability density at each time step and the expectation 
!! results as a function of time.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_time_evolution
!! write_expectation_values
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types
use quantum, only : analytic_width, analytic_position

implicit none

private
public :: read_input, write_time_evolution, write_expectation_values

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in several inputs, sets up the structure of the 
!! namelist file to be read. Then sets default values for the initial 
!! conditions in the case no file is enetered as an argument.
!!
!! In the other case, it then checks if the number of arguments is greater 
!! then zero, then calls for the argument number and checks if the namelist
!! file exists. 
!! If so, it reads each sections values in the namelist file and stops
!! the program if there is an error.
!! 
!! Then outputs the provided namelist values for use in the 
!! computation of the the wave function. 
!!----------------------------------------------------------------------
!! Input:
!!
!!----------------------------------------------------------------------
!! Output:
!! 
!! length                   real        the size of the box
!! n_points                 integer     the number of sample points in x
!! n_steps                  integer     the number of time steps
!! delta_t                  real        the size of the time step
!! width                    real        the width of the Gaussian wave function
!! center                   real        the center of the Gaussian wave function
!! k_oscillator             real        the oscillator parameter
!! potential_scenario       character   the scenario of the potential
!! time_file                character   file name for the expectation results as a function of time
!! density_file             character   file name for the results of the probability density
!!----------------------------------------------------------------------
subroutine read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator &
    , potential_scenario, time_file, density_file)
    implicit none
    real(dp), intent(out) :: length, delta_t, width, center, k_oscillator
    integer, intent(out) :: n_points, n_steps
    
    character(len=*) :: time_file, density_file 

    integer :: n_arguments, unit, ierror
    character(len = 1024) :: namelist_file, potential_scenario
    logical :: file_exists

    namelist /integration/ length, n_points, n_steps, delta_t
    namelist /wave_function/ width, center, potential_scenario
    namelist /oscillator/ k_oscillator
    namelist /output/ time_file, density_file

    ! Default values
    length = 5._dp
    n_points = 100 !100
    n_steps = 100 !100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 0._dp
    potential_scenario = 'zero' !'zero' 'harmonic_well'
    k_oscillator = 1.0_dp
    time_file = 'zero_time_results.dat'
    density_file = 'zero_density_results.dat'

    n_arguments = command_argument_count()

    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit = unit, file = namelist_file)
            read(unit, nml = integration, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading integration namelist'
                stop
            endif
            read(unit, nml = wave_function, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading wave_function namelist'
                stop
            endif
            read(unit, nml = oscillator, iostat = ierror)
            if (ierror /= 0) then
                print*, 'Error reading oscillator in namelist.'
                stop
            end if
            read(unit, nml = output, iostat = ierror)
            if (ierror /= 0) then
                print*, 'Error reading output in namelist.'
                stop
            end if
            close(unit)
        else
            print*, 'Argument, ', trim(namelist_file)
            print*, 'does not exist. Ending program'
            stop
        endif
    elseif(n_arguments /= 0) then
        print*, 'Incorrect number of arguments'
        print*, 'The program takes either 0 or 1 arguments'
        print*, 'See documentation on README.md for details'
        stop
    endif

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine writes to the density_file file the probability 
!! density at different times. The first column contains the sample 
!! points along the x axis.
!!
!! The successive lines contain the probability density at  
!! each time step.
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name                character   the name of the file wriitten to
!! x_points                 real        array containing the sample points inside the box
!! time_wave_functions      real        matrix containing the time evolution of the wave function
!!----------------------------------------------------------------------
subroutine write_time_evolution(file_name, x_points, time_wave_functions)
    implicit none

    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: x_points(:), time_wave_functions(:,:)
    integer :: unit, i, j, n

    n = size(x_points)

    open(newunit=unit,file=trim(file_name))
    write(unit,'(5a20)') 'x', 'initial_state'!, 'time_step_1.65', 'time_step_3.35', 'time_step_5'

    do i=1,n
        write(unit,'(101F14.7)') x_points(i), (time_wave_functions(i,j), j=1,n)
    enddo

    close(unit)

    print *, 'the time evolutions of the wave function were written in ', file_name
    
end subroutine write_time_evolution

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine writes to the time_file file the expectation 
!! values as a function time. The first COLUMN contains the times
!! at which the wave function was calculated
!!
!! The successive columns contains the expectation values 
!! (normalization, position, width, analytic width/position) at the 
!! respective times.
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name             character   the name of the file wriitten to
!! potential_scenario    character   the scenario of the potential
!! initial_width         real        the width of the Gaussian wave function
!! norm                  real        array containing the norm at each time step
!! position              real        array containing the numerical position at each time step
!! width                 real        array containing the numerical width at each time step
!! k                     real        the oscillator parameter
!! center                real        the center of the Gaussian wave function
!!----------------------------------------------------------------------
subroutine write_expectation_values(file_name, potential_scenario, initial_width, norm, position, width, k, center)
    implicit none
    character(len=*), intent(in) :: file_name, potential_scenario
    real(dp), intent(in) :: norm(:), position(:), width(:), initial_width, k, center
    real(dp) :: t, width_analytic, position_analytic
    integer :: unit, i, n

    n = size(norm)

    open(newunit=unit,file=trim(file_name))

    t = 0._dp
    ! Includes analytic_width if the scenario is zero potential
    if (potential_scenario .eq. 'zero') then
        write(unit,'(5a20)') 'time', 'norm', 'position', 'width', 'analytic_width'
        do i=1,n
            width_analytic = analytic_width(t, initial_width)
            write(unit, *) t, norm(i), position(i), width(i), width_analytic
            t = t + 0.05
        enddo
    ! Includes analytic_position if the scenario has harmonic potential
    else if (potential_scenario .eq. 'harmonic_well') then
        write(unit,'(5a20)') 'time', 'norm', 'position', 'width', 'analytic_position'
        do i=1,n
            position_analytic = analytic_position(t, k, center)
            write(unit, *) t, norm(i), position(i), width(i), position_analytic
            t = t + 0.05
        enddo
    end if

    close(unit)

    print *, 'the expectation values were written in ', file_name

    
end subroutine write_expectation_values

end module read_write
