!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The purpose of this module is to perform all of the calculations 
!! neccesary for evolving the time dependent Schrodinger equation using 
!! the Cranck-Nicolson method for different potentials inside a box. 
!!
!! This module also provides functions for calculating expected values 
!! such as (norm, width, position) numerically and analyticlly. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! sample_box
!! construct_initial_wavefunction
!! construct_time_evolution_matrix
!! evolve_wave_function
!! expectation_values
!!----------------------------------------------------------------------
!! Included functions:
!!
!! analytic_width
!! analytic_position
!-----------------------------------------------------------------------
module quantum

use types
use hamiltonian, only : construct_hamiltonian, harmonic_potential_energy, h_bar, mass
use linear_algebra, only : invert_matrix, matrix_mult, mult_matrix_vector

implicit none
private
public sample_box, construct_initial_wavefunction, evolve_wave_function, &
       construct_time_evolution_matrix, expectation_values, analytic_width, analytic_position

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine specifies sample points x within our sample box of 
!! length L.
!!----------------------------------------------------------------------
!! Input:
!!
!! L          real        length of the sample box
!! N          integer     number of sample points
!!----------------------------------------------------------------------
!! Output:
!!
!! x          real        array containing the sample points inside the box
!!----------------------------------------------------------------------
subroutine sample_box(L, n, x)
    implicit none
    real(dp), intent(in) :: L
    integer, intent(in) ::  n
    real(dp), intent(out) :: x(:)

    real(dp) :: dx
    integer :: i

    dx = (2*L)/(n-1)

    x = 0._dp
    do i=1,n
        x(i) = -L + dx*(i-1)
	enddo

end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: construct_initial_wavefunction
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine calculates the initial gaussian probabilty densities for 
!! each of the sample points inside the box into the first half of the array.
!! While the second half of the array contains zeros because that is the 
!! imaginary portion and the Gaussian is a purely real function.
!!----------------------------------------------------------------------
!! Input:
!!
!! x                real        array containing the sample points in the box
!! sigma            real        the width of the Gaussian wave function
!! x_0              real        the center of the Gaussian wave function
!!----------------------------------------------------------------------
!! Output:
!!
!! wave_function    real        array containing the initial wave function
!!----------------------------------------------------------------------
subroutine construct_initial_wavefunction(x, sigma, x_0, wave_function)
    implicit none
    real(dp), intent(in) :: x(:), sigma, x_0
    real(dp), intent(out) :: wave_function(:)

    wave_function = 0._dp
    wave_function = ((2*pi*sigma**2)**(-1.0/4.0))*exp(-((x-x_0)**2/(4*sigma**2)))

end subroutine construct_initial_wavefunction

!-----------------------------------------------------------------------
!! Subroutine: construct_time_evolution_matrix
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine constructs the 2n by 2n evolution matrix used to evolve the 
!! wave function using the Crank-Nicolson method.
!!
!! The subroutine breaks the evolution matrix into two seperate super
!! "super" (matricies within the matrix) matricies and assigns each quadrant 
!! of each one individually. 
!!
!! The first and fourth quadrants of each matrix are set as the identity matrix.
!! 
!! The second and fourth quadrants contain the hamiltonian matrix multipled by a
!! constant factor delta_t/2*h_bar. The two matricies have negative values
!! in opposite quadrants.
!!
!! The subroutine then inverses the matrix with the negative in the second 
!! quadrant. Then multiplies the two matricies together to create the full
!! evolution matrix.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector             real        array containing the sample points inside the box
!! k                    real        the oscillator parameter
!! delta_t              real        the size of the time step
!! n_points             integer     the number of sample points in x
!! potential_scenario   character   the scenario of the potential
!!----------------------------------------------------------------------
!! Output:
!!
!! evolution_matrix     real        2n by 2n matrix used for evolving the wave function
!!----------------------------------------------------------------------
subroutine construct_time_evolution_matrix(x_vector, k, potential_scenario, delta_t, n_points, length, evolution_matrix)
    implicit none
    real(dp), intent(in) :: delta_t, x_vector(:), k
    integer, intent(in) :: n_points
    character(len=1024), intent(in) :: potential_scenario
    real(dp), intent(out) :: evolution_matrix(:,:)

    real(dp) :: delta, length
    integer :: i, j, n
    real(dp), allocatable :: potential_diagonal(:), hamiltonian(:,:)
    real(dp), allocatable :: pos_super_matrix(:,:), neg_super_matrix(:,:), neg_super_matrix_inv(:,:)

    n = n_points

    delta = (2*length)/(n-1)

    ! Allocate Arrays
    ! Arrays for the hamiltonian
    if(allocated(potential_diagonal)) deallocate(potential_diagonal)
       allocate(potential_diagonal(1:n))
    if(allocated(hamiltonian)) deallocate(hamiltonian)
       allocate(hamiltonian(1:n,1:n))

    ! Allocate the super matricies
    if(allocated(pos_super_matrix)) deallocate(pos_super_matrix)
        allocate(pos_super_matrix(1:2*n,1:2*n))
    if(allocated(neg_super_matrix)) deallocate(neg_super_matrix)
        allocate(neg_super_matrix(1:2*n,1:2*n))
    if(allocated(neg_super_matrix_inv)) deallocate(neg_super_matrix)
        allocate(neg_super_matrix_inv(1:2*n,1:2*n))

    ! Set potential based on weather it is harmonic or zero
    if (potential_scenario .eq. 'harmonic_well') then
        call harmonic_potential_energy(x_vector, k, potential_diagonal)
    else
        potential_diagonal = 0
    end if

    call construct_hamiltonian(delta, potential_diagonal, hamiltonian)

    neg_super_matrix = 0._dp
    pos_super_matrix = 0._dp
    
    ! Quadrant 1 of the super matix is identity matrix
    do i = 1,n
        do j = 1,n
            if (i == j) then 
                pos_super_matrix(i,j) = 1._dp
                neg_super_matrix(i,j) = 1._dp
            end if
        end do 
    enddo

    ! Quadrant 4 of the super matix is identity matrix
    do i = n,2*n
        do j = n,2*n
            if (i == j) then 
                pos_super_matrix(i,j) = 1._dp
                neg_super_matrix(i,j) = 1._dp
            end if
        end do 
    enddo

    ! Quadrant 2 of the super matix 
    do i = 1,n
        do j = n+1,2*n
            pos_super_matrix(i,j) = (delta_t/(2*h_bar))*hamiltonian(i,j-n)
            neg_super_matrix(i,j) = (-delta_t/(2*h_bar))*hamiltonian(i,j-n)
        end do 
    enddo

    ! Quadrant 3 of the super matix 
    do i = n+1,2*n
        do j = 1,n
            pos_super_matrix(i,j) = (-delta_t/(2*h_bar))*hamiltonian(i-n,j)
            neg_super_matrix(i,j) = (delta_t/(2*h_bar))*hamiltonian(i-n,j)
        end do 
    enddo

    call invert_matrix(neg_super_matrix, neg_super_matrix_inv)
    call matrix_mult(pos_super_matrix, neg_super_matrix_inv, evolution_matrix)

end subroutine construct_time_evolution_matrix

!-----------------------------------------------------------------------
!! Subroutine: evolve_wave_function
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine is responsible for evolving the initial wave function
!! using the evolution matrix provided by the Crank-Nicolson method.
!!
!! It sets the initial wave function as the first column of the time evolution
!! matrix the multiplies the array by the evolution matrix to get a new
!! wave function.
!!
!! Then proceeds to use this new wave function to calculate the next evolution
!! of the wave and does so for each time step in evolution.
!!----------------------------------------------------------------------
!! Input:
!!
!! wave_function        real     2n array containing the initial wave function     
!! evolution_matrix     real     2n by 2n matrix used for evolving the wave function
!!----------------------------------------------------------------------
!! Output:
!!
!! time_wave_function    real    matrix containing the time evolution of the wave function
!!----------------------------------------------------------------------
subroutine evolve_wave_function(wave_function, evolution_matrix, time_wave_function)
    implicit none
    real(dp), intent(in) :: wave_function(:), evolution_matrix(:,:)
    real(dp), intent(out) :: time_wave_function(:,:)

    real(dp), allocatable :: input_vec(:), output_vec(:)
    integer :: i, j, n, n_steps, n_shape(1:2)

    n_shape = shape(time_wave_function)
    n = n_shape(1)
    n_steps = n_shape(2)

    ! Settting the initial Gaussian in the time evolution matrix
    do i=1,n
        time_wave_function(i,1) = wave_function(i)**2 + wave_function(i+n)**2
    enddo

    allocate(input_vec(1:2*n))
    allocate(output_vec(1:2*n))

    do i=2,n_steps
        if (i /= 2) then
            output_vec = 0._dp
            call mult_matrix_vector(input_vec, evolution_matrix, output_vec)
        else
            call mult_matrix_vector(wave_function, evolution_matrix, output_vec)
        end if 
    
        do j=1,n
            time_wave_function(j,i) = output_vec(j)**2 + output_vec(j+n)**2
        enddo

        input_vec = output_vec
    enddo
    
end subroutine evolve_wave_function

!-----------------------------------------------------------------------
!! Subroutine: expectation_values
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in the sample points in the box and the probabilty 
!! densities for each evolution of the wave function.
!!
!! Uses these values to numericallycalculate the expectation values (normalization, position
!! width) for each evolution of the wave function as a function of time.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_points              real        array containing the sample points inside the box
!! time_wave_function    real        matrix containing the time evolution of the wave function
!!----------------------------------------------------------------------
!! Output:
!!
!! norm                  real        array containing the norm at each time step
!! position              real        array containing the numerical position at each time step
!! width                 real        array containing the numerical width at each time step
!!----------------------------------------------------------------------
subroutine expectation_values(x_points, time_wave_function, norm, position, width)
    implicit none
    real(dp), intent(in) :: x_points(:), time_wave_function(:,:)
    real(dp), intent(out) :: norm(:), position(:), width(:)

    real(dp), allocatable :: x(:), x_squared(:)
    real(dp) :: sum_j
    integer :: i, j, n, n_shape(1:2)

    n_shape = shape(time_wave_function)

    ! Time steps + 1
    n = n_shape(2)
    
    ! Calculate normalization
    do i=1,n
        norm(i) = sum(time_wave_function(:,i))
    enddo

    ! Calculate numerical position
    do i=1,n
        sum_j = 0._dp
        do j=1,n_shape(1)
            sum_j = sum_j + x_points(j)*time_wave_function(j,i)
        enddo
        position(i) = sum_j/norm(i)
    enddo

    ! Calculate numerical width
    allocate(x_squared(1:n))
    do i=1,n
        sum_j = 0._dp
        do j=1,n_shape(1)
            sum_j = sum_j + x_points(j)**2*time_wave_function(j,i)
        enddo
        x_squared(i) = sum_j/norm(i)
    enddo

    width = sqrt(x_squared - position**2)


end subroutine expectation_values

!-----------------------------------------------------------------------
!! Function: analytic_width
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function calculates analytic time evolution of the width 
!!----------------------------------------------------------------------
!! Input:
!!
!! t        real        time step 
!! width    real        width of the initial Gaussian wave function
!!----------------------------------------------------------------------
!! Output:
!!
!! k        real        analytic time evolution of the width 
!!----------------------------------------------------------------------
real(dp) function analytic_width(t, width) result(k)
    implicit none
    real(dp), intent(in) :: t, width

    k = sqrt(width**2 + (h_bar**4*t**2)/(4*mass**2*width**2))

end function analytic_width

!-----------------------------------------------------------------------
!! Function: analytic_position
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function calculates position of the harmonic oscillator in the 
!! classical case.
!!----------------------------------------------------------------------
!! Input:
!!
!! t                real        time step 
!! k_oscillator     real        the oscillator parameter
!! x_0              real        center of the initial Gaussian wave function
!!----------------------------------------------------------------------
!! Output:
!!
!! k                real        classical position 
!!----------------------------------------------------------------------
real(dp) function analytic_position(t, k_oscillator, x_0) result(k)
    implicit none
    real(dp), intent(in) :: t, k_oscillator, x_0

    k = x_0*cos(sqrt(k_oscillator/mass)*t)

end function analytic_position


  
end module quantum