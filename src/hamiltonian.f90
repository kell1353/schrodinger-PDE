!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! The purpose of this module to evaluate the kinectic and potential 
!! energies and construct the hamiltonian. Which will then be used to
!! solve for the eignvalues and eigenfunctions of the wave equation.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! construct_hamiltonian
!! harmonic_potential_energy
!! WS_potential_energy
!! tri_diagonal_kinetic_energy
!!----------------------------------------------------------------------

module hamiltonian
use types
implicit none


real(dp), parameter :: h_bar = 1           !MeV/fm
real(dp), parameter :: a = .2              !fm (fermi)
integer, parameter :: mass = 1             !MeV/c^2 (let c=1)
integer, parameter :: potential_depth = 50 !MeV


private

public :: construct_hamiltonian, harmonic_potential_energy, h_bar, mass
contains


!-----------------------------------------------------------------------
!! Subroutine: construct_hamiltonian
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the delta, and specific potential energies and 
!! then calls a subroutine to calculate the kinetic energies.
!!
!! It then constructs the hamiltonian matrix using  
!! the kinetic and potential energy as the on and off diagonal of the matrix.
!!----------------------------------------------------------------------
!! Input:
!!
!! delta						real	Real delta constant for the kinetic energy
!! potential_diagonal 			real	Array containing the specific potential energy cvalues
!-----------------------------------------------------------------------
!! Output:
!!
!! hamiltonian 		            real	Matrix containing the hamiltonian constants
!-----------------------------------------------------------------------
subroutine construct_hamiltonian(delta, potential_diagonal, hamiltonian)
implicit none
real(dp), intent(in) :: delta, potential_diagonal(:)
real(dp), intent(out) :: hamiltonian(:,:)

real(dp), allocatable :: kinetic_diagonal(:), kinetic_off_diagonal(:), hamiltonian_diagonal(:), hamiltonian_off_diagonal(:)

integer :: i, j, n

n = size(potential_diagonal)

!Allocate arrays
allocate(kinetic_diagonal(1:n))
allocate(kinetic_off_diagonal(1:n-1))
allocate(hamiltonian_diagonal(1:n))
allocate(hamiltonian_off_diagonal(1:n-1))

call tri_diagonal_kinetic_energy(delta, kinetic_diagonal, kinetic_off_diagonal)

! kinetic diagonal and off-diagonal terms are constant
hamiltonian_diagonal = kinetic_diagonal + potential_diagonal
hamiltonian_off_diagonal = kinetic_off_diagonal

! Create Hamiltonian Matrix
hamiltonian = 0._dp
do i = 1,n 
    do j = 1,n
        ! Diagonal terms
        if (i == j) then 
            hamiltonian(i,j) = hamiltonian_diagonal(i)
        ! Off-diagonal terms
        else if (i == j + 1 .or. i == j - 1) then
            hamiltonian(i,j) = hamiltonian_off_diagonal(1)
        end if
    end do
end do

end subroutine construct_hamiltonian



!-----------------------------------------------------------------------
!! Subroutine: harmonic_potential_energy
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the samples points in the box and calculates
!! the harmonic well potential energy into the diagonal array.
!!----------------------------------------------------------------------
!! Input:
!!
!! x 							real	Array containing the sample points
!! k                            real    Oscillation constant
!-----------------------------------------------------------------------
!! Output:
!!
!! harmonic_potential_diagonal 	real	Array containing the harmonic potential 
!! 										energy values along the diagonal 
!-----------------------------------------------------------------------
subroutine harmonic_potential_energy(x, k, harmonic_potential_diagonal)
implicit none
real(dp), intent(in) :: x(:), k
real(dp), intent(out) :: harmonic_potential_diagonal(:)

integer :: i, n

n = size(x)

harmonic_potential_diagonal = 0._dp
do i=1,n
    harmonic_potential_diagonal(i) = (k/2)*(x(i)**2)
enddo

end subroutine harmonic_potential_energy




!-----------------------------------------------------------------------
!! Subroutine: WS_potential_energy
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the sample points in the box and calculates
!! the Woods-Saxon potential energy into the diagonal array.
!!----------------------------------------------------------------------
!! Input:
!!
!! R 						real	Radius of the Woods_Saxon potential
!! x 						real	Array containing the sample points
!-----------------------------------------------------------------------
!! Output:
!!
!! WS_potential_diagonal 	real	Array containing the Woods-Saxon potential 
!! 									energy values along the diagonal 
!-----------------------------------------------------------------------
subroutine WS_potential_energy(R, x_vector, WS_potential_diagonal)
implicit none
real(dp), intent(in) :: R, x_vector(:)
real(dp), intent(out) :: WS_potential_diagonal(:)

integer :: i, n

n = size(x_vector)

WS_potential_diagonal = 0._dp
do i=1,n
    WS_potential_diagonal(i) = (-potential_depth)/(1+exp(abs(x_vector(i))-R)/a)
enddo

end subroutine WS_potential_energy




!-----------------------------------------------------------------------
!! Subroutine: tri_diagonal_kinetic_energy
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the delta from the sample points and calculates
!! the kinetic energy into diagonal and off diagonal arrays.
!!----------------------------------------------------------------------
!! Input:
!!
!! delta					real	Real delta constant for the kinetic energy
!-----------------------------------------------------------------------
!! Output:
!!
!! kinetic_diagonal 		real	Array containing the kinetic energy diagonal constants
!! kinetic_off_diagonal 	real	Array containing the kinetic energy off diagonal constants
!-----------------------------------------------------------------------
subroutine tri_diagonal_kinetic_energy(delta, kinetic_diagonal, kinetic_off_diagonal)
implicit none
real(dp), intent(in) :: delta
real(dp), intent(out) :: kinetic_diagonal(:), kinetic_off_diagonal(:)

kinetic_diagonal = ((h_bar**2)/mass)/delta**2
kinetic_off_diagonal = ((-0.5*h_bar**2)/mass)/delta**2

end subroutine tri_diagonal_kinetic_energy


    
end module hamiltonian