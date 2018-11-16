module global_variables

! mathematical parameters
  complex(8),parameter :: zI =(0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)


! simulation box
  integer :: nx           ! number of grid points
  real(8) :: length_x     ! length of the simulation box
  real(8) :: dx           ! dx = length_x/nx
  real(8),allocatable :: xx(:)

! electronic system
  integer :: nsate
  real(8) :: lattice_a   ! period of spatial potentials

! common quantities
  real(8),allocatable :: v_ext(:), v_int(:,:)
  real(8),allocatable :: rho(:), v_H(:)

! GS quantities
  real(8),allocatable :: psi(:,:), rho_dm(:,:)

! RT quantities
  complex(8),allocatable :: zpsi(:,:), zrho_dm(:,:)


end module global_variables
!----------------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none

  call input
  call initialize

  call calc_ground_state

end program main
!----------------------------------------------------------------------------------------!
! set initial parameters
subroutine input
  use global_variables
  implicit none


  nstate = 3
  lattice_a = 5d0

  nx = nstate*32
  length_x = nstate*lattice_a
  dx = length_x/nx

end subroutine input
!----------------------------------------------------------------------------------------!
subroutine initialize
  use global_variables
  implicit none
  integer :: ix, jx
  real(8) :: dist

  allocate(xx(0:nx-1))
  do ix = 0, nx-1
    xx(ix) = dx*ix
  end do


  allocate(v_ext(0:nx-1), v_int(0:nx-1,0:nx-1))
  allocate(rho(0:nx-1), v_H(0:nx-1))
  allocate(psi(0:nx-1,nstate), rho_dm(0:nx-1,0:nx-1))

  allocate(zpsi(0:nx-1,nstate), zrho_dm(0:nx-1,0:nx-1))


! set potentials
! external potential
  do ix = 0, nx-1
    v_ext(ix) = cos(2d0*pi*xx(ix)/lattice_a)
  end do

! interacting potential
  do ix = 0, nx-1
    do jx = 0, nx-1
      dist = abs(xx(ix)-xx(jx))
      v_int(ix,jx) = cos(pi*dist/lattice_a)**2
    end do
  end do



end subroutine initialize
!----------------------------------------------------------------------------------------!
subroutine calc_ground_state
  use global_variables
  implicit none
end subroutine calc_ground_state
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!


