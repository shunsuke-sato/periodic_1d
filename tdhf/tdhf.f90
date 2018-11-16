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
  real(8),allocatable :: psi(:,:), rho_dm(:,:), v_F(:,:)

! RT quantities
  complex(8),allocatable :: zpsi(:,:), zrho_dm(:,:), zv_F(:,:)


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
  allocate(v_F(0:nx-1, 0:nx-1))

  allocate(zpsi(0:nx-1,nstate), zrho_dm(0:nx-1,0:nx-1))
  allocate(zv_F(0:nx-1, 0:nx-1))


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
  integer :: istate, jstate
  real(8) :: ss

! initialize wavefunction
  do istate = 1, nstate
    call random_number(psi(:,istate))
  end do

! Gram-Schmidt
  do istate = 1, nstate
    do jstate = 1, istate-1
      ss = sum(psi(:,istate)*psi(:,jstate))*dx
      psi(:,istate) = psi(:istate) - ss * psi(:,jstate)
    end do
    ss = sum(psi(:,istate)**2)*dx
    psi(:,istate) = psi(:,istate)/sqrt(ss)
  end do


  call update_hamiltonian('gs')


end subroutine calc_ground_state
!----------------------------------------------------------------------------------------!
subroutine update_hamiltonian(gs_rt)
  use global_variables
  implicit none
  character(2),intent(in) :: gs_rt
  integer :: istate
  integer :: ix, jx

  select case(gs_rt)
  case('gs','GS')

! one-body density
    rho = psi(:,1)**2
    do istate = 2, nstate
      rho = rho + psi(:,istate)**2
    end do
    rho = 2d0*rho
    
! one-body reduced density matrix
    rho_dm = 0d0
    do istate = 1, nstate
      do jx = 0, nx -1
        do ix = 0, nx -1
          rho_dm(ix,jx) = rho_dm(ix,jx) + psi(ix,istate)*psi(jx,istate)
        end do
      end do
    end do
    v_F = rho_dm*v_int

  case('rt','RT')

! one-body density
    rho = abs(zpsi(:,1))**2
    do istate = 2, nstate
      rho = rho + abs(zpsi(:,istate))**2
    end do
    rho = 2d0*rho
    
! one-body reduced density matrix
    zrho_dm = 0d0
    do istate = 1, nstate
      do jx = 0, nx -1
        do ix = 0, nx -1
          zrho_dm(ix,jx) = zrho_dm(ix,jx) + psi(ix,istate)*conjg(psi(jx,istate))
        end do
      end do
    end do
    zv_F = zrho_dm*v_int
  case default
    stop 'Error update_hamiltonian'
  end select

  do ix = 0, nx
    v_H(ix) = sum(v_int(:,ix)*rho(:))*dx
  end do

end subroutine update_hamiltonian
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!


