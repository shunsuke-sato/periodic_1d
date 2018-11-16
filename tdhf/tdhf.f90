module global_variables

! mathematical parameters
  complex(8),parameter :: zI =(0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! finite difference coefficients
  real(8),parameter :: clap0 = -5d0/2d0
  real(8),parameter :: clap1 =  4d0/3d0
  real(8),parameter :: clap2 = -1d0/12d0
  real(8),parameter :: cgra0 =  0d0
  real(8),parameter :: cgra1 =  2d0/3d0
  real(8),parameter :: cgra2 = -1d0/12d0


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
  integer :: iscf
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

  do iscf = 1, 300
    call cg_eigen(5)
  end do


end subroutine calc_ground_state
!----------------------------------------------------------------------------------------!
! CG eigen solver, PRB 68, 165337 (2003)
subroutine cg_eigen(ncg_in)
  use global_variables
  implicit none
  integer,intent(in) :: ncg_in
  real(8) :: psi_t(0:nx-1)
  real(8) :: hpsi_t(0:nx-1)
  real(8) :: xi(0:nx-1)
  real(8) :: phi_t(0:nx-1),phi_old(0:nx-1)
  real(8) :: lambda, xixi, xixi_old
  integer :: istate
  real(8) :: ss

  do istate = 1 ,nstate

    psi_t(:) = psi(:,istate)
    do jstate = 1, istate-1
      ss = sum(psi_t(:)*psi(:,jstate))*dx
      psi_t(:) = psi_t(: - ss * psi(:,jstate)
    end do
    ss = sum(psi_t(:)**2)*dx
    psi_t(:) = psi_t(:)/sqrt(ss)

! calc xi
    call hpsi(psi_t, hpsi_t)
    lambda = sum(psi_t*hpsi_t)*dx
    xi = lambda*psi_t - hpsi_t
    do jstate = 1, istate-1
      ss = sum(xi(:)*psi(:,jstate))*dx
      xi(:) = xi(:) - ss * psi(:,jstate)
    end do

    xixi = sum(xi**2)*dx
    xixi_old = xixi

    do icg = 0, ncg

      if(icg == 0)then
        gamma = 0d0
        phi_old = 0d0
      else
        gamma = xixi/xixi_old
        xixi_old = xixi
      end if

      phi_t = xi + gamma*phi_old
      do jstate = 1, istate-1
        ss = sum(phi_t(:)*psi(:,jstate))*dx
        phi_t(:) = phi_t(:) - ss * psi(:,jstate)
      end do
      ss = sum(phi_t**2)*dx
      phi_t = phi_t/sqrt(ss)

! implementing theta ....

    end do


  end do


end subroutine cg_eigen
  
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
    v_F = -rho_dm*v_int

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
    zv_F = -zrho_dm*v_int
  case default
    stop 'Error update_hamiltonian'
  end select

  do ix = 0, nx
    v_H(ix) = sum(v_int(:,ix)*rho(:))*dx
  end do

end subroutine update_hamiltonian
!----------------------------------------------------------------------------------------!
subroutine hpsi(psi_in, hpsi_out)
  use global_variables
  implicit none
  real(8),intent(in) :: psi_in(0:nx-1)
  real(8),intent(out) :: hpsi_out(0:nx-1)
  integer :: ix


! kinetic energy
  ix = 0
  hpsi_out(ix) = -0.5d0*(clap0*psi_in(ix) &
    +clap1*(psi_in(ix+1)+psi_in(ix-1+nx)) &
    +clap2*(psi_in(ix+2)+psi_in(ix-2+nx)))

  ix = 1
  hpsi_out(ix) = -0.5d0*(clap0*psi_in(ix) &
    +clap1*(psi_in(ix+1)+psi_in(ix-1)) &
    +clap2*(psi_in(ix+2)+psi_in(ix-2+nx)))
  do ix = 0+2, nx-1-2
    hpsi_out(ix) = -0.5d0*(clap0*psi_in(ix) &
      +clap1*(psi_in(ix+1)+psi_in(ix-1)) &
      +clap2*(psi_in(ix+2)+psi_in(ix-2)))
  end do
  ix = nx -1 -1
  hpsi_out(ix) = -0.5d0*(clap0*psi_in(ix) &
    +clap1*(psi_in(ix+1)+psi_in(ix-1)) &
    +clap2*(psi_in(ix+2-nx)+psi_in(ix-2)))

  ix = nx -1
    hpsi_out(ix) = -0.5d0*(clap0*psi_in(ix) &
      +clap1*(psi_in(ix+1-nx)+psi_in(ix-1)) &
      +clap2*(psi_in(ix+2-nx)+psi_in(ix-2)))


! local potential
  hpsi_out = hpsi_out + (v_ext + v_H)*psi_in

! nonlocal potential
! this part has to be optimized by BLAS routine
  hpsi_out = hpsi_out + matmul(v_F,psi_in)


end subroutine hpsi
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!


