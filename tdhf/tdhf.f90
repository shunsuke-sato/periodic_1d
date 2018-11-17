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
  integer :: nstate, nstate_occ
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


  nstate = 5
  nstate_occ = 3
  lattice_a = 5d0

  nx = 256*nstate_occ
  length_x = nstate_occ*lattice_a
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

  allocate(zpsi(0:nx-1,nstate_occ), zrho_dm(0:nx-1,0:nx-1))
  allocate(zv_F(0:nx-1, 0:nx-1))


! set potentials
! external potential
  do ix = 0, nx-1
!    v_ext(ix) = -cos(pi*xx(ix)/lattice_a)**2
    v_ext(ix) = 1.0d0*(cos(pi*xx(ix)/lattice_a)**2 + cos(pi*(xx(ix)/lattice_a-0.25d0))**4)
!    v_ext(ix) = 0d0
  end do

! interacting potential
  do ix = 0, nx-1
    do jx = 0, nx-1
      dist = abs(xx(ix)-xx(jx))
      v_int(ix,jx) = 0.1d0*cos(pi*dist/lattice_a)**2
!      v_int(ix,jx) = 0d0
    end do
  end do


end subroutine initialize
!----------------------------------------------------------------------------------------!
subroutine calc_ground_state
  use global_variables
  implicit none
  integer :: iscf,ix
  integer :: istate, jstate
  real(8) :: ss, diff_rho_dm
  real(8),parameter :: epsilon_rho_dm = 1d-6
  real(8), allocatable :: v_H_t(:), v_F_t(:,:), rho_dm_old(:,:)
  real(8),parameter :: mixing = 0.5d0

  allocate(v_H_t(0:nx-1), v_F_t(0:nx-1,0:nx-1))
  allocate(rho_dm_old(0:nx-1,0:nx-1))
  v_H_t = 0d0
  v_F_t = 0d0


! initialize wavefunction
  do istate = 1, nstate
!    call random_number(psi(:,istate))
    do ix = 0, nx-1
      psi(ix,istate) = cos(nstate*pi*xx(ix)/length_x)**2
    end do
  end do

! Gram-Schmidt
  do istate = 1, nstate
    do jstate = 1, istate-1
      ss = sum(psi(:,istate)*psi(:,jstate))*dx
      psi(:,istate) = psi(:,istate) - ss * psi(:,jstate)
    end do
    ss = sum(psi(:,istate)**2)*dx
    psi(:,istate) = psi(:,istate)/sqrt(ss)
  end do


  call update_hamiltonian('gs')
  rho_dm_old = rho_dm
  v_H = 0d0; v_F = 0d0

  do iscf = 1, 600
    call cg_eigen(20)
    v_H_t = v_H
    v_F_t = v_F
    call update_hamiltonian('gs')
    v_H = mixing*v_H + (1d0-mixing)*v_H_t
    v_F = mixing*v_F + (1d0-mixing)*v_F_t

    diff_rho_dm = sqrt(sum((rho_dm-rho_dm_old)**2)*dx**2)
!    write(*,*)"diff_rho_dm",iscf,diff_rho_dm
    if(diff_rho_dm < epsilon_rho_dm)exit
    rho_dm_old = rho_dm

  end do

  write(*,"(A,2x,I7,2x,e16.6e3)")"diff. rho_dm",iscf, diff_rho_dm

  open(20,file='gs_wfn_rho_v.out')
  write(20,"(A)")"# xx, rho, v_ext, v_H, psi"
  do ix = 0,nx-1
    write(20,"(999e26.16e3)")xx(ix),rho(ix),v_ext(ix),v_H(ix),psi(ix,:)
  end do
  close(20)


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
  real(8) :: lambda, xixi, xixi_old, gamma
  integer :: istate, jstate, icg
  real(8) :: ss, aa, bb, theta, res

  do istate = 1 ,nstate

    psi_t(:) = psi(:,istate)
    do jstate = 1, istate-1
      ss = sum(psi_t(:)*psi(:,jstate))*dx
      psi_t(:) = psi_t(:) - ss * psi(:,jstate)
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

    do icg = 0, ncg_in

      if(icg == 0)then
        gamma = 0d0
        phi_old = 0d0
      else
        gamma = xixi/xixi_old
        xixi_old = xixi
      end if

      phi_t = xi + gamma*phi_old
      phi_old = phi_t
      ss = sum(phi_t*psi_t)*dx
      phi_t(:) = phi_t(:) - ss * psi_t(:)
      ss = sum(phi_t**2)*dx
      phi_t = phi_t/sqrt(ss)

      bb = 2d0*sum(phi_t*hpsi_t)*dx
      call hpsi(phi_t,hpsi_t)
      aa = sum(phi_t*hpsi_t)*dx - lambda
      aa = - aa ! fix
      if(aa /= 0d0)then
        theta = 0.5d0*atan(bb/aa)
      else
        theta = 0.25d0*pi
      end if
      psi_t = cos(theta)*psi_t + sin(theta)*phi_t


      do jstate = 1, istate-1
        ss = sum(psi_t(:)*psi(:,jstate))*dx
        psi_t(:) = psi_t(:) - ss * psi(:,jstate)
      end do
      ss = sum(psi_t(:)**2)*dx
      psi_t(:) = psi_t(:)/sqrt(ss)
      if(icg == ncg_in)exit

! calc xi
      call hpsi(psi_t, hpsi_t)
      lambda = sum(psi_t*hpsi_t)*dx
      xi = lambda*psi_t - hpsi_t
      do jstate = 1, istate-1
        ss = sum(xi(:)*psi(:,jstate))*dx
        xi(:) = xi(:) - ss * psi(:,jstate)
      end do

      xixi = sum(xi**2)*dx

    end do
    psi(:,istate) = psi_t


  end do


  write(*,"(A)")"Eigenstates"
  do istate = 1, nstate
    call hpsi(psi(:,istate),hpsi_t)
    lambda = sum(psi(:,istate)*hpsi_t)*dx
    psi_t = hpsi_t - lambda*psi(:,istate)
    res = sum(psi_t**2)*dx
    write(*,"(A,2x,I4,2x,999e26.16e3)")"# orb.",istate,lambda,res

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
    do istate = 2, nstate_occ
      rho = rho + psi(:,istate)**2
    end do
    rho = 2d0*rho
    
! one-body reduced density matrix
    rho_dm = 0d0
    do istate = 1, nstate_occ
      do jx = 0, nx -1
        do ix = 0, nx -1
          rho_dm(ix,jx) = rho_dm(ix,jx) + psi(ix,istate)*psi(jx,istate)
        end do
      end do
    end do
    v_F = -rho_dm*v_int*dx

  case('rt','RT')

! one-body density
    rho = abs(zpsi(:,1))**2
    do istate = 2, nstate_occ
      rho = rho + abs(zpsi(:,istate))**2
    end do
    rho = 2d0*rho
    
! one-body reduced density matrix
    zrho_dm = 0d0
    do istate = 1, nstate_occ
      do jx = 0, nx -1
        do ix = 0, nx -1
          zrho_dm(ix,jx) = zrho_dm(ix,jx) + zpsi(ix,istate)*conjg(zpsi(jx,istate))
        end do
      end do
    end do
    zv_F = -zrho_dm*v_int*dx
  case default
    stop 'Error update_hamiltonian'
  end select

  do ix = 0, nx-1
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
  real(8) :: c0, c1, c2

  c0 = -0.5d0*clap0/dx**2
  c1 = -0.5d0*clap1/dx**2
  c2 = -0.5d0*clap2/dx**2


! kinetic energy
  ix = 0
  hpsi_out(ix) = c0*psi_in(ix) &
    +c1*(psi_in(ix+1)+psi_in(ix-1+nx)) &
    +c2*(psi_in(ix+2)+psi_in(ix-2+nx))

  ix = 1
  hpsi_out(ix) = c0*psi_in(ix) &
    +c1*(psi_in(ix+1)+psi_in(ix-1)) &
    +c2*(psi_in(ix+2)+psi_in(ix-2+nx))
  do ix = 0+2, nx-1-2
    hpsi_out(ix) = c0*psi_in(ix) &
      +c1*(psi_in(ix+1)+psi_in(ix-1)) &
      +c2*(psi_in(ix+2)+psi_in(ix-2))
  end do
  ix = nx -1 -1
  hpsi_out(ix) = c0*psi_in(ix) &
    +c1*(psi_in(ix+1)+psi_in(ix-1)) &
    +c2*(psi_in(ix+2-nx)+psi_in(ix-2))

  ix = nx -1
    hpsi_out(ix) = c0*psi_in(ix) &
      +c1*(psi_in(ix+1-nx)+psi_in(ix-1)) &
      +c2*(psi_in(ix+2-nx)+psi_in(ix-2))


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


