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

  integer,allocatable :: ixp1(:),ixm1(:),iyp1(:),iym1(:)
  integer,allocatable :: ixp2(:),ixm2(:),iyp2(:),iym2(:)


! simulation box
  integer :: nx,ny          ! number of grid points
  real(8) :: length_x,length_y     ! length of the simulation box
  real(8) :: dx,dy           ! dx = length_x/nx
  real(8),allocatable :: xx(:),yy(:)

! electronic system
  integer :: nstate, nstate_occ
  real(8) :: lattice_ax,lattice_ay   ! period of spatial potentials

! common quantities
  real(8),allocatable :: v_ext(:,:)
  real(8),allocatable :: rho(:,:), v_H(:,:)

! GS quantities
  real(8),allocatable :: psi(:,:,:)

! RT quantities
  complex(8),allocatable :: zpsi(:,:,:)

!
  integer :: nt
  real(8) :: dt, total_time
  real(8),allocatable :: jt(:,:),Act(:,:)

  real(8) :: E0, omega0, Tpulse0

end module global_variables
!----------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none

  call input
  call initialize
  
end program main
!----------------------------------------------------------------------------------!
subroutine input
  use global_variables
  implicit none


  nstate = 5
  nstate_occ = 3
  lattice_ax = 5d0
  lattice_ay = 5d0

  nx = 32
  ny = 32
  length_x = lattice_ax
  length_y = lattice_ay
  
  dx = length_x/nx
  dy = length_y/ny


  E0 = -1d-5
  omega0 = 1.17685180813d0 !0.358d0 !1d-2
  Tpulse0 = 40d0*2d0*pi/omega0

  dt = 0.01d0
  total_time = tpulse0*2d0
  nt = aint(total_time/dt)+1


end subroutine input
!----------------------------------------------------------------------------------!
subroutine initialize
  use global_variables
  implicit none
  integer :: ix,iy

  allocate(xx(0:nx-1),yy(0:ny-1))
  do ix = 0, nx-1
     xx(ix) = dx*ix
  end do
  do iy = 0, ny-1
     yy(ix) = dy*iy
  end do

  allocate(v_ext(0:nx-1,0:ny-1))
  allocate(rho(0:nx-1,0:ny-1))
  allocate(psi(0:nx-1,0:ny-1,nstate))

  allocate(zpsi(0:nx-1,0:ny-1,nstate))

  ! set potential
  do iy = 0,ny-1
     do ix = 0,nx-1
        v_ext(ix,iy) = -0.5d0*(cos(pi*xx(ix)/lattice_ax)**2&
                              *cos(pi*yy(iy)/lattice_ay)**2&
                              +sin(pi*xx(ix)/lattice_ax)**2&
                              *cos(pi*(yy(iy)-0.25d0)/lattice_ay)**2)
     end do
  end do

  
  

end subroutine initialize
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!

