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
  integer :: nstate, nstate_occ, nk, nkx, nky
  real(8) :: lattice_ax,lattice_ay   ! period of spatial potentials
  real(8),allocatable :: kx(:),ky(:),kx0(:),ky0(:)

! common quantities
  real(8),allocatable :: v_ext(:,:),v_int(:,:)
  real(8),allocatable :: rho(:,:), v_H(:,:)


! GS quantities
  real(8),allocatable :: psi(:,:,:,:)

! RT quantities
  real(8),allocatable :: occ(:,:)
  complex(8),allocatable :: zpsi(:,:,:,:)

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

  call calc_ground_state
  
end program main
!----------------------------------------------------------------------------------!
subroutine input
  use global_variables
  implicit none


  nstate = 5
  nstate_occ = 3
  nkx = 3
  nky = 3
  nk = nkx*nky
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
  integer :: ix,iy,ixt,iyt
  integer :: ikx, iky, ik

  allocate(xx(0:nx-1),yy(0:ny-1))
  do ix = 0, nx-1
     xx(ix) = dx*ix
  end do
  do iy = 0, ny-1
     yy(ix) = dy*iy
  end do


  allocate(kx(nk),ky(nk),kx0(nk),ky0(nk))
  ik = 0
  do ikx =1,nkx
     do iky =1,nky
       ik = ik + 1
        if(ikx-1 <= nkx/2)then
           kx0(ik) = 2d0*pi*dble(ikx-1)/dble(nkx)
        else
           kx0(ik) = 2d0*pi*dble(ikx-1-nkx)/dble(nkx)
        end if

        if(iky-1 <= nky/2)then
           ky0(ik) = 2d0*pi*dble(iky-1)/dble(nky)
        else
           ky0(ik) = 2d0*pi*dble(iky-1-nky)/dble(nky)
        end if        
        
     end do
  end do

  kx = kx0
  ky = ky0
  

  allocate(ixp1(0:nx-1),ixm1(0:nx-1),ixp2(0:nx-1),ixm2(0:nx-1))
  allocate(iyp1(0:ny-1),iym1(0:ny-1),iyp2(0:ny-1),iym2(0:ny-1))

  do ix = 0,nx-1
! xp1     
     ixt = mod(ix + 1 + nx,nx)
     ixp1(ix) = ixt
! xp2     
     ixt = mod(ix + 2 + nx,nx)
     ixp2(ix) = ixt
! xm1     
     ixt = mod(ix - 1 + nx,nx)
     ixm1(ix) = ixt
! xm2
     ixt = mod(ix - 2 + nx,nx)
     ixm2(ix) = ixt
  end do

  do iy = 0,ny-1
! yp1     
     iyt = mod(iy + 1 + ny,ny)
     iyp1(iy) = iyt
! yp2     
     iyt = mod(iy + 2 + ny,ny)
     iyp2(iy) = iyt
! ym1     
     iyt = mod(iy - 1 + ny,ny)
     iym1(iy) = iyt
! ym2
     iyt = mod(iy - 2 + ny,ny)
     iym2(iy) = iyt
  end do  
  

  allocate(v_ext(0:nx-1,0:ny-1))
  allocate(v_int(0:nx-1,0:ny-1))
  allocate(rho(0:nx-1,0:ny-1))
  allocate(psi(0:nx-1,0:ny-1,nstate,nk))

  allocate(zpsi(0:nx-1,0:ny-1,nstate,nk))
  allocate(occ(nstate,nk))
  occ = 0d0
  occ(1:nstate_occ,:) = 2d0/dble(nk)
  

  ! set potential
  do iy = 0,ny-1
     do ix = 0,nx-1
        v_ext(ix,iy) = -0.5d0*(cos(pi*xx(ix)/lattice_ax)**2&
                              *cos(pi*yy(iy)/lattice_ay)**2&
                              +sin(pi*xx(ix)/lattice_ax)**2&
                              *cos(pi*(yy(iy)-0.25d0)/lattice_ay)**2)
     end do
  end do

  ! set potential
  do iy = 0,ny-1
     do ix = 0,nx-1
        v_int(ix,iy) = 0.5d0*(cos(pi*xx(ix)/lattice_ax)**2&
                              *cos(pi*yy(iy)/lattice_ay)**2)

     end do
  end do

  
  

end subroutine initialize
!----------------------------------------------------------------------------------!
subroutine calc_ground_state
  use global_variables
  implicit none
  integer :: ix,iy,istate,ik
  integer :: jstate
  real(8) :: ss
  complex(8) :: zs
  real(8) :: rand_tmp(0:nx-1,0:ny-1)
  kx = kx0
  ky = ky0

  do ik = 1, nk
     do istate = 1, nstate
        call random_number(rand_tmp)
        zpsi(:,:,istate,ik) = rand_tmp
     end do
  end do


! Gram-Schmidt  
  do ik = 1, nk
     do istate = 1, nstate
        do jstate = istate + 1, nstate
           zs = sum(conjg(zpsi(:,:,jstate,ik))*zpsi(:,:,istate,ik))*dx*dy
           zpsi(:,:,istate,ik) = zpsi(:,:,istate,ik) &
                -zs*zpsi(:,:,jstate,ik)
        end do
        ss = sum(abs(zpsi(:,:,istate,ik))**2)*dx*dy
        zpsi(:,:,istate,ik) = zpsi(:,:,istate,ik)/sqrt(ss)
     end do
  end do

  call update_hamiltonian

  
end subroutine calc_ground_state
!----------------------------------------------------------------------------------!
subroutine cg_eigen(ncg_in)
  use global_variables
  implicit none
  integer,intent(in) :: ncg_in

  complex(8) :: psi_t(0:nx-1,0:ny-1)
  complex(8) :: hpsi_t(0:nx-1,0:ny-1)
  complex(8) :: xi(0:nx-1,0:ny-1)
  complex(8) :: phi_t(0:nx-1,0:ny-1),phi_old(0:nx-1,0:ny-1)
  real(8) :: lambda, xixi, xixi_old, gamma
  integer :: istate, jstate, icg, i
  real(8) :: ss, aa, bb, theta 
  complex(8) :: zs
  real(8) :: eps(nstate),res(nstate)
  real(8) :: kxy(2)
  integer :: ik

  do ik = 1, nk
    kxy(1) = kx0(ik); kxy(2) = ky0(ik)

    do istate = 1, nstate

      psi_t(:,:) = zpsi(:,:,istate,ik)
      do jstate = 1, istate -1
        zs = sum(conjg(zpsi(:,:,jstate,ik))*psi_t(:,:))*dx*dy
        psi_t(:,:) = psi_t(:,:) - zs*zpsi(:,:,jstate,ik)
      end do
      ss = sum(abs(psi_t)**2)*dx*dy
      psi_t = psi_t/sqrt(ss)

! calc xi
      call zhpsi(psi_t, hpsi_t, kxy)
      lambda = sum(conjg(psi_t)*hpsi_t)*dx*dy
      xi = lambda*psi_t - hpsi_t
      do jstate = 1, istate -1
        zs = sum(conjg(zpsi(:,:,jstate,ik))*xi(:,:))*dx*dy
        xi(:,:) = xi(:,:) - zs*zpsi(:,:,jstate,ik)
      end do
      
      xixi = sum(abs(xi)**2)*dx*dy
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
        zs = sum(conjg(psi_t)*phi_t)*dx*dy
        phi_t = phi_t - zs*psi_t
        ss = sum(abs(phi_t)**2)*dx*dy
        phi_t = phi_t/sqrt(ss)

        bb = 2d0*sum(conjg(phi_t)*hpsi_t)*dx*dy
        call zhpsi(phi_t,hpsi_t, kxy)
        aa = sum(conjg(phi_t)*hpsi_t)*dx*dy-lambda
        aa = -aa ! fix: there is a typo in the reference paper.

        if(aa /= 0d0)then
          theta = 0.5d0*atan(bb/aa)
        else
          theta = 0.25d0*pi
        end if
        psi_t = cos(theta)*psi_t + sin(theta)*phi_t

        do jstate = 1, istate -1
          zs = sum(conjg(zpsi(:,:,jstate,ik))*psi_t(:,:))*dx*dy
          psi_t(:,:) = psi_t(:,:) - zs*zpsi(:,:,jstate,ik)
        end do
        ss = sum(abs(psi_t)**2)*dx*dy
        psi_t = psi_t/sqrt(ss)

        if(icg == ncg_in)exit

! calc xi
        call zhpsi(psi_t,hpsi_t, kxy)
        lambda = sum(conjg(psi_t)*hpsi_t)*dx*dy
        xi = lambda*psi_t - hpsi_t        
        do jstate = 1, istate -1
          zs = sum(conjg(zpsi(:,:,jstate,ik))*xi(:,:))*dx*dy
          xi(:,:) = xi(:,:) - zs*zpsi(:,:,jstate,ik)
        end do

        xixi = sum(abs(xi)**2)*dx*dy
      end do
      zpsi(:,:,istate,ik) = psi_t

    end do
  end do

end subroutine cg_eigen
!----------------------------------------------------------------------------------!
subroutine zhpsi(zpsi_in,zhpsi_out,ktmp)
  use global_variables
  implicit none
  complex(8),intent(in) :: zpsi_in(0:nx-1,0:ny-1)
  complex(8),intent(out) :: zhpsi_out(0:nx-1,0:ny-1)
  real(8),intent(in) :: ktmp(2)
  integer :: ix,iy
  real(8) :: c0x,c1x,c2x,c0y,c1y,c2y
  real(8) :: g0x,g1x,g2x,g0y,g1y,g2y
  real(8) :: c0

  c0x = -0.5d0*clap0/dx**2
  c1x = -0.5d0*clap1/dx**2
  c2x = -0.5d0*clap2/dx**2

  c0y = -0.5d0*clap0/dy**2
  c1y = -0.5d0*clap1/dy**2
  c2y = -0.5d0*clap2/dy**2

  g1x = ktmp(1)*cgra1/dx
  g2x = ktmp(1)*cgra2/dx

  g1y = ktmp(2)*cgra1/dy
  g2y = ktmp(2)*cgra2/dy


  c0 = c0x+c0y-0.5d0*sum(ktmp(:)**2)

  do iy = 0,ny-1
    do ix = 0,nx-1

      zhpsi_out(ix,iy) = c0*zpsi_in(ix,iy) &
        +c1x*(zpsi_in(ixp1(ix),iy)+zpsi_in(ixm1(ix),iy)) &
        +c2x*(zpsi_in(ixp2(ix),iy)+zpsi_in(ixm2(ix),iy)) &
        +c1y*(zpsi_in(ix,iyp1(iy))+zpsi_in(ix,iym1(iy))) &
        +c2y*(zpsi_in(ix,iyp2(iy))+zpsi_in(ix,iym2(iy))) &
        -zI*(&
         g1x*(zpsi_in(ixp1(ix),iy)-zpsi_in(ixm1(ix),iy)) &
        +g2x*(zpsi_in(ixp2(ix),iy)-zpsi_in(ixm2(ix),iy)) &
        +g1y*(zpsi_in(ix,iyp1(iy))-zpsi_in(ix,iym1(iy))) &
        +g2y*(zpsi_in(ix,iyp2(iy))-zpsi_in(ix,iym2(iy))))

    end do
  end do

! adding external and hartree potentials
  zhpsi_out = zhpsi_out + (v_ext+v_H)*zpsi_in

! forck operator should be applied here
!!!!!!

end subroutine zhpsi
!----------------------------------------------------------------------------------!
subroutine update_hamiltonian
  use global_variables
  implicit none
  integer :: ik,istate
  integer :: ix,iy, ixt,iyt, ixd,iyd

  rho = 0d0
  do ik = 1,nk
    do istate = 1, nstate_occ
      rho = rho + occ(istate,ik)*abs(zpsi(:,:,istate,ik))**2
    end do
  end do

! calculate Hartree potential
! this part should be optimized, probably, by FFT.
  v_h = 0d0
  do iy = 0, ny-1
    do ix = 0, nx-1

      do iyt=0,ny-1
        do ixt=0,nx-1
          v_h(ix,iy) = v_h(ix,iy) &
            + v_int(abs(ixt-ix),abs(iyt-iy))*rho(ixt,iyt)

        end do
      end do
    end do
  end do

  v_h = v_h*dx*dy
  

end subroutine update_hamiltonian
!----------------------------------------------------------------------------------!

