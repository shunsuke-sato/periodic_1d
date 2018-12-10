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

  real(8) :: E0, omega0, Tpulse0, dir_pol(2)

end module global_variables
!----------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none

  call input
  call initialize

  call calc_ground_state

  call calc_time_propagation
  
end program main
!----------------------------------------------------------------------------------!
subroutine input
  use global_variables
  implicit none


  nstate = 5
  nstate_occ = 1
  nkx = 1
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
  dir_pol(1:2) = (/ 1d0, 0d0 /)

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
     yy(iy) = dy*iy
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
  allocate(v_h(0:nx-1,0:ny-1))
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
!      v_int(ix,iy) = 0d0
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
  integer :: iscf

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
        do jstate = 1, istate-1
           zs = sum(conjg(zpsi(:,:,jstate,ik))*zpsi(:,:,istate,ik))*dx*dy
           zpsi(:,:,istate,ik) = zpsi(:,:,istate,ik) &
                -zs*zpsi(:,:,jstate,ik)
        end do
        ss = sum(abs(zpsi(:,:,istate,ik))**2)*dx*dy
        zpsi(:,:,istate,ik) = zpsi(:,:,istate,ik)/sqrt(ss)
     end do
  end do

  call update_hamiltonian

  do iscf = 1, 600
    call cg_eigen(5)
    call update_hamiltonian

  end do

  
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


    do istate = 1, nstate
      call zhpsi(zpsi(:,:,istate,ik),hpsi_t, kxy)
      eps(istate) = sum(conjg(zpsi(:,:,istate,ik))*hpsi_t)*dx*dy
      psi_t = hpsi_t - eps(istate)*zpsi(:,:,istate,ik)
      res(istate) = sum(psi_t**2)*dx
    end do

    do i = 1,nstate
      do istate = 1,nstate-1
        jstate = istate + 1
        if(eps(istate)>eps(jstate))then
          
          psi_t = zpsi(:,:,istate,ik)
          lambda = eps(istate)
          ss = res(istate)
          
          zpsi(:,:,istate,ik) = zpsi(:,:,jstate,ik)
          eps(istate)= eps(jstate)
          res(istate) = res(jstate)

          zpsi(:,:,jstate,ik) = psi_t 
          eps(jstate) = lambda
          res(jstate) = ss
        end if
      end do
    end do
    


    do istate = 1,nstate
      write(*,"(A,2x,I4,2x,I4,2x,999e26.16e3)")"ik, # orb.",ik,istate,eps(istate),res(istate)
    end do
    write(*,*)
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


  c0 = c0x+c0y+0.5d0*sum(ktmp(:)**2)

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
subroutine calc_time_propagation
  use global_variables
  implicit none
  integer :: it, it_t
  real(8) :: jt_t(2)

  call update_hamiltonian
  call init_laser_field
  call calc_current(jt_t, 0)
  jt(:,0) = jt_t(:)

  do it = 0, nt
    write(*,*)'it=',it,nt
    call dt_evolve(it)
    call calc_current(jt_t, it+1)
    jt(:,it+1) = jt_t(:)

    if(mod(it,max(1,nt/100)) == 0 .or. it == nt)then
      open(30,file='Act_jt.out')
      do it_t = 0, nt+1
        write(30,"(999e26.16e3)")dt*it_t,Act(:,it_t),jt(:,it_t)
      end do
      close(30)
    end if

  end do

end subroutine calc_time_propagation
!----------------------------------------------------------------------------------!
subroutine init_laser_field
  use global_variables
  implicit none
  integer :: it
  real(8) :: a0, tt, x
  
  allocate(Act(2,-1:nt+1), jt(2,-1:nt+1))
  Act = 0d0
  jt = 0d0

  a0 = -E0/omega0

! impulse
!  Act(0:nt+1) = a0
!  return

! laser
  do it = 0, nt
    tt = dt*it
    x = tt-0.5d0*tpulse0
    if(abs(x)<0.5d0*tpulse0)then
      Act(:,it) = dir_pol(:)*a0*cos(pi*x/tpulse0)**2*sin(omega0*x)
    end if
  end do

end subroutine init_laser_field
!----------------------------------------------------------------------------------!
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8)  :: zpsi_t(0:nx-1,0:ny-1,nstate_occ,nk)
  real(8) :: Act_t(2)

  zpsi_t(:,:,1:nstate_occ,:) = zpsi(:,:,1:nstate_occ,:)

! predictor step
  Act_t(:) = Act(:,it)

  kx(:) = kx0(:) + Act_t(1)
  ky(:) = ky0(:) + Act_t(2)

  call propagator(dt*0.5d0)
  call update_hamiltonian

! corrector step
  Act_t(:) = 0.5d0*(Act(:,it) + Act(:,it+1))
  zpsi(:,:,1:nstate_occ,:) = zpsi_t(:,:,1:nstate_occ,:)

  kx(:) = kx0(:) + Act_t(1)
  ky(:) = ky0(:) + Act_t(2)

  call propagator(dt)
  call update_hamiltonian


end subroutine dt_evolve
!----------------------------------------------------------------------------------!
subroutine propagator(dt_t)
  use global_variables
  implicit none
  real(8),intent(in) :: dt_t
  integer,parameter :: nexp_Taylor = 4
  integer :: iexp
  complex(8) :: zpsi_t(0:nx-1,0:ny-1)
  complex(8) :: zhpsi_t(0:nx-1,0:ny-1)
  complex(8) :: zfact
  real(8) :: kxy(2)
  integer :: ik,istate


  do ik = 1, nk
    do istate = 1, nstate_occ

      kxy(1) = kx(ik); kxy(2) = ky(ik)

      zpsi_t(0:nx-1,0:ny-1) = zpsi(0:nx-1,0:ny-1,istate,ik)
      zfact = 1d0
      do iexp = 1, nexp_Taylor
        zfact = zfact*(-zi*dt_t)/iexp
        call zhpsi(zpsi_t,zhpsi_t,kxy)
        zpsi(0:nx-1,0:ny-1,istate,ik) = zpsi(0:nx-1,0:ny-1,istate,ik) &
          + zfact*zhpsi_t(0:nx-1,0:ny-1)

        if(iexp == nexp_Taylor)exit
        zpsi_t = zhpsi_t
      end do


    end do
  end do

end subroutine propagator
!----------------------------------------------------------------------------------!
subroutine calc_current(jt_t, it_t)
  use global_variables
  implicit none
  integer,intent(in)  :: it_t
  real(8),intent(out) :: jt_t(2)
  real(8) :: jt_k_t(2)
  real(8) :: kxy(2)
  integer :: ik, istate

  kx(:) = kx0(:) + Act(1,it_t) 
  ky(:) = ky0(:) + Act(2,it_t) 

  jt_t = 0d0
  do ik = 1, nk
    kxy(1) = kx(ik); kxy(2) = ky(ik)
    do istate = 1, nstate_occ
      
      call calc_current_k(zpsi(:,:,istate,ik),kxy, jt_k_t)
      jt_t = jt_t + occ(istate,ik)*jt_k_t

    end do
  end do

  
end subroutine calc_current
!----------------------------------------------------------------------------------!
subroutine calc_current_k(zpsi_t,kxy,jt_t)
  use global_variables
  implicit none
  complex(8),intent(in) :: zpsi_t(0:nx-1,0:ny-1)
  real(8),intent(in) :: kxy(2)
  real(8),intent(out) :: jt_t(2)
  real(8) :: g1x,g2x,g1y,g2y
  integer :: ix, iy
  complex(8) :: zs1, zs2


  g1x = cgra1/dx
  g2x = cgra2/dx

  g1y = cgra1/dy
  g2y = cgra2/dy
  
  jt_t = 0d0

  do iy = 0,ny-1
    do ix = 0,nx-1

      zs1 = -zi*(g1x*(zpsi_t(ixp1(ix),iy)-zpsi_t(ixm1(ix),iy)) &
                +g2x*(zpsi_t(ixp2(ix),iy)-zpsi_t(ixm2(ix),iy)))&
                +kxy(1)*zpsi_t(ix,iy)

      zs2 = -zi*(g1y*(zpsi_t(ix,iyp1(iy))-zpsi_t(ix,iym1(iy))) &
                +g2y*(zpsi_t(ix,iyp1(iy))-zpsi_t(ix,iym2(iy))))&
                +kxy(2)*zpsi_t(ix,iy)

      jt_t(1) = jt_t(1) + real(conjg(zpsi_t(ix,iy))*zs1)
      jt_t(2) = jt_t(2) + real(conjg(zpsi_t(ix,iy))*zs2)

    end do
  end do

  jt_t = jt_t *dx*dy

end subroutine calc_current_k
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!

