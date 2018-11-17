program main
  implicit none
  integer,parameter :: nt = 42712
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8) :: tt(0:nt), jt_p(0:nt), jt_m(0:nt)
  real(8) :: jt_2nd(0:nt), jt_ave
  real(8) :: window(0:nt)
  real(8) :: sigma = 2d0*pi/1.17685180813d0
  integer :: it,it_t
  real(8) :: dt, xx, f1


  open(20,file="Act_jt_p.out")
  do it = 0,nt
    read(20,*)tt(it),f1,jt_p(it)
  end do
  close(20)

  open(20,file="Act_jt_m.out")
  do it = 0,nt
    read(20,*)tt(it),f1,jt_m(it)
  end do
  close(20)

  dt = tt(1)-tt(0)
  jt_2nd = 0.5d0*(jt_p + jt_m)

  do it = 0,nt
    xx = dt*it
    window(it) = exp(-0.5d0*(xx/sigma)**2)/sqrt(2d0*pi*sigma**2)*dt
  end do


  open(20,file='jt_2nd.out')
  do it = 0,nt

    jt_ave = 0d0
    do it_t = 0,nt
      jt_ave = jt_ave + window(abs(it - it_t))*jt_2nd(it_t)
    end do

    write(20,"(999e26.16e3)")tt(it),jt_2nd(it),jt_ave


  end do
  close(20)




end program main
