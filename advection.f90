!=======================================================================
!   This program solves the advection equation with the Lax Method
!=======================================================================

!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
  integer, parameter :: nx=200
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=1, dx=xmax/float(nx)
  !   This is the speed of propagation of the advection equation
  !   (a in the notes)
  real, parameter :: speed = 1.5
  !
  !   This is a vector that contains u(x)
  real :: u(0:nx+1)
  !  
end module globals
!=======================================================================
!   main program
program advection
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $\Delta t$
  real, parameter :: tmax= 3.             ! maximumn integration time
  real, parameter :: dtprint=0.1          ! interval between outputs
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

  ! This subroutine generates the initial conditions
  call initflow(time, tprint, itprint)
  
  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

  ! output at tprint intervals
    if(time.ge.tprint) then
      write(*,*) time,tmax,dt
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $\Delta t$ allowed by the CFL criterium
    call timestep(dt)
    !
    ! Integrate u fom t to t+dt
    call tstep(dt,time) 
    ! time counter increases
    time=time+dt

  end do

  stop
end program advection
!=======================================================================
! generates initial condition
subroutine initflow(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here is
  ! f(x,t=0)=1.   for   0.4<x<0.6
  ! f(x,t=0)=0.55 for x=0.4 or x=0.6
  ! f(x,t=0)=0.1  elsewhere
  real, parameter :: u1=0.1, u2=1., x1=0.4, x2=0.6
  real :: x
  integer :: i

  !  fill the vector u
  do i=0,nx+1
    x=float(i)*dx   ! obtain the position $x_i$
    if( (x < x1 ) .or. (x > x2 ) ) then
      u(i)=u1
    else
      u(i)=u2
    end if
    !
    if( (x-0.5*dx <= x1).and.(x+0.5*dx >= x1) ) then
      u(i)=(u1+u2)/2.
    end if
    if( (x-0.5*dx <= x2).and.(x+0.5*dx >= x2) ) then
      u(i)=(u1+u2)/2.
    end if
  end do
  
  print*,u(1)
  
  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0
  
  return
end subroutine initflow
!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=20) file1
  integer :: i

  ! open output file
  write(file1,'(a,i2.2,a)') 'adv_lax-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and u
  do i=1,nx
    write(10,*) float(i)*dx,u(i)
  end do
  
  ! closes output file
  close(10)

  return
end subroutine output
!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  ! Courant number =0.9
  real, parameter :: Co=0.9
  real :: del
  integer :: i

  dt=Co*dx/abs(speed)

  return
end subroutine timestep
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(0:nx+1), f(0:nx+1)
  real :: dtx
  integer :: i

  !  obtain the fluxes
  call fluxes(nx,speed,u,f)

  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  dtx=dt/dx
  do i=1,nx
    up(i)=0.5*(u(i-1)+u(i+1)-dtx*(f(i+1)-f(i-1)))
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,up)

  ! copy the up to the u
  u(:)=up(:)

  return
end subroutine tstep
!=======================================================================
! Obtain the fluxes F
subroutine fluxes(nx,speed,u,f)
  implicit none
  integer, intent(in) :: nx
  real, intent(in) :: speed
  real, intent(in) :: u(0:nx+1)
  real, intent(out):: f(0:nx+1)
  integer :: i

  do i=0,nx+1
    f(i)=speed*u(i)
  end do

  return
end subroutine fluxes
!=======================================================================
! Set boundary conditions
subroutine boundaries(nx,u)
  implicit none
  integer, intent(in) :: nx
  real,    intent(out):: u(0:nx+1)
  integer :: i

  !   Periodic boundary conditions
  u(0 )=u(nx)
  u(nx+1)=u(1)
  
  return
end subroutine boundaries
!=======================================================================
