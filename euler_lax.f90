!=======================================================================
!   This program solves the Euler equations with the Lax Method
!=======================================================================

!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
  integer, parameter :: nx=100
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=1, dx=xmax/float(nx)
  real, parameter :: gamma=1.4
  !   This is a vector that contains u(x)
  real :: u(3,0:nx+1)
  !  
end module globals
!=======================================================================
!   main program
program euler_lax
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real           :: time, dt             !  t, $\Delta t$
  real, parameter :: tmax= .3             ! maximumn integration time
  real, parameter :: dtprint=0.01          ! interval between outputs
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

  ! This subroutine generates the initial conditions
  call initconds(time, tprint, itprint)
  
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
end program euler_lax
!=======================================================================
! generates initial condition
subroutine initconds(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here is the Sod tube test
  integer :: i
  real :: x

  !  fill the vector u
  do i=0,nx+1
    x=float(i)*dx   ! obtain the position $x_i$
    if (x < 0.5 ) then
      u(1,i)=1.0
      u(2,i)=0.0
      u(3,i)=1.0/(gamma-1)
    else
      u(1,i)=.125
      u(2,i)=0.0
      u(3,i)=0.1/(gamma-1)
    end if
    !
    if( (x-0.5*dx <= 0.5).and.(x+0.5*dx >= 0.5) ) then
      u(1,i)=1.125/2.
      u(2,i)=0.0
      u(3,i)=1.1/2./(gamma-1)
    end if
  end do
  
  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0
  
  return
end subroutine initconds
!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=20) file1
  integer :: i
  real :: rho,vx,P
  ! open output file
  write(file1,'(a,i2.2,a)') 'euler_lax-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and rho, u(=vx) and P
  do i=1,nx
    rho=u(1,i)
    vx=u(2,i)/rho
    P=(u(3,i)-0.5*rho*vx**2)*(gamma-1.)
    write(10,*) float(i)*dx,rho, vx, P
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
  real :: rho, vx , P, cs
  integer :: i

  dt=1E30
  do i=0,nx+1
    rho=u(1,i)
    vx=u(2,i)/rho
    P=(u(3,i)-0.5*rho*vx**2)*(gamma-1.)
    cs=sqrt(gamma*P/rho)
    dt=min( dt,Co*dx/(abs(vx)+cs) )
  end do

  return
end subroutine timestep
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(3,0:nx+1), f(3,0:nx+1)
  real :: dtx
  integer :: i

  !  obtain the fluxes
  call fluxes(nx,gamma,u,f)

  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  dtx=dt/dx
  do i=1,nx
    up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,up)

  ! copy the up to the u
  u(:,:)=up(:,:)

  return
end subroutine tstep
!=======================================================================
! Obtain the fluxes F
subroutine fluxes(nx,gamma,u,f)
  implicit none
  integer, intent(in) :: nx
  real, intent(in) :: gamma
  real, intent(in) :: u(3,0:nx+1)
  real, intent(out):: f(3,0:nx+1)
  integer :: i
  real :: rho, vx , P

  do i=0,nx+1
    rho=u(1,i)
    vx= u(2,i)/rho
    P= (u(3,i)-0.5*rho*vx**2)*(gamma-1.)

    f(1,i)=rho*vx
    f(2,i)=rho*vx**2+P
    f(3,i)=vx*(u(3,i)+P)
    
  end do

  return
end subroutine fluxes
!=======================================================================
! Set boundary conditions
subroutine boundaries(nx,u)
  implicit none
  integer, intent(in) :: nx
  real,    intent(out):: u(3,0:nx+1)
  integer :: i

  !   Periodic boundary conditions
  !u(:,0 )=u(:,nx)
  !u(:,nx+1)=u(:,1)

  !  open (outlow) boundary conditions
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)

  return
end subroutine boundaries
!=======================================================================
