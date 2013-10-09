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
  real, parameter :: speed = 1.5, gamma = 5.0/3.0
  !
  !   This is a vector that contains u(x)
  real :: u(3,0:nx+1)
  !  
end module globals
!=======================================================================
!   main program
program euler
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $\Delta t$
  real, parameter :: tmax= 0.2             ! maximumn integration time
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
end program euler
!=======================================================================
! generates initial condition
subroutine initflow(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here is
  
  real :: x
  integer :: i

  !  fill the vector u
  do i=0,nx+1
    x=float(i)*dx   ! obtain the position $x_i$

    ! u(1,i) = density
    ! u(2, i) = momentum
    ! u(3, i) = total energy density    
        
    if (x < 0.5) then
        u(1, i) = 1.0        
        u(3, i) = 1.0 / (gamma - 1)
    elseif (x.eq.0.5) then
        u(1, i) = (1.0 + 0.1)/2.0
        u(3, i) = (1.0 + 0.125) / (gamma - 1)
    else
        u(1, i) = 0.1
        u(3, i) = 0.125 / (gamma - 1)
    end if

    u(2, i) = 0.0

  end do
  
  print*,u(1,1)
  
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
    write(10,*) float(i)*dx, u(1, i), u(2, i), u(3, i)
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
  real :: del, const, min_u
  real :: velocities(0:nx+1)
  real :: coeffs(0:nx+1)
  integer :: i, minimum_i(1)
  
  const = (Co*dx)/(abs(speed))
  ! velocity = momentum / density
  velocities = abs(u(2, 0:nx+1) / u(1, 0:nx+1))

  coeffs = const * 1./velocities
  minimum_i = minloc(coeffs)
  dt=coeffs(minimum_i(1))

  return
end subroutine timestep
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(3, 0:nx+1), f(3, 0:nx+1)
  real :: dtx
  integer :: i

  !  obtain the fluxes
  call fluxes(nx,speed,u,f)

  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  dtx=dt/dx
  do i=1,nx
    up(1,i)=0.5*(u(1,i-1)+u(1,i+1)-dtx*(f(1,i+1)-f(1,i-1)))
    up(2,i)=0.5*(u(2,i-1)+u(2,i+1)-dtx*(f(2,i+1)-f(2,i-1)))
    up(3,i)=0.5*(u(3,i-1)+u(3,i+1)-dtx*(f(3,i+1)-f(3,i-1)))
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,up)

  ! copy the up to the u
  u(:,:) = up(:,:)

  return
end subroutine tstep
!=======================================================================
! Obtain the fluxes F
subroutine fluxes(nx,speed,u,f)
  implicit none
  integer, intent(in) :: nx
  real, intent(in) :: speed
  real, intent(in) :: u(3,0:nx+1)
  real, intent(out):: f(3,0:nx+1)
  integer :: i

  do i=0,nx+1
    f(1,i)=speed*u(1,i)
    f(2,i)=speed*u(2,i)
    f(3,i)=speed*u(3,i)
  end do

  return
end subroutine fluxes
!=======================================================================
! Set boundary conditions
subroutine boundaries(nx,u)
  implicit none
  integer, intent(in) :: nx
  real,    intent(out):: u(3, 0:nx+1)
  integer :: i

  !   Periodic boundary conditions: outflow
  u(1,0)=u(1,1)
  u(1,nx+1)=u(1,nx)
  u(2,0)=u(2,1)
  u(2,nx+1)=u(2,nx)
  u(3,0)=u(3,1)
  u(3,nx+1)=u(3,nx)
  
  return
end subroutine boundaries
!=======================================================================

