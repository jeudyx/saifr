!=======================================================================
!   This program solves the Euler equations with the Lax Method in 2D
!=======================================================================

!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
  integer, parameter :: nx=200, ny=100
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=0.5, ymax=0.25, dx=xmax/float(nx), dy=ymax/float(ny)
  real, parameter :: gamma=1.4, pi = 3.14159265359
  !   This is a vector that contains u(x)
  real :: u(4,0:nx+1,0:ny+1)
  !  
end module globals
!=======================================================================
!   main program
program euler_variable_jet
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real           :: time, dt             !  t, $\Delta t$
  real, parameter :: tmax= 1.5             ! maximumn integration time
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
end program euler_variable_jet
!=======================================================================
! generates initial condition
subroutine initconds(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here is the Sod tube test
  integer :: i, j
  real :: x, y

  !  fill the vector u
  do i=0,nx+1
    do j = 0, ny+1
        x=float(i)*dx   ! obtain the position $x_i$
        y=float(j)*dy   ! obtain the position $y_i$        
        u(1,i,j)=1.0    ! rho
        u(2,i,j)=0.0    ! vx
        u(3,i,j)=0.0    ! vy
        u(4,i,j)=0.01/(gamma-1.)
    end do
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
  character (len=20) filerho, fileE, fileVel
  integer :: i,j
  real :: rho,vx,P
  ! open output files
  write(filerho,'(a,i2.2,a)') 'rho-jet-',itprint,'.dat'
  open(unit=10,file=filerho,status='unknown')

  write(fileE,'(a,i2.2,a)') 'ener-jet-',itprint,'.dat'
  open(unit=11,file=fileE,status='unknown')

  write(fileVel,'(a,i2.2,a)') 'vel-jet-',itprint,'.dat'
  open(unit=12,file=fileVel,status='unknown')

  ! writes data
  do j=1,ny
    write(10,*) u(1,:,j)
    write(11,*) u(4,:,j)

    ! the velocity in different format
    do i=1,nx
      write(12,*) float(i)*dx,float(j)*dy,u(2,i,j)/u(1,i,j),u(3,i,j)/u(1,i,j)
    end do

  end do

  ! closes output files
  close(10)
  close(11)
  close(12)

return
end subroutine output
!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  ! Courant number =0.9
  real, parameter :: Co=0.7
  real :: rho, vx , vy, P, cs
  integer :: i, j

  dt=1E30
  do i=0,nx+1
    do j=0, ny+1
        rho=u(1,i,j)
        vx=u(2,i, j)/rho
        vy=u(3,i, j)/rho
        P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)
        cs=sqrt(gamma*P/rho)
        dt=min( dt,Co*dx/(abs(vx)+cs), dt,Co*dy/(abs(vy)+cs) )
    end do
  end do

  return
end subroutine timestep
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(4,0:nx+1,0:ny+1), f(4,0:nx+1,0:ny+1), g(4,0:nx+1,0:ny+1)
  real :: dtx, dty
  integer :: i, j

  !  obtain the fluxes
  call fluxes(nx,ny,gamma,u,f,g)

  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  dtx=dt/(2.*dx)
  dty=dt/(2.*dy)
  do i=1,nx
    do j=1, ny
        up(:,i,j)=0.25*(u(:,i-1,j)+u(:,i+1,j)+u(:,i,j-1)+u(:,i,j+1)) - dtx*(f(:,i+1,j)-f(:,i-1,j)) - dty*(g(:,i,j+1)-g(:,i,j-1))
    end do
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,ny,up,time,pi,gamma)

  ! copy the up to the u
  u(:,:,:)=up(:,:,:)

  return
end subroutine tstep
!=======================================================================
! Obtain the fluxes F
subroutine fluxes(nx,ny,gamma,u,f,g)
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: gamma
  real, intent(in) :: u(4,0:nx+1,0:ny+1)
  real, intent(out):: f(4,0:nx+1,0:ny+1), g(4,0:nx+1,0:ny+1)
  integer :: i, j
  real :: rho, vx, vy , P

  do i=0,nx+1
    do j=0,ny+1
        rho=u(1,i,j)
        vx= u(2,i,j)/rho
        vy= u(3,i,j)/rho
        P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)

        f(1,i,j)=rho*vx
        f(2,i,j)=rho*vx**2+P
        f(3,i,j)=rho*vx*vy
        f(4,i,j)=vx*(u(4,i,j)+P)

        g(1,i,j)=rho*vy
        g(2,i,j)=rho*vy*vx
        g(3,i,j)=rho*vy**2+P
        g(4,i,j)=vy*(u(4,i,j)+P)
    end do
  end do

  return
end subroutine fluxes
!=======================================================================
! Set boundary conditions
subroutine boundaries(nx,ny,u,t,pi,gamma)
  implicit none
  integer, intent(in) :: nx, ny
  real,    intent(in)::t, pi, gamma
  real,    intent(out):: u(4,0:nx+1,0:ny+1)
  integer :: i,j
  real :: vx, rho

  do j=0,ny+1      
      ! Open boundary condition at i = Nx, j from 0 to ny
      u(:,nx+1,j)=u(:,nx,j)
      ! Reflective boundary condition on i = 0 except in a region 10 cells wide
      if(j.ge.45 .and. j.le.54) then
        vx = 1.0 + (0.5 * sin((2.0 * pi * t) / 0.2))
        rho = 1.0
        u(1,0,j) = rho
        u(2,0,j) = vx
        u(3,0,j) = 0.0  ! vy
        u(4,0,j) = rho * (0.5 * vx**2) + (0.1 / (gamma - 1.0))
      else
        u(:,0,j) = -1. * u(:,1,j)
      endif
  end do
  
  ! Bottom and top
  do i=0,nx+1
      u(:,i,0)=u(:,i,1)
      u(:,i,ny+1)=u(:,i,ny)
  end do

  return
end subroutine boundaries
!=======================================================================

