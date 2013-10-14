!=======================================================================
!   This program solves the Euler equations with the Godunov Method
!   using the Rusanov Fluxes
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
  real :: u(3,0:nx+1), prim(3,0:nx+1)
  !
end module globals
!=======================================================================
!   main program
program euler_rk2_hll
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real           :: time, dt              !  t, $\Delta t$
  real, parameter :: tmax= .3             ! maximumn integration time
  real, parameter :: dtprint=0.01         ! interval between outputs
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

  ! This subroutine generates the initial conditions
  call initconds(time, tprint, itprint)

  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

    ! updates the primitives
    call u2prim(nx,gamma,u,prim)

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
    call tstep_rk2(dt,time)
    ! time counter increases
    time=time+dt

  end do

  stop
end program euler_rk2_hll
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
      u(3,i)=1.0/(gamma-1.)
    else
      u(1,i)=.125
      u(2,i)=0.0
      u(3,i)=0.1/(gamma-1.)
    end if
    !
    if( (x-0.5*dx <= 0.5).and.(x+0.5*dx >= 0.5) ) then
      u(1,i)=1.125/2.
      u(2,i)=0.0
      u(3,i)=1.1/2./(gamma-1.)
    end if
  end do

  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0

  return
end subroutine initconds
!=======================================================================
! computes the primitives as a function of the Us
subroutine u2prim(nx,gamma,u,prim)
  implicit none
  integer, intent(in) :: nx
  real, intent(in)    :: gamma
  real , intent(in)   :: u(3,0:nx+1)
  real , intent(out)  :: prim(3,0:nx+1)
  integer :: i
  real :: rho, vx

  do i=0,nx+1

    prim(1,i) = u(1,i)
    prim(2,i) = u(2,i)/u(1,i)
    prim(3,i) = (u(3,i)-0.5*prim(1,i)*prim(2,i)**2)*(gamma-1.)

  end do

  return
end subroutine u2prim
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
  write(file1,'(a,i2.2,a)') 'euler_god_hll-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and rho, u(=vx) and P
  do i=1,nx
    write(10,*) float(i)*dx,prim(:,i)
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
  real :: del, cs
  integer :: i

  del=1E30
  do i=0,nx+1
    cs=sqrt(gamma*prim(3,i)/prim(1,i))
    del=min( del,dx/(abs(prim(2,i))+cs) )
  end do
  dt=Co*del

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
!   call rusanov_fluxes(nx,gamma,prim,f)
  call hll_fluxes(nx,gamma,prim,f)

  !   Here is the upwind Godunov method
  dtx=dt/dx
  do i=1,nx
    up(:,i)=u(:,i)-dtx*(f(:,i)-f(:,i-1))
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,up)

  ! copy the up to the u
  u(:,:)=up(:,:)

  return
end subroutine tstep
!=======================================================================
! integration from t to t+dt with the method of Runge-Kutta
subroutine tstep_rk2(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(3,0:nx+1), u1(3,0:nx+1), f(3,0:nx+1)
  real :: dtx
  integer :: i

!== 1st STEP OF RK2 ===
!
  ! updates the primitives
  call u2prim(nx,gamma,u,prim)

  !  obtain the fluxes
  call hll_fluxes(nx,gamma,prim,f)

  !   Here is the upwind Godunov method
  dtx=dt/dx
  do i=1,nx
    up(:,i)=u(:,i)-dtx*(f(:,i)-f(:,i-1))
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(nx,up)

  ! copy the up to the u1
  u1(:,:)=up(:,:)

!== 2nd STEP OF RK2 ===
!
  ! updates the primitives for the intermediate integration step
  call u2prim(nx,gamma,u1,prim)

  !  obtain the fluxes
  call hll_fluxes(nx,gamma,prim,f)

  !   Here is the upwind Godunov method
  do i=1,nx
    up(:,i)=0.5*(u(:,i)+u1(:,i))-0.5*dtx*(f(:,i)-f(:,i-1))
  end do

  ! copy the up to the u
  u(:,:)=up(:,:)

  return
end subroutine tstep_rk2

!=======================================================================
!  computes the HLL fluxes on the entire domain
subroutine hll_fluxes(nx,gamma,prim,f)
  implicit none
  integer, intent(in) :: nx
  real,    intent(in) :: gamma
  real,    intent(in) :: prim(3,0:nx+1)
  real,    intent(out):: f(3,0:nx+1)
  integer :: i
  real :: priml(3), primr(3), ff(3), deltaL(3), deltaR(3), deltaMM(3,0:nx+1)

  do i = 1, nx
    !   we reconstruct the Left and Right states using minmod limited linear
    !   interpolation
    !
    deltaL(:) = prim(:,i  ) - prim(:,i-1)
    deltaR(:) = prim(:,i+1) - prim(:,i  )

    deltaMM(:,i) = (sign(0.5,deltaL(:)) + sign(0.5,deltaR(:))) &
                  * min(abs(deltaL(:)), abs(deltaR(:)))
  end do
  deltaMM(:,   0) = 0.0
  deltaMM(:,nx+1) = 0.0

  do i=0,nx

    !   these are the Left and Right states that enter
    !   the Riemann problem
    ! Ver pagina 7 de presentacion part3.pdf
    primL(:)= prim(:,i  ) + 0.5 * deltaMM(:,i  )
    primR(:)= prim(:,i+1) - 0.5 * deltaMM(:,i+1)

    call prim2hll(gamma, primL, primR, ff)
    f(:,i)=ff(:)

  end do

end subroutine hll_fluxes

!=======================================================================
!  computes the Rusanov fluxes on the entire domain
subroutine rusanov_fluxes(nx,gamma,prim,f)
  implicit none
  integer, intent(in) :: nx
  real,    intent(in) :: gamma
  real,    intent(in) :: prim(3,0:nx+1)
  real,    intent(out):: f(3,0:nx+1)
  integer :: i
  real :: priml(3), primr(3), ff(3)

  do i=0,nx

    !   these are the Left and Right states that enter
    !   the Riemann problem
    primL(:)= prim(:,i  )
    primR(:)= prim(:,i+1)

    call prim2rus(gamma, primL, primR, ff)
    f(:,i)=ff(:)

  end do

end subroutine rusanov_fluxes
!=======================================================================
! Obtain the Rusanov fluxes
subroutine prim2rus(gamma,primL,primR,ff)
  implicit none
  real, intent(in) :: gamma, primL(3), primR(3)
  real, intent(out):: ff(3)
  real :: fL(3), fR(3),ur(3), ul(3)
  real :: lambda, vmax

  lambda=max(abs(primL(2)) + sqrt(gamma*primL(3)/primL(1)) &
           , abs(primR(2)) + sqrt(gamma*primR(3)/primR(1)) )

  call eulerfluxes(gamma,primL,fL)
  call eulerfluxes(gamma,primR,fR)
  call prim2u(gamma, primL, uL)
  call prim2u(gamma, primR, uR)

  ff(:)=0.5*( fl(:)+fr(:) -lambda*( ur(:)-ul(:) )  )


return
end subroutine prim2rus
!=======================================================================
! Obtain the HLL fluxes
subroutine prim2hll(gamma,primL,primR,ff)
  implicit none
  real, intent(in) :: gamma, primL(3), primR(3)
  real, intent(out):: ff(3)
  real :: fL(3), fR(3),ur(3), ul(3)
  real :: lambdaL, lambdaR, sL, sR

  lambdaL = sqrt(gamma*primL(3)/primL(1))
  lambdaR = sqrt(gamma*primR(3)/primR(1))
  sL = min(primL(2) - lambdaL, primR(2) - lambdaR)
  sR = max(primL(2) + lambdaL, primR(2) + lambdaR)

  call eulerfluxes(gamma,primL,fL)
  call eulerfluxes(gamma,primR,fR)
  call prim2u(gamma, primL, uL)
  call prim2u(gamma, primR, uR)

  if (sL >= 0.0) then
    ! Supersonic, only one side flux
    ff(:) = fL(:)
  else if (sR <= 0.0) then
    ! Supersonic, only one side flux
    ff(:) = fR(:)
  else
    ! Average weighting the fluxes
    ff(:) = (sR * fL(:) - sL * fR(:) + sR * sL * ( uR(:) - uL(:) )) / (sR - sL)
  endif

return
end subroutine prim2hll
!=======================================================================
!  computes the euler fluxes, one cell
subroutine eulerfluxes(gamma,pp,ff)
  implicit none
  real, intent(in) :: gamma, pp(3)
  real, intent(out):: ff(3)

  ff(1)=pp(1)*pp(2)
  ff(2)=pp(1)*pp(2)**2+pp(3)
  ff(3)=pp(2)*(0.5*pp(1)*pp(2)**2+gamma*pp(3)/(gamma-1.) )

  return
end subroutine eulerfluxes

!=======================================================================
! computes the primitives as a function of the Us, only in one cell
subroutine prim2u(gamma,pp,uu)
  implicit none
  real, intent(in)    :: gamma
  real , intent(in)   :: pp(3)
  real , intent(out)  :: uu(3)

    uu(1) = pp(1)
    uu(2) = pp(1)*pp(2)
    uu(3) = 0.5*pp(1)*pp(2)**2 +pp(3)/(gamma-1.)

  return
end subroutine prim2u
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

  !  open boundary conditions
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)

  return
end subroutine boundaries
!=======================================================================
