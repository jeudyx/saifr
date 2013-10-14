program prueba
    real :: u(3,0:10)
    real :: z(3,0:10)
    real :: otro(0:10)
    real :: velocities(0:10)
    integer :: i, minimum_i(1)
    real :: dt
    real :: prim(3)
    i = 0
    
    do while (i.le.10)
        u(1,i) = 1
        i = i + 1
    end do

    i = 0

    do while (i.le.10)
        u(2,i) = 2
        i = i + 1
    end do
    
    i = 0

    do while (i.le.10)
        u(3,i) = -3
        i = i + 1
    end do

    !otro(:) = abs(u(3, 0:10) / u(2, 0:10))
    !otro(5) = -1
    !print*, otro(minloc(otro)-1)
    !print*, otro    

    velocities = abs(u(2, 0:nx+1) / u(1, 0:nx+1))
    minimum_i = minloc(velocities)
    dt = velocities(minimum_i(1)-1)

    prim(:) = 10
    print*, prim

end program prueba
