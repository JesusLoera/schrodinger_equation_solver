! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 10/02/23

! En este programa resolvemos ED lineales de la forma 

!       y''(x) + P(x)y'(x) + Q(x)y(x) = R(x)  ;  y(xo)=yo,  y(xn)=yn
 
        program main
 
            implicit none
            real :: xo, xn, yo, yn, dx, s
            integer :: i, N
            real, allocatable :: X(:), Y(:), Z(:)

            ! Condiciones en la frontera
            xo = 0.0 ; xn = 5
            yo = 0.0 ; yn = 50

            ! Parametro libre de disparo
            s = 1

            ! Subintervalo de la malla de integración
            dx = 0.001
            N = (xn - xo)/dx + 1

            allocate(X(1:N), Y(1:N), Z(1:N))

            call linear_shooting(xo, xn, yo, yn, dx, X, Y, Z, N, s)

            ! write array to a file
            open(10, file='linear_shooting.dat')
            do i=1,N
                write(10,*) X(i), Y(i), Z(i)
            end do
            close(10)

 
        end program main

        ! Parametro de la Ed
        real function P(x)
            implicit none
            real :: x
            P = 0.0
            return
        end function

        ! Parametro de la ED
        real function Q(x)
            implicit none
            real :: x
            Q = 0.0
            return
        end function

        ! Parametro de la ED
        real function R(x)
            implicit none
            real :: x
            R = -9.8
            return
        end function

        ! Función para Runge-Kutta
        real function f(z, y, x)
            real z, y, x
            real P, Q, R
            f = R(x)-P(x)*z-Q(x)*y 
            return
        end function

        ! Función para Runge-Kutta
        real function g(z, y, x)
            real z, y, x
            real P, Q
            f = -P(x)*z-Q(x)*y 
            return
        end function

        ! Método de Runge-Kutta para resolver edo ivp
        subroutine rk4o2_f(xo, xn, yo, zo, dx, X, Y, Z, N)
            implicit none
            integer :: i
            integer :: N
            real :: xo, xn, yo, zo, dx
            real, dimension(N) :: X, Y, Z
            real :: f, k1, k2, k3, k4, i1, i2, i3, i4

            ! Definimos la malla de integración
            N = (xn - xo)/dx + 1
            do i = 0, N-1
                X(i+1) = xo + i*dx
                Y(i+1) = 0
                Z(i+1) = 0
            end do
            Y(1) = yo
            Z(1) = zo

            do i = 1, N-1
                k1 = Z(i)
                i1 = f(Z(i), Y(i), X(i))
                k2 = Z(i) + (0.5)*i1*dx
                i2 = f(Z(i)+0.5*i1*dx, Y(i)+0.5*k1*dx, X(i)+0.5*dx)
                k3 = Z(i) + (0.5)*i2*dx
                i3 = f(Z(i)+0.5*i2*dx, Y(i)+0.5*k2*dx, X(i)+0.5*dx)
                k4 = Z(i) + i3*dx
                i4 = f(Z(i)+i3*dx, Y(i)+k3*dx, X(i)+dx)
                Y(i+1) = Y(i) + (1.0/6)*dx*(k1+2*k2+2*k3+k4)
                Z(i+1) = Z(i) + (1.0/6)*dx*(i1+2*i2+2*i3+i4)
            end do
        end subroutine rk4o2_f

        ! Método de Runge-Kutta para resolver edo ivp
        subroutine rk4o2_g(xo, xn, yo, zo, dx, X, Y, Z, N)
            implicit none
            integer :: i
            integer :: N
            real :: xo, xn, yo, zo, dx
            real, dimension(N) :: X, Y, Z
            real :: g, k1, k2, k3, k4, i1, i2, i3, i4

            ! Definimos la malla de integración
            N = (xn - xo)/dx + 1
            do i = 0, N-1
                X(i+1) = xo + i*dx
                Y(i+1) = 0
                Z(i+1) = 0
            end do
            Y(1) = yo
            Z(1) = zo

            do i = 1, N-1
                k1 = Z(i)
                i1 = g(Z(i), Y(i), X(i))
                k2 = Z(i) + (0.5)*i1*dx
                i2 = g(Z(i)+0.5*i1*dx, Y(i)+0.5*k1*dx, X(i)+0.5*dx)
                k3 = Z(i) + (0.5)*i2*dx
                i3 = g(Z(i)+0.5*i2*dx, Y(i)+0.5*k2*dx, X(i)+0.5*dx)
                k4 = Z(i) + i3*dx
                i4 = g(Z(i)+i3*dx, Y(i)+k3*dx, X(i)+dx)
                Y(i+1) = Y(i) + (1.0/6)*dx*(k1+2*k2+2*k3+k4)
                Z(i+1) = Z(i) + (1.0/6)*dx*(i1+2*i2+2*i3+i4)
            end do
        end subroutine rk4o2_g

        ! Linear Shooting Method para resvolver ed bvp
        subroutine linear_shooting(xo, xn, yo, yn, dx, X, Y, Z, N, s)
            integer :: N, i
            real :: xo, xn, yo, dx, s
            real, dimension(N) :: X, Y, Z
            real, dimension(N) :: Ux, Uy, Uz
            real, dimension(N) :: Vx, Vy, Vz

            write(*,*) "Iniciando metodo de disparo lineal"

            call rk4o2_f(xo, xn, yo, 0.0, dx, Ux, Uy, Uz, N)
            call rk4o2_g(xo, xn, 0.0, s, dx, Vx, Vy, Vz, N)

            if ( Vy(N) .ne. 0 ) then
                do i = 1, N
                    Y(i) = Uy(i) + (yn - Uy(N))/(Vy(N))*Vy(i)
                    Z(i) = Uz(i) + (yn - Uy(N))/(Vy(N))*Vz(i)
                    X(i) = xo + (i-1)*dx
                end do
                write(*,*)"El parametro inicial de disparo fue correcto"
            else 
                write(*,*) "Pruebe con otro valor del parametro de"
                write(*,*) "disparo libre"
            end if

        end subroutine