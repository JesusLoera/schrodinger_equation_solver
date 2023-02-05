! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 03/02/23

! En este programa se implementa el método de Runge-Kutta de
! cuarto orden para resolver ecuaciones diferenciales del tipo

!     y''(x) = f(y'(x),y(x),x) ; y(0)=yo y y'(0)=zo

! Usando los siguientes cambios de variable y'(x)=z(X), y''(x)=z'(x)

!    z'(x) = f(z(x),y(x),x)   y   y'(x)=z(x)  ;   y(0)=yo y z(0)=zo

! Este algoritmo es más preciso que el algoritmo de Euler al basarse en
! la derivada central en lugar de la derivada hacia delante.
        
        real function f(z, y, x)
            implicit none
            real :: x,y,z
            f = -9.8
            return
        end  

        subroutine rk4o2(X, Y, Z, N, dx)
            implicit none
            integer :: N, i
            real, dimension(N) :: X, Y, Z
            real :: dx, f, k1, k2, k3, k4, i1, i2, i3, i4
            do i = 1, N-1
                k1 = dx*Z(i)
                i1 = dx*f(Z(i), Y(i), X(i))
                k2 = dx*(Z(i) + (0.5)*i1)
                i2 = dx*f(Z(i)+0.5*i1, Y(i)+0.5*k1, X(i)+0.5*dx)
                k3 = dx*(Z(i) + (0.5)*i2)
                i3 = dx*f(Z(i)+0.5*i2, Y(i)+0.5*k2, X(i)+0.5*dx)
                k4 = dx*(Z(i) + i3)
                i4 = dx*f(Z(i)+i3, Y(i)+k3, X(i)+dx)
                Y(i+1) = Y(i) + (1.0/6.0)*(k1+2*k2+2*k3+k4)
                Z(i+1) = Z(i) + (1.0/6.0)*(i1+2*i2+2*i3+i4)
            end do
        end subroutine

        program main

            implicit none
            integer :: i, N
            real :: xmin, xmax, dx, xo, yo, zo
            real, allocatable :: X(:), Y(:), Z(:)

            ! Definimos la malla de integración
            dx = 0.001 ; xmin = 0.0 ; xmax = 5.0
            N = (xmax - xmin)/dx + 1

            ! Condicion inicial
            xo = xmin  ;  yo = 0  ;  zo = 34.5

            ! Definimos el tamaño del arreglo
            allocate(X(1:N), Y(1:N), Z(1:N))
            do i = 0, N-1
                X(i+1) = xmin + i*dx
                Y(i+1) = 0
                Z(i+1) = 0
            end do
            Y(1)=yo  ; Z(1)=zo

            ! Llamamos el método de Euler
            call rk4o2(X, Y, Z, N, dx)

            ! write array to a file
            open(10, file='rk4o2.dat')
            do i=1,N
                write(10,*) X(i), Y(i), Z(i)
            end do
            close(10)

            ! Despejamos los arreglos de la memoria
            deallocate(X,Y,Z)

        end program main
