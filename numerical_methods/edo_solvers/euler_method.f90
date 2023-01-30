! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 30/01/23
 
        real function f(x, y)
            implicit none
            real :: x,y
            f = x*y
            return
        end 

        subroutine euler_method(X, Y, N, dx)
            implicit none
            integer :: N, i
            real, dimension(N) :: X, Y
            real :: dx, f
            do i = 1, N-1
                Y(i+1) = Y(i) + dx*f(X(i),Y(i))
            end do
        end subroutine

        program main
 
            implicit none
            integer :: i, N
            real :: xmin, xmax, dx, xo, yo
            real, allocatable :: X(:), Y(:)

            ! Definimos la malla de integración
            dx = 0.001
            xmin = 0.0
            xmax = 1.0
            N = (xmax - xmin)/dx + 1

            ! Condicion inicial
            xo = xmin
            yo = 2

            ! Definimos el tamaño del arreglo
            allocate(X(1:N), Y(1:N))
            do i = 0, N-1
                X(i+1) = xmin + i*dx
                Y(i+1) = 0
            end do
            Y(1) = yo

            ! Llamamos el método de Euler
            call euler_method(X, Y, N, dx)

            ! write array to a file
            open(10, file='resultados.txt')
            do i=1,N
                write(10, '(f8.6,A4,f8.6)')X(i), "    " ,Y(i)
            end do
            close(10)

            ! Despejamos los arreglos de la memoria
            deallocate(X,Y)

 
        end program main
