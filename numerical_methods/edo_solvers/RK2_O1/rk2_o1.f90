! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 01/02/23

! En este programa se implementa el método de leap-frog o del salto de 
! rana para resolver una ecuación diferencial de orden uno del tipo

!           y'(x) = f(y(X),x)

! Este algoritmo es más preciso que el algoritmo de Euler al basarse en
! la derivada central en lugar de la derivada hacia delante.
 
        real function f(x, y)
        implicit none
        real :: x,y
        f = x*y
        return
    end  

    subroutine rk2_o1(X, Y, N, dx)
        implicit none
        integer :: N, i
        real, dimension(N) :: X, Y
        real :: dx, f, k1, k2
        ! El segundo elemento inicial se calcula con la formula del
        ! método de euler y después los otros N-2 con leap-frog
        do i = 1, N-1
            k1 = f(X(i),Y(i))
            k2 = f(X(i) + 0.5*dx , Y(i) + 0.5*k1*dx )
            Y(i+1) = Y(i) + k2*dx
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
        call rk2_o1(X, Y, N, dx)

        ! write array to a file
        open(10, file='rk2_o1.dat')
        do i=1,N
            write(10, '(f8.6,A4,f8.6)')X(i), "    " ,Y(i)
        end do
        close(10)

        ! Despejamos los arreglos de la memoria
        deallocate(X,Y)

    end program main
