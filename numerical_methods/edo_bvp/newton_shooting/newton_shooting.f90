! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 04/02/23

! En este programa se implementa el método de shooting para resolver
! edo sujetas a bvp del tipo

!     y''(x) = f(y'(x),y(x),x) ; y(xo)=a , y(xn)=b

! Convertimos la edo bvp en una edo ivp y tratamos de adivinar
! la condición inicial que satisfaga las condiciones en la frontera

!     u_E''(x) = f(u_E'(x),u_E(x),x) ; u_E(xo)=a , u'_E(xo)=E

! Consideramos la diferencia entre dos soluciones en la frontera

!     g(E) = u_E(xn) - y(xn) = u_E(xn) - b 

! Se adivina E y se resuelve la ivp hasta satisfacer g(E) = 0

        ! Definimos la edo a resolver
        real(8) function f(z, y, x)
        implicit none
        real(8) :: x,y,z
        f = -2*z - y
        return
    end  

    real(8) function k(z, y, x)
        implicit none
        real(8) :: x,y,z
        k = -2*y-z
        return
    end  

    ! Definimos el programa principal    
    program main

        implicit none
        integer :: i, N, step
        real(8) :: xmin, xmax, dx, xo, yo, yn
        real(8) :: eo, root
        real(8), allocatable :: X(:), Y(:), Z(:), Exact(:)

        ! Definimos la malla de integración
        dx = 0.0001 ; xmin = 0.0 ; xmax = 1
        N = (xmax - xmin)/dx + 1

        ! Definimos la estimación incial del shooting
        eo = 1000

        ! Condicion frontera
        xo = xmin  ;  yo = 1   ;  yn = 3

        ! Definimos el tamaño del arreglo
        allocate(X(1:N), Y(1:N), Z(1:N), Exact(1:N))
        do i = 0, N-1
            X(i+1) = xmin + i*dx
            Y(i+1) = 0
            Z(i+1) = 0
        end do
        Y(1) = yo
        Z(1) = eo

        do i=1, N
            Exact(i) = exp(-X(i))+(3*(2.71828)-1)*X(i)*exp(-X(i))
        end do
        
        ! Llamamos al método de bisection shooting
        call newton_shooting(X,Y,Z,N,dx,eo,yn,root,step)
        write(*,*) "Se encontro una raiz en ", step," iteraciones"
        write(*,*) "La raiz encontrada es x = ", root

        ! Solución numérica a un archivo
        open(10, file='newton_shooting.dat')
        do i=1,N
            write(10,*) X(i), Y(i), Z(i)
        end do
        close(10)

        ! Solución exacta en un archivo
        open(12, file='exact_sol.dat')
        do i=1,N
            write(12,*) X(i), Exact(i)
        end do
        close(12)

        ! Despejamos los arreglos de la memoria
        deallocate(X,Y,Z, Exact)

    end program main

    ! Método de Runge-Kutta para resolver edo ivp

    subroutine rk4o2_f(X, Y, Z, N, dx)
        implicit none
        integer :: N, i
        real(8), dimension(N) :: X, Y, Z
        real(8) :: dx, f, k1, k2, k3, k4, i1, i2, i3, i4
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

    subroutine rk4o2_k(X, Y, Z, N, dx)
        implicit none
        integer :: N, i
        real(8), dimension(N) :: X, Y, Z
        real(8) :: dx, k, k1, k2, k3, k4, i1, i2, i3, i4
        do i = 1, N-1
            k1 = Z(i)
            i1 = k(Z(i), Y(i), X(i))
            k2 = Z(i) + (0.5)*i1*dx
            i2 = k(Z(i)+0.5*i1*dx, Y(i)+0.5*k1*dx, X(i)+0.5*dx)
            k3 = Z(i) + (0.5)*i2*dx
            i3 = k(Z(i)+0.5*i2*dx, Y(i)+0.5*k2*dx, X(i)+0.5*dx)
            k4 = Z(i) + i3*dx
            i4 = k(Z(i)+i3*dx, Y(i)+k3*dx, X(i)+dx)
            Y(i+1) = Y(i) + (1.0/6)*dx*(k1+2*k2+2*k3+k4)
            Z(i+1) = Z(i) + (1.0/6)*dx*(i1+2*i2+2*i3+i4)
        end do
    end subroutine rk4o2_k

    ! Método de bisección para hallar el valor inicial
    ! que satisfaga las condiciones frontera exigidas.

    subroutine newton_shooting(X,Y,Z,N,dx,Eo,Fn,root,step)
        implicit none
        real(8) :: Eo, E, root, tol, Fn, dx
        real(8), dimension(N) :: X, Y, Z
        integer :: steps, step, i, N
        real(8), allocatable :: Vx(:), Vy(:), Vz(:)
        ! definimos una tolerancia para la convergencia
        tol = 0.000001
        ! definimos un número iteraciones maxima para la convergencia
        steps = 10000

        ! Definimos el tamaño del arreglo
        allocate(Vx(1:N), Vy(1:N), Vz(1:N))
        do i = 1, N
            Vx(i) = X(i)
            Vy(i) = 0
            Vz(i) = 0
        end do
        Vz(1) = 1

        open(11, file='energy_shooting.dat')

        do i = 1, steps

            Z(1) = Eo
            call rk4o2_f(X, Y, Z, N, dx)
            call rk4o2_k(Vx, Vy, Vz, N, dx)

            E = Eo - (Y(N)-Fn)/(Vy(N))

            write(11,*) E
            if (  abs(E-Eo) .le. tol ) then
                step = i
                root = E 
                exit
            end if

            Eo = E

        end do
      
        if ( i == steps ) then
            write(*,*) "No se encontro raiz que satisfaga la"
            write(*,*) "tolerancia", tol, " en ", steps, " pasos" 
        end if
        close(11)

    end subroutine newton_shooting
