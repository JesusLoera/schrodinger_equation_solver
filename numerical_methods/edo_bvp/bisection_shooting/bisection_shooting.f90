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
            f = -9.8
            return
        end  

        ! Definimos la función de convergencia
        real(8) function g(Un, Fn)
            implicit none 
            real(8) Un, Fn
            g = Un - Fn
            return
        end

        ! Definimos el programa principal    
        program main

            implicit none
            integer :: i, N, step
            real(8) :: xmin, xmax, dx, xo, yo, yn
            real(8) :: emax, emin, root
            real(8), allocatable :: X(:), Y(:), Z(:)

            ! Definimos la malla de integración
            dx = 0.0001 ; xmin = 0.0 ; xmax = 5.0
            N = (xmax - xmin)/dx + 1

            ! Definimos el intervalo de busqueda para adivinar
            ! la condición inicial que define el bvp
            emin = 0.0 ; emax = 50.0

            ! Condicion frontera
            xo = xmin  ;  yo = 0   ;  yn = 50

            ! Definimos el tamaño del arreglo
            allocate(X(1:N), Y(1:N), Z(1:N))
            do i = 0, N-1
                X(i+1) = xmin + i*dx
                Y(i+1) = 0
                Z(i+1) = 0
            end do
            Y(1) = yo
            
            ! Llamamos al método de bisection shooting
            call bisect_shooting(X,Y,Z,N,dx,emin,emax,yn,root,step)
            write(*,*) "Se encontro una raiz en ", step," iteraciones"
            write(*,*) "La raiz encontrada es x = ", root

            ! write array to a file
            open(10, file='bisect_shooting.dat')
            do i=1,N
                write(10,*) X(i), Y(i), Z(i)
            end do
            close(10)

            ! Despejamos los arreglos de la memoria
            deallocate(X,Y,Z)

        end program main

        ! Método de Runge-Kutta para resolver edo ivp

        subroutine rk4o2(X, Y, Z, N, dx)
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
        end subroutine rk4o2

        ! Método de bisección para hallar el valor inicial
        ! que satisfaga las condiciones frontera exigidas.

        subroutine bisect_shooting(X,Y,Z,N,dx,a,b,Fn,root,step)
            implicit none
            real(8) :: a, b, E, root, g, tol, Fn, dx, Un, dE, dy, Ua
            real(8), dimension(N) :: X, Y, Z
            integer :: steps, step, i, N
            ! definimos una tolerancia para la convergencia
            tol = 0.000001
            ! definimos un número iteraciones maxima para la convergencia
            steps = 10000

            open(11, file='energy_shooting.dat')

            do i = 1, steps
                E = 0.5*(a+b)
                dE = b-a
                Z(1) = a 
                call rk4o2(X, Y, Z, N, dx)
                Ua = Y(N) 
                Z(1) = E
                call rk4o2(X, Y, Z, N, dx)
                Un = Y(N)
                dy = Un - Fn 
                write(11,*) E
                if (  abs(dE) .le. tol ) then
                    step = i
                    root = E 
                    exit
                end if

                if ( g(Ua, Fn)*g(Un, Fn)>0 ) then
                    a = E
                else if ( g(Ua, Fn)*g(Un, Fn)<0 ) then
                    b = E
                else 
                    if ( g(Ua, Fn)==0 ) then
                        write(*,*) "Hay una raiz en x = ", a
                        step = i
                        root = a
                        exit
                    end if
                    if ( g(Un, Fn)==0 ) then
                        write(*,*) "Hay una raiz en x = ", E
                        step = i
                        root = E
                        exit
                    end if
                end if
            end do
            if ( i == steps ) then
                write(*,*) "No se encontro raiz que satisfaga la"
                write(*,*) "tolerancia", tol, " en ", steps, " pasos" 
            end if
            close(11)
    
        end subroutine bisect_shooting
