! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 04/02/23
 
        real function f(x)
            implicit none 
            real :: x
            f = sqrt(x+1)*( (cos(x/2.0))**3 )
            return
        end  

        subroutine bisection(a, b, root, step)
            implicit none
            real :: a, b, c, root, dx, f, tol
            integer :: steps, step, i
            ! definimos una tolerancia para la convergencia
            tol = 0.0001 
            ! definimos un número iteraciones maxima para la convergencia
            steps = 1000
            do i = 1, steps
                c = 0.5*(a+b)
                dx = b-a 

                if ( dx .le. tol ) then
                    step = i
                    root = c 
                    exit
                end if

                if ( f(a)*f(c)>0 ) then
                    a = c
                else if ( f(a)*f(c)<0 ) then
                    b = c
                else 
                    if ( f(a)==0 ) then
                        write(*,*) "Hay una raiz en x = ", a
                        step = i
                        root = a
                        exit
                    end if
                    if ( f(c)==0 ) then
                        write(*,*) "Hay una raiz en x = ", c 
                        step = i
                        root = c 
                        exit
                    end if
                end if
            end do
            if ( i == steps ) then
                write(*,*) "No se encontro raiz que satisfaga la"
                write(*,*) "tolerancia", tol, " en ", steps, " pasos" 
            end if
    
        end subroutine bisection


        program main
 
            implicit none
            real :: a, b, root
            integer :: step

            ! intervalo de inspección
            a = 0.2 ; b = 4.0

            call bisection(a, b, root, step)

            write(*,*) "Se encontro una raiz en ", step," iteraciones"
            write(*,*) "La raiz enocntrada es x = ", root

        end program main