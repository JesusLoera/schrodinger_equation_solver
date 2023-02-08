! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 29/02/23

        real(8) function f(x)
            implicit none 
            real(8) :: x
            f = sqrt(x+1)*( (cos(x/2.0))**3 )
            return
        end  

        real(8) function derivada(x)
            implicit none 
            real(8) :: f, x, h
            h = 0.00001
            derivada = (f(x+h) - f(x-h))/(2.0*h)
            return
        end

        subroutine newton_raphson(xo, root, N)
            implicit none
            real(8), intent(in) :: xo
            real(8), intent(out) ::  root
            integer, intent(out) ::  N
            real(8) :: tolerancia, f, derivada, x, dx, xi
            integer :: steps, i
            logical :: rootFound = .false.

            tolerancia = 1.d-8
            steps = 1000
            xi = xo

            do i = 1, steps
                x = xi - f(xi)/derivada(xi)
                dx = abs(x-xi)
                if ( dx < tolerancia ) then
                    rootFound = .true.
                    root = x
                    N = i
                    exit
                end if
                xi = x
            end do

            if ( rootFound ) then
                write(*,*) "Se encontro una raiz en x=", x
                write(*,*) "Numero de pasos para encontrar la raiz: ",N
            else 
                write(*,*) "El metodo de Newton-Raphson no convergio"                
            end if
        end subroutine newton_raphson

        program main
            implicit none
            real(8) :: xo, root
            integer :: N
            xo = 1
            call newton_raphson(xo, root, N)
        end program main