! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 17/02/23

! Este programa resulve la ecn. 1D de Schrodinger
! con V = V(x), Phi(xo) = 0 y Phi(xn) = 0.

! Se resuelve en unidades atómicas
! m_e = 1, hbar = 1, q_e = 1

        real function V(x)
            real :: x, omega
            real :: me
            me = 1.0 ; omega = 1.0
            V = (0.5)*(me*omega**2)*(x**2)
            return
        end function
 
        program main
            implicit none
            real :: Emin, Emax, En
            real :: xo, xn, dx
            real :: yo, s
            integer :: n
            ! Valores en la frontera
            xo = -5.0 ; xn = 5.0 ; dx = 0.001
            yo = 0.0
            ! Parametro libre de disparo
            s = 1.0e-5
            ! Nivel de energía a buscar
            n = 3
            ! Región [Emin, Emax] para buscar En
            Emin = -10 ; Emax = 10

            call BEM(Emin, Emax, n, xo, xn, yo, s, dx, En)

        end program main

        subroutine integrate_SE(E, X, Y, dx, N_size)
            real :: E, V, dx, Vi
            real :: me, hbar
            integer :: N_size, i
            real, dimension(N_size) :: X, Y
            ! Unidades atómicas    
            hbar = 1.0 ; me = 1.0
            ! Integración de la ED de Schrodinger
            do i = 3, N_size
                Vi = V(X(i-1))
                Y(i)=2*(1+me*(dx**2)/(hbar**2)*(Vi-E))*Y(i-1)-Y(i-2)
            end do
        end subroutine

        subroutine BEM(Emin, Emax, n, xo, xn, yo, s, dx, En)

            implicit none
            real :: Emin, Emax, E, dE, En
            real :: xo, xn, dx
            real :: yo, s, tol1, tol2
            real, allocatable :: X(:), Y(:)
            real, allocatable :: X_aux(:), Y_aux(:)
            real :: A, int_simpson
            integer :: n, nodes, N_size, i, count_nodes, iter
            logical :: search_energy

            ! Definimos las tolerancias
            tol1 = 0.0001; tol2 = 0.0001

            ! Definimos el tamaño del arreglo
            N_size = (xn - xo)/dx + 1
            allocate(X(1:N_size), Y(1:N_size))
            allocate(X_aux(1:N_size), Y_aux(1:N_size))

            do i = 1, N_size
                X(i) = xo+(i-1)*dx
                Y(i) = 0
                X_aux(i) = X(i) ; Y_aux(i) = Y(i)
            end do
            Y(1) = yo ; Y(2) = s
            Y_aux(1) = yo ; Y_aux(2) = s
            
            ! Variable para identificar cuando hallemos En
            search_energy = .true.

            ! Limitamos el número de pasos para hallar En
            iter = 0
            do while (search_energy .and. (iter < 100000))

                iter = iter + 1

                do i = 1, 100000
                    E = 0.5*(Emin + Emax)
                    call integrate_SE(E, X, Y, dx, N_size)
                    nodes = count_nodes(Y, N_size)
                    if ( nodes < n ) then
                        Emin = E 
                    else if (nodes > n) then
                        Emax = E
                    else if (nodes .eq. n) then
                        exit
                    end if
                end do

                call integrate_SE(Emin,X_aux,Y_aux,dx,N_size)

                if ((Y(N_size)*Y_aux(N_size))>0) then
                    Emin = E 
                else if ((Y(N_size)*Y_aux(N_size))<0) then
                    Emax = E
                end if

                dE = abs(Emax-Emin)
                if ((dE.le.tol1) .or. (abs(Y(N_size)).le.tol2)) then
                    En = E
                    search_energy = .false.
                end if

            end do

            write(*,*) iter, E

            ! Normalización de la función de onda

            A = int_simpson(Y, dx, N_size)
            A = sqrt(A)

            do i = 1, N_size
                Y(i) = Y(i)/A
            end do

            ! write array to a file
            open(10, file='wave_function.dat')
            do i=1,N_size
                write(10,*) X(i), Y(i)
            end do
            close(10)

            ! write array to a file
            open(11, file='square_wave_function.dat')
            do i=1,N_size
                write(11,*) X(i), Y(i)**2
            end do
            close(11)

        end subroutine

        real function int_simpson(Y, dx, N_size)
            real ::  dx, integral, sum2, sum4
            real, dimension(N_size) :: Y 
            integer :: i
            sum4 = 0.0 ; sum2 = 0.0
            do i = 2, N_size-1, 2
                sum4 = sum4 + Y(i)**2
            end do
            do i = 3, N_size-2, 2
                sum2 = sum2 + Y(i)**2
            end do
            integral = Y(1)**2 + Y(N_size)**2 + 2*sum2 + 4*sum4
            integral = (dx/3.0)*(integral)
            int_simpson = integral
            return
        end function

        integer function count_nodes(Y, N_size)
            integer :: N_size, nodos, i
            real, dimension(N_size) :: Y
            nodos = 0
            do i = 2, N_size - 1
                if ( Y(i)*Y(i+1) < 0 ) then
                    nodos = nodos + 1
                end if
            end do
            count_nodes = nodos
            return
        end function