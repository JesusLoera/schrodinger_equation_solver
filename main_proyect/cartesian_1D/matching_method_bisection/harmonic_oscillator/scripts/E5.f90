! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 23/02/23

! Este programa resulve la ecn. 1D de Schrodinger
! con V = V(x), Phi(xo) = 0 y Phi(xn) = 0.

! Se resuelve en unidades atómicas
! m_e = 1, hbar = 1, q_e = 1

        real(8) function V(x)
            real(8) :: x, omega
            real(8) :: me
            me = 1.0 ; omega = 1.0
            V = (0.5)*(me*omega**2)*(x**2)
            return
        end function


        program main
            implicit none
            real(8) :: Emin, Emax, En
            real(8) :: xo, xn, dx
            real(8) :: yo, yn, s1, s2
            integer :: n
            ! Valores en la frontera
            xo = -5.0 ; xn = 5.0 ; dx = 0.001
            yo = 0.0  ; yn = 0.0
            ! Parametro libre de disparo
            s1 = 1.0e-5   ;   s2 = 1.0e-5
            ! Nivel de energía a buscar
            n=5
            ! Región [Emin, Emax] para buscar En
            Emin = -10.0 ; Emax = 10.0

            call MMB(Emin, Emax, n, xo, xn, yo, yn, s1, s2, dx, En)

        end program main


        subroutine MMB(Emin, Emax, n, xo, xn, yo, yn, s1, s2, dx, En)

            implicit none
            real(8) :: Emin, Emax, E, dE, En
            real(8) :: xo, xn, dx, xm, xm_aux
            real(8) :: yo, yn, s1, s2, tol1, tol2, coef, g1, g2, g3
            real(8), allocatable :: X(:), Y(:)
            real(8), allocatable :: Yi(:), Yb(:), Y_aux(:), Y_aux2(:)
            real(8) :: A, int_simpson
            integer :: n, nodes, N_size, i, j
            integer :: count_nodes, iter, xm_index
            logical :: search_energy

            ! Definimos las tolerancias
            tol1 = 0.0000001; tol2 = 0.00001

            ! Definimos el tamaño del arreglo
            N_size = (xn - xo)/dx + 1
            allocate(X(1:N_size), Y(1:N_size))
            allocate(Y_aux(1:N_size))
            allocate(Y_aux2(1:N_size))
            allocate(Yi(1:N_size), Yb(1:N_size))
            do i = 1, N_size
                X(i) = xo+(i-1)*dx
            end do

            ! Variable para identificar cuando hallemos En
            search_energy = .true.

            ! Limitamos el número de pasos para hallar En
            iter = 0
            do while (search_energy .and. (iter < 100000))
                iter = iter + 1
                do i = 1, 100000
                    ! Bisectamos el intervalo de energía
                    E = 0.5*(Emin + Emax)
                    ! Hallamos xm y su indice
                    xm = sqrt(2*E)
                    xm_index = 1  ;  xm_aux = xo
                    do while (xm_aux < xm)
                        xm_aux = xm_aux + dx
                        xm_index = xm_index + 1
                    end do
                    ! Integramos hacia delante y hacia atrás
                    Yi(1) = yo ; Yi(2) = s1
                    Yb(N_size) = yn ; Yb(N_size-1) = s2
                    call integrate_in(E, X, Yi, dx, N_size, xm_index)
                    call integrate_bw(E, X, Yb, dx, N_size, xm_index)
                    ! Emparejamos las soluciones de Yi y Yb para E
                    coef = Yi(xm_index)/Yb(xm_index)
                    do j = N_size, xm_index, -1
                        Yb(j) = coef*Yb(j)
                    end do
                    ! Definimos el arreglo Y formado por Yi y Yb
                    do j = 1, N_size
                        if ( j <= xm_index ) then
                            Y(j) = Yi(j)
                        else
                            Y(j) = Yb(j)
                        end if
                    end do
                    ! Contamos los nodos
                    nodes = count_nodes(Y, N_size)
                    if ( nodes < n ) then
                        Emin = E 
                    else if (nodes > n) then
                        Emax = E
                    else if (nodes .eq. n) then
                        exit
                    end if
                end do
                ! Verificamos si encontramos el número de nodos
                ! correcto antes del número máximo de pasos
                if( j .ge. 100000 ) then
                    write(*,*) "No se ha encontrado una solucion con"
                    write(*,*) "el numero correcto de nodos ..."
                    write(*,*) " "
                    write(*,*) "Amplie el rango de la energia o el "
                    write(*,*) "numero de pasos para hallar los nodos."
                    exit
                end if
                ! Calculamos Y_aux(Emin)
                Yi(1) = yo ; Yi(2) = s1
                Yb(N_size) = yn ; Yb(N_size-1) = s2
                call integrate_in(Emin, X, Yi, dx, N_size, xm_index)
                call integrate_bw(Emin, X, Yb, dx, N_size, xm_index)
                ! Emparejamos las soluciones de Yi y Yb para Emin
                coef = Yi(xm_index)/Yb(xm_index)
                do j = N_size, xm_index, -1
                    Yb(j) = coef*Yb(j)
                end do
                ! Definimos el arreglo Y_aux formado por Yi y Yb
                do j = 1, N_size
                    if ( j <= xm_index ) then
                        Y_aux(j) = Yi(j)
                    else
                        Y_aux(j) = Yb(j)
                    end if
                end do
                ! Calculamos Y_aux2(Emin)
                Yi(1) = yo ; Yi(2) = s1
                Yb(N_size) = yn ; Yb(N_size-1) = s2
                call integrate_in(Emax, X, Yi, dx, N_size, xm_index)
                call integrate_bw(Emax, X, Yb, dx, N_size, xm_index)
                ! Emparejamos las soluciones de Yi y Yb para E
                coef = Yi(xm_index)/Yb(xm_index)
                do j = N_size, xm_index, -1
                    Yb(j) = coef*Yb(j)
                end do
                ! Definimos el arreglo Y_aux2 formado por Yi y Yb
                do j = 1, N_size
                    if ( j <= xm_index ) then
                        Y_aux2(j) = Yi(j)
                    else
                        Y_aux2(j) = Yb(j)
                    end if
                end do
                ! Bisectamos la energía
                g1 = Y(xm_index+1) + Y(xm_index-1) - 2*Y(xm_index)
                g2 = Y_aux(xm_index+1) + Y_aux(xm_index-1)
                g2 = g2 - 2*Y_aux(xm_index)
                !g3 = Y_aux2(xm_index+1) + Y_aux2(xm_index-1)
                !g3 = g3 - 2*Y_aux2(xm_index)
                if ( (g1*g2)>0 ) then
                    Emin = E 
                else if ( (g1*g2)<0 ) then
                    Emax = E
                end if
                ! Evaluamos el criterio de convergencia
                dE = abs(Emax-Emin)

                if ((dE.le.tol1) .or. (abs(g1).le.tol2)) then
                    if ( dE.le.tol1 ) then
                        write(*,*) "Hubo convergencia en dE"
                    else if (abs(g1).le.tol2) then
                        write(*,*) "Hubo convergencia en G(E)"
                    end if
                    En = E
                    search_energy = .false.
                end if
            end do
            ! Normalización de la función de onda
            A = int_simpson(Y, dx, N_size)
            A = sqrt(A)
            do i = 1, N_size
                Y(i) = Y(i)/A
            end do
            ! Imprimimos información útil
            write(*,*) "E = ", En
            write(*,"(A7, F10.7)") "G(E) = ", g1
            write(*,*) "Iteraciones = ", iter
            ! Escribimos la solución en un archivo de texto
            open(10, file='../data/wf5.dat')
            write(10,"(A12, I2, A3, F10.4)") "# Eigenvalor", n," = ",E  
            write(10,*) "# X      Psi(x)     |Psi(x)|^2"
            do i=1,N_size
                write(10,*) X(i), Y(i), Y(i)**2
            end do
            close(10)

        end subroutine



        subroutine integrate_in(E, X, Y, dx, N_size, xm_index)
            real(8) :: E, V, dx, Vi
            real(8) :: me, hbar
            integer :: N_size, i, xm_index
            real(8), dimension(N_size) :: X, Y
            ! Unidades atómicas    
            hbar = 1.0 ; me = 1.0
            ! Integración de la ED de Schrodinger
            do i = 3, xm_index
                Vi = V(X(i-1))
                Y(i)=2*(1+me*(dx**2)/(hbar**2)*(Vi-E))*Y(i-1)-Y(i-2)
            end do
        end subroutine



        subroutine integrate_bw(E, X, Y, dx, N_size, xm_index)
            real(8) :: E, V, dx, Vi
            real(8) :: me, hbar
            integer :: N_size, i, xm_index
            real(8), dimension(N_size) :: X, Y
            ! Unidades atómicas    
            hbar = 1.0 ; me = 1.0
            ! Integración de la ED de Schrodinger
            do i = N_size-2, xm_index, -1
                Vi = V(X(i))
                Y(i)=2*(1+me*(dx**2)/(hbar**2)*(Vi-E))*Y(i+1)-Y(i+2)
            end do
        end subroutine



        real(8) function int_simpson(Y, dx, N_size)
            real(8) ::  dx, integral, sum2, sum4
            real(8), dimension(N_size) :: Y 
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
            real(8), dimension(N_size) :: Y
            nodos = 0
            do i = 2, N_size - 1
                if ( Y(i)*Y(i+1) < 0 ) then
                    nodos = nodos + 1
                end if
            end do
            count_nodes = nodos
            return
        end function