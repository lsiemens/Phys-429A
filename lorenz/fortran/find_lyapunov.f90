!
! find_lyapunov.f90
! Calculate lyapunov exponents for the lorenz attractor
! solutions calculated using rungekutta
!

program main
    use constants, only : rp
    implicit none

    real(rp) :: tmax, perturbation_size, x_0, y_0, z_0, normA, normB, normC, A(3), B(3), C(3)
    real(rp) :: h, invh, invperturbation_size, L1, L2, L3, L1err, L2err, L3err, dim, dimerr
    real(rp), allocatable :: exponents(:, :)
    real(rp) :: lorenz_solution(3), perturbed1(3), perturbed2(3), perturbed3(3)
    integer :: steps, renorm, i
!    integer :: uid

    call load_settings(steps, tmax, perturbation_size, renorm, x_0, y_0, z_0)
    call initalize_data(exponents, steps / renorm)
!    uid = save_solution_open("solution.txt")
    lorenz_solution = 0.0_rp
    perturbed1 = 0.0_rp
    perturbed2 = 0.0_rp
    perturbed3 = 0.0_rp
    lorenz_solution(:) = [x_0, y_0, z_0]
    perturbed1(:) = lorenz_solution(:) + [1.0_rp, 0.0_rp, 0.0_rp]*perturbation_size
    perturbed2(:) = lorenz_solution(:) + [0.0_rp, 1.0_rp, 0.0_rp]*perturbation_size
    perturbed3(:) = lorenz_solution(:) + [0.0_rp, 0.0_rp, 1.0_rp]*perturbation_size
    h = tmax/real(steps, rp)
    invh = 1.0_rp/h
    invperturbation_size = 1.0_rp/perturbation_size
    do i = 1, steps - 1
!        lorenz_solution(:) = runge_kutta4(lorenz_solution(:), h, lorenz)
!        perturbed1(:) = runge_kutta4(perturbed1(:), h, lorenz)
!        perturbed2(:) = runge_kutta4(perturbed2(:), h, lorenz)
!        perturbed3(:) = runge_kutta4(perturbed3(:), h, lorenz)

        lorenz_solution(:) = runge_kutta4(lorenz_solution(:), h, simplified_lorenz)
        perturbed1(:) = runge_kutta4(perturbed1(:), h, simplified_lorenz)
        perturbed2(:) = runge_kutta4(perturbed2(:), h, simplified_lorenz)
        perturbed3(:) = runge_kutta4(perturbed3(:), h, simplified_lorenz)

!        lorenz_solution(:) = runge_kutta4(lorenz_solution(:), h, general_polynomial)
!        perturbed1(:) = runge_kutta4(perturbed1(:), h, general_polynomial)
!        perturbed2(:) = runge_kutta4(perturbed2(:), h, general_polynomial)
!        perturbed3(:) = runge_kutta4(perturbed3(:), h, general_polynomial)

!        call save_data_point(uid, lorenz_solution)
        if (mod(i + 1, renorm) == 0) then
            A = perturbed1(:) - lorenz_solution(:)
            B = perturbed2(:) - lorenz_solution(:)
            C = perturbed3(:) - lorenz_solution(:)

            normA = norm(A)
            A = A/normA
            B = orthogonalize(A, B)
            normB = norm(B)
            B = B/normB
            C = orthogonalize(A, C)
            C = orthogonalize(B, C)
            normC = norm(C)
            C = C/normC

            perturbed1(:) = lorenz_solution(:) + A*perturbation_size
            perturbed2(:) = lorenz_solution(:) + B*perturbation_size
            perturbed3(:) = lorenz_solution(:) + C*perturbation_size
            exponents(:, (i + 1)/renorm) = [log(normA*invperturbation_size), log(normB*invperturbation_size), log(normC*invperturbation_size)]*invh
        end if
    end do
!    call save_solution_close(uid)
    exponents = exponents/real(renorm, rp)
    L1 = mean(exponents(1, :))
    L1err = stderr(exponents(1, :))
    L2 = mean(exponents(2, :))
    L2err = stderr(exponents(2, :))
    L3 = mean(exponents(3, :))
    L3err = stderr(exponents(3, :))
    call sort3(L1, L2, L3, L1err, L2err, L3err)
    call kaplan_yorke_dimension([L1, L2, L3], [L1err, L2err, L3err], dim, dimerr)
    print *, "L1:", L1, "+-", L1err
    print *, "L2:", L2, "+-", L2err
    print *, "L3:", L3, "+-", L3err
    
    print *, "KY dimension:", dim, "+-", dimerr

    print *, "trace:", L1 + L2 + L3, "+-", L1err + L2err + L3err

    contains
        pure function general_polynomial(xn)
            real(rp), intent(in) :: xn(3)
            real(rp) :: general_polynomial(3)
            real(rp) :: R

            R = 1.498_rp

            include "ode_polynomial.f90"
!            general_polynomial = 0.0_rp + xn(1) + R
            return
        end function general_polynomial

        pure function simplified_lorenz(xn)
            real(rp), intent(in) :: xn(3)
            real(rp) :: simplified_lorenz(3)
            
            real(rp) :: R

            R = 1.498_rp
            simplified_lorenz = [xn(2) - xn(1), -2.07_rp*xn(1)*xn(3), 1.41_rp*xn(1)*xn(2) - 1.42_rp*R]
            return
        end function simplified_lorenz

        pure function lorenz(xn)
            real(rp), intent(in) :: xn(3)
            real(rp) :: lorenz(3)
            
            real(rp) :: sigma, rho, beta

            sigma = 10.0_rp
            rho = 28.0_rp
            beta = 8.0_rp/3.0_rp
            lorenz = [sigma*(xn(2) - xn(1)), xn(1)*(rho - xn(3)) - xn(2), xn(1)*xn(2) - beta*xn(3)]
            return
        end function lorenz
    
        ! calculate and return x_{n+1}
        ! x_{n+1} = x_{n} + h*(k1 + 2*k2 + 2*k3 + k4)/6
        pure function runge_kutta4(xn, h, F)
            real(rp), intent(in) :: xn(3), h
            real(rp) :: runge_kutta4(3)
            
            real(rp) :: k1(3), k2(3), k3(3), k4(3)
            
            interface
                pure function F(vin)
                    use constants, only : rp
                    implicit none
                    real(rp), intent(in) :: vin(3)
                    real(rp) :: F(3)
                end function F
            end interface
            
            k1 = F(xn)
            k2 = F(xn + 0.5_rp*h*k1)
            k3 = F(xn + 0.5_rp*h*k2)
            k4 = F(xn + h*k3)
            runge_kutta4 = xn + h*(k1 + k4 + 2.0_rp*(k2 + k3))/6.0_rp
            return
        end function runge_kutta4
        
        pure function orthogonalize(V, A)
            real(rp), intent(in) :: V(3), A(3)
            real(rp) :: orthogonalize(3)
            
            orthogonalize = A - sum(A*V)*V
            return
        end function orthogonalize
        
        pure function norm(A)
            real(rp), intent(in) :: A(3)
            real(rp) :: norm
            norm = sqrt(sum(A*A))
            return
        end function norm

        function stderr(A)
            real(rp), intent(in) :: A(:)
            real(rp) :: stderr
            if (size(A) < 2) then
                print *, "You must have atleast two data pints to find the standard error."
                stop
            end if
            stderr = sqrt((mean(A*A) - mean(A)**2)/real(size(A) - 1, rp))
            return
        end function stderr

        pure function mean(A)
            real(rp), intent(in) :: A(:)
            real(rp) :: mean
            
            mean = sum(A)/max(1, size(A))
            return
        end function mean

        subroutine kaplan_yorke_dimension(L, Lerr, dim, dimerr)
            real(rp), intent(in) :: L(3), Lerr(3)
            real(rp), intent(out) :: dim, dimerr
            
            real(rp) :: temp
            integer :: i, j
            
            j = 0
            do i = 1, 3
                temp = sum(L(:i))
                if (temp > 0.0_rp) then
                    j = i
                end if
            end do
            
            if (j == 3) then
                print *, "Error: sum of all lyapunov exponents is positive"
                stop
            end if
             
            dim = real(j, rp) + sum(L(:j))/abs(L(j + 1))
            dimerr = sqrt(sum(Lerr(:j)**2/L(j + 1)**2) &
                           + (sum(L(:j))*Lerr(j + 1)/L(j + 1)**2)**2)
            return
        end subroutine kaplan_yorke_dimension
        
        subroutine sort3(x, y, z, xerr, yerr, zerr)
            real(rp), intent(inout) :: x, y, z, xerr, yerr, zerr
            
            call reorder(x, y, xerr, yerr)
            call reorder(y, z, yerr, zerr)
            call reorder(x, y, xerr, yerr)
            ! x >= y >= z
            return
        end subroutine sort3
        
        subroutine reorder(x1, y1, x2, y2)
            real(rp), intent(inout) :: x1, y1, x2, y2
            real(rp) :: temp1, temp2
            
            if (x1 >= y1) then
                return
            else
                temp1 = x1
                temp2 = x2
                
                x1 = y1
                x2 = y2
                
                y1 = temp1
                y2 = temp2
                return
            endif
        end subroutine reorder
        
        subroutine initalize_data(data, length)
            integer, intent(in) :: length
            real(rp), intent(out), allocatable :: data(:, :)

            integer :: alloc_error

            allocate(data(3, length), stat = alloc_error)
            if (alloc_error .ne. 0) then
                print *, "Error: data allocation failed!"
                stop
            end if
            data = 0.0_rp ! initalize to zero
            return
        end subroutine initalize_data
        
        subroutine load_settings(steps, tmax, perturbation, renorm, x_0,  y_0, z_0)
            integer, intent(out) :: steps, renorm
            real(rp), intent(out) :: tmax, perturbation, x_0, y_0, z_0

            integer :: uid
            
            open(newunit=uid, file="settings", status="old", action="read")
            read(uid, *) steps, tmax, perturbation, renorm
            read(uid, *) x_0, y_0, z_0
            close(uid)
            return
        end subroutine load_settings
        
        function save_solution_open(fname)
            character(len=*), intent(in) :: fname
            
            integer :: save_solution_open
            
            open(newunit=save_solution_open, file=fname, status="replace", action="write")
            return
        end function save_solution_open

        subroutine save_data_point(uid, data)
            real(rp), intent(in) :: data(3)
            integer, intent(in) :: uid
            
            write(uid, *) data
            return
        end subroutine save_data_point

        subroutine save_solution_close(uid)
            integer, intent(in) :: uid
            
            close(uid)
            return
        end subroutine save_solution_close
end program main
