!    enter the ode equation to be used. I must be saved as general_polynomial(3)
!    and defined interms of xn(3) where x = xn(1), y = xn(2) and z = xn(3)
!    example: simplified lorenz polynomial
!    general_polynomial = [xn(2) - xn(1), -xn(1)*xn(3), xn(1)*xn(2) - 1.0_rp]
!

    general_polynomial = [-.5_rp -10*xn(1) + 11*xn(1) - 0.1_rp*xn(3) + &
    0.17_rp*xn(1)**2-0.11_rp*xn(1)*xn(2), -xn(1)*xn(3), xn(1)*xn(2) - 1.0_rp]
