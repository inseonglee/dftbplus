!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks_reksfon
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_message, only : error
  use dftbp_math_eigensolver, only : heev
  use dftbp_reks_rekscommon, only : getFactor
  use dftbp_reks_reksvar, only : TReksCalc, reksTypes

  implicit none

  private
  public :: optimizeFons

  !> Parameter to distinguish between two asymptotic regimes
  !> in the behavior of the derivative of the REKS coefficient, fx
  real(dp) :: Threshold = 5.0_dp

  !> Convergence for NR solver
  real(dp) :: ConvergeLimit = 1.0E-10_dp
!  real(dp) :: ConvergeLimit = 1.0E-8_dp

  contains

  !> Optimize the fractional occupation numbers (FONs) in REKS
  subroutine optimizeFons(this, tConverged)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Has the calculation converged?
    logical, intent(in) :: tConverged

    real(dp) :: enLtot(6)
    real(dp) :: x1, x2, tmp_hess
    logical :: optFON1, optFON2, optFON3

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)

      call getFONs22_(x1, this%hess, this%enLtot, this%delta, this%FonMaxIter, this%Plevel)
      ! FONs(1,1) = n_a, FONs(2,1) = n_b
      this%FONs(1,1) = 2.0_dp * x1
      this%FONs(2,1) = 2.0_dp - this%FONs(1,1)

    case (reksTypes%ssr44)

      if (.not. tConverged) then
        optFON1 = .true.
        optFON2 = .false.
        optFON3 = .false.
        if (this%Efunction == 3) then
          optFON2 = .true.
        else if (this%Efunction == 4 .or. this%Efunction == 5) then
          optFON2 = .true.
          optFON3 = .true.
        end if
      else
        optFON1 = .false.
        optFON2 = .false.
        optFON3 = .false.
        if (this%Efunction == 2 .or. this%Efunction == 6) then
          optFON2 = .true.
          optFON3 = .true.
        else if (this%Efunction == 3) then
          optFON3 = .true.
        end if
      end if

      if (optFON1) then

        ! For PPS state
        call getFONs44_(x1, x2, this%enLtot, this%delta, this%FonMaxIter, this%Plevel)
        ! FONs(1,1) = n_a, FONs(2,1) = n_b, FONs(3,1) = n_c, FONs(4,1) = n_d
        this%FONs(1,1) = 2.0_dp * x1
        this%FONs(2,1) = 2.0_dp * x2
        this%FONs(3,1) = 2.0_dp - this%FONs(2,1)
        this%FONs(4,1) = 2.0_dp - this%FONs(1,1)

      end if

      if (optFON2) then

        ! For OSS1 state
        enLtot(1) = this%enLtot(9); enLtot(2) = this%enLtot(13)
        enLtot(3) = this%enLtot(5); enLtot(4) = this%enLtot(6)
        enLtot(5) = this%enLtot(7); enLtot(6) = this%enLtot(8)
        ! FONs(1,2) = n'_a, FONs(4,2) = n'_d
        call getFONs22_(x1, tmp_hess, enLtot, this%delta, this%FonMaxIter, this%Plevel)
        this%FONs(1,2) = 2.0_dp * x1
        this%FONs(4,2) = 2.0_dp - this%FONs(1,2)

        ! For OSS2 state
        enLtot(1) = this%enLtot(5); enLtot(2) = this%enLtot(15)
        enLtot(3) = this%enLtot(9); enLtot(4) = this%enLtot(10)
        enLtot(5) = this%enLtot(11); enLtot(6) = this%enLtot(12)
        ! FONs(2,2) = n'_b, FONs(3,2) = n'_c
        call getFONs22_(x2, tmp_hess, enLtot, this%delta, this%FonMaxIter, this%Plevel)
        this%FONs(2,2) = 2.0_dp * x2
        this%FONs(3,2) = 2.0_dp - this%FONs(2,2)

      end if

      if (optFON3) then

        ! For OSS3 state
        enLtot(1) = this%enLtot(17); enLtot(2) = this%enLtot(25)
        enLtot(3) = this%enLtot(21); enLtot(4) = this%enLtot(22)
        enLtot(5) = this%enLtot(23); enLtot(6) = this%enLtot(24)
        ! FONs(1,3) = m'_a, FONs(3,3) = m'_c
        call getFONs22_(x1, tmp_hess, enLtot, this%delta, this%FonMaxIter, this%Plevel)
        this%FONs(1,3) = 2.0_dp * x1
        this%FONs(3,3) = 2.0_dp - this%FONs(1,3)

        ! For OSS4 state
        enLtot(1) = this%enLtot(21); enLtot(2) = this%enLtot(27)
        enLtot(3) = this%enLtot(17); enLtot(4) = this%enLtot(18)
        enLtot(5) = this%enLtot(19); enLtot(6) = this%enLtot(20)
        ! FONs(2,3) = m'_b, FONs(4,3) = m'_d
        call getFONs22_(x2, tmp_hess, enLtot, this%delta, this%FonMaxIter, this%Plevel)
        this%FONs(2,3) = 2.0_dp * x2
        this%FONs(4,3) = 2.0_dp - this%FONs(2,3)

      end if

    end select

  end subroutine optimizeFons


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Optimize FONs in REKS(2,2) case with Newton-Raphson method
  subroutine getFONs22_(x_a, hess0, enLtot, delta, maxIter, opt)

    !> converged x1 (= n_a/2)
    real(dp), intent(out) :: x_a

    !> converged Hessian of FONs
    real(dp), intent(out) :: hess0

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Maximum iteration used in FON optimization
    integer, intent(in) :: maxIter

    !> Print level in standard output file
    integer, intent(in) :: opt

    real(dp) :: Const, ConDeno
    real(dp) :: xUpLim, xDownLim
    real(dp) :: root, x0, x1, grad, hess
    real(dp) :: fac, eps
    integer :: iter

    ! Calculate Const in equation 12.c
    ! Reference : JCP, 147, 034113 (2017) and its supporting information
    ConDeno = enLtot(5) - enLtot(3)
    ! ConDeno should be negative
    ! In general, E2 - E1 > 0 due to the MO swap (always, n_a > n_b)
    ! So, negative ConDeno generates negative Const
    if (ConDeno > 0.0_dp) ConDeno = -ConDeno
    Const = (enLtot(2) - enLtot(1)) / ConDeno

    ! Set up the starting value for x, x0
    !
    ! if |Const| > threshold solve the equation
    ! abs[1-2*x] = (2*(2 + d)/(1 + d)) * (4*x*(1 - x))**(-d/(2*(1 + d)))
    ! In the equation abs[1-2*x] is actually abs[Const]
    !
    ! else solve the equation:
    ! Const = (-2)*(2*x - 1)*(4*x*(1 - x))**(-1/2)
    !
    if (abs(Const) > Threshold) then
      fac = (2.0_dp*(2.0_dp+delta)/((1.0_dp+delta)*abs(Const))) &
          & **(2.0_dp*(1.0_dp+delta)/delta)
      if (Const > 0.0_dp) then
        ! x ~= 0 case
        root = 0.5_dp * (1.0_dp - sqrt(1.0_dp - fac))
      else
        ! x ~= 1 case
        root = 0.5_dp * (1.0_dp + sqrt(1.0_dp - fac))
      end if
    else
      ! x ~= near 1/2 case
      root = 0.5_dp * (1.0_dp - Const/sqrt(4.0_dp+Const**2))
    end if

    xDownLim = 1.0E-14_dp
    xUpLim = 1.0_dp - xDownLim
    x0 = root
    eps = 0.0_dp

    NRsolver: do iter = 1, maxIter

      ! Update x1 value
      x1 = x0 + eps
      ! Move to the inside of limit ( 0 < x < 1 )
      if (x1 > xUpLim) x1 = xUpLim
      if (x1 < xDownLim) x1 = xDownLim
      ! Calculate gradient and hessian of f(x) using x = n_a/2
      grad = dfdx(x1, delta)
      hess = d2fdx2(x1, delta)
      ! Update eps value
      eps = (Const - grad) / hess
      if (opt >= 2) then
        write(stdOut,'(2x,a,1x,i4,4x,a,F18.14,1x,a,F18.14)') &
            & 'NR solver: Iteration', iter, 'X =', x1, 'Eps =', eps
      end if

      ! Convergence check
      if (abs(eps) > ConvergeLimit) then
        x0 = x1
        if (opt >= 2 .and. iter == maxIter) then
          write(stdOut,'(2x,a,i4,a)') &
              & 'Warning! Maximum number of iterations (', maxIter, &
              & ') is exceeded in NR solver'
        end if
      else
        if (opt >= 2) then
          write(stdOut,'(2x,a,1x,i4,1x,a)') &
              & 'Convergence reached in NR solver after', iter, 'iterations'
        end if
        exit NRsolver
      end if

    end do NRsolver

    ! Converged x and hessian value
    x_a = x1
    hess0 = hess

  end subroutine getFONs22_


  !> Optimize FONs (n_a and n_b) in REKS(4,4) case with Newton-Raphson method
  subroutine getFONs44_(x_a, x_b, enLtot, delta, maxIter, opt)

    !> converged x1 (= n_a/2) and x2 (= n_b/2)
    real(dp), intent(out) :: x_a, x_b

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Maximum iteration used in FON optimization
    integer, intent(in) :: maxIter

    !> Print level in standard output file
    integer, intent(in) :: opt

    real(dp) :: eps(2), x0(2), x1(2)
    real(dp) :: dE(2), grad(2), hess(2,2), hess_inv(2,2)
    real(dp) :: H(4,4), eval(4)
    real(dp) :: xUpLim, xDownLim
    real(dp) :: norm, det, factor, golden, e0, e1
    integer :: ii, iter, miter

    xDownLim = 1.0E-14_dp
    xUpLim = 1.0_dp - xDownLim

    ! K_ad = energy difference with ad transition
    dE(1) = enLtot(7) - enLtot(5)
    ! K_bc = energy difference with bc transition
    dE(2) = enLtot(11) - enLtot(9)

    ! H = | E_aabb  K_bc   K_ad     0    |
    !     |  K_bc  E_aacc   0      K_ad  |
    !     |  K_ad    0    E_bbdd   K_bc  |
    !     |   0     K_ad   K_bc   E_ccdd |
    ! Obtain initial guess for FONs from diagonalization of H matrix
    H(:,:) = 0.0_dp
    ! energy for four microstates
    do ii = 1, 4
      H(ii,ii) = enLtot(ii)
    end do
    ! K_ad
    H(1,3) = -dE(1); H(3,1) = -dE(1)
    H(2,4) = -dE(1); H(4,2) = -dE(1)
    ! K_bc
    H(1,2) = -dE(2); H(2,1) = -dE(2)
    H(3,4) = -dE(2); H(4,3) = -dE(2)

    eval(:) = 0.0_dp
    call heev(H, eval, 'U', 'V')

    ! Now H is eigenvectors of the matrix; first n_a/2 and n_b/2
    x1(1) = H(1,1)**2 + H(1,2)**2
    x1(2) = H(1,1)**2 + H(1,3)**2

    ! TODO : If close to one, set it a bit apart
!    if (dabs(1.0_dp - x1(1)) <= 0.1_dp) x1(1) = 0.9_dp
!    if (dabs(1.0_dp - x1(2)) <= 0.1_dp) x1(2) = 0.9_dp

    if (dabs(dE(1)) < 1.0E-10_dp) x1(1) = xUpLim
    if (dabs(dE(2)) < 1.0E-10_dp) x1(2) = xUpLim

    e0 = getEpps(x1(1), x1(2), enLtot, delta)

    x0(:) = x1(:)
    eps(:) = 0.0_dp
    NRsolver: do iter = 1, maxIter

      ! Calculate gradient of the energy of PPS state
      grad(1) = x1(2) * enLtot(1) + (1.0_dp - x1(2)) * enLtot(2) - &
          & x1(2) * enLtot(3) - (1.0_dp - x1(2)) * enLtot(4) + dE(1) * dfdx(x1(1), delta)
      grad(2) = x1(1) * enLtot(1) + (1.0_dp - x1(1)) * enLtot(3) - &
          & x1(1) * enLtot(2) - (1.0_dp - x1(2)) * enLtot(4) + dE(2) * dfdx(x1(2), delta)

      ! Calculate hessian of the energy of PPS state
      hess(1,1) = dE(1) * d2fdx2(x1(1), delta)
      hess(2,2) = dE(2) * d2fdx2(x1(2), delta)
      hess(1,2) = enLtot(1) - enLtot(2) - enLtot(3) + enLtot(4)
      hess(2,1) = hess(1,2)
      ! Get inverse hessian
      hess_inv(1,1) = hess(2,2)
      hess_inv(2,2) = hess(1,1)
      hess_inv(1,2) = -hess(1,2)
      hess_inv(2,1) = hess_inv(1,2)

      ! Calculate the increment for x1 and x2
      ! eps = inverse hessian matrix dot gradient vector
      if (dabs(dE(1)) < 1.0E-8_dp .or. dabs(dE(2)) < 1.0E-8_dp) then
        ! Check if one of the occupations is fixed
        if (opt >= 2) then
          write(stdOut,'(2x,a,1x,i4,1x,a)') &
              & 'abs(K_ad) or abs(K_bc) is too small, one of the occupations is fixed!'
        end if
        eps(:) = 0.0_dp
        if (dabs(hess_inv(2,2)) > 1.0E-10_dp) then
          eps(1) = -grad(1) / hess_inv(2,2)
        end if
        if (dabs(hess_inv(1,1)) > 1.0E-10_dp) then
          eps(2) = -grad(2) / hess_inv(1,1)
        end if
        ! Inverse the sign if K_ad or K_bc is positive, it must be negative
        if (dE(1) > 0.0) then
          eps(1) = -eps(1)
        end if
        if (dE(2) > 0.0) then
          eps(2) = -eps(2)
        end if
      else
        ! General case
        det = hess(1,1) * hess(2,2) - hess(1,2)**2
        eps(1) = -grad(1) * hess_inv(1,1) - grad(2) * hess_inv(1,2)
        eps(2) = -grad(1) * hess_inv(2,1) - grad(2) * hess_inv(2,2)
        if (det > 1.0E-10_dp) then
          eps(1) = eps(1) / det
          eps(2) = eps(2) / det
        end if
      end if 
      norm = dsqrt(eps(1)**2 + eps(2)**2)

      factor = 1.0_dp
      golden = (1.0_dp + dsqrt(5.0_dp)) / 2.0_dp
      line_search: do miter = 1, maxIter

        ! Update x1 and x2
        x1(:) = x0(:) + eps(:) * factor

        ! Move to the inside of limit ( 0 < x < 1 )
        if (x1(1) > xUpLim) x1(1) = xUpLim
        if (x1(1) < xDownLim) x1(1) = xDownLim
        if (x1(2) > xUpLim) x1(2) = xUpLim
        if (x1(2) < xDownLim) x1(2) = xDownLim

        e1 = getEpps(x1(1), x1(2), enLtot, delta) 

        if (e1 < e0) then
          exit line_search
        end if
        if (dabs(factor*norm) < 1.0E-10_dp) then
          exit line_search
        end if

        ! Update the factor
        factor = factor / (golden * golden)
        ! Speed up the micro iteration
        if (miter > maxIter / 2) then
          factor = factor / (golden * golden)
        end if

      end do line_search

      e0 = e1
      eps(:) = x1(:) - x0(:)
      norm = dsqrt(eps(1)**2 + eps(2)**2)

      if (opt >= 2) then
        write(stdOut,'(2x,a,1x,i4,4x,a,F18.14,4x,a,F18.14,1x,a,F18.14,1x,a,F15.8)') &
            & 'NR solver: Iteration', iter, 'x1 =', x1(1), 'x2 =', x1(2), &
            & 'Eps =', norm, 'En =', e1
      end if

      ! Convergence check
      if (norm > ConvergeLimit) then
        x0(:) = x1(:)
        if (opt >= 2 .and. iter == maxIter) then
          write(stdOut,'(2x,a,i4,a)') &
              & 'Warning! Maximum number of iterations (', maxIter, &
              & ') is exceeded in NR solver of {n}'
        end if
      else
        if (opt >= 2) then
          write(stdOut,'(2x,a,1x,i4,1x,a)') &
              & 'Convergence reached in NR solver after', iter, 'iterations'
        end if
        exit NRsolver
      end if

    end do NRsolver

    ! Converged x1 and x2
    x_a = x1(1)
    x_b = x1(2)

  end subroutine getFONs44_


  !> Calculate gradient of f(x)
  function dfdx(x, delta) result(grad)

    !> x (= n/2)
    real(dp), intent(in) :: x

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    real(dp) :: y, fac1, fac2, fac3
    real(dp) :: grad

    ! Decide y = 4*x*(1-x) where x = n/2
    y = 4.0_dp * x * (1.0_dp - x)
    ! TODO : Is this setting helpful for convergence?
!    if (x * (1.0_dp-x) < 1.0E-10_dp) y = 4.0E-10_dp

    fac1 = 2.0_dp * (1.0_dp - 2.0_dp * x) / (1.0_dp + delta)
    fac2 = -0.5_dp * (y + delta) / (1.0_dp + delta)
    fac3 = 2.0_dp + delta - y - y * dlog(y)

    grad = fac1 * fac3 * y**fac2

  end function dfdx


  !> Calculate hessian of f(x)
  function d2fdx2(x, delta) result (hess)

    !> x (= n/2)
    real(dp), intent(in) :: x

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    real(dp) :: y, fac1, fac2, fac3, fac4, fac5, fac6
    real(dp) :: hess

    ! Decide y = 4*x*(1-x) where x = n/2
    y = 4.0_dp * x * (1.0_dp - x)
    ! TODO : Is this setting helpful for convergence?
!    if (x * (1.0_dp-x) < 1.0E-10_dp) y = 4.0E-10_dp

    fac1 = 4.0_dp / (1.0_dp + delta)**2
    fac2 = -1.0_dp * (2.0_dp + 3.0_dp * delta + y) / (2.0_dp + 2.0_dp * delta)
    fac4 = -1.0_dp * delta**2 - y * ( 8.0_dp + (-8.0_dp + y) * y)
    fac5 = delta * ( -2.0_dp + 5.0_dp * (-1.0_dp + y) * y)
    fac6 = 4.0_dp + delta * (2.0_dp-3.0_dp * y) + y * (-7.0_dp + 2.0_dp * y) + &
        & (-1.0_dp + y) * y * dlog(y)
    fac3 = fac4 + fac5 - y * dlog(y) * fac6

    hess = fac1 * fac3 * y**fac2

  end function d2fdx2


  !> Calculate the energy of PPS state using current FONs
  function getEpps(x1, x2, enLtot, delta) result (en)

    !> converged x1 (= n_a/2) and x2 (= n_b/2)
    real(dp), intent(in) :: x1, x2

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    real(dp) :: n_a, n_b, n_c, n_d, fac1, fac2
    real(dp) :: en

    n_a = 2.0_dp * x1; n_b = 2.0_dp * x2
    n_c = 2.0_dp - n_b; n_d = 2.0_dp - n_a

    fac1 = getFactor(n_a, n_d, delta)
    fac2 = getFactor(n_b, n_c, delta)

    en = 0.25_dp * (n_a * n_b * enLtot(1) + n_a * n_c * enLtot(2) &
        & + n_b * n_d * enLtot(3) + n_c * n_d * enLtot(4)) &
        & + fac1 * (enLtot(5) + enLtot(6) - enLtot(7) - enLtot(8)) &
        & + fac2 * (enLtot(9) + enLtot(10) - enLtot(11) - enLtot(12))

  end function getEpps


end module dftbp_reks_reksfon
