!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'error.fypp'

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks_reksproperty
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_reks_rekscommon, only : getTwoIndices, qm2udL, ud2qmL, matAO2MO,&
      & matMO2AO, assignFilling, assignIndex
  use dftbp_reks_reksio, only : printRelaxedFONs, printRelaxedFONsL, printUnrelaxedFONs
  use dftbp_reks_reksvar, only : reksTypes

  implicit none

  private
  public :: getUnrelaxedDensMatAndTdp, getRelaxedDensMat, getRelaxedDensMatL
  public :: getTdpParameters, buildTdpVectors, getSsrCoefDeriv, TDPshift
  public :: getDipoleIntegral, getDipoleMomentMatrix, getReksOsc

  contains

  !> Calculate unrelaxed density and transition density for target
  !> SA-REKS or SSR state (or L-th state)
  subroutine getUnrelaxedDensMatAndTdp(eigenvecs, overSqr, rhoL, FONs, &
      & eigvecsSSR, Lpaired, Nc, Na, rstate, Lstate, reksAlg, tSSR, tTDP, &
      & unrelRhoSqr, unrelTdm, densityMatrix, errStatus)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Dense density matrix for target microstate
    real(dp), intent(in) :: rhoL(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Calculate transition dipole moments
    logical, intent(in) :: tTDP

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(out) :: unrelRhoSqr(:,:)

    !> unrelaxed transition density matrix between SSR or SA-REKS states
    real(dp), allocatable, intent(inout) :: unrelTdm(:,:,:)

    !> Holds density matrix settings and pointers
    type(TDensityMatrix), intent(inout) :: densityMatrix

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: rhoX(:,:,:)
    real(dp), allocatable :: rhoXdel(:,:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpFilling(:,:)

    integer :: nOrb, nstates, nstHalf
    integer :: ii, ist, jst, kst, lst, ia, ib

    nOrb = size(eigenvecs,dim=1)
    nstates = size(eigvecsSSR,dim=1)
    nstHalf = nstates * (nstates - 1) / 2

    if (tSSR) then
      allocate(rhoX(nOrb,nOrb,nstates))
    else
      allocate(rhoX(nOrb,nOrb,1))
    end if
    if (Lstate == 0) then
      allocate(rhoXdel(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpFilling(nOrb,nstates))

    ! Core orbitals fillings
    tmpFilling(:,:) = 0.0_dp
    do ii = 1, Nc
      tmpFilling(ii,:) = 2.0_dp
    end do
    ! Active orbitals fillings
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getActiveFilling22_(FONs, Nc, tmpFilling)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Unrelaxed density matrix for SA-REKS or L-th state
    rhoX(:,:,:) = 0.0_dp
    if (tSSR) then
      do ist = 1, nstates
        call densityMatrix%getDensityMatrix(rhoX(:,:,ist), eigenvecs, tmpFilling(:,ist), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call adjointLowerTriangle(rhoX(:,:,ist))
      end do
    else
      if (Lstate == 0) then
        call densityMatrix%getDensityMatrix(rhoX(:,:,1), eigenvecs, tmpFilling(:,rstate), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call adjointLowerTriangle(rhoX(:,:,1))
      else
        ! find proper index for up+down in rhoSqrL
        rhoX(:,:,1) = rhoL
      end if
    end if

    ! Unrelaxed transition density matrix between SA-REKS states
    if (Lstate == 0) then
      rhoXdel(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getUnrelaxedTDM22_(eigenvecs, FONs, Nc, nstates, rhoXdel)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! Final unrelaxed density matrix for target state
    if (tSSR) then
      ! unrelRhoSqr is unrelaxed density matrix for target SSR state
      kst = 0
      unrelRhoSqr(:,:) = 0.0_dp
      do ist = 1, nstates
        do jst = ist, nstates
          if (ist == jst) then
            unrelRhoSqr(:,:) = unrelRhoSqr + eigvecsSSR(ist,rstate)**2 * rhoX(:,:,ist)
          else
            kst = kst + 1
            unrelRhoSqr(:,:) = unrelRhoSqr + 2.0_dp * eigvecsSSR(ist,rstate) * &
                & eigvecsSSR(jst,rstate) * rhoXdel(:,:,kst)
          end if
        end do
      end do
    else
      ! unrelRhoSqr is unrelaxed density matrix for target SA-REKS or L-th state
      unrelRhoSqr(:,:) = rhoX(:,:,1)
    end if

    ! Final unrelaxed transition density matrix between states
    if (tTDP .and. Lstate == 0) then
      if (tSSR) then
        ! unrelTdm is unrelaxed transition density matrix between SSR states
        unrelTdm(:,:,:) = 0.0_dp
        do lst = 1, nstHalf

          call getTwoIndices(nstates, lst, ia, ib, 1)

          kst = 0
          do ist = 1, nstates
            do jst = ist, nstates
              if (ist == jst) then
                unrelTdm(:,:,lst) = unrelTdm(:,:,lst) + rhoX(:,:,ist) &
                    & * eigvecsSSR(ist,ia) * eigvecsSSR(ist,ib)
              else
                kst = kst + 1
                ! <PPS|OSS> = <OSS|PPS>, etc
                unrelTdm(:,:,lst) = unrelTdm(:,:,lst) + rhoXdel(:,:,kst) &
                    & * ( eigvecsSSR(ist,ia) * eigvecsSSR(jst,ib) &
                    & + eigvecsSSR(jst,ia) * eigvecsSSR(ist,ib) )
              end if
            end do
          end do

        end do
      else
        ! unrelTdm is unrelaxed transition density matrix between SA-REKS states
        unrelTdm(:,:,:) = rhoXdel
      end if
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, unrelRhoSqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printUnrelaxedFONs(tmpRho, rstate, Lstate, Nc, Na, tSSR)

  end subroutine getUnrelaxedDensMatAndTdp


  !> Calculate relaxed density for target SA-REKS or SSR state
  subroutine getRelaxedDensMat(eigenvecs, overSqr, unrelRhoSqr, ZT, omega, &
      & FONs, eigvecsSSR, SAweight, Rab, G1, Nc, Na, rstate, reksAlg, &
      & tSSR, tNAC, relRhoSqr)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(in) :: unrelRhoSqr(:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> state-interaction term used in SSR gradients
    real(dp), allocatable, intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> relaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(out) :: relRhoSqr(:,:)

    real(dp), allocatable :: resRho(:,:)
    real(dp), allocatable :: resTdm(:,:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    real(dp) :: fp, fq
    integer :: Nv, superN, nOrb, nstates, nstHalf
    integer :: ii, ist, jst, kst, pq, p, q

    superN = size(ZT,dim=1)
    nOrb = size(eigenvecs,dim=1)
    Nv = nOrb - Nc - Na
    nstates = size(eigvecsSSR,dim=1)
    nstHalf = nstates * (nstates - 1) / 2

    allocate(resRho(nOrb,nOrb))
    if (tSSR) then
      allocate(resTdm(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! a part of transition density matrix originating from the
    ! response of the orbital occupation numbers
    if (tSSR) then
      resTdm(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getResponseTDM22_(eigenvecs, FONs, SAweight, Rab, &
            & G1, Nc, resTdm)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! a part of relaxed density matrix originating from the
    ! response of the orbital occupation numbers with XT
    resRho(:,:) = 0.0_dp

    if (tNAC) then
      ii = rstate
    else
      ii = 1
    end if

    tmpRho(:,:) = 0.0_dp
    do pq = 1, superN
      ! assign index p and q from pq
      call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)
      ! assign average filling for pth orbital
      call assignFilling(FONs, SAweight, Nc, p, reksAlg, fp)
      ! assign average filling for qth orbital
      call assignFilling(FONs, SAweight, Nc, q, reksAlg, fq)
      tmpRho(p,q) = (fp - fq) * ZT(pq,ii)
    end do
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs, transB='T')
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getResponseDM22_(eigenvecs, ZT(:,ii), tmpRho, omega, &
          & SAweight, G1, Nc, resRho)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Final relaxed density matrix for target state
    relRhoSqr(:,:) = 0.0_dp
    if (tSSR) then
      kst = 0
      ! relRhoSqr is relaxed density matrix for target SSR state
      relRhoSqr(:,:) = unrelRhoSqr - 2.0_dp * resRho
      do ist = 1, nstates
        do jst = ist + 1, nstates
          kst = kst + 1
          relRhoSqr(:,:) = relRhoSqr + 2.0_dp * eigvecsSSR(ist,rstate) * &
              & eigvecsSSR(jst,rstate) * resTdm(:,:,kst)
        end do
      end do
    else
      ! relRhoSqr is relaxed density matrix for target SA-REKS state
      relRhoSqr(:,:) = unrelRhoSqr - 2.0_dp * resRho
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, relRhosqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printRelaxedFONs(tmpRho, rstate, Nc, Na, tSSR)

  end subroutine getRelaxedDensMat


  !> Calculate relaxed density for target L-th microstate
  subroutine getRelaxedDensMatL(eigenvecs, rhoSqrL, overSqr, weight, &
      & SAweight, unrelRhoSqr, RmatL, ZT, omega, weightIL, G1, orderRmatL, &
      & Lpaired, Nc, Na, Lstate, reksAlg, relRhoSqr)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(inout) :: rhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> unrelaxed density matrix for target L-th state
    real(dp), intent(in) :: unrelRhoSqr(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> relaxed density matrix for target L-th state
    real(dp), intent(out) :: relRhoSqr(:,:)

    real(dp), allocatable :: resRhoL(:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: nOrb, tmpL, iL, Lmax

    nOrb = size(rhoSqrL,dim=1)
    Lmax = size(rhoSqrL,dim=4)

    allocate(resRhoL(nOrb,nOrb))
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! rhoSqrL has (my_ud) component
    call qm2udL(rhoSqrL, Lpaired)

    ! resRhoL is response part for L-th state
    resRhoL(:,:) = 0.0_dp

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getResponseDML22_(rhoSqrL, SAweight, ZT, omega, &
          & weightIL, G1, Lpaired, resRhoL)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! rhoSqrL has (my_qm) component
    call ud2qmL(rhoSqrL, Lpaired)

    do iL = 1, Lmax
      ! find proper index for RmatL
      tmpL = orderRmatL(iL)
      resRhoL(:,:) = resRhoL - 2.0_dp * RmatL(:,:,tmpL,1) * weight(iL)
    end do

    ! relRhoSqr is relaxed density matrix for L-th state
    relRhoSqr(:,:) = unrelRhoSqr + resRhoL

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, relRhosqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printRelaxedFONsL(tmpRho, Lstate, Nc, Na)

  end subroutine getRelaxedDensMatL


  !> Compute necessary parameters for gradients of transidion dipole
  subroutine getTdpParameters(eigenvecs, eigvecsSSR, rhoSqrL, FONs, weightIL,&
      & Nc, Lstate, reksAlg, Lpaired, tSSR, tranOcc, preTdp, unrelDp, unrelTdp,&
      & dipoleInt, densityMatrix, errStatus)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(in) :: rhoSqrL(:,:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> transition occupation matrix for transition dipole between SSR or SA-REKS states
    real(dp), intent(out) :: tranOcc(:,:,:)

    !> prefactor before the derivatives of FONs with dipole integral
    real(dp), intent(out) :: preTdp(:,:)

    !> unrelaxed dipole moment for SA-REKS state
    real(dp), allocatable, intent(inout) :: unrelDp(:,:)

    !> unrelaxed transition dipole moment between SA-REKS states
    real(dp), allocatable, intent(inout) :: unrelTdp(:,:)

    !> dipole integral in DFTB formalism
    real(dp), intent(inout) :: dipoleInt(:,:,:)

    !> Holds density matrix settings and pointers
    type(TDensityMatrix), intent(inout) :: densityMatrix

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: tmpFilling(:,:)
    real(dp), allocatable :: tmpMat(:,:,:)
    real(dp), allocatable :: matX(:,:,:)
    real(dp), allocatable :: matDel(:,:,:)

    integer :: nOrb, nstates, nstHalf
    integer :: ii, ist, jst, kst, lst, iOrb, ia, ib

    nOrb = size(eigenvecs,dim=1)
    nstates = size(eigvecsSSR,dim=1)
    nstHalf = nstates * (nstates - 1) / 2

    allocate(tmpFilling(nOrb,nstates))
    allocate(tmpMat(nOrb,nOrb,nstHalf))
    if (tSSR) then
      allocate(matX(nOrb,nOrb,nstates))
    end if
    if (Lstate == 0) then
      allocate(matDel(nOrb,nOrb,nstHalf))
    end if

    ! Common variables

    ! Core orbitals fillings
    tmpFilling(:,:) = 0.0_dp
    do iOrb = 1, Nc
      tmpFilling(iOrb,:) = 2.0_dp
    end do
    ! Active orbitals fillings
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getActiveFilling22_(FONs, Nc, tmpFilling)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Parameters for first contribution to TDP gradient

    ! Convert average fillings for SA-REKS state to matrix
    if (tSSR) then
      matX(:,:,:) = 0.0_dp
      do ist = 1, nstates
        do iOrb = 1, nOrb
          matX(iOrb,iOrb,ist) = tmpFilling(iOrb,ist)
        end do
      end do
    end if

    ! Calculate transition occupation matrix between SA-REKS states
    if (Lstate == 0) then
      matDel(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getTranOccMat22_(FONs, Nc, nstates, matDel)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! Calculate final average occupations for transition dipole
    if (tSSR) then
      ! tranOcc is transition occupation matrix between SSR states
      tranOcc(:,:,:) = 0.0_dp
      do lst = 1, nstHalf

        call getTwoIndices(nstates, lst, ia, ib, 1)

        kst = 0
        do ist = 1, nstates
          do jst = ist, nstates
            if (ist == jst) then
              tranOcc(:,:,lst) = tranOcc(:,:,lst) + matX(:,:,ist) &
                  & * eigvecsSSR(ist,ia) * eigvecsSSR(ist,ib)
            else
              kst = kst + 1
              ! <PPS|OSS> = <OSS|PPS>, etc
              tranOcc(:,:,lst) = tranOcc(:,:,lst) + matDel(:,:,kst) &
                  & * ( eigvecsSSR(ist,ia) * eigvecsSSR(jst,ib) &
                  & + eigvecsSSR(jst,ia) * eigvecsSSR(ist,ib) )
            end if
          end do
        end do

      end do
    else
      ! tranOcc is transition occupation matrix between SA-REKS states
      tranOcc(:,:,:) = matDel
    end if

    ! Parameters for second contribution to TDP gradient

    ! Density matrix contribution for SA-REKS state
    if (tSSR) then
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getDensFon22_(rhoSqrL, weightIL, Lpaired, matX)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! Density matrix contribution for transitions between SA-REKS states
    if (Lstate == 0) then
      matDel(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getTranFon22_(eigenvecs, FONs, Nc, nstates, matDel)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! Calculate final density matrix contributions
    if (tSSR) then
      tmpMat(:,:,:) = 0.0_dp
      do lst = 1, nstHalf

        call getTwoIndices(nstates, lst, ia, ib, 1)

        kst = 0
        do ist = 1, nstates
          do jst = ist, nstates
            if (ist == jst) then
              tmpMat(:,:,lst) = tmpMat(:,:,lst) + matX(:,:,ist) &
                  & * eigvecsSSR(ist,ia) * eigvecsSSR(ist,ib)
            else
              kst = kst + 1
              ! <PPS|OSS> = <OSS|PPS>, etc
              tmpMat(:,:,lst) = tmpMat(:,:,lst) + matDel(:,:,kst) &
                  & * ( eigvecsSSR(ist,ia) * eigvecsSSR(jst,ib) &
                  & + eigvecsSSR(jst,ia) * eigvecsSSR(ist,ib) )
            end if
          end do
        end do

      end do
    else
      tmpMat(:,:,:) = matDel
    end if

    ! Prefactors from final density matrices with dipole integral
    preTdp(:,:) = 0.0_dp
    do lst = 1, nstHalf
      do ii = 1, 3
        preTdp(ii,lst) = sum(tmpMat(:,:,lst)*dipoleInt(:,:,ii))
      end do
    end do

    ! Parameters for third contribution to TDP gradient

    ! Calculate dipole integral multiplied by SA-REKS or SI density matrix
    if (tSSR) then

      ! Unrelaxed density matrix for SA-REKS state
      matX(:,:,:) = 0.0_dp
      do ist = 1, nstates
        call densityMatrix%getDensityMatrix(matX(:,:,ist), eigenvecs(:,:,1),&
            & tmpFilling(:,ist), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call adjointLowerTriangle(matX(:,:,ist))
      end do

      ! Unrelaxed transition density matrix between SA-REKS states
      if (Lstate == 0) then
        matDel(:,:,:) = 0.0_dp
        select case (reksAlg)
        case (reksTypes%noReks)
        case (reksTypes%ssr22)
          call getUnrelaxedTDM22_(eigenvecs(:,:,1), FONs, Nc, nstates, matDel)
        case (reksTypes%ssr44)
          call error("SSR(4,4) is not implemented yet")
        end select
      end if

      ! Calculate unrelaxed dipole moment
      do ist = 1, nstates
        call getDipoleMomentMatrix(matX(:,:,ist), dipoleInt, unrelDp(:,ist))
      end do

      ! Calculate unrelaxed transition dipole moment
      do lst = 1, nstHalf
        call getDipoleMomentMatrix(matDel(:,:,lst), dipoleInt, unrelTdp(:,lst))
      end do

    end if

    ! Transform dipole moment integrals in MO basis
    do ii = 1, 3
      call matAO2MO(dipoleInt(:,:,ii), eigenvecs(:,:,1))
    end do

  end subroutine getTdpParameters


  !> Calculate X^T vector for transition dipoles between SA-REKS or SSR states
  subroutine buildTdpVectors(dipoleInt, tranOcc, preTdp, omega, G1, Nc, Na, &
      & reksAlg, XTtdp, symTdpVec)

    !> dipole integral in DFTB formalism
    real(dp), intent(in) :: dipoleInt(:,:,:)

    !> transition occupation matrix for transition dipole between SSR or SA-REKS states
    real(dp), intent(in) :: tranOcc(:,:,:)

    !> prefactor before the derivatives of FONs with dipole integral
    real(dp), intent(in) :: preTdp(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> transition dipole vector used for TDP gradient
    real(dp), intent(inout) :: XTtdp(:,:,:)

    !> symmetric part of transition dipole vector for gradients of TDP
    real(dp), intent(inout) :: symTdpVec(:,:,:,:)

    real(dp), allocatable :: tmpMat1(:,:)
    real(dp), allocatable :: tmpMat2(:,:)

    real(dp) :: tmp1, tmp2
    integer :: nOrb, superN, nstHalf, Nv
    integer :: ii, ij, i, j, ist

    nOrb = size(tranOcc,dim=1)
    superN = size(XTtdp,dim=1)
    nstHalf = size(XTtdp,dim=3)
    Nv = nOrb - Nc - Na

    allocate(tmpMat1(nOrb,nOrb))
    allocate(tmpMat2(nOrb,nOrb))

    XTtdp(:,:,:) = 0.0_dp
    do ist = 1, nstHalf
      do ii = 1, 3

        call gemm(tmpMat1, dipoleInt(:,:,ii), tranOcc(:,:,ist))
        call gemm(tmpMat2, tranOcc(:,:,ist), dipoleInt(:,:,ii))

        do ij = 1, superN
          ! assign index i and j from ij
          call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
          tmp1 = 0.5_dp * (tmpMat1(i,j) - tmpMat1(j,i))
          tmp2 = 0.5_dp * (tmpMat2(i,j) - tmpMat2(j,i))
          XTtdp(ij,ii,ist) = -(tmp1 - tmp2)
        end do

        symTdpVec(:,:,ii,ist) = 0.25_dp * (tmpMat1 + transpose(tmpMat1) &
            & + tmpMat2 + transpose(tmpMat2))

      end do
    end do

    do ist = 1, nstHalf
      do ii = 1, 3
        do ij = 1, superN
          XTtdp(ij,ii,ist) = XTtdp(ij,ii,ist) - G1*omega(ij)*preTdp(ii,ist)
        end do
      end do
    end do

  end subroutine buildTdpVectors


  ! TODO : Currently, only 2 state problem can be solved analytically
  subroutine getSsrCoefDeriv(energy, eigvecsSSR, SAgrad, SIgrad, eigvecsSSRderiv)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> gradient of SA-REKS state
    real(dp), intent(in) :: SAgrad(:,:,:)

    !> gradient of state-interaction term
    real(dp), intent(in) :: SIgrad(:,:,:)

    !> derivatives of eigenvectors from SA-REKS state
    real(dp), intent(out) :: eigvecsSSRderiv(:,:,:,:)

    real(dp), allocatable :: tmpEn(:,:)
    real(dp), allocatable :: tmpHam(:,:)
    real(dp), allocatable :: gVec22(:,:)
    real(dp), allocatable :: hVec22(:,:)
    real(dp), allocatable :: phaseDeriv(:,:)
    real(dp), allocatable :: tmpDeriv(:,:)

    real(dp) :: tmpValue1, tmpValue2
    integer :: ist, nstates, nAtom

    nstates = size(eigvecsSSR,dim=1)
    nAtom = size(SAgrad,dim=2)

    allocate(tmpEn(nstates,nstates))
    allocate(tmpHam(nstates,nstates))
    if (nstates == 2) then
      allocate(gVec22(3,nAtom))
      allocate(hVec22(3,nAtom))
      allocate(phaseDeriv(3,nAtom))
      allocate(tmpDeriv(3,nAtom))
    end if

    ! Obtain diabatic Hamiltonian matrix from energy and eigvecsSSR; H * U = U * E
    tmpEn(:,:) = 0.0_dp
    do ist = 1, nstates
      tmpEn(ist,ist) = energy(ist)
    end do
    tmpHam(:,:) = matmul(eigvecsSSR,matmul(tmpEn,transpose(eigvecsSSR)))

    ! Initialize the derivatives for SSR coefficients
    eigvecsSSRderiv(:,:,:,:) = 0.0_dp

    if (nstates == 2) then

      ! For 2 state SSR case, determinant becomes zero due to symmetry
      ! Thus, we directly calculate the root of quadratic equation
      ! derived from the determinant equation

      gVec22(:,:) = 0.0_dp
      gVec22(:,:) = SAgrad(:,:,1) - SAgrad(:,:,2)

      hVec22(:,:) = 0.0_dp
      hVec22(:,:) = SIgrad(:,:,1)

      phaseDeriv(:,:) = 0.0_dp

      ! first contribution
      tmpValue1 = 0.5_dp / tmpHam(1,2)
      tmpValue2 = 0.5_dp * (tmpHam(1,1) - tmpHam(2,2)) / tmpHam(1,2)**2

      tmpDeriv(:,:) = 0.0_dp
      tmpDeriv(:,:) = (gVec22 * tmpValue1 - tmpValue2 * hVec22)

      phaseDeriv(:,:) = phaseDeriv + eigvecsSSR(1,1)**2 * tmpDeriv

      ! second contribution - 1
      tmpValue1 = 2.0_dp * (tmpHam(1,1) - tmpHam(2,2))
      tmpValue2 = 8.0_dp * tmpHam(1,2)

      tmpDeriv(:,:) = 0.0_dp
      tmpDeriv(:,:) = tmpValue1 * gVec22 + tmpValue2 * hVec22

      tmpValue1 = 0.25_dp / tmpHam(1,2)
      tmpValue2 = 1.0_dp / sqrt( (tmpHam(1,1) - tmpHam(2,2))**2 + 4.0_dp * tmpHam(1,2)**2 )

      phaseDeriv(:,:) = phaseDeriv + eigvecsSSR(1,1)**2 * tmpValue1 * tmpValue2 * tmpDeriv

      ! second contribution - 2
      tmpValue1 = sqrt( (tmpHam(1,1) - tmpHam(2,2))**2 + 4.0_dp * tmpHam(1,2)**2 )
      tmpValue2 = -0.5_dp / tmpHam(1,2)**2

      tmpDeriv(:,:) = 0.0_dp
      tmpDeriv(:,:) = tmpValue1 * tmpValue2 * hVec22

      phaseDeriv(:,:) = phaseDeriv + eigvecsSSR(1,1)**2 * tmpDeriv

      eigvecsSSRderiv(:,:,1,1) = - eigvecsSSR(1,2) * phaseDeriv
      eigvecsSSRderiv(:,:,1,2) = eigvecsSSR(1,1) * phaseDeriv
      eigvecsSSRderiv(:,:,2,1) = - eigvecsSSR(1,1) * phaseDeriv
      eigvecsSSRderiv(:,:,2,2) = - eigvecsSSR(1,2) * phaseDeriv

    end if

  end subroutine getSsrCoefDeriv


  !> Calculate TDP gradient except RT shift and SSR state coefficients terms
  subroutine TDPshift(iSquare, eigenvecs, over, coord0, mOrb, Qmat, symTdpVec, gradL,&
      & preTdp, Sderiv, unrelTdm, ZTtdp, SAweight, omega, weightIL, G1, getDenseAO,&
      & getAtomIndex, ii, grad)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> auxiliary matrix in AO basis related to transition dipole term
    real(dp), intent(inout) :: Qmat(:,:)

    !> symmetric part of transition dipole vector for gradients of TDP
    real(dp), intent(inout) :: symTdpVec(:,:)

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(in) :: gradL(:,:,:)

    !> prefactor before the derivatives of FONs with dipole integral
    real(dp), intent(in) :: preTdp(:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> unrelaxed transition density matrix between SSR or SA-REKS states
    real(dp), intent(in) :: unrelTdm(:,:)

    !> solution of A * Z = X equation with X is XTtdp
    real(dp), intent(in) :: ZTtdp(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Current index for x, y, z axis
    integer, intent(in) :: ii

    !> all gradient components except R*T shift term
    real(dp), intent(out) :: grad(:,:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: Rderiv(:)

    real(dp) :: tmpValue, tmpR, derivTmp
    integer :: iL, Lmax, nOrb, nAtom, sparseSize
    integer :: mu, nu, iAtom, iAtom1, iAtom2, l

    Lmax = size(weightIL,dim=1)
    nOrb = size(Sderiv,dim=1)
    nAtom = size(gradL,dim=2)
    sparseSize = size(getDenseAO,dim=1)

    allocate(tmpMat(nOrb,nOrb))
    allocate(Rderiv(sparseSize))

    call matMO2AO(Qmat, eigenvecs(:,:,1))
    call matMO2AO(symTdpVec, eigenvecs(:,:,1))

    tmpValue = sum(ZTtdp(:)*omega(:))
    ! energy derivative term originating from orbital response
    grad(:,:) = 0.0_dp
    do iL = 1, Lmax
      grad(:,:) = grad - gradL(:,:,iL) * &
          & SAweight(1)*G1*weightIL(iL)*tmpValue
    end do
    ! Q * S_deriv term originating from orbital response
    call shiftQSgrad_(Qmat + transpose(Qmat), Sderiv, -1.0_dp, iSquare, mOrb, grad)

    ! symmetric TDP * S_deriv term originating from orbital response
    call shiftQSgrad_(symTdpVec, Sderiv, 2.0_dp, iSquare, mOrb, grad)

    ! energy derivative terms originating from FON derivatives
    do iL = 1, Lmax
      grad(:,:) = grad + gradL(:,:,iL) * G1*weightIL(iL)*preTdp(ii)
    end do

    ! unrelaxed transition density * overlap * coordinates term; overlap derivative
    tmpMat(:,:) = 0.0_dp
    do mu = 1, nOrb
      do nu = 1, nOrb
        ! find proper atom index
        iAtom1 = getAtomIndex(mu)
        iAtom2 = getAtomIndex(nu)
        tmpMat(mu,nu) = 0.5_dp * unrelTdm(mu,nu) * (coord0(ii,iAtom1) + coord0(ii,iAtom2))
      end do
    end do
    call shiftQSgrad_(tmpMat + transpose(tmpMat), Sderiv, 1.0_dp, iSquare, mOrb, grad)

    ! unrelaxed transition density * overlap * coordinates term; coordinate derivative
    tmpMat(:,:) = unrelTdm + transpose(unrelTdm)
    do mu = 1, nOrb
      tmpMat(mu,mu) = unrelTdm(mu,mu)
    end do
    do iAtom = 1, nAtom

      Rderiv(:) = 0.0_dp
      do l = 1, sparseSize

        ! set the AO indices with respect to sparsity
        mu = getDenseAO(l,1)
        nu = getDenseAO(l,2)
        ! find proper atom index
        iAtom1 = getAtomIndex(mu)
        iAtom2 = getAtomIndex(nu)

        if (mu <= nu) then
          if (abs(over(l)) >= epsilon(1.0_dp)) then

            tmpR = 0.0_dp
            ! mu in alpha
            if (iAtom1 == iAtom) then
              tmpR = tmpR + 1.0_dp
            end if
            ! nu in alpha
            if (iAtom2 == iAtom) then
              tmpR = tmpR + 1.0_dp
            end if

            ! calculate the contribution of position derivatives to Rderiv
            Rderiv(l) = Rderiv(l) + 0.5_dp*tmpR*tmpMat(mu,nu)

          end if
        end if

      end do
      derivTmp = sum(over(:)*Rderiv(:))
      grad(ii,iAtom) = grad(ii,iAtom) + derivTmp

    end do

    contains

      ! TODO : This routine is temporary, must be modified!
      !> Calculate Q*S gradient contribution
      subroutine shiftQSgrad_(Qmat, Sderiv, fac, iSquare, mOrb, deriv)

        !> auxiliary matrix in AO basis related to SA-REKS term
        real(dp), intent(in) :: Qmat(:,:)

        !> Dense overlap derivative in AO basis
        real(dp), intent(in) :: Sderiv(:,:,:)

        !> factor of Q*S term, for SI term it becomes -1 otherwise 1
        real(dp), intent(in) :: fac

        !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
        integer, intent(in) :: iSquare(:)

        !> Max. nr. of orbitals for any species
        integer, intent(in) :: mOrb

        !> gradient contribution from tr(Q*S)
        real(dp), intent(inout) :: deriv(:,:)

        real(dp), allocatable :: sqrQtmp(:,:)
        real(dp), allocatable :: sqrStmp(:,:,:)

        real(dp) :: derivTmp(3)
        integer :: iAtom1, nOrb1, iAtom2, nOrb2, ii
        integer :: nAtom, nAO1, nAO2, iOrb1, iOrb2

        nAtom = size(deriv,dim=2)

        allocate(sqrQtmp(mOrb,mOrb))
        allocate(sqrStmp(mOrb,mOrb,3))

        do iAtom1 = 1, nAtom
          nOrb1 = iSquare(iAtom1+1) - iSquare(iAtom1)
          nAO1 = iSquare(iAtom1) - 1
          do iAtom2 = iAtom1, nAtom
            nOrb2 = iSquare(iAtom2+1) - iSquare(iAtom2)
            nAO2 = iSquare(iAtom2) - 1
            if (iAtom1 /= iAtom2) then

              sqrQtmp(:,:) = 0.0_dp
              sqrStmp(:,:,:) = 0.0_dp
              do iOrb1 = 1, nOrb1
                do iOrb2 = 1, nOrb2
                  sqrQtmp(iOrb2,iOrb1) = Qmat(nAO2+iOrb2,nAO1+iOrb1)
                  do ii = 1, 3
                    sqrStmp(iOrb2,iOrb1,ii) = Sderiv(nAO2+iOrb2,nAO1+iOrb1,ii)
                  end do
                end do
              end do

              derivTmp(:) = 0.0_dp
              do ii = 1, 3
                derivTmp(ii) = sum(sqrQtmp(1:nOrb2,1:nOrb1)*sqrStmp(1:nOrb2,1:nOrb1,ii))
              end do

              ! forces from atom 1 on atom 2f and 2f onto 1
              deriv(:,iAtom1) = deriv(:,iAtom1) + fac*derivTmp(:)
              deriv(:,iAtom2) = deriv(:,iAtom2) - fac*derivTmp(:)

            end if
          end do
        end do

      end subroutine shiftQSgrad_

  end subroutine TDPshift


  !> Calculate dipole integral in DFTB formalism
  subroutine getDipoleIntegral(coord0, over, getAtomIndex, dipoleInt)

    !> central cell coordinates of atoms
    real(dp), intent(in) :: coord0(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: over(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> dipole integral in DFTB formalism
    real(dp), intent(out) :: dipoleInt(:,:,:)

    real(dp), allocatable :: R(:,:)

    integer :: ii, mu, nu, nOrb

    nOrb = size(over,dim=1)

    allocate(R(nOrb,3))

    R(:,:) = 0.0_dp
    do mu = 1, nOrb
      ii = getAtomIndex(mu)
      R(mu,:) = coord0(:,ii)
    end do

    ! TODO : Only Mulliken term is considered for dipole integral at the moment even though
    !        the onsite integrals are included during SCC iteration
    dipoleInt(:,:,:) = 0.0_dp
    do ii = 1, 3
      do mu = 1, nOrb
        do nu = 1, nOrb
          ! A negative sign must be included due to the electron charge (e = -1)
          dipoleInt(mu,nu,ii) = -0.5_dp * over(mu,nu) * (R(mu,ii) + R(nu,ii))
        end do
      end do
    end do

  end subroutine getDipoleIntegral


  !> Calculate dipole moment using dipole integral and density matrix
  subroutine getDipoleMomentMatrix(Pmat, dipoleInt, dipole)

    !> density matrix related to dipole
    real(dp), intent(in) :: Pmat(:,:)

    !> dipole integral in DFTB formalism
    real(dp), intent(in) :: dipoleInt(:,:,:)

    !> resulting dipole moment
    real(dp), intent(out) :: dipole(:)

    integer :: ii

    dipole(:) = 0.0_dp
    do ii = 1, 3
      dipole(ii) = -sum(dipoleInt(:,:,ii)*Pmat(:,:))
    end do

  end subroutine getDipoleMomentMatrix


  !> get the oscillator strength between the states
  subroutine getReksOsc(tdp, energy)

    !> transition dipole moment between states
    real(dp), intent(in) :: tdp(:,:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    real(dp) :: osc
    integer :: ia, ib, ist, nstates, nstHalf

    nstates = size(energy,dim=1)
    nstHalf = size(tdp,dim=2)

    write(stdOut,*)
    write(stdOut,'(A)') " Oscillator Strength (au)"
    do ist = 1, nstHalf

      call getTwoIndices(nstates, ist, ia, ib, 1)

      osc = 2.0_dp / 3.0_dp * (energy(ib) - energy(ia)) * sum(tdp(:,ist)**2)

      write(stdOut,'(A4,I1,A6,I1,A5)',advance="no") " ( S", ia - 1, " <-> S", ib - 1, " ) : "
      write(stdOut,'(1(f12.6))') osc

    end do

  end subroutine getReksOsc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate filling information for active space in (2,2) case
  subroutine getActiveFilling22_(FONs, Nc, tmpFilling)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> temporary average filling for SA-REKS state
    real(dp), intent(inout) :: tmpFilling(:,:)

    real(dp) :: n_a, n_b
    integer :: a, b, nstates

    nstates = size(tmpFilling,dim=2)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    ! PPS state fillings
    tmpFilling(a,1) = n_a
    tmpFilling(b,1) = n_b
    ! OSS state fillings
    tmpFilling(a,2) = 1.0_dp
    tmpFilling(b,2) = 1.0_dp
    if (nstates == 3) then
      ! DES state fillings
      tmpFilling(a,3) = n_b
      tmpFilling(b,3) = n_a
    end if

  end subroutine getActiveFilling22_


  !> Calculate unrelaxed transition density between SA-REKS states in (2,2) case
  subroutine getUnrelaxedTDM22_(eigenvecs, FONs, Nc, nstates, rhoXdel)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> unrelaxed transition density matrix between SA-REKS states
    real(dp), intent(inout) :: rhoXdel(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        rhoXdel(nu,mu,1) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
            & (sqrt(n_a) - sqrt(n_b))
      end do
    end do
    if (nstates == 3) then
      do mu = 1, nOrb
        do nu = 1, nOrb
          rhoXdel(nu,mu,3) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
              & (sqrt(n_a) + sqrt(n_b))
        end do
      end do
    end if

  end subroutine getUnrelaxedTDM22_


  !> Calculate response part of relaxed density for transition density contribution
  subroutine getResponseTDM22_(eigenvecs, FONs, SAweight, Rab, &
      & G1, Nc, resTdm)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> a part of transition density matrix originating from the
    !> response of the orbital occupation numbers
    real(dp), intent(inout) :: resTdm(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b, nstHalf

    nstHalf = size(resTdm,dim=3)
    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1); n_b = FONs(2,1)
    a = Nc + 1; b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        resTdm(nu,mu,1) = resTdm(nu,mu,1) + SAweight(1) * &
            & ( (n_a-1.0_dp)*sqrt(n_a) - (n_b-1.0_dp)*sqrt(n_b) ) * &
            & eigenvecs(mu,a) * eigenvecs(nu,b)
        resTdm(nu,mu,1) = resTdm(nu,mu,1) - 2.0_dp * G1 * &
            & Rab(1,2) * ( eigenvecs(mu,a) * eigenvecs(nu,a) - &
            & eigenvecs(mu,b) * eigenvecs(nu,b) )
      end do
    end do
    if (nstHalf == 3) then
      do mu = 1, nOrb
        do nu = 1, nOrb
          resTdm(nu,mu,3) = resTdm(nu,mu,3) + SAweight(1) * &
              & ( (n_a-1.0_dp)*sqrt(n_a) + (n_b-1.0_dp)*sqrt(n_b) ) * &
              & eigenvecs(mu,a) * eigenvecs(nu,b)
          resTdm(nu,mu,3) = resTdm(nu,mu,3) - 2.0_dp * G1 * &
              & Rab(2,3) * ( eigenvecs(mu,a) * eigenvecs(nu,a) - &
              & eigenvecs(mu,b) * eigenvecs(nu,b) )
        end do
      end do
    end if

  end subroutine getResponseTDM22_


  !> Calculate response part of relaxed density for density contribution
  subroutine getResponseDM22_(eigenvecs, ZT, tmpZ, omega, &
      & SAweight, G1, Nc, resRho)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> temporary matrix including ZT in MO basis
    real(dp), intent(in) :: tmpZ(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> a part of relaxed density matrix originating from the
    !> response of the orbital occupation numbers with XT
    real(dp), intent(out) :: resRho(:,:)

    real(dp) :: tmpValue
    integer :: a, b, mu, nu, nOrb

    nOrb = size(eigenvecs,dim=1)

    a = Nc + 1
    b = Nc + 2

    tmpValue = sum(ZT(:)*omega(:))
    do mu = 1, nOrb
      do nu = 1, nOrb
        resRho(nu,mu) = tmpZ(mu,nu) - SAweight(1) * G1 * tmpValue * &
            & (eigenvecs(mu,a)*eigenvecs(nu,a) - eigenvecs(mu,b)*eigenvecs(nu,b))
      end do
    end do

  end subroutine getResponseDM22_


  !> Calculate response part of L-th relaxed density for density contribution
  subroutine getResponseDML22_(rhoSqrL, SAweight, ZT, omega, &
      & weightIL, G1, Lpaired, resRhoL)

    !> Dense density matrix for each microstate
    real(dp), intent(in) :: rhoSqrL(:,:,:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> response part of relaxed density matrix for target L-th state
    real(dp), intent(inout) :: resRhoL(:,:)

    real(dp) :: tmpValue
    integer :: tmpL, iL, Lmax

    Lmax = size(rhoSqrL,dim=4)

    tmpValue = sum(ZT(:,1)*omega(:))
    do iL = 1, Lmax

      ! find proper index for down spin in rhoSqrL
      if (iL <= Lpaired) then
        tmpL = iL
      else
        if (mod(iL,2) == 1) then
          tmpL = iL + 1
        else
          tmpL = iL - 1
        end if
      end if
      resRhoL(:,:) = resRhoL + SAweight(1) * G1 * tmpValue * &
          & weightIL(iL) * (rhoSqrL(:,:,1,iL) + rhoSqrL(:,:,1,tmpL))

    end do

  end subroutine getResponseDML22_


  !> Calculate transition occupation matrix between SA-REKS states in (2,2) case
  subroutine getTranOccMat22_(FONs, Nc, nstates, matDel)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> temporary transition occupation matrix between SA-REKS states
    real(dp), intent(inout) :: matDel(:,:,:)

    real(dp) :: n_a, n_b
    integer :: a, b

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    matDel(a,b,1) = sqrt(n_a) - sqrt(n_b)
    if (nstates == 3) then
      matDel(a,b,3) = sqrt(n_a) + sqrt(n_b)
    end if

  end subroutine getTranOccMat22_


  !> Calculate density matrix contribution from FON derivative of weighting factors in (2,2) case
  subroutine getDensFon22_(rhoSqrL, weightIL, Lpaired, matX)

    !> Dense density matrix for each microstate
    real(dp), intent(in) :: rhoSqrL(:,:,:,:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> temporary density matrix contribution with FON derivatives
    real(dp), intent(inout) :: matX(:,:,:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: fac(:)
    integer :: ist, iL, nOrb, nstates, Lmax

    nOrb = size(rhoSqrL,dim=1)
    nstates = size(matX,dim=3)
    Lmax = size(weightIL,dim=1)

    allocate(tmpMat(nOrb,nOrb))
    allocate(fac(nstates))

    ! Loop for contributions of SA-REKS(2,2) states
    do ist = 1, nstates
      if (ist == 1) then
        fac(ist) = 1.0_dp
      else if (ist == 2) then
        ! The weighting factors for OSS states are fixed numbers,
        ! so they do not contribute the TDP gradient
        fac(ist) = 0.0_dp
      else if (ist == 3) then
        fac(ist) = -1.0_dp
      end if
    end do

    ! Derivative of n_a generates I_L contribution
    tmpMat(:,:) = 0.0_dp
    do iL = 1, Lmax
      if (iL <= Lpaired) then
        tmpMat(:,:) = tmpMat + rhoSqrL(:,:,1,iL) * weightIL(iL)
      else
        if (mod(iL,2) == 1) then
          tmpMat(:,:) = tmpMat + rhoSqrL(:,:,1,iL) * weightIL(iL)
        else
          tmpMat(:,:) = tmpMat + rhoSqrL(:,:,1,iL-1) * weightIL(iL)
        end if
      end if
    end do

    matX(:,:,:) = 0.0_dp
    do ist = 1, nstates
      matX(:,:,ist) = fac(ist) * tmpMat
    end do

  end subroutine getDensFon22_


  !> Calculate density matrix contribution from FON derivative of SI terms in (2,2) case
  subroutine getTranFon22_(eigenvecs, FONs, Nc, nstates, matDel)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> temporary transition density matrix contribution with FON derivatives
    real(dp), intent(inout) :: matDel(:,:,:)

    real(dp) :: n_a, n_b
    integer :: a, b, mu, nu, nOrb

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        matDel(mu,nu,1) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
          & (1.0_dp/sqrt(n_a) + 1.0_dp/sqrt(n_b))
        if (nstates == 3) then
          matDel(mu,nu,3) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
            & (1.0_dp/sqrt(n_a) - 1.0_dp/sqrt(n_b))
        end if
      end do
    end do

  end subroutine getTranFon22_


end module dftbp_reks_reksproperty
