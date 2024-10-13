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
module dftbp_reks_reksen
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : globalTimers, TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_elecsolvers_elecsolvers, only: TElectronicSolver
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_reks_rekscommon, only : getFactor, getTwoIndices, matAO2MO
  use dftbp_reks_reksio, only : printReksSSRInfo
  use dftbp_reks_reksvar, only : TReksCalc, reksTypes
  use dftbp_type_densedescr, only : TDenseDescr

  implicit none

  private
  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: adjustEigenval, solveSecularEqn
  public :: setReksTargetEnergy

  contains

  !> Construct L, spin dependent microstates from identical KS orbitals
  subroutine constructMicrostates(this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getFillingL22_(this%Nc, this%fillingL)
    case (reksTypes%ssr44)
      call getFillingL44_(this%Nc, this%fillingL)
    end select

  end subroutine constructMicrostates


  !> Calculate the weight of each microstate for current cycle, C_L
  subroutine calcWeights(this, tConverged)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Has the calculation converged?
    logical, intent(in) :: tConverged

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getWeightL22_(this%FONs, this%delta, this%SAweight, this%weightL, this%weight)
    case (reksTypes%ssr44)
      call getWeightL44_(this%FONs, this%delta, this%SAweight, this%Efunction, &
          & this%tAllStates, tConverged, this%weightL, this%weight)
    end select

  end subroutine calcWeights


  !> Swap the active orbitals for feasible occupation in REKS
  subroutine activeOrbSwap(this, eigenvecs)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call MOswap22_(eigenvecs, this%SAweight, this%FONs, this%Efunction, this%Nc)
    case (reksTypes%ssr44)
      call MOswap44_(eigenvecs, this%SAweight, this%FONs, this%Efunction, this%Nc)
    end select

  end subroutine activeOrbSwap


  !> Calculate filling for minimzed state with optimized FONs
  subroutine getFilling(this, filling)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> occupations (level)
    real(dp), intent(out) :: filling(:)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getFilling22_(filling, this%SAweight, this%FONs, this%Efunction, this%Nc)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getFilling


  !> Calculate the energy of SA-REKS states and averaged state
  subroutine calcSaReksEnergy(this, energy)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Energy terms in the system
    type(TEnergies), intent(inout) :: energy

    integer :: ist

    ! Compute the energy contributions for target SA-REKS state
    ! electronic energy = nonSCC + scc + spin + 3rd + fock + onsite
    energy%EnonSCC = sum(this%weightL(this%rstate,:)*this%enLnonSCC(:))
    energy%Escc = sum(this%weightL(this%rstate,:)*this%enLscc(:))
    energy%Espin = sum(this%weightL(this%rstate,:)*this%enLspin(:))
    if (this%t3rd) then
      energy%e3rd = sum(this%weightL(this%rstate,:)*this%enL3rd(:))
    end if
    if (this%isOnsite) then
      energy%eOnSite = sum(this%weightL(this%rstate,:)*this%enLonSite(:))
    end if
    if (this%isHybridXc) then
      energy%Efock = sum(this%weightL(this%rstate,:)*this%enLfock(:))
    end if
    if (this%isRS_OnsCorr) then
      energy%EfockOnSite = sum(this%weightL(this%rstate,:)*this%enLfockOnSite(:))
    end if
    if (this%isDispersion) then
      energy%Edisp = sum(this%weightL(this%rstate,:)*this%enLdisp(:))
    end if

    energy%Eelec = energy%EnonSCC + energy%Escc + energy%Espin + &
        & energy%e3rd + energy%eOnSite + energy%Efock + energy%EfockOnSite
    energy%Etotal = energy%Eelec + energy%Erep + energy%Edisp

    ! Compute the total energy for SA-REKS states
    do ist = 1, this%nstates
      this%energy(ist) = sum(this%weightL(ist,:)*this%enLtot(:))
    end do

!    if (abs(energy%Etotal - this%energy(this%rstate)) >= epsilon(1.0_dp)) then
    if (abs(energy%Etotal - this%energy(this%rstate)) >= 1.0e-8_dp) then
      call error("Wrong energy contribution for target SA-REKS state")
    end if

    ! In this step Eavg becomes the energy of averaged state
    ! From this energy we can check the variational principle
    energy%Eavg = 0.0_dp
    do ist = 1, this%SAstates
      energy%Eavg = energy%Eavg + this%SAweight(ist) * this%energy(ist)
    end do

  end subroutine calcSaReksEnergy


  !> Make pseudo-fock operator with Hamiltonian of each microstate
  !> and diagonalize the fock matrix
  subroutine getFockandDiag(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, &
      & electronicSolver, eigen, this)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: orbFON(:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: ii, nOrb

    nOrb = size(this%fockFc,dim=1)

    allocate(orbFON(nOrb))
    allocate(tmpMat(nOrb,nOrb))

    call getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, this%hamSqrL, this%hamSpL, this%weight, &
        & this%fillingL, this%Nc, this%Na, this%Lpaired, this%isHybridXc, &
        & this%isRS_OnsCorr, orbFON, this%fockFc, this%fockFa)

    call matAO2MO(this%fockFc, eigenvecs(:,:,1))
    do ii = 1, this%Na
      call matAO2MO(this%fockFa(:,:,ii), eigenvecs(:,:,1))
    end do

    call getPseudoFock_(this%fockFc, this%fockFa, orbFON, this%Nc, this%Na, this%fock)

    call levelShifting_(this%fock, this%shift, this%Nc, this%Na)

    ! Diagonalize the pesudo-Fock matrix
    tmpMat(:,:) = this%fock

    eigen(:,1,1) = 0.0_dp
    call heev(tmpMat, eigen(:,1,1), 'U', 'V')
    this%eigvecsFock(:,:) = tmpMat

  end subroutine getFockandDiag


  !> guess new eigenvectors from Fock eigenvectors
  subroutine guessNewEigvecs(eigenvecs, eigvecsFock)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> eigenvectors from pesudo-fock matrix
    real(dp), intent(in) :: eigvecsFock(:,:)

    real(dp), allocatable :: tmpVec(:,:)
    integer :: nOrb

    nOrb = size(eigvecsFock,dim=1)

    allocate(tmpVec(nOrb,nOrb))

    tmpVec(:,:) = 0.0_dp
    call gemm(tmpVec, eigenvecs, eigvecsFock)
    eigenvecs(:,:) = tmpVec

  end subroutine guessNewEigvecs


  !> adjust the eigenvalues (eliminate shift values)
  subroutine adjustEigenval(this, eigen)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> eigenvalues
    real(dp), intent(inout) :: eigen(:,:,:)

    integer :: nOrb, ind, ii

    nOrb = size(eigen,dim=1)

    do ii = this%Nc + 1, this%Nc + this%Na
      ind = ii - this%Nc
      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * this%shift
    end do

    do ii = this%Nc + this%Na + 1, nOrb
      ind = this%Na + 1
      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * this%shift
    end do

  end subroutine adjustEigenval


  !> Solve secular equation with coupling element between SA-REKS states
  subroutine solveSecularEqn(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, electronicSolver, &
      & eigenvecs, this)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: Wab(:,:)
    real(dp), allocatable :: StateCoup(:,:)
    real(dp), allocatable :: tmpState(:,:)
    real(dp), allocatable :: tmpEigen(:)
    real(dp), allocatable :: tmpEn(:)

    integer :: ist, jst, nActPair

    nActPair = this%Na * (this%Na - 1) / 2

    allocate(Wab(nActPair,2))
    allocate(StateCoup(this%nstates,this%nstates))
    allocate(tmpState(this%nstates,this%nstates))
    allocate(tmpEigen(this%nstates))
    allocate(tmpEn(this%nstates))

    call getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, eigenvecs(:,:,1), this%hamSqrL, &
        & this%hamSpL, this%weight, this%fillingL, this%Nc, this%Na, &
        & this%Lpaired, this%isHybridXc, this%isRS_OnsCorr, Wab)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getStateCoup22_(Wab, this%FONs, StateCoup)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! diagonalize the state energies
    ! obtain SSR energies & state-interaction term
    tmpEigen(:) = 0.0_dp

    tmpState(:,:) = 0.0_dp
    do ist = 1, this%nstates
      do jst = 1, this%nstates
        if (ist == jst) then
          tmpState(ist,jst) = this%energy(ist)
        else
          tmpState(ist,jst) = StateCoup(ist,jst)
        end if
      end do
    end do

    ! save state energies to print information
    tmpEn(:) = this%energy
    if (this%tSSR) then
      call heev(tmpState, tmpEigen, 'U', 'V')
      this%eigvecsSSR(:,:) = tmpState
      this%energy(:) = tmpEigen
    end if

    ! print state energies and couplings
    call printReksSSRInfo(this, Wab, tmpEn, StateCoup)

  end subroutine solveSecularEqn


  !> Set correct final energy values for target state or microstate
  subroutine setReksTargetEnergy(this, energy, cellVol, pressure)

    !> data type for REKS
    type(TReksCalc), intent(in) :: this

    !> Energy terms in the system
    type(TEnergies), intent(inout) :: energy

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> External pressure
    real(dp), intent(in) :: pressure

    ! get correct energy values
    if (this%Lstate == 0) then

      ! get energy contributions for target state
      energy%Etotal = this%energy(this%rstate)
      if (this%nstates > 1) then
        energy%Eexcited = this%energy(this%rstate) - this%energy(1)
      else
        energy%Eexcited = 0.0_dp
      end if

    else

      ! get energy contributions for target microstate
      energy%EnonSCC = this%enLnonSCC(this%Lstate)
      energy%ESCC = this%enLSCC(this%Lstate)
      energy%Espin = this%enLspin(this%Lstate)
      if (this%t3rd) then
        energy%e3rd = this%enL3rd(this%Lstate)
      end if
      if (this%isOnsite) then
        energy%eOnSite = this%enLonsite(this%Lstate)
      end if
      if (this%isHybridXc) then
        energy%Efock = this%enLfock(this%Lstate)
      end if
      if (this%isRS_OnsCorr) then
        energy%EfockOnSite = this%enLfockOnsite(this%Lstate)
      end if
      if (this%isDispersion) then
        energy%Edisp = this%enLdisp(this%Lstate)
      end if

      energy%Eelec = energy%EnonSCC + energy%Escc + energy%Espin + &
          & energy%e3rd + energy%eOnSite + energy%Efock + energy%EfockOnsite
      energy%Etotal = energy%Eelec + energy%Erep + energy%Edisp
      energy%Eexcited = 0.0_dp

!      if (abs(energy%Etotal - this%enLtot(this%Lstate)) > epsilon(1.0_dp)) then
      if (abs(energy%Etotal - this%enLtot(this%Lstate)) > 1.0e-8_dp) then
        call error("Wrong energy contribution for target microstate")
      end if

    end if

    ! REKS is not affected by filling, so TS becomes 0
    energy%EMermin = energy%Etotal
    ! extrapolated to 0 K
    energy%Ezero = energy%Etotal
    energy%EGibbs = energy%EMermin + cellVol * pressure
    energy%EForceRelated = energy%EGibbs

  end subroutine setReksTargetEnergy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate filling of each microstate in REKS(2,2)
  subroutine getFillingL22_(Nc, fillingL)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Filling for each microstate
    real(dp), intent(out) :: fillingL(:,:,:)

    integer :: iL, iSpin, ii, nSpin, Lmax

    nSpin = size(fillingL,dim=2)
    Lmax = size(fillingL,dim=3)

    fillingL(:,:,:) = 0.0_dp

    ! Filling of core orbitals
    do iL = 1, Lmax
      do iSpin = 1, nSpin
        do ii = 1, Nc
          fillingL(ii,iSpin,iL) = 1.0_dp
        end do
      end do
    end do

    ! Filling of active orbitals for REKS(2,2) case
    ! 1 = aa'
    fillingL(Nc+1,1,1) = 1.0_dp; fillingL(Nc+1,2,1) = 1.0_dp
    ! 2 = bb'
    fillingL(Nc+2,1,2) = 1.0_dp; fillingL(Nc+2,2,2) = 1.0_dp
    ! 3 = ab'
    fillingL(Nc+1,1,3) = 1.0_dp; fillingL(Nc+2,2,3) = 1.0_dp
    ! 4 = a'b
    fillingL(Nc+2,1,4) = 1.0_dp; fillingL(Nc+1,2,4) = 1.0_dp
    ! 5 = ab
    fillingL(Nc+1,1,5) = 1.0_dp; fillingL(Nc+2,1,5) = 1.0_dp
    ! 6 = a'b'
    fillingL(Nc+1,2,6) = 1.0_dp; fillingL(Nc+2,2,6) = 1.0_dp

  end subroutine getFillingL22_


  !> Calculate filling of each microstate in REKS(4,4)
  subroutine getFillingL44_(Nc, fillingL)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Filling for each microstate
    real(dp), intent(out) :: fillingL(:,:,:)

    integer :: iL, iSpin, ii, nSpin, Lmax

    nSpin = size(fillingL,dim=2)
    Lmax = size(fillingL,dim=3)

    fillingL(:,:,:) = 0.0_dp

    ! Filling of core orbitals
    do iL = 1, Lmax
      do iSpin = 1, nSpin
        do ii = 1, Nc
          fillingL(ii,iSpin,iL) = 1.0_dp
        end do
      end do
    end do

    ! Filling of active orbitals for REKS(4,4) case
    ! 1 = aa'bb'
    fillingL(Nc+1,1,1) = 1.0_dp; fillingL(Nc+1,2,1) = 1.0_dp
    fillingL(Nc+2,1,1) = 1.0_dp; fillingL(Nc+2,2,1) = 1.0_dp
    ! 2 = aa'cc' 
    fillingL(Nc+1,1,2) = 1.0_dp; fillingL(Nc+1,2,2) = 1.0_dp
    fillingL(Nc+3,1,2) = 1.0_dp; fillingL(Nc+3,2,2) = 1.0_dp
    ! 3 = bb'dd'
    fillingL(Nc+2,1,3) = 1.0_dp; fillingL(Nc+2,2,3) = 1.0_dp
    fillingL(Nc+4,1,3) = 1.0_dp; fillingL(Nc+4,2,3) = 1.0_dp
    ! 4 = cc'dd'
    fillingL(Nc+3,1,4) = 1.0_dp; fillingL(Nc+3,2,4) = 1.0_dp
    fillingL(Nc+4,1,4) = 1.0_dp; fillingL(Nc+4,2,4) = 1.0_dp
    ! 5 = abb'd'
    fillingL(Nc+1,1,5) = 1.0_dp; fillingL(Nc+2,1,5) = 1.0_dp
    fillingL(Nc+2,2,5) = 1.0_dp; fillingL(Nc+4,2,5) = 1.0_dp
    ! 6 = a'bb'd
    fillingL(Nc+1,2,6) = 1.0_dp; fillingL(Nc+2,1,6) = 1.0_dp
    fillingL(Nc+2,2,6) = 1.0_dp; fillingL(Nc+4,1,6) = 1.0_dp
    ! 7 = abb'd
    fillingL(Nc+1,1,7) = 1.0_dp; fillingL(Nc+2,1,7) = 1.0_dp
    fillingL(Nc+2,2,7) = 1.0_dp; fillingL(Nc+4,1,7) = 1.0_dp
    ! 8 = a'bb'd'
    fillingL(Nc+1,2,8) = 1.0_dp; fillingL(Nc+2,1,8) = 1.0_dp
    fillingL(Nc+2,2,8) = 1.0_dp; fillingL(Nc+4,2,8) = 1.0_dp
    ! 9 = aa'bc'
    fillingL(Nc+1,1,9) = 1.0_dp; fillingL(Nc+1,2,9) = 1.0_dp
    fillingL(Nc+2,1,9) = 1.0_dp; fillingL(Nc+3,2,9) = 1.0_dp
    ! 10 = aa'b'c
    fillingL(Nc+1,1,10) = 1.0_dp; fillingL(Nc+1,2,10) = 1.0_dp
    fillingL(Nc+2,2,10) = 1.0_dp; fillingL(Nc+3,1,10) = 1.0_dp
    ! 11 = aa'bc
    fillingL(Nc+1,1,11) = 1.0_dp; fillingL(Nc+1,2,11) = 1.0_dp
    fillingL(Nc+2,1,11) = 1.0_dp; fillingL(Nc+3,1,11) = 1.0_dp
    ! 12 = aa'b'c'
    fillingL(Nc+1,1,12) = 1.0_dp; fillingL(Nc+1,2,12) = 1.0_dp
    fillingL(Nc+2,2,12) = 1.0_dp; fillingL(Nc+3,2,12) = 1.0_dp
    ! 13 = bc'dd'
    fillingL(Nc+2,1,13) = 1.0_dp; fillingL(Nc+3,2,13) = 1.0_dp
    fillingL(Nc+4,1,13) = 1.0_dp; fillingL(Nc+4,2,13) = 1.0_dp
    ! 14 = b'cdd'
    fillingL(Nc+2,2,14) = 1.0_dp; fillingL(Nc+3,1,14) = 1.0_dp
    fillingL(Nc+4,1,14) = 1.0_dp; fillingL(Nc+4,2,14) = 1.0_dp
    ! 15 = acc'd'
    fillingL(Nc+1,1,15) = 1.0_dp; fillingL(Nc+3,1,15) = 1.0_dp
    fillingL(Nc+3,2,15) = 1.0_dp; fillingL(Nc+4,2,15) = 1.0_dp
    ! 16 = a'cc'd
    fillingL(Nc+1,2,16) = 1.0_dp; fillingL(Nc+3,1,16) = 1.0_dp
    fillingL(Nc+3,2,16) = 1.0_dp; fillingL(Nc+4,1,16) = 1.0_dp
    ! 17 = aa'bd'
    fillingL(Nc+1,1,17) = 1.0_dp; fillingL(Nc+1,2,17) = 1.0_dp
    fillingL(Nc+2,1,17) = 1.0_dp; fillingL(Nc+4,2,17) = 1.0_dp
    ! 18 = aa'b'd
    fillingL(Nc+1,1,18) = 1.0_dp; fillingL(Nc+1,2,18) = 1.0_dp
    fillingL(Nc+2,2,18) = 1.0_dp; fillingL(Nc+4,1,18) = 1.0_dp
    ! 19 = aa'bd
    fillingL(Nc+1,1,19) = 1.0_dp; fillingL(Nc+1,2,19) = 1.0_dp
    fillingL(Nc+2,1,19) = 1.0_dp; fillingL(Nc+4,1,19) = 1.0_dp
    ! 20 = aa'b'd'
    fillingL(Nc+1,1,20) = 1.0_dp; fillingL(Nc+1,2,20) = 1.0_dp
    fillingL(Nc+2,2,20) = 1.0_dp; fillingL(Nc+4,2,20) = 1.0_dp
    ! 21 = abb'c'
    fillingL(Nc+1,1,21) = 1.0_dp; fillingL(Nc+2,1,21) = 1.0_dp
    fillingL(Nc+2,2,21) = 1.0_dp; fillingL(Nc+3,2,21) = 1.0_dp
    ! 22 = a'bb'c
    fillingL(Nc+1,2,22) = 1.0_dp; fillingL(Nc+2,1,22) = 1.0_dp
    fillingL(Nc+2,2,22) = 1.0_dp; fillingL(Nc+3,1,22) = 1.0_dp
    ! 23 = abb'c
    fillingL(Nc+1,1,23) = 1.0_dp; fillingL(Nc+2,1,23) = 1.0_dp
    fillingL(Nc+2,2,23) = 1.0_dp; fillingL(Nc+3,1,23) = 1.0_dp
    ! 24 = a'bb'c'
    fillingL(Nc+1,2,24) = 1.0_dp; fillingL(Nc+2,1,24) = 1.0_dp
    fillingL(Nc+2,2,24) = 1.0_dp; fillingL(Nc+3,2,24) = 1.0_dp
    ! 25 = bcc'd'
    fillingL(Nc+2,1,25) = 1.0_dp; fillingL(Nc+3,1,25) = 1.0_dp
    fillingL(Nc+3,2,25) = 1.0_dp; fillingL(Nc+4,2,25) = 1.0_dp
    ! 26 = b'cc'd
    fillingL(Nc+2,2,26) = 1.0_dp; fillingL(Nc+3,1,26) = 1.0_dp
    fillingL(Nc+3,2,26) = 1.0_dp; fillingL(Nc+4,1,26) = 1.0_dp
    ! 27 = ac'dd'
    fillingL(Nc+1,1,27) = 1.0_dp; fillingL(Nc+3,2,27) = 1.0_dp
    fillingL(Nc+4,1,27) = 1.0_dp; fillingL(Nc+4,2,27) = 1.0_dp
    ! 28 = a'cdd'
    fillingL(Nc+1,2,28) = 1.0_dp; fillingL(Nc+3,1,28) = 1.0_dp
    fillingL(Nc+4,1,28) = 1.0_dp; fillingL(Nc+4,2,28) = 1.0_dp
    ! 29 = abc'd'
    fillingL(Nc+1,1,29) = 1.0_dp; fillingL(Nc+2,1,29) = 1.0_dp
    fillingL(Nc+3,2,29) = 1.0_dp; fillingL(Nc+4,2,29) = 1.0_dp
    ! 30 = a'b'cd
    fillingL(Nc+1,2,30) = 1.0_dp; fillingL(Nc+2,2,30) = 1.0_dp
    fillingL(Nc+3,1,30) = 1.0_dp; fillingL(Nc+4,1,30) = 1.0_dp
    ! 31 = ab'cd'
    fillingL(Nc+1,1,31) = 1.0_dp; fillingL(Nc+2,2,31) = 1.0_dp
    fillingL(Nc+3,1,31) = 1.0_dp; fillingL(Nc+4,2,31) = 1.0_dp
    ! 32 = a'bc'd
    fillingL(Nc+1,2,32) = 1.0_dp; fillingL(Nc+2,1,32) = 1.0_dp
    fillingL(Nc+3,2,32) = 1.0_dp; fillingL(Nc+4,1,32) = 1.0_dp
    ! 33 = ab'c'd
    fillingL(Nc+1,1,33) = 1.0_dp; fillingL(Nc+2,2,33) = 1.0_dp
    fillingL(Nc+3,2,33) = 1.0_dp; fillingL(Nc+4,1,33) = 1.0_dp
    ! 34 = a'bcd'
    fillingL(Nc+1,2,34) = 1.0_dp; fillingL(Nc+2,1,34) = 1.0_dp
    fillingL(Nc+3,1,34) = 1.0_dp; fillingL(Nc+4,2,34) = 1.0_dp
    ! 35 = abcd
    fillingL(Nc+1,1,35) = 1.0_dp; fillingL(Nc+2,1,35) = 1.0_dp
    fillingL(Nc+3,1,35) = 1.0_dp; fillingL(Nc+4,1,35) = 1.0_dp
    ! 36 = a'b'c'd'
    fillingL(Nc+1,2,36) = 1.0_dp; fillingL(Nc+2,2,36) = 1.0_dp
    fillingL(Nc+3,2,36) = 1.0_dp; fillingL(Nc+4,2,36) = 1.0_dp

  end subroutine getFillingL44_


  !> Make (2e,2o) weights, C_L used in SA-REKS
  subroutine getWeightL22_(FONs, delta, SAweight, weightL, weight)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Weight for each microstate per state
    real(dp), intent(out) :: weightL(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(out) :: weight(:)

    real(dp) :: n_a, n_b, fac
    integer :: iL, Lmax, ist, SAstates, nstates

    Lmax = size(weightL,dim=2)
    SAstates = size(SAweight,dim=1)
    nstates = size(weightL,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    fac = getFactor(n_a, n_b, delta)

    weightL(1,1) = 0.5_dp*n_a
    weightL(1,2) = 0.5_dp*n_b
    weightL(1,3) = fac
    weightL(1,4) = fac
    weightL(1,5) = -fac
    weightL(1,6) = -fac

    if (nstates >= 2) then
      weightL(2,1) = 0.0_dp
      weightL(2,2) = 0.0_dp
      weightL(2,3) = 1.0_dp
      weightL(2,4) = 1.0_dp
      weightL(2,5) = -0.5_dp
      weightL(2,6) = -0.5_dp
    end if

    if (nstates >= 3) then
      weightL(3,1) = 0.5_dp*n_b
      weightL(3,2) = 0.5_dp*n_a
      weightL(3,3) = -fac
      weightL(3,4) = -fac
      weightL(3,5) = fac
      weightL(3,6) = fac
    end if

    ! Decide which state will be optimized; single-state REKS or SA-REKS
    ! Efunction = 1 -> PPS state is optimized
    ! Efunction = 2 -> (PPS+OSS)/2 state is optimized
    weight(:) = 0.0_dp
    do iL = 1, Lmax
      do ist = 1, SAstates
        weight(iL) = weight(iL) + SAweight(ist)*weightL(ist,iL)
      end do
    end do

  end subroutine getWeightL22_


  !> Make (4e,4o) weights, C_L used in SA-REKS
  subroutine getWeightL44_(FONs, delta, SAweight, Efunction, tAllStates, tConverged,&
      & weightL, weight)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Decide the energy states in SA-REKS
    logical, intent(in) :: tAllStates

    !> Has the calculation converged?
    logical, intent(in) :: tConverged

    !> Weight for each microstate per state
    real(dp), intent(out) :: weightL(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(out) :: weight(:)

    real(dp) :: n_a, n_b, n_c, n_d
    real(dp) :: np_a, np_b, np_c, np_d
    real(dp) :: mp_a, mp_b, mp_c, mp_d
    real(dp) :: fac1, fac2
    integer :: iL, Lmax, ist, SAstates
    integer :: indPPS, indOSS1, indOSS2
    integer :: indOSS3, indOSS4, indDOSS
    integer :: indDSPS, indDES1, indDES2

    Lmax = size(weightL,dim=2)
    SAstates = size(SAweight,dim=1)

    ! Reset the indices for SA-REKS(4,4) states
    indPPS = 0; indOSS1 = 0; indOSS2 = 0
    indOSS3 = 0; indOSS4 = 0; indDOSS = 0
    indDSPS = 0; indDES1 = 0; indDES2 = 0

    if (tAllStates .and. tConverged) then
      ! Calculate weight of 9 states
      indPPS = 1; indOSS1 = 2; indOSS2 = 3
      indOSS3 = 4; indOSS4 = 5; indDOSS = 6
      indDSPS = 7; indDES1 = 8; indDES2 = 9
    else
      ! Calculate weight of states used for state-averaging
      ! PPS state
      indPPS = 1
      if (Efunction == 2) then
        ! DSPS state
        indDSPS = 2
      else if (Efunction == 3 .or. Efunction == 4) then
        ! OSS1, OSS2 state
        indOSS1 = 2; indOSS2 = 3
        if (Efunction == 4) then
          ! OSS3, OSS4 state
          indOSS3 = 4; indOSS4 = 5
        end if
      end if
    end if

    n_a = FONs(1,1); n_b = FONs(2,1); n_c = FONs(3,1); n_d = FONs(4,1)
    np_a = FONs(1,2); np_b = FONs(2,2); np_c = FONs(3,2); np_d = FONs(4,2)
    mp_a = FONs(1,3); mp_b = FONs(2,3); mp_c = FONs(3,3); mp_d = FONs(4,3)

    if (indPPS > 0) then
      weightL(indPPS,:) = 0.0_dp
      fac1 = getFactor(n_a, n_d, delta)
      fac2 = getFactor(n_b, n_c, delta)
      weightL(indPPS,1) = 0.25_dp * n_a * n_b
      weightL(indPPS,2) = 0.25_dp * n_a * n_c
      weightL(indPPS,3) = 0.25_dp * n_b * n_d
      weightL(indPPS,4) = 0.25_dp * n_c * n_d
      weightL(indPPS,5) = fac1; weightL(indPPS,6) = fac1
      weightL(indPPS,7) = -fac1; weightL(indPPS,8) = -fac1
      weightL(indPPS,9) = fac2; weightL(indPPS,10) = fac2
      weightL(indPPS,11) = -fac2; weightL(indPPS,12) = -fac2
    end if

    if (indOSS1 > 0) then
      weightL(indOSS1,:) = 0.0_dp
      fac1 = getFactor(np_a, np_d, delta)
      weightL(indOSS1,5) = fac1; weightL(indOSS1,6) = fac1
      weightL(indOSS1,7) = -fac1; weightL(indOSS1,8) = -fac1
      weightL(indOSS1,9)  = 0.25_dp * np_a + 0.5_dp 
      weightL(indOSS1,10) = 0.25_dp * np_a + 0.5_dp
      weightL(indOSS1,11) = -0.5_dp
      weightL(indOSS1,12) = -0.5_dp
      weightL(indOSS1,13) = 0.25_dp * np_d
      weightL(indOSS1,14) = 0.25_dp * np_d
    end if

    if (indOSS2 > 0) then
      weightL(indOSS2,:) = 0.0_dp
      fac2 = getFactor(np_b, np_c, delta)
      weightL(indOSS2,5) = 0.25_dp * np_b + 0.5_dp 
      weightL(indOSS2,6) = 0.25_dp * np_b + 0.5_dp
      weightL(indOSS2,7) = -0.5_dp
      weightL(indOSS2,8) = -0.5_dp
      weightL(indOSS2,9)  = fac2; weightL(indOSS2,10) = fac2
      weightL(indOSS2,11) = -fac2; weightL(indOSS2,12) = -fac2
      weightL(indOSS2,15) = 0.25_dp * np_c
      weightL(indOSS2,16) = 0.25_dp * np_c
    end if

    if (indOSS3 > 0) then
      weightL(indOSS3,:) = 0.0_dp
      fac1 = getFactor(mp_a, mp_c, delta)
      weightL(indOSS3,17) = 0.25_dp * mp_a + 0.5_dp 
      weightL(indOSS3,18) = 0.25_dp * mp_a + 0.5_dp
      weightL(indOSS3,19) = -0.5_dp
      weightL(indOSS3,20) = -0.5_dp
      weightL(indOSS3,21) = fac1; weightL(indOSS3,22) = fac1
      weightL(indOSS3,23) = -fac1; weightL(indOSS3,24) = -fac1
      weightL(indOSS3,25) = 0.25_dp * mp_c
      weightL(indOSS3,26) = 0.25_dp * mp_c
    end if

    if (indOSS4 > 0) then
      weightL(indOSS4,:) = 0.0_dp
      fac2 = getFactor(mp_b, mp_d, delta)
      weightL(indOSS4,17) = fac2; weightL(indOSS4,18) = fac2
      weightL(indOSS4,19) = -fac2; weightL(indOSS4,20) = -fac2
      weightL(indOSS4,21) = 0.25_dp * mp_b + 0.5_dp 
      weightL(indOSS4,22) = 0.25_dp * mp_b + 0.5_dp
      weightL(indOSS4,23) = -0.5_dp
      weightL(indOSS4,24) = -0.5_dp
      weightL(indOSS4,27) = 0.25_dp * mp_d
      weightL(indOSS4,28) = 0.25_dp * mp_d
    end if

    if (indDOSS > 0) then
      weightL(indDOSS,:) = 0.0_dp
      weightL(indDOSS,29) = 0.5_dp; weightL(indDOSS,30) = 0.5_dp
      weightL(indDOSS,31) = 0.5_dp; weightL(indDOSS,32) = 0.5_dp
      weightL(indDOSS,33) = -0.25_dp; weightL(indDOSS,34) = -0.25_dp
      weightL(indDOSS,35) = -0.25_dp; weightL(indDOSS,36) = -0.25_dp
    end if

    if (indDSPS > 0) then
      weightL(indDSPS,:) = 0.0_dp
      weightL(indDSPS,33) = 0.75_dp; weightL(indDSPS,34) = 0.75_dp
      weightL(indDSPS,35) = -0.25_dp; weightL(indDSPS,36) = -0.25_dp
    end if

    if (indDES1 > 0) then
      weightL(indDES1,:) = 0.0_dp
      fac1 = getFactor(n_a, n_d, delta)
      fac2 = getFactor(n_b, n_c, delta)
      weightL(indDES1,1) = 0.25_dp * n_a * n_c
      weightL(indDES1,2) = 0.25_dp * n_a * n_b
      weightL(indDES1,3) = 0.25_dp * n_c * n_d
      weightL(indDES1,4) = 0.25_dp * n_b * n_d
      weightL(indDES1,5) = fac1; weightL(indDES1,6) = fac1
      weightL(indDES1,7) = -fac1; weightL(indDES1,8) = -fac1
      weightL(indDES1,9) = -fac2; weightL(indDES1,10) = -fac2
      weightL(indDES1,11) = fac2; weightL(indDES1,12) = fac2
    end if

    if (indDES2 > 0) then
      weightL(indDES2,:) = 0.0_dp
      fac1 = getFactor(n_a, n_d, delta)
      fac2 = getFactor(n_b, n_c, delta)
      weightL(indDES2,1) = 0.25_dp * n_d * n_b
      weightL(indDES2,2) = 0.25_dp * n_d * n_c
      weightL(indDES2,3) = 0.25_dp * n_b * n_a
      weightL(indDES2,4) = 0.25_dp * n_c * n_a
      weightL(indDES2,5) = -fac1; weightL(indDES2,6) = -fac1
      weightL(indDES2,7) = fac1; weightL(indDES2,8) = fac1
      weightL(indDES2,9) = fac2; weightL(indDES2,10) = fac2
      weightL(indDES2,11) = -fac2; weightL(indDES2,12) = -fac2
    end if

    ! Decide which state will be optimized; single-state REKS or SA-REKS
    ! Efunction = 1 -> PPS state is optimized
    ! Efunction = 2 -> (PPS+DSPS)/2 state is optimized
    ! Efunction = 3 -> (PPS+OSS1+OSS2)/3 state is optimized
    ! Efunction = 4 -> (PPS+OSS1+OSS2+OSS3+OSS4)/5 state is optimized
    weight(:) = 0.0_dp
    do iL = 1, Lmax
      do ist = 1, SAstates
        weight(iL) = weight(iL) + SAweight(ist)*weightL(ist,iL)
      end do
    end do

  end subroutine getWeightL44_


  !> Swap active orbitals when fa < fb in REKS(2,2) case
  subroutine MOswap22_(eigenvecs, SAweight, FONs, Efunction, Nc)

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    real(dp), allocatable :: tmpMO(:)

    real(dp) :: n_a, n_b, fa, fb
    integer :: nOrb

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    allocate(tmpMO(nOrb))

    if (Efunction == 1) then
      ! REKS charge
      fa = n_a * 0.5_dp
      fb = n_b * 0.5_dp
    else if (Efunction == 2) then
      ! 2SA-REKS charge
      fa = (n_a*SAweight(1) + SAweight(2)) * 0.5_dp
      fb = (n_b*SAweight(1) + SAweight(2)) * 0.5_dp
    end if

    if (fa < fb) then
      write(stdOut,'(A6,F9.6,A20,I4,A8,I4,A8)') " fa = ", fa, &
          & ", MO swap between a(", Nc+1, ") and b(", Nc+2, ") occurs"
      tmpMO(:) = eigenvecs(:,Nc+1)
      eigenvecs(:,Nc+1) = eigenvecs(:,Nc+2)
      eigenvecs(:,Nc+2) = tmpMO
    end if

  end subroutine MOswap22_


  !> Swap active orbitals when fa < fd or fb < fc in REKS(4,4) case
  subroutine MOswap44_(eigenvecs, SAweight, FONs, Efunction, Nc)

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    real(dp), allocatable :: tmpMO(:)

    real(dp) :: n_a, n_b, n_c, n_d
    real(dp) :: np_a, np_b, np_c, np_d
    real(dp) :: mp_a, mp_b, mp_c, mp_d
    real(dp) :: fa, fb, fc, fd
    integer :: nOrb

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1); n_b = FONs(2,1); n_c = FONs(3,1); n_d = FONs(4,1)
    np_a = FONs(1,2); np_b = FONs(2,2); np_c = FONs(3,2); np_d = FONs(4,2)
    mp_a = FONs(1,3); mp_b = FONs(2,3); mp_c = FONs(3,3); mp_d = FONs(4,3)

    allocate(tmpMO(nOrb))

    if (Efunction == 1) then
      ! REKS charge
      fa = n_a * 0.5_dp
      fb = n_b * 0.5_dp
      fc = n_c * 0.5_dp
      fd = n_d * 0.5_dp
    else if (Efunction == 2) then
      ! 2SA-REKS charge
      fa = (n_a*SAweight(1) + SAweight(2)) * 0.5_dp
      fb = (n_b*SAweight(1) + SAweight(2)) * 0.5_dp
      fc = (n_c*SAweight(1) + SAweight(2)) * 0.5_dp
      fd = (n_d*SAweight(1) + SAweight(2)) * 0.5_dp
    else if (Efunction == 3) then
      ! 3SA-REKS charge
      fa = (n_a*SAweight(1) + np_a*SAweight(2) + SAweight(3)) * 0.5_dp
      fb = (n_b*SAweight(1) + SAweight(2) + np_b*SAweight(3)) * 0.5_dp
      fc = (n_c*SAweight(1) + SAweight(2) + np_c*SAweight(3)) * 0.5_dp
      fd = (n_d*SAweight(1) + np_d*SAweight(2) + SAweight(3)) * 0.5_dp
    else if (Efunction == 4) then
      ! 5SA-REKS charge
      fa = (n_a*SAweight(1) + np_a*SAweight(2) + SAweight(3) &
          &+ mp_a*SAweight(4) + SAweight(5)) * 0.5_dp
      fb = (n_b*SAweight(1) + SAweight(2) + np_b*SAweight(3) &
          &+ SAweight(4) + mp_b*SAweight(5)) * 0.5_dp
      fc = (n_c*SAweight(1) + SAweight(2) + np_c*SAweight(3) &
          &+ mp_c*SAweight(4) + SAweight(5)) * 0.5_dp
      fd = (n_d*SAweight(1) + np_d*SAweight(2) + SAweight(3) &
          &+ SAweight(5) + mp_d*SAweight(5)) * 0.5_dp
    end if

    if (fa < fd) then
      write(stdOut,'(A6,F9.6,A20,I4,A8,I4,A8)') " fa = ", fa, &
          & ", MO swap between a(", Nc+1, ") and d(", Nc+4, ") occurs"
      tmpMO(:) = eigenvecs(:,Nc+1)
      eigenvecs(:,Nc+1) = eigenvecs(:,Nc+4)
      eigenvecs(:,Nc+4) = tmpMO
    end if

    if (fb < fc) then
      write(stdOut,'(A6,F9.6,A20,I4,A8,I4,A8)') " fb = ", fb, &
          & ", MO swap between b(", Nc+2, ") and c(", Nc+3, ") occurs"
      tmpMO(:) = eigenvecs(:,Nc+2)
      eigenvecs(:,Nc+2) = eigenvecs(:,Nc+3)
      eigenvecs(:,Nc+3) = tmpMO
    end if

  end subroutine MOswap44_


  !> Calculate filling for minimzed state with optimized FONs in REKS(2,2)
  subroutine getFilling22_(filling, SAweight, FONs, Efunction, Nc)

    !> occupations (level)
    real(dp), intent(out) :: filling(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    real(dp) :: n_a, n_b
    integer :: ii

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    filling(:) = 0.0_dp
    do ii = 1, Nc
      filling(ii) = 2.0_dp
    end do
    if (Efunction == 1) then
      ! REKS charge
      filling(Nc+1) = n_a
      filling(Nc+2) = n_b
    else if (Efunction == 2) then
      ! SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + SAweight(2)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2)
    end if

  end subroutine getFilling22_


  !> Calculate Fc and Fa from Hamiltonian of each microstate
  subroutine getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, hamSqrL, hamSpL, weight, fillingL, &
      & Nc, Na, Lpaired, isHybridXc, isRS_OnsCorr, orbFON, Fc, Fa)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> state-averaged occupation numbers
    real(dp), intent(inout) :: orbFON(:)

    !> dense fock matrix for core orbitals
    real(dp), intent(out) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(out) :: Fa(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: isHybridXc

    !> Whether to run onsite correction with range-separated functional
    logical, intent(in) :: isRS_OnsCorr

    real(dp), allocatable :: tmpHam(:,:)

    integer :: iL, Lmax, nOrb

    nOrb = size(Fc,dim=1)
    Lmax = size(weight,dim=1)

    if (.not. (isHybridXc .or. isRS_OnsCorr)) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    call fockFON_(fillingL, weight, orbFON)

    Fc(:,:) = 0.0_dp
    Fa(:,:,:) = 0.0_dp
    do iL = 1, Lmax

      if (.not. (isHybridXc .or. isRS_OnsCorr)) then
        tmpHam(:,:) = 0.0_dp
        ! convert from sparse to dense for hamSpL in AO basis
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call adjointLowerTriangle(tmpHam)
      end if

      ! compute the Fock operator with core, a, b orbitals in AO basis
      if (isHybridXc .or. isRS_OnsCorr) then
        call fockFcAO_(hamSqrL(:,:,1,iL), weight, Lpaired, iL, Fc)
        call fockFaAO_(hamSqrL(:,:,1,iL), weight, fillingL, orbFON, &
            & Nc, Na, Lpaired, iL, Fa)
      else
        call fockFcAO_(tmpHam, weight, Lpaired, iL, Fc)
        call fockFaAO_(tmpHam, weight, fillingL, orbFON, &
            & Nc, Na, Lpaired, iL, Fa)
      end if

    end do

  end subroutine getFockFcFa_


  !> Calculate pseudo-fock matrix from Fc and Fa
  subroutine getPseudoFock_(Fc, Fa, orbFON, Nc, Na, fock)

    !> dense pseudo-fock matrix
    real(dp), intent(out) :: fock(:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> state-averaged occupation numbers
    real(dp), intent(in) :: orbFON(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    real(dp) :: res
    integer :: ii, jj, ind1, ind2, nOrb

    nOrb = size(fock,dim=1)

    fock(:,:) = 0.0_dp
    do ii = 1, Nc
      do jj = ii, Nc
        fock(jj,ii) = Fc(ii,jj)
      end do
      do jj = Nc + 1, Nc + Na
        ind1 = jj - Nc
        call fockFijMO_(res, Fc(ii,jj), Fa(ii,jj,ind1), &
            & orbFON(ii), orbFON(jj))
        fock(jj,ii) = res
      end do
      do jj = Nc + Na + 1, nOrb
        fock(jj,ii) = Fc(ii,jj)
      end do
    end do

    do jj = Nc + Na + 1, nOrb
      do ii = Nc + 1, Nc + Na
        ind1 = ii - Nc
        call fockFijMO_(res, Fc(jj,ii), Fa(jj,ii,ind1), &
            & orbFON(jj), orbFON(ii))
        fock(jj,ii) = res
      end do
      do ii = jj, nOrb
        fock(ii,jj) = Fc(jj,ii)
      end do
    end do

    do ii = Nc + 1, Nc + Na
      ind1 = ii - Nc
      do jj = Nc + 1, Nc + Na
        ind2 = jj - Nc
        if (ii == jj) then
          fock(jj,ii) = Fa(ii,jj,ind1)
        else
          call fockFijMO_(res, Fa(ii,jj,ind1), Fa(ii,jj,ind2), &
              & orbFON(ii), orbFON(jj))
          fock(jj,ii) = res
        end if
      end do
    end do

    call adjointLowerTriangle(fock)

  end subroutine getPseudoFock_


  !> Avoid changing the order of MOs
  !> Required number of cycles increases as the number of shift increases
  subroutine levelShifting_(fock, shift, Nc, Na)

    !> dense pseudo-fock matrix
    real(dp), intent(inout) :: fock(:,:)

    !> Shift value in SCC cycle
    real(dp), intent(in) :: shift

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    integer :: nOrb, ind, ii

    nOrb = size(fock,dim=1)

    do ii = Nc + 1, Nc + Na
      ind = ii - Nc
      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
    end do

    do ii = Nc + Na + 1, nOrb
      ind = Na + 1
      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
    end do

  end subroutine levelShifting_


  !> Calculate state-averaged FONs
  subroutine fockFON_(fillingL, weight, orbFON)

    !> state-averaged occupation numbers
    real(dp), intent(out) :: orbFON(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    integer :: Lmax, iL

    Lmax = size(weight,dim=1)

    orbFON(:) = 0.0_dp
    do iL = 1, Lmax
      orbFON(:) = orbFON(:) + 0.5_dp * weight(iL) * &
          & ( fillingL(:,1,iL) + fillingL(:,2,iL) )
    end do

  end subroutine fockFON_


  !> Calculate fock matrix for core orbitals in AO basis
  subroutine fockFcAO_(hamSqr, weight, Lpaired, iL, Fc)

    !> dense fock matrix for core orbitals
    real(dp), intent(inout) :: Fc(:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current index in loop L
    integer, intent(in) :: iL

    if (iL <= Lpaired) then
      Fc(:,:) = Fc + 0.5_dp * hamSqr * &
          & ( weight(iL) + weight(iL) )
    else
      if (mod(iL,2) == 1) then
        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
            & ( weight(iL) + weight(iL+1) )
      else
        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
            & ( weight(iL) + weight(iL-1) )
      end if
    end if

  end subroutine fockFcAO_


  !> Calculate fock matrix for active orbitals in AO basis
  subroutine fockFaAO_(hamSqr, weight, fillingL, orbFON, Nc, Na, &
      & Lpaired, iL, Fa)

    !> dense fock matrix for active orbitals
    real(dp), intent(inout) :: Fa(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> state-averaged occupation numbers
    real(dp), intent(in) :: orbFON(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current index in loop L
    integer, intent(in) :: iL

    integer :: ind, ind_a

    do ind = 1, Na
      ind_a = Nc + ind
      if (iL <= Lpaired) then
        Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
            & ( weight(iL) + weight(iL) ) * hamSqr / orbFON(ind_a)
      else
        if (mod(iL,2) == 1) then
          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
              & ( weight(iL) + weight(iL+1) ) * hamSqr / orbFON(ind_a)
        else
          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
              & ( weight(iL) + weight(iL-1) ) * hamSqr / orbFON(ind_a)
        end if
      end if
    end do

  end subroutine fockFaAO_


  !> Calculate pseudo-fock off-diagonal element in MO basis
  subroutine fockFijMO_(res, fock_i, fock_j, f_i, f_j)

    !> temporary pseudo-fock value
    real(dp), intent(out) :: res

    !> temporary Fc or Fa values
    real(dp), intent(in) :: fock_i, fock_j

    !> temporary orbFON values
    real(dp), intent(in) :: f_i, f_j

    real(dp) :: eps = 1.0E-3_dp

    res = 0.0_dp
    if (abs(f_j-f_i) .LT. eps) then
      if (f_j >= f_i) then
        res = ( f_j*fock_j - f_i*fock_i )
      else
        res = -( f_j*fock_j - f_i*fock_i )
      end if
    else
      res = ( f_j*fock_j - f_i*fock_i ) / (f_j - f_i)
    end if

  end subroutine fockFijMO_


  !> Calculate converged Lagrangian values
  subroutine getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, hamSqrL, hamSpL, weight, &
      & fillingL, Nc, Na, Lpaired, isHybridXc, isRS_OnsCorr, Wab)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: isHybridXc

    !> Whether to run onsite correction with range-separated functional
    logical, intent(in) :: isRS_OnsCorr

    !> converged Lagrangian values within active space
    real(dp), intent(out) :: Wab(:,:)

    real(dp), allocatable :: tmpHam(:,:)
    real(dp), allocatable :: tmpHamL(:,:,:)

    integer :: nOrb, iL, Lmax
    integer :: ia, ib, ist, nActPair

    nOrb = size(eigenvecs,dim=1)
    Lmax = size(fillingL,dim=3)
    nActPair = Na * (Na - 1) / 2

    if (.not. (isHybridXc .or. isRS_OnsCorr)) then
      allocate(tmpHam(nOrb,nOrb))
    end if
    allocate(tmpHamL(nActPair,1,Lmax))

    tmpHamL(:,:,:) = 0.0_dp
    do ist = 1, nActPair

      call getTwoIndices(Na, ist, ia, ib, 1)

      do iL = 1, Lmax

        if (isHybridXc .or. isRS_OnsCorr) then
          ! convert hamSqrL from AO basis to MO basis
          ! hamSqrL has (my_ud) component
          if (ist == 1) then
            call matAO2MO(hamSqrL(:,:,1,iL), eigenvecs)
          end if
          tmpHamL(ist,1,iL) = hamSqrL(Nc+ia,Nc+ib,1,iL)
        else
          tmpHam(:,:) = 0.0_dp
          ! convert from sparse to dense for hamSpL in AO basis
          ! hamSpL has (my_ud) component
          call env%globalTimer%startTimer(globalTimers%sparseToDense)
          call unpackHS(tmpHam, hamSpL(:,1,iL), &
              & neighbourList%iNeighbour, nNeighbourSK, &
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
          call env%globalTimer%stopTimer(globalTimers%sparseToDense)
          call adjointLowerTriangle(tmpHam)
          ! convert tmpHam from AO basis to MO basis
          call matAO2MO(tmpHam, eigenvecs)
          ! save F_{L,ab}^{\sigma} in MO basis
          tmpHamL(ist,1,iL) = tmpHam(Nc+ia,Nc+ib)
        end if

      end do

      ! calculate the Lagrangian eps_{ab} and state-interaction term
      Wab(ist,1) = 0.0_dp
      Wab(ist,2) = 0.0_dp
      do iL = 1, Lmax
        if (iL <= Lpaired) then
          Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
          Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
        else
          if (mod(iL,2) == 1) then
            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
          else
            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
          end if
        end if
      end do

    end do

  end subroutine getLagrangians_


  !> calculate state-interaction terms between SA-REKS states in (2,2) case
  subroutine getStateCoup22_(Wab, FONs, StateCoup)

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(out) :: StateCoup(:,:)

    real(dp) :: n_a, n_b
    integer :: nstates

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    nstates = size(StateCoup,dim=1)

    StateCoup(:,:) = 0.0_dp
    StateCoup(1,2) = sqrt(n_a) * Wab(1,1) - sqrt(n_b) * Wab(1,1)
    StateCoup(2,1) = StateCoup(1,2)
    if (nstates == 3) then
      StateCoup(2,3) = sqrt(n_a) * Wab(1,1) + sqrt(n_b) * Wab(1,1)
      StateCoup(3,2) = StateCoup(2,3)
    end if

  end subroutine getStateCoup22_


end module dftbp_reks_reksen
