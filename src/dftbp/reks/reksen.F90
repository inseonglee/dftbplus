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
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_elecsolvers_elecsolvers, only: TElectronicSolver
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_reks_rekscommon, only : getFactor, getTwoIndices, matAO2MO,&
      & getFullLongRangePars
  use dftbp_reks_reksio, only : printReksSSRInfo
  use dftbp_reks_reksvar, only : TReksCalc, reksTypes
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_orbitals, only : TOrbitals

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
      call getFilling44_(filling, this%SAweight, this%FONs, this%Efunction, this%Nc)
    end select

  end subroutine getFilling


  !> Calculate the energy of SA-REKS states and averaged state
  subroutine calcSaReksEnergy(this, energy, tConverged)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Energy terms in the system
    type(TEnergies), intent(inout) :: energy

    !> Has the calculation converged?
    logical, intent(in) :: tConverged

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

    if (.not. tConverged) then
      ! In this step Eavg becomes the energy of averaged state
      ! From this energy we can check the variational principle
      energy%Eavg = 0.0_dp
      do ist = 1, this%SAstates
        energy%Eavg = energy%Eavg + this%SAweight(ist) * this%energy(ist)
      end do
    end if

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
  subroutine solveSecularEqn(env, denseDesc, hybridXc, orb, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, species, electronicSolver, &
      & eigenvecs, spinW, onSiteElements, this)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> spin constants
    real(dp), intent(in) :: spinW(:,:,:)

    !> Correction to energy from on-site matrix elements
    real(dp), intent(in) :: onSiteElements(:,:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: Wab(:,:)
    real(dp), allocatable :: ERI(:,:)
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
      allocate(ERI(nActPair,nActPair))
      call getERI_(env, denseDesc, hybridXc, orb, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, species, eigenvecs(:,:,1), spinW, &
          & onSiteElements, this%hamSqrL, this%hamSpL, this%overSqr, this%weight, &
          & this%fillingL, this%getAtomIndex, this%Nc, this%Na, this%Lpaired, &
          & this%isOnsite, this%isHybridXc, this%isRS_OnsCorr, ERI)
      call getStateCoup44_(Wab, ERI, this%enLtot, this%FONs, this%Efunction, &
          & this%tAllStates, StateCoup)
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
    call printReksSSRInfo(this, Wab, ERI, tmpEn, StateCoup)

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
        weight(iL) = weight(iL) + SAweight(ist) * weightL(ist,iL)
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
    if (.not. tConverged) then
      weight(:) = 0.0_dp
      do iL = 1, Lmax
        do ist = 1, SAstates
          weight(iL) = weight(iL) + SAweight(ist) * weightL(ist,iL)
        end do
      end do
    end if

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
      ! 2SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + SAweight(2)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2)
    end if

  end subroutine getFilling22_


  !> Calculate filling for minimzed state with optimized FONs in REKS(4,4)
  subroutine getFilling44_(filling, SAweight, FONs, Efunction, Nc)

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

    real(dp) :: n_a, n_b, n_c, n_d
    real(dp) :: np_a, np_b, np_c, np_d
    real(dp) :: mp_a, mp_b, mp_c, mp_d
    real(dp) :: fa, fb, fc, fd
    integer :: ii

    n_a = FONs(1,1); n_b = FONs(2,1); n_c = FONs(3,1); n_d = FONs(4,1)
    np_a = FONs(1,2); np_b = FONs(2,2); np_c = FONs(3,2); np_d = FONs(4,2)
    mp_a = FONs(1,3); mp_b = FONs(2,3); mp_c = FONs(3,3); mp_d = FONs(4,3)

    filling(:) = 0.0_dp
    do ii = 1, Nc
      filling(ii) = 2.0_dp
    end do

    if (Efunction == 1) then
      ! REKS charge
      filling(Nc+1) = n_a
      filling(Nc+2) = n_b
      filling(Nc+3) = n_c
      filling(Nc+4) = n_d
    else if (Efunction == 2) then
      ! 2SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + SAweight(2)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2)
      filling(Nc+3) = n_c*SAweight(1) + SAweight(2)
      filling(Nc+4) = n_d*SAweight(1) + SAweight(2)
    else if (Efunction == 3) then
      ! 3SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + np_a*SAweight(2) + SAweight(3)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2) + np_b*SAweight(3)
      filling(Nc+3) = n_c*SAweight(1) + SAweight(2) + np_c*SAweight(3)
      filling(Nc+4) = n_d*SAweight(1) + np_d*SAweight(2) + SAweight(3)
    else if (Efunction == 4) then
      ! 5SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + np_a*SAweight(2) + SAweight(3) &
          &+ mp_a*SAweight(4) + SAweight(5)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2) + np_b*SAweight(3) &
          &+ SAweight(4) + mp_b*SAweight(5)
      filling(Nc+3) = n_c*SAweight(1) + SAweight(2) + np_c*SAweight(3) &
          &+ mp_c*SAweight(4) + SAweight(5)
      filling(Nc+4) = n_d*SAweight(1) + np_d*SAweight(2) + SAweight(3) &
          &+ SAweight(4) + mp_d*SAweight(5)
    end if

  end subroutine getFilling44_


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


  !> Calculate electron repulsion integrals (ERI) for 3 or 4 indices
  subroutine getERI_(env, denseDesc, hybridXc, orb, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, species, eigenvecs, spinW, onSiteElements, &
      & hamSqrL, hamSpL, overSqr, weight, fillingL, getAtomIndex, Nc, Na, Lpaired, &
      & isOnsite, isHybridXc, isRS_OnsCorr, ERI)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:)

    !> spin constants
    real(dp), intent(in) :: spinW(:,:,:)

    !> Correction to energy from on-site matrix elements
    real(dp), intent(in) :: onSiteElements(:,:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Are on-site corrections being used?
    logical, intent(in) :: isOnsite

    !> Whether to run a range separated calculation
    logical, intent(in) :: isHybridXc

    !> Whether to run onsite correction with range-separated functional
    logical, intent(in) :: isRS_OnsCorr

    !> electron repulsion integrals within active space (4,4); ab ac ad bc bd cd
    real(dp), intent(inout) :: ERI(:,:)

    type(TScc), allocatable :: sccCalc ! never allocated

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: auxRmat(:,:)
    real(dp), allocatable :: tmpHamL(:,:,:)

    real(dp), allocatable :: auxRhalf_fr(:,:)
    real(dp), allocatable :: auxZhalf_fr(:,:)
    real(dp), allocatable :: tmpHxc_fr(:)

    real(dp), allocatable :: auxRhalf_lr(:,:)
    real(dp), allocatable :: auxZhalf_lr(:,:)
    real(dp), allocatable :: tmpHxc_lr(:)

    real(dp), allocatable :: SpinAO(:,:)
    real(dp), allocatable :: LrGammaAO(:,:)
    real(dp), allocatable :: OnsiteAO(:,:,:,:)
    real(dp), allocatable :: LrOnsiteAO(:,:,:)

    ! Following three variables are not used, just defined
    real(dp), allocatable :: GammaAO(:,:)
    real(dp), allocatable :: GammaDeriv(:,:,:)
    real(dp), allocatable :: LrGammaDeriv(:,:,:)

    real(dp) :: tmp22!, tmp11
    integer :: nOrbHalf, nOrb, nSpin, nActPair, Lmax
    integer :: ii, jj, mu, nu, tau, gam
    integer :: ia, ib, ist, ja, jb, jst, kst, iL

    real(dp) :: tmpS1, tmpS2, tmpS3, tmpS4
    real(dp) :: tmpL1, tmpL2, tmpL3, tmpL4

    real(dp) :: tmpO1s, tmpO2s, tmpO3s, tmpO4s
    real(dp) :: tmpO1d, tmpO2d, tmpO3d, tmpO4d
    real(dp) :: tmpvalue1s, tmpvalue1d
    integer :: kap

    nSpin = size(fillingL,dim=2)
    Lmax = size(fillingL,dim=3)
    nOrb = size(eigenvecs,dim=1)
    nOrbHalf = nOrb * (nOrb + 1) / 2
    nActPair = Na * (Na - 1) / 2

    allocate(tmpMat(nOrb,nOrb))
    allocate(auxRmat(nOrb,nOrb))
    if (.not. (isHybridXc .or. isRS_OnsCorr)) then
      allocate(tmpHamL(nOrb,nOrb,Lmax))
    end if

    allocate(auxRhalf_fr(nOrbHalf,3))
    allocate(auxZhalf_fr(nOrbHalf,3))
    allocate(tmpHxc_fr(nOrbHalf))

    if (isHybridXc .or. isRS_OnsCorr) then
      allocate(auxRhalf_lr(nOrbHalf,3))
      allocate(auxZhalf_lr(nOrbHalf,3))
      allocate(tmpHxc_lr(nOrbHalf))
    end if

    ! Here we use temporary variables (e.g. SpinAO) instead of original
    ! defined variable (e.g. this%SpinAO) because they can be used for
    ! gradient calculation
    allocate(SpinAO(nOrb,nOrb))
    if (isOnsite) then
      allocate(OnsiteAO(nOrb,nOrb,nSpin,2))
    end if
    if (isHybridXc) then
      allocate(LrGammaAO(nOrb,nOrb))
    end if
    if (isRS_OnsCorr) then
      allocate(LrOnsiteAO(nOrb,nOrb,2))
    end if

    ! get spinW, lr-gamma, on-site constants
    call getFullLongRangePars(env, sccCalc, hybridXc, orb, species,&
        & neighbourList%iNeighbour, img2CentCell, denseDesc%iAtomStart,&
        & spinW, onSiteElements, getAtomIndex, isOnsite, isHybridXc,&
        & isRS_OnsCorr, GammaAO, GammaDeriv, SpinAO, OnsiteAO,&
        & LrGammaAO, LrGammaDeriv, LrOnsiteAO, optionERI=.true.)

    ! Calculate auxiliary R matrices with half dense form
    kst = 1
    auxRhalf_fr(:,:) = 0.0_dp
    if (isHybridXc .or. isRS_OnsCorr) then
      auxRhalf_lr(:,:) = 0.0_dp
    end if

    do ist = 1, nActPair
      call getTwoIndices(Na, ist, ia, ib, 1)
      do jst = 1, nActPair
        call getTwoIndices(Na, jst, ja, jb, 1)

        if (ist + jst == nActPair + 1) then
          if (ist < jst) then
            cycle
          end if

          ! kst = 1; ist = 4(bc), jst = 3(ad)
          ! kst = 2; ist = 5(bd), jst = 2(ac)
          ! kst = 3; ist = 6(cd), jst = 1(ab)
          auxRmat(:,:) = 0.0_dp
          do mu = 1, nOrb
            do nu = 1, nOrb
              auxRmat(mu,nu) = eigenvecs(mu,Nc+ja) * eigenvecs(nu,Nc+ib)
            end do
          end do

          tmpMat(:,:) = auxRmat(:,:) + transpose(auxRmat(:,:))
          do mu = 1, nOrb
            tmpMat(mu,mu) = auxRmat(mu,mu)
          end do
          do ii = 1, nOrbHalf
            call getTwoIndices(nOrb, ii, mu, nu, 2)
            auxRhalf_fr(ii,kst) = tmpMat(mu,nu)
          end do

          if (isHybridXc .or. isRS_OnsCorr) then

            auxRmat(:,:) = 0.0_dp
            do mu = 1, nOrb
              do nu = 1, nOrb
                auxRmat(mu,nu) = eigenvecs(mu,Nc+ia) * eigenvecs(nu,Nc+ib)
              end do
            end do

            tmpMat(:,:) = auxRmat(:,:) + transpose(auxRmat(:,:))
            do mu = 1, nOrb
              tmpMat(mu,mu) = auxRmat(mu,mu)
            end do
            do ii = 1, nOrbHalf
              call getTwoIndices(nOrb, ii, mu, nu, 2)
              auxRhalf_lr(ii,kst) = tmpMat(mu,nu)
            end do

          end if

          kst = kst + 1
        end if

      end do
    end do
    deallocate(auxRmat)

    ! Calculate auxiliary Z matrices from Hxc kernel and eigenvectors
    ! Note that the full-range and long-range parts are treated separately
    ! and the full-range gamma does not contribute the ERI term

    ! zeroing for ZmatL
    auxZhalf_fr(:,:) = 0.0_dp
    if (isHybridXc .or. isRS_OnsCorr) then
      auxZhalf_lr(:,:) = 0.0_dp
    end if
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mu,nu,tau,gam,tmp22, &
!$OMP& tmpHxc_fr,tmpHxc_lr,tmpS1,tmpS2,tmpS3,tmpS4,tmpvalue1s, &
!$OMP& tmpvalue1d,tmpL1,tmpL2,tmpL3,tmpL4,tmpO1s,tmpO1d,tmpO2s, &
!$OMP& tmpO2d,tmpO3s,tmpO3d,tmpO4s,tmpO4d) SCHEDULE(RUNTIME)
    do jj = 1, nOrbHalf

      call getTwoIndices(nOrb, jj, tau, gam, 2)

      tmpHxc_fr(:) = 0.0_dp
      if (isHybridXc .or. isRS_OnsCorr) then
        tmpHxc_lr(:) = 0.0_dp
      end if
      do ii = 1, nOrbHalf

        !call getTwoIndices(nOrb, jj, tau, gam, 2)
        tmp22 = ( real(2.0_dp*nOrb+3.0_dp, dp) - sqrt( (2.0_dp*nOrb+ &
            & 3.0_dp)**2.0_dp - 8.0_dp*(nOrb+ii) ) )/2.0_dp
        mu = int( real(tmp22, dp) )
        nu = mu**2/2 - mu/2 - nOrb*mu + nOrb + ii

        ! spin term
        tmpS1 = SpinAO(mu,tau)
        tmpS2 = SpinAO(nu,tau)
        tmpS3 = SpinAO(mu,gam)
        tmpS4 = SpinAO(nu,gam)

        tmpHxc_fr(ii) = tmpHxc_fr(ii) - 0.5_dp * overSqr(mu,nu) * &
            & overSqr(tau,gam) * (tmpS1+tmpS2+tmpS3+tmpS4)

        ! LC term
        if (isHybridXc) then

          tmpL1 = LrGammaAO(mu,tau)
          tmpL2 = LrGammaAO(mu,gam)
          tmpL3 = LrGammaAO(nu,tau)
          tmpL4 = LrGammaAO(nu,gam)
          tmpvalue1s = 0.25_dp * overSqr(mu,nu) * &
              & overSqr(tau,gam) * (tmpL1+tmpL2+tmpL3+tmpL4)

          tmpHxc_lr(ii) = tmpHxc_lr(ii) + tmpvalue1s

        end if

        ! full-range onsite terms
        if (isOnsite) then

          tmpO1s = OnsiteAO(mu,tau,1,1)
          tmpO2s = OnsiteAO(nu,gam,1,1)
          tmpO3s = OnsiteAO(mu,gam,1,1)
          tmpO4s = OnsiteAO(nu,tau,1,1)

          tmpvalue1s = overSqr(mu,gam) * overSqr(nu,tau) * (tmpO1s + tmpO2s) &
              & + overSqr(mu,tau) * overSqr(nu,gam) * (tmpO3s + tmpO4s)

          tmpO1d = OnsiteAO(mu,tau,2,1)
          tmpO2d = OnsiteAO(nu,gam,2,1)
          tmpO3d = OnsiteAO(mu,gam,2,1)
          tmpO4d = OnsiteAO(nu,tau,2,1)

          tmpvalue1d = overSqr(mu,gam) * overSqr(nu,tau) * (tmpO1d + tmpO2d) &
              & + overSqr(mu,tau) * overSqr(nu,gam) * (tmpO3d + tmpO4d)

          tmpHxc_fr(ii) = tmpHxc_fr(ii) + 0.25_dp * (tmpvalue1d - tmpvalue1s)

          tmpvalue1s = 0.0_dp
          tmpvalue1d = 0.0_dp
          if (mu == tau) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,nu) * overSqr(kap,gam) &
                  & * OnsiteAO(mu,kap,1,1)
              tmpvalue1d = tmpvalue1d + overSqr(kap,nu) * overSqr(kap,gam) &
                  & * OnsiteAO(mu,kap,2,1)
            end do
          end if
          if (mu == gam) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,nu) * overSqr(kap,tau) &
                  & * OnsiteAO(mu,kap,1,1)
              tmpvalue1d = tmpvalue1d + overSqr(kap,nu) * overSqr(kap,tau) &
                  & * OnsiteAO(mu,kap,2,1)
            end do
          end if
          if (nu == tau) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,mu) * overSqr(kap,gam) &
                  & * OnsiteAO(nu,kap,1,1)
              tmpvalue1d = tmpvalue1d + overSqr(kap,mu) * overSqr(kap,gam) &
                  & * OnsiteAO(nu,kap,2,1)
            end do
          end if
          if (nu == gam) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,mu) * overSqr(kap,tau) &
                  & * OnsiteAO(nu,kap,1,1)
              tmpvalue1d = tmpvalue1d + overSqr(kap,mu) * overSqr(kap,tau) &
                  & * OnsiteAO(nu,kap,2,1)
            end do
          end if

          tmpHxc_fr(ii) = tmpHxc_fr(ii) + 0.25_dp * (tmpvalue1d - tmpvalue1s)

          tmpO1s = OnsiteAO(mu,tau,1,2)
          tmpO2s = OnsiteAO(nu,tau,1,2)
          tmpO3s = OnsiteAO(mu,gam,1,2)
          tmpO4s = OnsiteAO(nu,gam,1,2)

          tmpO1d = OnsiteAO(mu,tau,2,2)
          tmpO2d = OnsiteAO(nu,tau,2,2)
          tmpO3d = OnsiteAO(mu,gam,2,2)
          tmpO4d = OnsiteAO(nu,gam,2,2)

          tmpvalue1s = overSqr(mu,nu) * overSqr(tau,gam) * &
              & (tmpO1s + tmpO2s + tmpO3s + tmpO4s)
          tmpvalue1d = overSqr(mu,nu) * overSqr(tau,gam) * &
              & (tmpO1d + tmpO2d + tmpO3d + tmpO4d)

          tmpHxc_fr(ii) = tmpHxc_fr(ii) - (tmpvalue1d - tmpvalue1s)

        end if

        ! long-range onsite terms
        if (isRS_OnsCorr) then

          tmpO1s = LrOnsiteAO(mu,tau,1)
          tmpO2s = LrOnsiteAO(nu,gam,1)
          tmpO3s = LrOnsiteAO(mu,gam,1)
          tmpO4s = LrOnsiteAO(nu,tau,1)

          tmpvalue1s = overSqr(mu,gam) * overSqr(nu,tau) * (tmpO1s + tmpO2s) &
              & + overSqr(mu,tau) * overSqr(nu,gam) * (tmpO3s + tmpO4s)

          tmpHxc_lr(ii) = tmpHxc_lr(ii) + 0.25_dp * tmpvalue1s

          tmpvalue1s = 0.0_dp
          if (mu == tau) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,nu) * overSqr(kap,gam) &
                  & * LrOnsiteAO(mu,kap,1)
            end do
          end if
          if (mu == gam) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,nu) * overSqr(kap,tau) &
                  & * LrOnsiteAO(mu,kap,1)
            end do
          end if
          if (nu == tau) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,mu) * overSqr(kap,gam) &
                  & * LrOnsiteAO(nu,kap,1)
            end do
          end if
          if (nu == gam) then
            do kap = 1, nOrb
              tmpvalue1s = tmpvalue1s + overSqr(kap,mu) * overSqr(kap,tau) &
                  & * LrOnsiteAO(nu,kap,1)
            end do
          end if

          tmpHxc_lr(ii) = tmpHxc_lr(ii) + 0.25_dp * tmpvalue1s

          tmpO1s = LrOnsiteAO(mu,tau,2)
          tmpO2s = LrOnsiteAO(mu,gam,2)
          tmpO3s = LrOnsiteAO(nu,tau,2)
          tmpO4s = LrOnsiteAO(nu,gam,2)

          tmpvalue1s = overSqr(mu,nu) * overSqr(tau,gam) * &
              & (tmpO1s + tmpO2s + tmpO3s + tmpO4s)

          tmpHxc_lr(ii) = tmpHxc_lr(ii) - tmpvalue1s

        end if

      end do

      ! calculate the ZmatL
      do kst = 1, 3
        auxZhalf_fr(jj,kst) = sum(auxRhalf_fr(:,kst)*tmpHxc_fr)
        if (isHybridXc .or. isRS_OnsCorr) then
          auxZhalf_lr(jj,kst) = sum(auxRhalf_lr(:,kst)*tmpHxc_lr)
        end if
      end do

    end do
!$OMP END PARALLEL DO

    ! Calculate ERIs from auxiliary Z matrices and eigenvectors
    kst = 1
    ERI(:,:) = 0.0_dp
    do ist = 1, nActPair
      call getTwoIndices(Na, ist, ia, ib, 1)
      do jst = 1, nActPair
        call getTwoIndices(Na, jst, ja, jb, 1)

        if (ist + jst == 7) then
          ! This condition is 4eri
          if (ist < jst) then
            cycle
          end if
        else
          ! This condition is 3eri
          cycle
        end if

        do jj = 1, nOrbHalf
          call getTwoIndices(nOrb, jj, tau, gam, 2)

          ERI(ist,jst) = ERI(ist,jst) + auxZhalf_fr(jj,kst) * eigenvecs(tau,Nc+ia) * eigenvecs(gam,Nc+jb)
          if (isHybridXc .or. isRS_OnsCorr) then
            ERI(ist,jst) = ERI(ist,jst) + auxZhalf_lr(jj,kst) * eigenvecs(tau,Nc+ja) * eigenvecs(gam,Nc+jb)
          end if
          if (tau /= gam) then
            ERI(ist,jst) = ERI(ist,jst) + auxZhalf_fr(jj,kst) * eigenvecs(tau,Nc+jb) * eigenvecs(gam,Nc+ia)
            if (isHybridXc .or. isRS_OnsCorr) then
              ERI(ist,jst) = ERI(ist,jst) + auxZhalf_lr(jj,kst) * eigenvecs(tau,Nc+jb) * eigenvecs(gam,Nc+ja)
            end if
          end if
        end do
        kst = kst + 1

      enddo
    enddo

    ! 3ERIs can be computed from Hamiltonians for each microstate
    if (.not. (isHybridXc .or. isRS_OnsCorr)) then
      do iL = 1, Lmax
        tmpMat(:,:) = 0.0_dp
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpMat, hamSpL(:,1,iL), &
            & neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call adjointLowerTriangle(tmpMat)
        ! convert the hamiltonians from AO basis to MO basis
        call matAO2MO(tmpMat, eigenvecs)
        tmpHamL(:,:,iL) = tmpMat(:,:)
      end do
      ERI(2,1) = tmpHamL(Nc+3, Nc+2, 22) - tmpHamL(Nc+3, Nc+2, 23)
      ERI(3,1) = tmpHamL(Nc+4, Nc+2, 6) - tmpHamL(Nc+4, Nc+2, 7)
      ERI(4,1) = tmpHamL(Nc+3, Nc+1, 10) - tmpHamL(Nc+3, Nc+1, 11)
      ERI(5,1) = tmpHamL(Nc+4, Nc+1, 18) - tmpHamL(Nc+4, Nc+1, 19)
      ERI(3,2) = tmpHamL(Nc+3, Nc+4, 22) - tmpHamL(Nc+3, Nc+4, 23)
!      ERI(3,2) = tmpHamL(Nc+4, Nc+3, 6) - tmpHamL(Nc+4, Nc+3, 7)
      ERI(4,2) = tmpHamL(Nc+1, Nc+2, 21) - tmpHamL(Nc+1, Nc+2, 23)
!      ERI(4,2) = tmpHamL(Nc+2, Nc+1, 9) - tmpHamL(Nc+2, Nc+1, 11)
      ERI(6,2) = tmpHamL(Nc+1, Nc+4, 21) - tmpHamL(Nc+1, Nc+4, 23)
      ERI(5,3) = tmpHamL(Nc+1, Nc+2, 5) - tmpHamL(Nc+1, Nc+2, 7)
!      ERI(5,3) = tmpHamL(Nc+2, Nc+1, 17) - tmpHamL(Nc+2, Nc+1, 19)
      ERI(6,3) = tmpHamL(Nc+1, Nc+3, 5) - tmpHamL(Nc+1, Nc+3, 7)
      ERI(5,4) = tmpHamL(Nc+3, Nc+4, 10) - tmpHamL(Nc+3, Nc+4, 11)
!      ERI(5,4) = tmpHamL(Nc+4, Nc+3, 18) - tmpHamL(Nc+4, Nc+3, 19)
      ERI(6,4) = tmpHamL(Nc+2, Nc+4, 9) - tmpHamL(Nc+2, Nc+4, 11)
      ERI(6,5) = tmpHamL(Nc+2, Nc+3, 17) - tmpHamL(Nc+2, Nc+3, 19)
    else
      ERI(2,1) = hamSqrL(Nc+3, Nc+2, 1, 22) - hamSqrL(Nc+3, Nc+2, 1, 23)
      ERI(3,1) = hamSqrL(Nc+4, Nc+2, 1, 6) - hamSqrL(Nc+4, Nc+2, 1, 7)
      ERI(4,1) = hamSqrL(Nc+3, Nc+1, 1, 10) - hamSqrL(Nc+3, Nc+1, 1, 11)
      ERI(5,1) = hamSqrL(Nc+4, Nc+1, 1, 18) - hamSqrL(Nc+4, Nc+1, 1, 19)
      ERI(3,2) = hamSqrL(Nc+3, Nc+4, 1, 22) - hamSqrL(Nc+3, Nc+4, 1, 23)
!      ERI(3,2) = hamSqrL(Nc+4, Nc+3, 1, 6) - hamSqrL(Nc+4, Nc+3, 1, 7)
      ERI(4,2) = hamSqrL(Nc+1, Nc+2, 1, 21) - hamSqrL(Nc+1, Nc+2, 1, 23)
!      ERI(4,2) = hamSqrL(Nc+2, Nc+1, 1, 9) - hamSqrL(Nc+2, Nc+1, 1, 11)
      ERI(6,2) = hamSqrL(Nc+1, Nc+4, 1, 21) - hamSqrL(Nc+1, Nc+4, 1, 23)
      ERI(5,3) = hamSqrL(Nc+1, Nc+2, 1, 5) - hamSqrL(Nc+1, Nc+2, 1, 7)
!      ERI(5,3) = hamSqrL(Nc+2, Nc+1, 1, 17) - hamSqrL(Nc+2, Nc+1, 1, 19)
      ERI(6,3) = hamSqrL(Nc+1, Nc+3, 1, 5) - hamSqrL(Nc+1, Nc+3, 1, 7)
      ERI(5,4) = hamSqrL(Nc+3, Nc+4, 1, 10) - hamSqrL(Nc+3, Nc+4, 1, 11)
!      ERI(5,4) = hamSqrL(Nc+4, Nc+3, 1, 18) - hamSqrL(Nc+4, Nc+3, 1, 19)
      ERI(6,4) = hamSqrL(Nc+2, Nc+4, 1, 9) - hamSqrL(Nc+2, Nc+4, 1, 11)
      ERI(6,5) = hamSqrL(Nc+2, Nc+3, 1, 17) - hamSqrL(Nc+2, Nc+3, 1, 19)
    end if

    call adjointLowerTriangle(ERI)

  end subroutine getERI_


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

    ! Initialize SI terms
    StateCoup(:,:) = 0.0_dp

    ! <PPS|H|OSS>
    StateCoup(1,2) = sqrt(n_a) * Wab(1,1) - sqrt(n_b) * Wab(1,1)
    StateCoup(2,1) = StateCoup(1,2)
    if (nstates == 3) then
      ! <OSS|H|DES>
      StateCoup(2,3) = sqrt(n_a) * Wab(1,1) + sqrt(n_b) * Wab(1,1)
      StateCoup(3,2) = StateCoup(2,3)
    end if

  end subroutine getStateCoup22_


  !> calculate state-interaction terms between SA-REKS states in (4,4) case
  subroutine getStateCoup44_(Wab, ERI, enLtot, FONs, Efunction, tAllStates, StateCoup)

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> electron repulsion integrals within active space (4,4); ab ac ad bc bd cd
    real(dp), intent(in) :: ERI(:,:)

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Decide the energy states in SA-REKS
    logical, intent(in) :: tAllStates

    !> state-interaction term between SA-REKS states
    real(dp), intent(inout) :: StateCoup(:,:)

    real(dp) :: n_a, n_b, n_c, n_d
    real(dp) :: np_a, np_b, np_c, np_d
    real(dp) :: mp_a, mp_b, mp_c, mp_d
    integer :: indPPS, indOSS1, indOSS2
    integer :: indOSS3, indOSS4, indDOSS
    integer :: indDSPS, indDES1, indDES2

    indPPS = 0; indOSS1 = 0; indOSS2 = 0
    indOSS3 = 0; indOSS4 = 0; indDOSS = 0
    indDSPS = 0; indDES1 = 0; indDES2 = 0

    if (tAllStates) then
      indPPS = 1; indOSS1 = 2; indOSS2 = 3
      indOSS3 = 4; indOSS4 = 5; indDOSS = 6
      indDSPS = 7; indDES1 = 8; indDES2 = 9
    else
      indPPS = 1
      if (Efunction == 2) then
        indDSPS = 2
      else if (Efunction == 3 .or. Efunction == 4) then
        indOSS1 = 2; indOSS2 = 3
        if (Efunction == 4) then
          indOSS3 = 4; indOSS4 = 5
        end if
      end if
    end if

    n_a = FONs(1,1); n_b = FONs(2,1); n_c = FONs(3,1); n_d = FONs(4,1)
    np_a = FONs(1,2); np_b = FONs(2,2); np_c = FONs(3,2); np_d = FONs(4,2)
    mp_a = FONs(1,3); mp_b = FONs(2,3); mp_c = FONs(3,3); mp_d = FONs(4,3)

    ! Initialize SI terms
    StateCoup(:,:) = 0.0_dp

    if (indPPS > 0) then

      if (indOSS1 > 0) then
        ! <PPS|H|OSS1>
        StateCoup(indOSS1,indPPS) = 0.5_dp * (dsqrt(n_a*np_a) + dsqrt(n_d*np_d)) &
            & * (dsqrt(n_b) - dsqrt(n_c)) * Wab(4,1)
      end if

      if (indOSS2 > 0) then
        ! <PPS|H|OSS2>
        StateCoup(indOSS2,indPPS) = 0.5_dp * (dsqrt(n_b*np_b) + dsqrt(n_c*np_c)) &
            & * (dsqrt(n_a) - dsqrt(n_d)) * Wab(3,1)
      end if

      if (indOSS3 > 0) then
        ! <PPS|H|OSS3>
        StateCoup(indOSS3,indPPS) = 0.5_dp * dsqrt(mp_a) * ( dsqrt(n_a*n_b) * Wab(5,1) &
            & - dsqrt(n_a*n_c) * ERI(4,6) + dsqrt(n_b*n_d) * ERI(1,3) ) &
            & - 0.5_dp * dsqrt(mp_c) * ( dsqrt(n_c*n_d) * Wab(5,1) &
            & - dsqrt(n_a*n_c) * ERI(1,3) + dsqrt(n_b*n_d) * ERI(4,6) )
      end if

      if (indOSS4 > 0) then
        ! <PPS|H|OSS4>
        StateCoup(indOSS4,indPPS) = 0.5_dp * dsqrt(mp_b) * ( dsqrt(n_a*n_b) * Wab(2,1) &
            & + dsqrt(n_a*n_c) * ERI(1,4) - dsqrt(n_b*n_d) * ERI(3,6) ) &
            & - 0.5_dp * dsqrt(mp_d) * ( dsqrt(n_c*n_d) * Wab(2,1) &
            & + dsqrt(n_a*n_c) * ERI(3,6) - dsqrt(n_b*n_d) * ERI(1,4) )
      end if

      if (indDOSS > 0) then
        ! <PPS|H|DOSS>
        StateCoup(indDOSS,indPPS) = 0.5_dp * (dsqrt(n_a*n_b) + dsqrt(n_c*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(2,5)) &
            & - 0.5_dp * (dsqrt(n_a*n_c) + dsqrt(n_b*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(1,6))
      end if

      if (indDSPS > 0) then
        ! <PPS|H|DSPS>
        StateCoup(indDSPS,indPPS) = 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_b) + dsqrt(n_c*n_d)) * ERI(2,5) &
            & + 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_c) + dsqrt(n_b*n_d)) * ERI(1,6)
      end if

    end if

    if (indOSS1 > 0) then

      if (indOSS2 > 0) then
        ! <OSS1|H|OSS2>
        StateCoup(indOSS2,indOSS1) = 0.5_dp * (dsqrt(np_a*np_b) + dsqrt(np_c*np_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(1,6)) &
            & - 0.5_dp * (dsqrt(np_a*np_c) + dsqrt(np_b*np_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(2,5))
      end if

      if (indOSS3 > 0) then
        ! <OSS1|H|OSS3>
        StateCoup(indOSS3,indOSS1) = 0.5_dp * dsqrt(mp_a) * ( dsqrt(np_a) * Wab(6,1) &
            & + dsqrt(np_a) * ERI(4,5) + dsqrt(np_d) * ERI(2,3) ) &
            & - 0.5_dp * dsqrt(mp_c) * ( dsqrt(np_d) * Wab(6,1) &
            & - dsqrt(np_d) * ERI(4,5) - dsqrt(np_a) * ERI(2,3) )
      end if

      if (indOSS4 > 0) then
        ! <OSS1|H|OSS4>
        StateCoup(indOSS4,indOSS1) = - 0.5_dp * dsqrt(mp_b) * ( dsqrt(np_a) * Wab(1,1) &
            & - dsqrt(np_a) * ERI(2,4) - dsqrt(np_d) * ERI(3,5) ) &
            & + 0.5_dp * dsqrt(mp_d) * ( dsqrt(np_d) * Wab(1,1) &
            & + dsqrt(np_d) * ERI(2,4) + dsqrt(np_a) * ERI(3,5) )
      end if

      if (indDOSS > 0) then
        ! <OSS1|H|DOSS>
        StateCoup(indDOSS,indOSS1) = (dsqrt(np_a) - dsqrt(np_d)) * Wab(3,1)
      end if

      if (indDSPS > 0) then
        ! <OSS1|H|DSPS>
        StateCoup(indDSPS,indOSS1) = 1.0_dp / dsqrt(3.0_dp) * (dsqrt(np_a) + dsqrt(np_d)) * (ERI(2,6) - ERI(1,5))
      end if

      if (indDES1 > 0) then
        ! <OSS1|H|DES1>
        StateCoup(indDES1,indOSS1) = 0.5_dp * (dsqrt(n_a*np_a) + dsqrt(n_d*np_d)) &
            & * (dsqrt(n_b) + dsqrt(n_c)) * Wab(4,1)
      end if

      if (indDES2 > 0) then
        ! <OSS1|H|DES2>
        StateCoup(indDES2,indOSS1) = 0.5_dp * (dsqrt(n_d*np_a) - dsqrt(n_a*np_d)) &
            & * (dsqrt(n_b) - dsqrt(n_c)) * Wab(4,1)
      end if

    end if

    if (indOSS2 > 0) then

      if (indOSS3 > 0) then
        ! <OSS2|H|OSS3>
        StateCoup(indOSS3,indOSS2) = - 0.5_dp * dsqrt(mp_a) * ( dsqrt(np_b) * Wab(1,1) &
            & - dsqrt(np_b) * ERI(3,5) - dsqrt(np_c) * ERI(2,4) ) &
            & + 0.5_dp * dsqrt(mp_c) * ( dsqrt(np_c) * Wab(1,1) &
            & + dsqrt(np_c) * ERI(3,5) + dsqrt(np_b) * ERI(2,4) )
      end if

      if (indOSS4 > 0) then
        ! <OSS2|H|OSS4>
        StateCoup(indOSS4,indOSS2) = 0.5_dp * dsqrt(mp_b) * ( dsqrt(np_b) * Wab(6,1) &
            & + dsqrt(np_b) * ERI(2,3) + dsqrt(np_c) * ERI(4,5) ) &
            & - 0.5_dp * dsqrt(mp_d) * ( dsqrt(np_c) * Wab(6,1) &
            & - dsqrt(np_c) * ERI(2,3) - dsqrt(np_b) * ERI(4,5) )
      end if

      if (indDOSS > 0) then
        ! <OSS2|H|DOSS>
        StateCoup(indDOSS,indOSS2) = (dsqrt(np_b) - dsqrt(np_c)) * Wab(4,1)
      end if

      if (indDSPS > 0) then
        ! <OSS2|H|DSPS>
        StateCoup(indDSPS,indOSS2) = 1.0_dp / dsqrt(3.0_dp) * (dsqrt(np_b) + dsqrt(np_c)) & 
            & * (ERI(5,6) - ERI(1,2))
      end if

      if (indDES1 > 0) then
        ! <OSS2|H|DES1>
        StateCoup(indDES1,indOSS2) = 0.5_dp * (dsqrt(n_c*np_b) - dsqrt(n_b*np_c)) &
            & * (dsqrt(n_a) - dsqrt(n_d)) * Wab(3,1)
      end if

      if (indDES2 > 0) then
        ! <OSS2|H|DES2>
        StateCoup(indDES2,indOSS2) = 0.5_dp * (dsqrt(n_b*np_b) + dsqrt(n_c*np_c)) &
            & * (dsqrt(n_a) + dsqrt(n_d)) * Wab(3,1)
      end if

    end if

    if (indOSS3 > 0) then

      if (indOSS4 > 0) then
        ! <OSS3|H|OSS4>
        StateCoup(indOSS4,indOSS3) = 0.5_dp * (dsqrt(mp_a*mp_b) + dsqrt(mp_c*mp_d)) &
            & * (2.0_dp * ERI(2,5) - ERI(1,6)) &
            & - 0.5_dp * (dsqrt(mp_a*mp_d) + dsqrt(mp_b*mp_c)) &
            & * (2.0_dp * ERI(2,5) - ERI(3,4))
      end if

      if (indDOSS > 0) then
        ! <OSS3|H|DOSS>
        StateCoup(indDOSS,indOSS3) = - (dsqrt(mp_a) - dsqrt(mp_c)) * Wab(2,1) &
            & + (dsqrt(mp_a) + dsqrt(mp_c)) * (ERI(3,6) - ERI(1,4))
      end if

      if (indDSPS > 0) then
        ! <OSS3|H|DSPS>
        StateCoup(indDSPS,indOSS3) = 0.5_dp * dsqrt(3.0_dp) * (dsqrt(mp_a) - dsqrt(mp_c)) * Wab(2,1) &
            & + 1.0_dp / dsqrt(3.0_dp) * (dsqrt(mp_a) + dsqrt(mp_c)) * (ERI(3,6) - ERI(1,4))
      end if

      if (indDES1 > 0) then
        ! <OSS3|H|DES1>
        StateCoup(indDES1,indOSS3) = 0.5_dp * dsqrt(mp_a) * ( dsqrt(n_a*n_c) * Wab(5,1) &
            & + dsqrt(n_a*n_b) * ERI(4,6) + dsqrt(n_c*n_d) * ERI(1,3) ) &
            & + 0.5_dp * dsqrt(mp_c) * ( dsqrt(n_b*n_d) * Wab(5,1) &
            & - dsqrt(n_a*n_b) * ERI(1,3) - dsqrt(n_c*n_d) * ERI(4,6) )
      end if

      if (indDES2 > 0) then
        ! <OSS3|H|DES2>
        StateCoup(indDES2,indOSS3) = 0.5_dp * dsqrt(mp_a) * ( dsqrt(n_b*n_d) * Wab(5,1) &
            & - dsqrt(n_a*n_b) * ERI(1,3) - dsqrt(n_c*n_d) * ERI(4,6) ) &
            & + 0.5_dp * dsqrt(mp_c) * ( dsqrt(n_a*n_c) * Wab(5,1) &
            & + dsqrt(n_a*n_b) * ERI(4,6) + dsqrt(n_c*n_d) * ERI(1,3) )
      end if

    end if

    if (indOSS4 > 0) then

      if (indDOSS > 0) then
        ! <OSS4|H|DOSS>
        StateCoup(indDOSS,indOSS4) = - (dsqrt(mp_b) - dsqrt(mp_d)) * Wab(5,1) &
            & + (dsqrt(mp_b) + dsqrt(mp_d)) * (ERI(4,6) - ERI(1,3))
      end if

      if (indDSPS > 0) then
        ! <OSS4|H|DSPS>
        StateCoup(indDSPS,indOSS4) = 0.5_dp * dsqrt(3.0_dp) * (dsqrt(mp_b) - dsqrt(mp_d)) * Wab(5,1) &
            & + 1.0_dp / dsqrt(3.0_dp) * (dsqrt(mp_b) + dsqrt(mp_d)) * (ERI(4,6) - ERI(1,3))
      end if

      if (indDES1 > 0) then
        ! <OSS4|H|DES1>
        StateCoup(indDES1,indOSS4) = 0.5_dp * dsqrt(mp_b) * ( dsqrt(n_a*n_c) * Wab(2,1) &
            & - dsqrt(n_a*n_b) * ERI(1,4) - dsqrt(n_c*n_d) * ERI(3,6) ) &
            & + 0.5_dp * dsqrt(mp_d) * ( dsqrt(n_b*n_d) * Wab(2,1) &
            & + dsqrt(n_a*n_b) * ERI(3,6) + dsqrt(n_c*n_d) * ERI(1,4) )
      end if

      if (indDES2 > 0) then
        ! <OSS4|H|DES2>
        StateCoup(indDES2,indOSS4) = 0.5_dp * dsqrt(mp_b) * ( dsqrt(n_b*n_d) * Wab(2,1) &
            & + dsqrt(n_a*n_b) * ERI(3,6) + dsqrt(n_c*n_d) * ERI(1,4) ) &
            & + 0.5_dp * dsqrt(mp_d) * ( dsqrt(n_a*n_c) * Wab(2,1) &
            & - dsqrt(n_a*n_b) * ERI(1,4) - dsqrt(n_c*n_d) * ERI(3,6) )
      end if

    end if

    if (indDOSS > 0) then

      if (indDSPS > 0) then
        ! <DOSS|H|DSPS>
        StateCoup(indDSPS,indDOSS) = 0.25_dp * dsqrt(3.0_dp) * (enLtot(31) + enLtot(32) - enLtot(29) - enLtot(30))
      end if

      if (indDES1 > 0) then
        ! <DOSS|H|DES1>
        StateCoup(indDES1,indDOSS) = 0.5_dp * (dsqrt(n_a*n_c) - dsqrt(n_b*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(2,5)) &
            & + 0.5_dp * (dsqrt(n_a*n_b) - dsqrt(n_c*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(1,6))
      end if

      if (indDES2 > 0) then
        ! <DOSS|H|DES2>
        StateCoup(indDES2,indDOSS) = - 0.5_dp * (dsqrt(n_a*n_c) - dsqrt(n_b*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(2,5)) &
            & + 0.5_dp * (dsqrt(n_a*n_b) - dsqrt(n_c*n_d)) &
            & * (2.0_dp * ERI(3,4) - ERI(1,6))
      end if

    end if

    if (indDSPS > 0) then

      if (indDES1 > 0) then
        ! <DSPS|H|DES1>
        StateCoup(indDES1,indDSPS) = 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_c) - dsqrt(n_b*n_d)) * ERI(2,5) &
            & - 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_b) - dsqrt(n_c*n_d)) * ERI(1,6)
      end if

      if (indDES2 > 0) then
        ! <DSPS|H|DES2>
        StateCoup(indDES2,indDSPS) = - 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_c) - dsqrt(n_b*n_d)) * ERI(2,5) &
            & - 0.5_dp * dsqrt(3.0_dp) * (dsqrt(n_a*n_b) - dsqrt(n_c*n_d)) * ERI(1,6)
      end if

    end if

    call adjointLowerTriangle(StateCoup)

  end subroutine getStateCoup44_


end module dftbp_reks_reksen
