!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the process of converting from a Hamiltonian and overlap in sparse form
!> into eigenvectors
module dftbp_dftbplus_eigenvects
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolver, electronicSolverTypes
  use dftbp_extlibs_elsiiface, only : elsi_write_mat_complex, elsi_finalize_rw, elsi_ev_complex,&
      & elsi_ev_real, elsi_write_mat_real
  use dftbp_io_message, only : error, cleanShutdown
  use dftbp_math_eigensolver, only : hegv, hegvd, gvr
#:if WITH_GPU
  use dftbp_math_eigensolver, only : gpu_gvd
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_, scalafx_phegv, scalafx_phegvd, scalafx_phegvr,&
      & scalafx_psygv, scalafx_psygvd, scalafx_psygvr
#:endif
  implicit none
  
  private
  public :: diagDenseMtx
#:if WITH_SCALAPACK
  public :: diagDenseMtxBlacs
#:endif


  interface diagDenseMtx
    module procedure diagDenseRealMtx
    module procedure diagDenseComplexMtx
  end interface diagDenseMtx

#:if WITH_SCALAPACK

  interface diagDenseMtxBlacs
    module procedure diagDenseRealMtxBlacs
    module procedure diagDenseCplxMtxBlacs
  end interface diagDenseMtxBlacs

#:endif

contains

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenvectors and values,
  !> as well as often the Cholesky factorized overlap matrix (due to a side effect of lapack)
  subroutine diagDenseRealMtx(env, electronicSolver, jobz, HSqrReal, SSqrReal, eigen)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' with vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> Eigen values.
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))
    @:ASSERT(all(shape(HSqrReal) == shape(SSqrReal)))
    @:ASSERT(size(HSqrReal, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call gvr(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%magma_gvd)
  #:if WITH_GPU
      call gpu_gvd(env%gpu%nGpu, HSqrReal, SSqrReal, eigen, 'L', jobz)
  #:else
      call error("This binary is compiled without GPU support")
  #:endif
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseRealMtx


  !> Diagonalizes a sparse represented Hamiltonian and overlap with k-points to give the
  !> eigenvectors and values, as well as often the Cholesky factorized overlap matrix (due to a side
  !> effect of lapack)
  subroutine diagDenseComplexMtx(env, electronicSolver, jobz, HSqrCplx, SSqrCplx, eigen)

    !> Environment
    type(TEnvironment), intent(in) :: env
    
    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    !> Large square matrix for the overlap workspace, overwritten with the Cholesky factorized form.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> The eigenvalues of the matrices
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%magma_gvd)
  #:if WITH_GPU
      call gpu_gvd(env%gpu%nGpu, HSqrCplx, SSqrCplx, eigen, 'L', jobz)
  #:else
      call error("This binary is compiled without GPU support")
  #:endif
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseComplexMtx


#:if WITH_SCALAPACK

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseRealMtxBlacs(electronicSolver, iCholesky, jobz, desc, HSqr, SSqr, eigenVals,&
      & eigenVecs)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Index of the overlap matrix where Cholesky should be stored.
    integer, intent(in) :: iCholesky

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix with hamiltonian
    real(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    real(dp), intent(out) :: eigenVecs(:,:)

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    if (electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%getCholesky(iCholesky, SSqr)
    end if

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call scalafx_psygv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%divideandconquer)
      call scalafx_psygvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & allocfix=.true., skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%relativelyrobust)
      call scalafx_psygvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%elpa)
      if (electronicSolver%elsi%tWriteHS) then
        call elsi_write_mat_real(electronicSolver%elsi%rwHandle, "ELSI_Hreal.bin", HSqr)
        call elsi_write_mat_real(electronicSolver%elsi%rwHandle, "ELSI_Sreal.bin", SSqr)
        call elsi_finalize_rw(electronicSolver%elsi%rwHandle)
        call cleanShutdown("Finished dense matrix write")
      end if
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_real(electronicSolver%elsi%handle, HSqr, SSqr, eigenVals, eigenVecs)

    case default
      call error('Unknown eigensolver')
    end select

    if (.not. electronicSolver%hasCholesky(1)) then
      call electronicSolver%storeCholesky(1, SSqr)
    end if

  end subroutine diagDenseRealMtxBlacs



  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseCplxMtxBlacs(electronicSolver, iCholesky, jobz, desc, HSqr, SSqr, eigenVals,&
      & eigenVecs)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Index of the overlap matrix where Cholesky should be stored.
    integer, intent(in) :: iCholesky

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix for the resulting eigenValsvectors
    complex(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    complex(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    complex(dp), intent(out) :: eigenVecs(:,:)

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    if (electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%getCholesky(iCholesky, SSqr)
    end if

    select case(electronicSolver%iSolver)

    case(electronicSolverTypes%QR)
      call scalafx_phegv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz, &
          & skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%divideandconquer)
      call scalafx_phegvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & allocfix=.true., skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%relativelyrobust)
      call scalafx_phegvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%elpa)
      if (electronicSolver%elsi%tWriteHS) then
        call elsi_write_mat_complex(electronicSolver%elsi%rwHandle, "ELSI_Hcmplx.bin", HSqr)
        call elsi_write_mat_complex(electronicSolver%elsi%rwHandle, "ELSI_Scmplx.bin", SSqr)
        call elsi_finalize_rw(electronicSolver%elsi%rwHandle)
        call cleanShutdown("Finished dense matrix write")
      end if
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_complex(electronicSolver%elsi%handle, HSqr, SSqr, eigenVals, eigenVecs)

    case default
      call error('Unknown eigensolver')
    end select

    if (.not. electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%storeCholesky(iCholesky, SSqr)
    end if

  end subroutine diagDenseCplxMtxBlacs

#:endif

end module dftbp_dftbplus_eigenvects
