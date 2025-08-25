subroutine DiagH0TAPW(N, ns, is, ELoc, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species
    use tbpar, only : U
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    real(dp), intent(out) :: ELoc(N)
    real(dp), intent(in) :: KLoc(3), cell(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    complex(dp), allocatable :: resid(:), v(:,:), workd(:), d(:), workev(:)
    complex(dp), allocatable :: workl(:)
    complex(dp), allocatable :: z(:,:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iparam(:), ipntr(:)
    logical, allocatable :: select(:)
    integer :: max_iter, nn, iter
    logical :: useShift
    double precision :: shift

    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)
    complex(dp) :: sigma

    ! PARDISO variables
    !integer, parameter :: mtype = 13  ! Complex unsymmetric matrix
    !integer :: pt(64), iparm(64), maxfct, mnum, phase, error, msglvl
    !integer :: nnn, nrhs
    !complex(dp), allocatable :: a(:)
    !integer, allocatable :: ia(:), ja(:)
    !complex(dp) :: ddum
    !integer :: idum

    INTEGER :: nnn, nnz
    !complex(dp) :: x(N) ! Pardiso needs X(N) even if sol is returned in b

    !complex(dp) :: b(N)
!C.. Internal solver memory pointer for 64-bit architectures
    INTEGER*8 pt(64)
!C.. Internal solver memory pointer for 32-bit architectures
!C.. INTEGER*4 pt(64)

!C.. This is OK in both cases
!    TYPE(MKL_PARDISO_HANDLE) pt(64)
!C.. All other variables
    INTEGER maxfct, mnum, mtype, phase, nrhs, msglvl
    INTEGER iparm(64)
    !allocate(row_ptr(N+1), col_ind(nnz_temp), values(nnz_temp))
    INTEGER, allocatable :: ia(:) ! row_ptr
    INTEGER, allocatable :: ja(:) ! col_ind
    INTEGER idum(1)
    complex(dp) ::  ddum(1)
    integer error
    complex(dp), allocatable :: a(:) ! values
    
    DATA nrhs /1/, maxfct /1/, mnum /1/

    complex(dp), allocatable :: EVectors(:,:)
    logical :: saveRitz, symmetric
    integer :: nconv
    complex(dp), allocatable :: ax(:)
    double precision, allocatable :: rd(:,:)
    double precision :: dznrm2, dlapy2

    ! === TAPW variables ===
    integer :: NG, Nlabel, M
    real(dp), allocatable :: x(:), y(:)
    integer, allocatable :: label(:)
    real(dp), allocatable :: Gx(:), Gy(:)
    complex(dp), allocatable :: X(:,:), Hproj(:,:)
    real(dp), allocatable :: eigvals(:)
    complex(dp), allocatable :: ZWorkLoc(:)
    real(dp), allocatable :: DWorkLoc(:)
    integer :: lwork

    
    ! Parameters
    !integer, parameter :: dp = c_double
    !integer, parameter :: int_kind = c_int
    
    ! Variables
    !type(C_PTR) :: mkl_handle
    !!integer(int_kind) :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !!integer(int_kind), allocatable :: row_ptr(:), col_ind(:)
    !integer :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !integer, allocatable :: row_ptr(:), col_ind(:)
    !real(dp), allocatable :: values(:), b(:), x(:), temp_value
    !type(SPARSE_MATRIX_DESCR) :: descr


    !integer :: nev, ncv, lworkl 
 
    nev=neig 
    ncv=nev*10
    lworkl= 3*NCV**2 + 5*NCV

    max_iter = 10000
    allocate(resid(N), v(N, ncv), workd(3*N), workl(lworkl), rwork(ncv), d(nev+1), iparam(11), ipntr(14), select(ncv), rand_real(N), rand_imag(N))
    allocate(z(N,nev))
    allocate(workev(2*ncv))
    allocate(ax(N))
    allocate(rd(ncv,3))

    ! Ensure correct size of nev and ncv
    if ( (nev < 1) .or. (nev >= ncv) .or. (ncv > N) ) then
        print *, 'Error: invalid parameters - nev=', nev, 'ncv=', ncv, 'N=', N
        stop
    end if

    ! Initialize arrays
    !v = 0.0d0
    !workd = 0.0d0
    !workl = 0.0d0
    !rwork = 0.0d0
    !d = 0.0d0

    bmat = 'I'

    iparam(1) = 1
    iparam(3) = max_iter
    call MIO_InputParameter('Diag.SparseSetShift',shift,0.0_dp)
    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    call MIO_InputParameter('Diag.SparseSaveRitz',saveRitz,.false.)
    if (useShift) then
       which = 'LM'
       iparam(7) = 3
       sigma = cmplx(shift, 0.0_dp)
       print*, "We will find eigenvalues close to ", shift
    else
       which = 'SM'
       sigma = (0.0_dp, 0.0_dp)
       iparam(7) = 1
    end if
    ! Initialize ARPACK parameters

    call MIO_InputParameter('Diag.SparseSetTol',tol,0.1_dp)

    ido = 0
    info = 0

    !iparam(4) = 1
    !iparam(7) = 1

    nn = N

    ! Initialize the starting vector resid with random values
    !call random_number(rand_real)
    !call random_number(rand_imag)
    !resid = cmplx(rand_real, rand_imag)
    !resid = 0

    ! Check resid for initial state
    if (size(resid) /= N) then
        print *, 'Error: resid size mismatch: ', size(resid), ' expected: ', N
        stop
    end if

     ! Debug print for resid
    if (any(resid /= resid)) then
        print *, 'Error: resid contains NaN values initially.'
        stop
    end if

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints
    !print *, 'Initial resid norm: ', resid_norm
    !print *, 'Initial iparam: ', iparam
    !print *, 'nn, nev, ncv: ', nn, nev, ncv

    ! Create CSR sparse matrix storage
    print *, "initialize the sparse matrix and put it in csr format"
    call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values,sigma)
    print *, "done"
    !call test_sparse_matvec(nn, row_ptr, col_ind, values)
    symmetric = is_structurally_symmetric(values, row_ptr, col_ind, n)

    print*, "is it symmetric?", symmetric

    ! Debug prints for CSR matrix
    if (any(values /= values)) then
        print *, 'Error: values contains NaN values after initialization.'
        stop
    end if

    
    if (maxval(abs(values)) > 1e10) then
        print *, 'Warning: values contains extremely large values.'
    end if

    if (any(col_ind < 1 .or. col_ind > N)) then
        print *, 'Error: col_ind contains invalid indices after initialization.'
        stop
    end if

    ! === TAPW settings (set these appropriately) ===
    NG = 20         ! for example
    Nlabel = 2      ! typically 2 (e.g., sublattice A/B)
    
    M = NG * Nlabel

    real(dp) :: k_ref(2)
    k_ref = (/ 4.0_dp*pi/(3.0_dp*a), 0.0_dp /)   ! set appropriately

    integer :: NGrange
    NGrange = 4   ! e.g., use all G such that -4 ≤ n1,n2 ≤ +4
    
    !call generate_G_list_from_rcell(rcell, NGrange, Gx, Gy, NG)
    call generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG)

    Nlabel = 2
    M = NG * Nlabel
    
    allocate(X(N, M))
    call build_X(X, x, y, label, Gx, Gy, N, NG, Nlabel)


    allocate(Hproj(M, M))
    call transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, X, Hproj)

    !call zheev('V', 'U', M, Hproj, M, eigvals, work, lwork, rwork, info)
    !call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
    allocate(eigvals(M))
    lwork = 2*M
    allocate(ZWorkLoc(2*M), DWorkLoc(3*M))
    
    call ZHEEV('V','U', M, Hproj, M, eigvals, ZWorkLoc, lwork, DWorkLoc, info)
    
    if (info /= 0) then
       print *, "Diagonalization failed: ZHEEV info =", info
       stop
    end if

! eigvals now contains eigenvalues of projected H




end subroutine DiagH0TAPW
