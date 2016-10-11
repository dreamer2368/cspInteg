MODULE eigen

  IMPLICIT NONE

CONTAINS

  SUBROUTINE eigenDec(nvars, dfx, V, U, eigv)

    !% input/output variables
    integer, intent(in)  :: nvars
    real*8,  intent(in)  :: dfx(nvars,nvars)
    real*8,  intent(out) :: V(nvars,nvars)
    real*8,  intent(out) :: U(nvars,nvars)
    real*8,  intent(out) :: eigv(nvars)
    !% local variables
    !  lapack
    integer              :: astat
    integer              :: info
    integer              :: lwork
    character            :: jobu = 'N'
    character            :: jobv = 'V'
    real*8               :: evr(nvars)
    real*8               :: evi(nvars)
    real*8               :: evv(nvars,2)
    real*8               :: rlam(nvars)
    real*8,  allocatable :: work(:)
    !  ordering
    integer              :: i,j,k
    integer              :: imv
    real*8               :: x,y,z
    real*8               :: temp
    real*8               :: norm
    !  inverse
    integer              :: liwork
    integer              :: ipiv(nvars)
    real*8,  allocatable :: iwork(:)
    
    
    !=======================================================!
    ! allocation
    if(allocated(work)) deallocate(work)

    lwork  = 4*nvars
    liwork = nvars
    
    allocate(&
         work(lwork),   &
         iwork(liwork), &
         stat=astat)
    if (astat .ne. 0) stop

    !=======================================================!
    ! eigendecomposition of dfx
    call dgeev(jobu, jobv, nvars, dfx, nvars, evr, evi, &
         U, nvars, V, nvars, work, lwork, info)

    evv(:,1) = evr
    evv(:,2) = evi
    rlam     = abs(evr)

    !=======================================================!
    ! order eigenvalues and eigenvectors    
    do i = 1,nvars-1
       x   = rlam(i)
       imv = i
       do j = i,nvars
          y = rlam(j)
          if(y .lt. x) then
             x   = y
             imv = j
          end if
       end do
       if(i .ne. imv) then
          do k = 1,2
             temp       = evv(i,k)
             evv(i,k)   = evv(imv,k)
             evv(imv,k) = temp
          end do
          temp      = rlam(i)
          rlam(i)   = rlam(imv)
          rlam(imv) = temp
          do k = 1,nvars
             temp     = V(k,i)
             V(k,i)   = V(k,imv)
             V(k,imv) = temp
          end do
       end if  
    end do

    !=======================================================!
    ! real canonical form
    do i = 1,nvars-1
       x = evv(i,2)
       y = evv(i+1,2)
       if((x .eq. 0.0d0) .or. (y .eq. 0.0d0)) then
          cycle
       end if
       z = abs(abs(x)-abs(y)/abs(x))
       if(z .gt. 1.0d-07) then
          cycle
       end if
       do j = 1,nvars
          x = V(j,i)
          y = V(j,i+1)
          V(j,i)   = x-y
          V(j,i+1) = x+y
       end do
    end do

    !=======================================================!
    ! normalize and get left eigenvectors
    do i = 1,nvars
       norm   = sqrt( sum( V(:,i) * V(:,i) ) )
       norm   = 1.0d0 / norm
       V(:,i) = norm * V(:,i)
       U(:,i) = V(:,i)
    end do

    call dgetrf(nvars, nvars, U, nvars, ipiv, info)
    call dgetri(nvars, U, nvars, ipiv, iwork, liwork, info)

    eigv = evv(:,1)
 
  END SUBROUTINE eigenDec
  
  SUBROUTINE schur(mode, nslow, nfast, nvars, dfx, qr, ql, qrs, qrf, qls, qlf, eig, eigv, eigs, eigf)

    ! input/output variables
    integer,   intent(in)            :: mode
    integer,   intent(in)            :: nslow
    integer,   intent(in)            :: nfast
    integer,   intent(in)            :: nvars
    real*8,    intent(inout)         :: dfx(nvars,nvars)
    real*8,    intent(out), optional :: qr(nvars,nvars)
    real*8,    intent(out), optional :: ql(nvars,nvars)
    real*8,    intent(out), optional :: qrs(nvars,nslow)
    real*8,    intent(out), optional :: qrf(nvars,nfast)
    real*8,    intent(out), optional :: qls(nslow,nvars)
    real*8,    intent(out), optional :: qlf(nfast,nvars)
    real*8,    intent(out), optional :: eig(nvars,nvars)
    real*8,    intent(out), optional :: eigv(nvars)
    real*8,    intent(out), optional :: eigs(nslow,nslow)
    real*8,    intent(out), optional :: eigf(nfast,nfast)
    ! local variables
    character                        :: jobvs = 'v'
    character                        :: sort  = 'n'
    integer                          :: i, astat
    integer                          :: info
    integer                          :: ifst
    integer                          :: ilst
    integer                          :: sdim = 0
    integer                          :: lwork1
    real*8                           :: jac(nvars,nvars)
    real*8                           :: qz(nvars,nvars)
    real*8                           :: qzt(nvars,nvars)
    real*8                           :: zrs(nvars,nslow)
    real*8                           :: zrf(nvars,nfast)
    real*8                           :: zls(nslow,nvars)
    real*8                           :: zlf(nfast,nvars)
    real*8                           :: wr(nvars)
    real*8                           :: wi(nvars)
    real*8                           :: work2(nvars)
    real*8,    allocatable           :: work1(:)
    logical                          :: sel, bwork(nvars)

    real*8 :: qdq(nvars,nvars)

800 format(4e16.8)

    !=======================================================!
    ! allocation
    if(allocated(work1)) deallocate(work1)

    lwork1 = 4*nvars
    
    allocate(work1(lwork1), stat=astat)
    if (astat .ne. 0) stop

    !=======================================================!
    ! schur decomposition of dfx
    jac = dfx
    call dgees(jobvs, sort, sel, nvars, jac, nvars, &
         sdim, wr, wi, qz, nvars, work1, lwork1, bwork, info)
    
    !=======================================================!
    ! order schur decomposition
    ilst = 1
100 ifst = ilst

    do i = ilst+1,nvars
       if (abs(jac(i,i)) .lt. abs(jac(ifst,ifst))) then
          ifst = i
       end if
    end do

    if(ifst .ne. ilst) then
       call dtrexc(jobvs, nvars, jac, nvars, qz, nvars, &
            ifst, ilst, work2, info)
    end if

    ilst = ilst + 1

    if (ilst < nvars) then
       goto 100
    end if

    !=======================================================!
    ! output
    select case (mode)
    case (0)

       if(present(qr)) then
          qr = qz
       end if

       if(present(ql)) then
          ql = qzt
       end if

       if(present(eig)) then
          eig = jac
       end if

       if(present(eigv)) then
          eigv = 0.0d0
          do i = 1,nvars
             eigv(i) = jac(i,i)
          end do
       end if

    case (1)

       ! invariant decomposition
       call split_schur(nslow=nslow, nfast=nfast, nvars=nvars,  &
            jac=jac, qz=qz, zrs=zrs, zrf=zrf, zls=zls, zlf=zlf, &
            info=info)
    
       if(present(qrs)) then
          qrs = zrs
       end if
       
       if(present(qrf)) then
          qrf = zrf
       end if

       if(present(qls)) then
          qls = zls
       end if

       if(present(qlf)) then
          qlf = zlf
       end if

       if(present(eigs)) then
          eigs = jac(1:nslow,1:nslow)
       end if

       if(present(eigf)) then
          eigf = jac(nslow+1:nvars,nslow+1:nvars)
       end if

       if(present(eigv)) then
          eigv = 0.0d0
          do i = 1,nvars
             eigv(i) = jac(i,i)
          end do
       end if

    end select

700 continue

    !=======================================================!
    ! deallocation
    if(allocated(work1)) deallocate(work1)
    
  END SUBROUTINE schur

  SUBROUTINE split_schur(nslow, nfast, nvars, jac, qz, &
                         myzr, myzl, zrs, zrf, zls, zlf, info)

  ! input/output variables
    integer, intent(in)              :: nslow
    integer, intent(in)              :: nfast
    integer, intent(in)              :: nvars
    real*8,  intent(inout)           :: jac(nvars,nvars)
    real*8,  intent(in)              :: qz(nvars,nvars)
    real*8,  intent(out), optional   :: myzr(nvars,nvars)
    real*8,  intent(out), optional   :: myzl(nvars,nvars)
    real*8,  intent(out), optional   :: zrs(nvars,nslow)
    real*8,  intent(out), optional   :: zrf(nvars,nfast)
    real*8,  intent(out), optional   :: zls(nslow,nvars)
    real*8,  intent(out), optional   :: zlf(nfast,nvars)
    integer, intent(out)             :: info
    ! local variables
    character                        :: trans = 'n'
    integer                          :: i, j, astat
    integer                          :: flag
    integer                          :: isgn = -1     
    real*8                           :: scale = 1.0d0
    real*8                           :: mat(nvars,nvars)
    real*8                           :: jss(nslow,nslow)
    real*8                           :: jff(nfast,nfast)
    real*8                           :: jsf(nslow,nfast)
    real*8                           :: zr(nvars,nvars)
    real*8                           :: zl(nvars,nvars)
    real*8                           :: szr(nvars,nslow)
    real*8                           :: fzr(nvars,nfast)
    real*8                           :: szl(nslow,nvars)
    real*8                           :: fzl(nfast,nvars)
    
600 format(2f10.4)

    !=======================================================!
    ! solve sylvester equation
    mat = 0.0d0
    forall (i = 1:nvars) mat(i,i) = 1.0d0
    
    jss =  jac(1:nslow,1:nslow)
    jff =  jac(nslow+1:nvars,nslow+1:nvars)
    jsf =  jac(1:nslow,nslow+1:nvars)
    jsf = -jsf
    
   call dtrsyl(trans, trans, isgn, nslow, nfast, jss, &
        nslow, jff, nfast, jsf, nslow, scale, info)

   !=======================================================!
   ! find invariant decomposition
   ! info = 0: dtrsyl successful
   ! info = 1: dtrsyl failed, no zr/zl decomposition
   select case (info)
   case (0)
      flag = 0
100   select case (flag)
         case (0)
            mat(1:nslow,1+nslow:nvars) =  jsf
            zr = matmul(qz, mat)
            flag = 1
            goto 100
         case (1)
            mat(1:nslow,1+nslow:nvars) = -jsf
            zl = matmul(mat, transpose(qz))
         end select
         szr = zr(1:nvars,1:nslow)
         fzr = zr(1:nvars,nslow+1:nvars)
         szl = zl(1:nslow,1:nvars)
         fzl = zl(nslow+1:nvars,1:nvars)        
    case (1)
       zr = qz
       zl = transpose(qz)
       szr = zr(1:nvars,1:nslow)
       fzr = zr(1:nvars,nslow+1:nvars)
       szl = zl(1:nslow,1:nvars)
       fzl = zl(nslow+1:nvars,1:nvars)
    end select

    continue

    jac(1:nslow,nslow+1:nvars) = 0.0d0

    if(present(myzr)) then
       myzr = zr
    end if

    if(present(myzl)) then
       myzl = zl
    end if

    if(present(zrs)) then
       zrs = szr
    end if

    if(present(zrf)) then
       zrf = fzr
    end if

    if(present(zls)) then
       zls = szl
    end if

    if(present(zlf)) then
       zlf = fzl
    end if
    
  END SUBROUTINE split_schur

  
!!!
!!!
  SUBROUTINE fronorm_schur(nslow, nfast, nvars, dfx, nsf, zr, zl, fronorm)

    ! input/output variables
    integer, intent(in)  :: nslow
    integer, intent(in)  :: nfast
    integer, intent(in)  :: nvars
    real*8,  intent(in)  :: dfx(nvars,nvars)
    real*8,  intent(in)  :: nsf(nvars,nvars)
    real*8,  intent(in)  :: zr(nvars,nvars)
    real*8,  intent(in)  :: zl(nvars,nvars)
    real*8,  intent(out) :: fronorm
    ! local variables
    integer              :: astat
    real*8,  allocatable :: &
         zreig(:,:),        &
         zreigzl(:,:),      &
         work(:)
    real*8, external     :: dlange

    !=======================================================!
    ! allocation
    if(allocated(work))    deallocate(work)
    if(allocated(zreig))   deallocate(zreig)
    if(allocated(zreigzl)) deallocate(zreigzl)

    allocate(work(nvars),      &
         zreig(nvars,nvars),   &
         zreigzl(nvars,nvars), &
         stat = astat)
    if (astat .ne. 0) stop

    !=======================================================!
    ! jacobian reconstruction
    zreig   = matmul(zr, nsf)
    zreigzl = matmul(zreig, zl) 
    zreigzl = zreigzl - dfx

    !=======================================================!
    ! frobenius norm
    fronorm = dlange('F', nvars, nvars, zreigzl, nvars, work)

    !=======================================================!
    ! deallocation
    deallocate(zreig, zreigzl, work)

  END SUBROUTINE fronorm_schur
  

!!!
!!!
  SUBROUTINE colspace(m, n, r, mat, col)

    IMPLICIT NONE
    ! input/output variables
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    integer, intent(in)  :: r
    integer, intent(in)  :: mat(m,n)
    real*8,  intent(out) :: col(m,r)
    ! local variables
    integer              :: astat
    integer              :: pivot
    integer              :: i, j, t
    integer              :: npiv
    integer, allocatable :: ipiv(:)
    real*8,  allocatable :: tmat(:,:)
    real*8,  allocatable :: trow(:)
    real*8,  allocatable :: temp(:,:)
   
    npiv = min(m,n)

    !=======================================================!
    ! allocation   
    if(allocated(ipiv)) deallocate(ipiv)
    if(allocated(tmat)) deallocate(tmat)
    if(allocated(trow)) deallocate(trow)
    if(allocated(temp)) deallocate(temp)
    
    allocate(&
         ipiv(npiv), &
         trow(m),    &
         tmat(n,m),  &
         temp(m,n),  &
         stat=astat)
    if(astat .ne. 0) stop

    !=======================================================!
    ! lu decomposition of mat-transpose
    tmat = real(transpose(mat),8)
    
    call dgetrf(n, m, tmat, n, ipiv, astat)

    forall(i = 1:n, j = 1:m, i .gt. j) tmat(i,j) = 0

    !=======================================================!
    ! row reduction of u
    call rref(n, m, tmat)

    !=======================================================!
    ! transpose and remove zero-cols   
    temp = transpose(tmat)

    col  = temp(1:m,1:r)

    !=======================================================!
    ! deallocation 
    if(allocated(ipiv)) deallocate(ipiv)
    if(allocated(tmat)) deallocate(tmat)
    if(allocated(trow)) deallocate(trow)
    if(allocated(temp)) deallocate(temp)
    
  END SUBROUTINE colspace


!!!
!!!
  SUBROUTINE pinv(m, n, r, mat, inv)

    IMPLICIT NONE
    ! input/ output variables
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    integer, intent(in)  :: r
    real*8,  intent(in)  :: mat(m,n)
    real*8,  intent(out) :: inv(n,m)
    ! local variables
    character            :: jobu = 'a', jobvt = 'a'
    integer              :: astat, i
    integer              :: nrhs
    integer              :: lsval
    integer              :: lwork
    real*8,  allocatable :: a(:,:)
    real*8,  allocatable :: u(:,:) 
    real*8,  allocatable :: sval(:)
    real*8,  allocatable :: vt(:,:)
    real*8,  allocatable :: s(:,:)
    real*8,  allocatable :: work(:)
    
    !=======================================================!
    ! allocation
    lsval = min(m,n)
    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    if(allocated(a))    deallocate(a)
    if(allocated(u))    deallocate(u)
    if(allocated(sval)) deallocate(sval)
    if(allocated(vt))   deallocate(vt)
    if(allocated(s))    deallocate(s)
    if(allocated(work)) deallocate(work)

    allocate(&
         a(m,n),         &
         u(m,m),         &
         sval(lsval),    &
         vt(n,n),        &
         s(n,m),         &
         work(lwork),    &
         stat=astat)
    if(astat .ne. 0) stop

    !=======================================================!
    ! solve the linear least squares problem
    a = mat
    
    call dgesvd(jobu, jobvt, m, n, a, m, sval, &
         u, m, vt, n, work, lwork, astat)

    sval = 1/sval

    s = 0
    forall(i = 1:r) s(i,i) = sval(i)

    s   = matmul(s,transpose(u))
    
    inv = matmul(transpose(vt),s)

    !=======================================================!
    ! deallocation
    if(allocated(u))    deallocate(u)
    if(allocated(sval)) deallocate(sval)
    if(allocated(vt))   deallocate(vt)
    if(allocated(work)) deallocate(work)
    
  END SUBROUTINE pinv


!!!
!!!
  SUBROUTINE rref(m, n, mat)

    IMPLICIT NONE
    ! input/output variables
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    real*8,  intent(inout) :: mat(m,n)
    ! local variables
    integer :: pivot
    integer :: i, j
    real*8, allocatable :: trow(:)

    !=======================================================!
    ! allocation
    allocate(trow(n))

    !=======================================================!
    ! elimination
    pivot = 1
    do j = 1, m
       if(n .le. pivot) exit
       i = j
       do while(mat(i,pivot) .eq. 0)
          i = i + 1
          if(m .eq. i) then
             i = j
             pivot = pivot + 1
             if(n .eq. pivot) return
          end if
       end do
       trow = mat(i,:)
       mat(i,:) = mat(j,:)
       mat(j,:) = trow
       mat(j,:) = mat(j,:)/mat(j,pivot)
       do i = 1, m
          if (i .ne. j) then
             mat(i,:) = mat(i,:) - mat(j,:) * mat(i,pivot)
          end if
       end do
       pivot = pivot + 1
    end do

    !=======================================================!
    ! deallocation
    deallocate(trow)
    
  END SUBROUTINE rref

  SUBROUTINE nullspace(m, n, r, k, mat, ker)
      
    ! input/output variables
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    integer, intent(in)  :: r
    integer, intent(in)  :: k
    real*8,  intent(in)  :: mat(m,n)
    real*8,  intent(out) :: ker(n,k)
    ! local variables
    character            :: jobu = 'a', jobvt = 'a'
    integer              :: astat, i, info
    integer              :: lsval
    integer              :: lwork
    real*8               :: a(m,n)
    real*8               :: u(m,m)
    real*8               :: vt(n,n)
    real*8,  allocatable :: sval(:)
    real*8,  allocatable :: work(:)
    
    !=======================================================!
    ! allocation
    lsval = min(m,n)
    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    if(allocated(sval)) deallocate(sval)
    if(allocated(work)) deallocate(work)

    allocate(&
         sval(lsval),    &
         work(lwork),    &
         stat=astat)
    if(astat .ne. 0) stop

    !=======================================================!
    ! solve the linear least squares problem
    a = mat

    info = 0
    call dgesvd(jobu, jobvt, m, n, a, m, sval, &
         u, m, vt, n, work, lwork, info)
    if(info .ne. 0) stop

    vt  = transpose(vt)
    ker = vt(:,r+1:n)

    !=======================================================!
    ! deallocation
    if(allocated(sval)) deallocate(sval)
    if(allocated(work)) deallocate(work)

  END SUBROUTINE nullspace

  SUBROUTINE conjugategradient(n, mat, b, x, ifail)

    ! input/output variables
    integer, intent(in)  :: n
    real*8,  intent(in)  :: mat(n,n)
    real*8,  intent(in)  :: b(n)
    real*8,  intent(out) :: x(n)
    integer, intent(out) :: ifail
    ! local variables
    integer              :: i, info
    integer              :: k, kmax
    integer              :: ipiv(n)
    real*8               :: alpha
    real*8               :: beta
    real*8               :: r(n)
    real*8               :: rp(n)
    real*8               :: y(n)
    real*8               :: yp(n)
    real*8               :: p(n)
    real*8               :: pcn(n,n)
    real*8               :: err
    real*8               :: tol

    !=======================================================!
    ! initialize
    kmax = 50
    tol  = 1
    tol  = epsilon(tol)

    pcn  =  0.0
    forall(i = 1:n) pcn(i,i) = 1/mat(i,i)

    x    =  0.0
    r    =  matmul(mat,x) - b
    y    =  matmul(pcn,r)
    rp   =  r
    yp   =  y
    p    = -y
    err  =  sqrt(sum(r*r))
    k    =  0
        
    !=======================================================!
    ! iterate
    do while((err > tol) .and. (k .le. kmax)) 

       alpha = sum(r*y)/sum(p*matmul(mat,p))
       x     = x + alpha*p
       r     = r + alpha*matmul(mat,p)
       y     = matmul(pcn,r)
       beta  = sum(r*y)/sum(rp*yp)
       p     = beta*p - y
       rp    = r
       yp    = y
       err   = sqrt(sum(r*r))
       k     = k + 1

    end do

    
    if(err > tol) then
       ifail = 1
       write(*,'(A23,i3,A11)') "PCG failed after ", k," iterations"
       write(*,'(A6,e14.4)') "err = ", err
    else
       ifail = 0
       write(*,'(A23,i3,A11)') "PCG converged within ", k," iterations"
       write(*,'(A6,e14.4)') "err = ", err
    end if
    
  END SUBROUTINE conjugategradient

END MODULE eigen
