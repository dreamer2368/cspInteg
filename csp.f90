MODULE csp
  
  USE chem
  USE eigen

  IMPLICIT NONE

  ! global variables
  integer                   :: ncons         ! number of cons variables
  integer                   :: nslow         ! number of slow variables
  integer                   :: nfast         ! number of fast variables
  integer                   :: nvars         ! number of variables

CONTAINS
  
  SUBROUTINE driver(dt, to, tf, atol, rtol, Ti, pi, xi)

    ! input/output variables
    real*8,  intent(inout)               :: dt     ! time  step
    real*8,  intent(in)                  :: to     ! init. time
    real*8,  intent(in)                  :: tf     ! final time
    real*8,  intent(in)                  :: atol   ! abs.  tol.
    real*8,  intent(in)                  :: rtol   ! rel.  tol.
    real*8,  intent(in)                  :: Ti     ! init. temperature
    real*8,  intent(in)                  :: pi     ! init. pressure
    real*8,  intent(inout), dimension(:) :: xi     ! init. condition
    ! local variables
    integer                              :: vio
    integer                              :: ver
    integer                              :: rkord
    real*8                               :: my_dt
    real*8                               :: y(nspec)
    real*8                               :: statevec(nspec+1)
        
    !=======================================================!
    ! initialize 
    my_dt = dt

    call initChem(Ti, pi, xi, y)
    
    statevec(1:nspec) = y
    statevec(nspec+1) = Ti
    
    nvars = nspec + 1
    ncons = nelem
    nslow = nvars-nelem-1
    nfast = nvars-nslow

    write(*,*)
    write(*,'(A10,i2)') "nvars = ",nvars
    write(*,'(A10,i2)') "ncons = ",ncons
    write(*,'(A10,i2)') "nslow = ",nslow
    write(*,'(A10,i2)') "nfast = ",nfast
    write(*,*)
    
    !=======================================================!
    ! integrate
    vio   = 1
    ver   = 1
    rkord = 4
    call cspIntegration(vio, ver, rkord, my_dt, to, tf, atol, rtol, statevec)

  END SUBROUTINE driver
  
  SUBROUTINE cspIntegration(vio, ver, rkord, dt, to, tf, atol, rtol, statevec)
    
    !% input/output variables
    integer, intent(in)    :: vio
    integer, intent(in)    :: ver
    integer, intent(in)    :: rkord
    real*8,  intent(inout) :: dt
    real*8,  intent(in)    :: to
    real*8,  intent(in)    :: tf
    real*8,  intent(in)    :: atol
    real*8,  intent(in)    :: rtol
    real*8,  intent(inout) :: statevec(nvars)
    !% local variables
    !  file name
    character(len=99)      :: runtemp
    character(len=99)      :: runphi
    character(len=99)      :: filename1
    character(len=99)      :: filename2
    character(len=99)      :: filename3
    character(len=99)      :: filename4
    !  integrator:
    integer                :: astat
    integer                :: i
    real*8                 :: tn
    real*8                 :: rhs(nvars)
    real*8                 :: jac(nvars,nvars)
    real*8                 :: qr(nvars,nvars)
    real*8                 :: ql(nvars,nvars)
    real*8                 :: pts(nvars)
    real*8                 :: Ti
    real*8                 :: ST
    real*8                 :: STmax
    real*8                 :: tign
    
    !=======================================================!
    ! IO
    if(vio .eq. 1) then
       write(runtemp,'(f14.0)') T
       runtemp   = runtemp(10:13)
       filename1 = trim(adjustl("outs/csp.sanDiego.")) // trim(adjustl(runtemp))
       filename1 = trim(adjustl(filename1))       // trim(adjustl("K.temp.bin"))

       filename2 = trim(adjustl("outs/csp.sanDiego.")) // trim(adjustl(runtemp))
       filename2 = trim(adjustl(filename2))       // trim(adjustl("K.nslow.bin"))
       
       open(unit=150,           &
            file=filename1,     &
            form="unformatted", &
            access="stream")

       open(unit=250,           &
            file=filename2,     &
            form="unformatted", &
            access="stream")
       
    end if
    
    !=======================================================!
    ! initialize
    Ti = T
    tn = to
    i  = 0

    if(vio .eq. 1) then
       write(150) tn, y, T
       write(250) nslow
    end if
    
    !=======================================================!
    ! ildm integration
100 continue
    
    ! subspaces
    call chemf(nvars, tn, statevec, rhs)
    call chemj(nvars, tn, statevec, 1, 1, jac, nvars)
    call cspSubspaces(atol, rtol, statevec, rhs, jac, qr, ql, dt)

    ! explicit time-integration
    call rungeKutta(rkord, tn, dt, qr, ql, statevec)

    y = statevec(1:nspec)
    T = statevec(nspec+1)
    
    ! advance
    tn   = tn + dt
    i    = i  + 1
    
    ! output
    if(vio .eq. 1) then
       write(150) tn, y, T
       write(250) nslow
    end if
    
    if(ver .eq. 1) then
       write(*,'(e14.4,2i4,12e14.6)') tn, nslow, nfast, dt, T
    end if
    
    if(tn .le. tf) goto 100
200 continue

    if(vio .eq. 1) then
       close(150)
       close(250)
    end if

    write(*,*)
    write(*,'(A20,e14.4)') "init  temp: ",Ti
    write(*,'(A20,e14.4)') "final temp:",T
    write(*,'(A20,i8)') "model time steps: ",i
    write(*,*)
    
  END SUBROUTINE cspIntegration

  SUBROUTINE cspSubspaces(atol, rtol, statevec, rhs, jac, qr, ql, dt)

    !% input/output variables
    real*8, intent(in)  :: atol
    real*8, intent(in)  :: rtol
    real*8, intent(in)  :: statevec(nvars)
    real*8, intent(in)  :: rhs(nvars)
    real*8, intent(in)  :: jac(nvars)
    real*8, intent(out) :: qr(nvars,nvars)
    real*8, intent(out) :: ql(nvars,nvars)
    real*8, intent(out) :: dt
    !% local variables
    integer             :: ns,nf,nc
    integer             :: i,j,k,l,m,n
    integer             :: expm
    integer             :: flag
    real*8              :: dstate(nvars)
    real*8              :: err(nvars)
    real*8              :: eig(nvars)
    real*8              :: tau(nvars)

    !=======================================================!
    ! eigendecomposition & time scales
    call eigenDec(nvars, jac, qr, ql, eig)
    tau = 1.0d0/abs(eig)
    
    !=======================================================!
    ! partition of bases
    nc  = ncons
    err = rtol * abs(statevec) + atol

    do i = 1,nvars-nc-1

       ns    = ncons+i
       nf    = nvars-ns
       expm  = 0
       flag  = 0
       dstate = -matmul(qr(:,ncons+1+i:nvars), &
            matmul(ql(ncons+1+i:nvars,:), rhs) / eig(ncons+1+i:nvars) )

       do j = 1,nf
          if(eig(nvars+1-j) .gt. 1.0d-10) then
             expm = 1
             exit
          end if
       end do

       do k = 1,nvars
          if(abs(dstate(k)) .gt. err(k)) then
             flag = 1
             exit
          end if
       end do

       if((expm .eq. 0) .and. (flag .eq. 0)) then
          exit
       end if

    end do

    nslow = ns-ncons
    nfast = nf+ncons
    dt    = max(tau(nslow+ncons) * 1.0d-01, 1.0d-09)
    
  END SUBROUTINE cspSubspaces

  SUBROUTINE cspDiagnostics(qr, ql, pts)

    !% input/output variables
    real*8, intent(in)  :: qr(nvars,nvars)
    real*8, intent(in)  :: ql(nvars,nvars)
    real*8, intent(out) :: pts(nvars)
    !% local variables
    integer             :: i
    real*8              :: pss(nvars,nvars)

    pss = matmul( qr(:,1:nslow+ncons), ql(1:nslow+ncons,:) )
    forall(i = 1:nvars) pts(i) = pss(i,i)
    
  END SUBROUTINE cspDiagnostics
  
  SUBROUTINE rungeKutta(rkord, tn, dt, qr, ql, statevec)

    !% input/output variables
    integer, intent(in)    :: rkord
    real*8,  intent(in)    :: tn
    real*8,  intent(in)    :: dt
    real*8,  intent(in)    :: qr(nvars,nvars)
    real*8,  intent(in)    :: ql(nvars,nvars)
    real*8,  intent(inout) :: statevec(nvars)
    !% local variables
    !  integrator:
    real*8,  parameter     :: OneSixth = 1.0d0/6.0d0
    real*8,  parameter     :: OneThird = 1.0d0/3.0d0
    real*8,  parameter     :: OneHalf  = 1.0d0/2.0d0
    real*8                 :: u(nslow+ncons)
    real*8                 :: k1(nslow+ncons)
    real*8                 :: k2(nslow+ncons)
    real*8                 :: k3(nslow+ncons)
    real*8                 :: k4(nslow+ncons)
    real*8                 :: statep(nvars)
    real*8                 :: qrs(nvars,nslow+ncons)
    real*8                 :: qls(nslow+ncons,nvars)

    !=======================================================!
    ! partitioned bases
    qrs = qr(:,1:nslow+ncons)
    qls = ql(1:nslow+ncons,:)
    
    !=======================================================!
    ! RK stages
    u = 0.0d0

    call slowf(tn, statevec, u, qrs, qls, k1)
    k1 = dt*k1

    u  = OneHalf * k1
    call slowf(tn, statevec, u, qrs, qls, k2)
    k2 = dt*k2

    select case (rkord)
       
    case (2)
       
       u  = k2
       statep = statevec + matmul(qrs, u)
       
    case (4)

       u  = OneHalf * k2
       call slowf(tn, statevec, u, qrs, qls, k3)
       k3 = dt*k3
       
       u  = k3
       call slowf(tn, statevec, u, qrs, qls, k4)
       k4 = dt*k4

       u  = OneSixth * k1 + OneThird * k2 + OneThird * k3 + OneSixth * k4
       statep = statevec + matmul(qrs, u)

    end select

    call cspcc(statep, statevec)
    
  END SUBROUTINE rungeKutta

  SUBROUTINE slowf(t, statevec, phi, qrs, qls, fs)

    !% input/output variables
    real*8,  intent(in)  :: t
    real*8,  intent(in)  :: statevec(nvars)
    real*8,  intent(in)  :: phi(nslow+ncons)
    real*8,  intent(in)  :: qrs(nvars,nslow+ncons)
    real*8,  intent(in)  :: qls(nslow+ncons,nvars)
    real*8,  intent(out) :: fs(nslow+ncons)
    !% local variables
    integer              :: i
    real*8               :: x(nvars)
    real*8               :: fx(nvars)
    
    !=======================================================!
    x  = statevec + matmul(qrs, phi)
    call chemf(nvars, t, x, fx)
    fs = matmul(qls, fx)
    
  END SUBROUTINE slowf

  SUBROUTINE cspcc(stateo, statem)

    !% input/output variables
    real*8,  intent(in)  :: stateo(nvars)
    real*8,  intent(out) :: statem(nvars)
    !  iterative process
    real*8               :: state(nvars)           ! state
    real*8               :: dstate(nvars)          ! increment
    real*8               :: rhs(nvars)             ! chemistry source
    real*8               :: jac(nvars,nvars)       ! chemistry jacobian
    real*8               :: gstate(nfast-ncons)    ! sim constraint
    real*8               :: qr(nvars,nvars)        ! right eigenvectors
    real*8               :: ql(nvars,nvars)        ! left  eigenvectors
    real*8               :: qrf(nvars,nfast-ncons) ! right fast eigenvectors
    real*8               :: qlf(nfast-ncons,nvars) ! left  fast eigenvectors
    real*8               :: eig(nvars)             ! eigenvalues
    !  linear solver
    integer              :: info
    integer              :: ifail
    integer              :: ipiv(nfast-ncons)

    !=======================================================!
    ! csp correction
    state = stateo

    call chemf(nvars, 0.0, state, rhs)
    call chemj(nvars, 0.0, state, 1, 1, jac, nvars)
    call eigenDec(nvars, jac, qr, ql, eig)

    qrf    = qr(:,nslow+ncons+1:nvars)
    qlf    = ql(nslow+ncons+1:nvars,:)
    gstate = matmul(qlf, rhs) / eig(nslow+ncons+1:nvars)
    dstate = -matmul(qrf, gstate)
    statem = stateo + dstate

  END SUBROUTINE cspcc

 END MODULE csp
