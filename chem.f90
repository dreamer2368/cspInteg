MODULE chem

  USE eigen
  
  !global variables
  logical             :: not_init_chem
  integer             :: nspec
  integer             :: nelem
  integer             :: nreac
  real*8, allocatable :: mw(:)
  real*8, allocatable :: y(:)
  real*8              :: p
  real*8              :: h
  real*8              :: T
  real*8              :: rho
  
CONTAINS

  SUBROUTINE newGas(ns, ne, nr)

    ! input/output variables
    integer, intent(out) :: ns
    integer, intent(out) :: ne
    integer, intent(out) :: nr
    ! local variables
    integer              :: astat, s, r, e


    !=======================================================!
    ! no. species, elements and reactions   
    call nSpecies(nspec)
    call nReactions(nreac)
    call nElements(nelem)

    if(allocated(mw)) deallocate(mw)
    if(allocated(y))  deallocate(y)
    allocate(mw(nspec), y(nspec), stat=astat)
    if(astat .ne. 0) stop

    call getMolecularWeights(mw)

    ns = nspec
    ne = nelem
    nr = nreac

  END SUBROUTINE newGas

  SUBROUTINE initChem(Ti, pi, xi, yi)

    ! input/output variables
    real*8, intent(in)  :: Ti        ! init temperature
    real*8, intent(in)  :: pi        ! init pressure
    real*8, intent(in)  :: xi(nspec) ! init mole fractions
    real*8, intent(out) :: yi(nspec) ! init mass fractions
    ! local variables
    integer             :: astat, i
    
    !=======================================================!
    ! flag is on
    not_init_chem = .true.
    
    !=======================================================!
    ! thermo state
    call getMassFractions(xi, yi)
    call getEnthalpyMass(Ti, yi, h)
    call getDensity(Ti, pi, yi, rho)

    T = Ti
    p = pi
    y = yi

    write(*,*)
    write(*,'(A30)') "Initial Thermodynamic State    "
    write(*,*)
    write(*,'(A6,e14.4,A10)')  "T   = ", T,             " [K]"
    write(*,'(A6,e14.4,A10)')  "p   = ", p/1.01325d+05, " [atm]"
    write(*,'(A6,e14.4,A10)')  "h   = ", h,             " [J/kg]"
    write(*,'(A6,e14.4,A10)')  "rho = ", rho,           " [kg/m^3]"
    write(*,*)
    write(*,'(A6,9e14.4,A10)') "y   = ", yi
    write(*,*)

    !=======================================================!
    ! turn flag off
    not_init_chem = .false.

    !=======================================================!
    ! deallocation
    
  END SUBROUTINE initChem

  SUBROUTINE chemf(n, tn, statevec, rhs)

    ! input/output variables
    integer, intent(in)  :: n       
    real*8,  intent(in)  :: tn      
    real*8,  intent(in)  :: statevec(n) ! mass fractions & temperature
    real*8,  intent(out) :: rhs(n)      ! chemistry source
    ! local variables
    
    !=======================================================!
    ! mixture state
    call getSource(p, statevec, rhs)
    
  END SUBROUTINE chemf

  SUBROUTINE chemj(n, tn, statevec, ml, mu, jac, lddfy)

    ! input/output variables
    integer, intent(in)  :: n
    integer, intent(in)  :: ml
    integer, intent(in)  :: mu
    integer, intent(in)  :: lddfy
    real*8,  intent(in)  :: tn
    real*8,  intent(in)  :: statevec(n)
    real*8,  intent(out) :: jac(lddfy,n)

    integer :: i
    
    !=======================================================!
    ! jacobian
    call getJacobian(p, statevec, jac)
    
  END SUBROUTINE chemj

END MODULE chem
