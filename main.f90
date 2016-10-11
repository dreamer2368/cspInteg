PROGRAM main

  USE csp
  
  IMPLICIT NONE
  ! local variables
  character(len=4)     :: Tin
  character(len=4)     :: pin
  integer              :: astat
  integer              :: ns, ne, nr
  real*8               :: dt
  real*8               :: to
  real*8               :: tf
  real*8               :: atol
  real*8               :: rtol
  real*8               :: pi
  real*8               :: Ti
  real*8               :: Ta
  real*8               :: phi
  real*8               :: nux
  real*8,  allocatable :: xi(:)

  !=======================================================!
  ! create new gas & allocate initial state
  call newGas(ns, ne, nr)

  if(allocated(xi))    deallocate(xi)
  allocate(xi(ns), stat=astat)
  if(astat .ne. 0) stop

  !=======================================================!
  ! Initialize
  call getarg(1, Tin)
  call getarg(2, pin)
  read(Tin,*) Ti
  read(pin,*) pi

  ! stoichiometry
  phi = 1.0d0

  ! initial conditions
  pi  = pi * 1.01325d+05
  xi  = 0.0d0 ! 1.0d-15

  !=======================================================!
  ! Hydrogen-Oxygen Combustion
  ! nux:
  nux   = 5.0d-01
  ! H2: 
  xi(1) = 3.0d-01
  ! O2: 
  xi(3) = nux * xi(1) / phi
  ! N2:
  xi(9) = 1.0d0 - (1.0d0 + nux / phi) * xi(1)
  
  dt   = 1.0d-12
  to   = 1.0d-12
  tf   = 1.0d-03
  rtol = 1.0d-03
  atol = 1.0d-10
  
  !=======================================================!
  ! integrate
  call driver(dt, to, tf, atol, rtol, Ti, pi, xi)

  !=======================================================!
  ! wrap up
  if(allocated(xi))    deallocate(xi)

END PROGRAM main
