! Include precision definitions
#include "defs.H90"

!! Module containing functions required for program.
MODULE test_sk_cl_mod

USE sampling
USE tree
USE rand

IMPLICIT NONE

! Set constants.
REAL(rk), PARAMETER :: pi = 3.1415926535
REAL(rk) :: min_inf = -1.0 * HUGE(pi)

INTEGER :: i

! Set model parameters.
INTEGER,  PARAMETER :: dim = 2
REAL(rk), PARAMETER :: sigma = 0.1
REAL(rk), PARAMETER :: hs = 0.5

CONTAINS

  ! Single Gaussian Likelihood function.
  REAL(rk) FUNCTION like1(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,n
    REAL(rk) :: r2
    r2 = 0.0
    n = size(par)
    
    ! Find n-d radius.
    DO i=1,n
       r2 = r2 + (par(i)**2)
    END DO
    
    ! Return log likelihood.
    like1 = -1.0 * r2 / (2 * sigma**2)
    
  END FUNCTION like1
  
  
  
  ! Two-Gaussian Likelihood function.
  REAL(rk) FUNCTION like2(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,n
    REAL(rk) :: r,r1,r2,e,e1,e2
    r = 0.0
    n = size(par)

    ! Find n-d radius.
    DO i=2,n
       r = r + (par(i)**2)
    END DO
    
    ! Return log likelihood.
    r1 = r + (par(1) - hs)**2
    r2 = r + (par(1) + hs)**2
    IF(r1 > 5 * r2) THEN
       like2 = -1.0 * r2 / (2 * sigma**2)
    ELSE IF(r2 > 5 * r1) THEN
       like2 = -1.0 * r1 / (2 * sigma**2)
    ELSE
       e1 = exp(-1.0 * r1 / (2 * sigma**2))
       e2 = exp(-1.0 * r2 / (2 * sigma**2))
       e = e1 + e2
       like2 = log(e)
       
    END IF
  END FUNCTION like2


  LOGICAL FUNCTION isinf(a) 
    
    REAL(rk), INTENT(IN) :: a
    
    REAL(rk) :: mi, pi
    INTEGER :: z = 0
    
    mi = -1.0 / z
    pi =  1.0 / z

    IF(a .EQ. mi .OR. a .EQ. pi) THEN
       isinf = .TRUE.
    ELSE
       isinf = .FALSE.
    END IF

  END FUNCTION isinf


  ! Three-Gaussian Likelihood function.
  REAL(rk) FUNCTION like3(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,n
    REAL(rk) :: r,r1,r2,r3,e1,e2,e3,e
    r = 0.0
    n = size(par)

    ! Find n-d radius.
    DO i=3,n
       r = r + (par(i)**2)
    END DO
    
    ! Return log likelihood.
    r1 = r + (par(1) - hs)**2 + (par(2))**2
    r2 = r + (par(1) + 0.5 * hs)**2 + (par(2) - (sqrt(3.0) * hs / 2.0))**2
    r3 = r + (par(1) + 0.5 * hs)**2 + (par(2) + (sqrt(3.0) * hs / 2.0))**2

    e1 = exp(-1.0 * r1 / (2 * sigma**2))
    e2 = exp(-1.0 * r2 / (2 * sigma**2))
    e3 = exp(-1.0 * r3 / (2 * sigma**2))
    e = e1 + e2 + e3
    like3 = log(e)
    IF(isnan(like3) .OR. isinf(like3)) THEN
       PRINT *,r1,r2,r3,e1,e2,e3
    END IF 
       
  END FUNCTION like3
  
  
  ! Uniform prior in the Unit sphere,
  REAL(rk) FUNCTION prior(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,n
    REAL(rk) :: r2
    r2 = 0.0
    n = size(par)
    
    ! Find n-d radius.
    DO i=1,n
       r2 = r2 + (par(i)**2)
    END DO
    
    ! Check if in unit n-sphere.
    IF(r2 > 1.0) THEN
       prior = min_inf
    ELSE
       prior = log(fact(n/2) / pi**(n/2))
    END IF

  END FUNCTION prior
  
  
  ! Factorial function.
  INTEGER FUNCTION fact(n) 
    
    INTEGER :: i,n
    fact = 1
    IF(n < 1) THEN 
       PRINT *,"n must be a positive integer."
    END IF
    
    DO i=1,n
       fact = fact * i
    END DO
    
  END FUNCTION fact
  
  

  ! Evidence value (analytic).
  REAL(rk) FUNCTION ev_analytic()
    
    ev_analytic = log(fact(dim / 2)*((2*sigma**2)**(dim/2)))
    
  END FUNCTION ev_analytic
  
  
  
  ! Spherical method for initial point selection.
  SUBROUTINE initial_sphere(tree, n, dpars, L, X)
    
    ! The point tree to add to.
    TYPE(pnode), POINTER :: tree
    ! Number of points to add.
    INTEGER, INTENT(IN) :: n
    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars
    
    INTERFACE
       ! The Likelihood function
       REAL(rk) FUNCTION L(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION L
       
       ! The Prior function
       REAL(rk) FUNCTION X(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION X
         
    END INTERFACE
    
    
    INTEGER :: i
    REAL(rk), DIMENSION(:), ALLOCATABLE :: point
    
    ALLOCATE(point(dpars))
    
    DO i=1,n
       
       CALL rng_sphere_nd(point, dpars)
       
       CALL insert_node(tree, X(point), L(point), point)
       
    END DO
      
  END SUBROUTINE initial_sphere
  
END MODULE test_sk_cl_mod


!! Calculate the example given in Skilling's Paper, modified to give
!! two clusters seperated by 2*hs.
PROGRAM test_sk_cl

  USE test_sk_cl_mod
  USE nested_cluster
  USE ellipsoid
  USE metropolis

  IMPLICIT NONE

  REAL(rk), DIMENSION(2) :: ev
  TYPE(pset), POINTER :: cs

  REAL(rk) :: en
  en = 1.0


  ! Set up sampling.
  CALL set_npoints(200)
  CALL set_nrepl(10000)
  CALL set_niter(40)
  CALL set_nevid(30)
  CALL set_volfrac(real(1e-2,rk))
  CALL set_volred(real(0.4,rk))
  CALL set_clustsep(real(0.4,rk))

  PRINT *,"=== Single Gaussian ==="
  PRINT *,"Analytic:",ev_analytic()
  CALL init_ellipsoid(dim, en, like1, prior)

  ! Sample
  CALL sample_points(dim, like1, prior, &
       initial_sphere, replace_ellipsoid, cs)

  PRINT *,"Uniform"
  ev = evidence(cs)
  PRINT *,"Evidence",ev(1),"+/-",ev(2)
  
  ! Print sampling efficiency.
  CALL print_sampling_eff()

  
  PRINT *,"=== Two Gaussians ==="
  PRINT *,"Analytic:",ev_analytic()+log(2.0)
  ! CALL set_cluster(0)
  CALL init_ellipsoid(dim, en, like2, prior)

  ! Sample
  CALL sample_points(dim, like2, prior, &
       initial_sphere, replace_ellipsoid, cs)

  PRINT *,"Uniform"
  ev = evidence(cs)
  PRINT *,"Evidence",ev(1),"+/-",ev(2)
  
  ! Print sampling efficiency.
  CALL print_sampling_eff()

  PRINT *,"Outputting point tree...."
  CALL print_pset(cs, "pset_2c.dat")
  PRINT *,"Done"


  PRINT *,"=== Three Gaussians ==="
  PRINT *,"Analytic:",ev_analytic()+log(3.0)
  CALL init_ellipsoid(dim, en, like3, prior)

  ! Sample
  CALL sample_points(dim, like3, prior, &
       initial_sphere, replace_ellipsoid, cs)

  PRINT *,"Uniform"
  ev = evidence(cs)
  PRINT *,"Evidence",ev(1),"+/-",ev(2)
  
  ! Print sampling efficiency.
  CALL print_sampling_eff()

  PRINT *,"Outputting point tree...."
  CALL print_pset(cs, "pset_3c.dat")
  PRINT *,"Done"

    
END PROGRAM test_sk_cl
