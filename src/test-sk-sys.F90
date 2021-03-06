! Include precision definitions
#include "defs.H90"

! Module containing required functions.
MODULE test_sk_sys_mod

USE sampling
USE tree
USE rand

IMPLICIT NONE


  ! Set constants.
  REAL(rk), PARAMETER :: pi = 3.1415926535
  REAL(rk) :: min_inf = -1.0 * HUGE(pi)

  INTEGER :: i

  ! Set model parameters.
  INTEGER,  PARAMETER :: dim = 8
  REAL(rk), PARAMETER :: sigma = 0.05

CONTAINS

  ! Gaussian Likelihood function.
  REAL(rk) FUNCTION like(par)
    
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
    like = -1.0 * r2 / (2 * sigma**2)
    
  END FUNCTION like

  
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
  
END MODULE test_sk_sys_mod


!! Calculate the example given in Skilling's Paper.
PROGRAM test_sk_sys

  USE nested_cluster
  USE ellipsoid
  USE metropolis
  USE test_sk_sys_mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: nl = 100

  REAL(rk), DIMENSION(nl,2) :: ev
  TYPE(pset), POINTER :: cs
  REAL(rk), DIMENSION(2) :: e
  !INTEGER, DIMENSION(4) :: np
  INTEGER :: np
  INTEGER :: is,ip

  REAL(rk), DIMENSION(4) :: en
  np = 30

  en = (/1.0, 1.25, 1.5, 1.75/) 

  !np = (/20, 30, 40, 60/)

  ! Set up sampling.
  CALL set_nrepl(50000)
  CALL set_niter(5)
  CALL set_nevid(30)
  CALL set_volfrac(real(1e-2,rk))
  CALL set_cluster(0)

  PRINT *,"=== Single Gaussian ==="
  PRINT *,"Analytic:",ev_analytic()
  
  DO ip=1,size(en)

     CALL set_npoints(np)
     CALL init_ellipsoid(dim, en(i), like, prior)
  
     e = 0.0
     DO is=1,nl
        
        NULLIFY(cs)
        
        ! Sample
        CALL sample_points(dim, like, prior, &
             initial_sphere, replace_ellipsoid, cs)

        PRINT *,evidence(cs)
        
        e = e + evidence(cs)

     END DO
     
     PRINT "(i5, f20.10, 2f20.10)",dim,en(ip),(e(1) / nl) - ev_analytic(),e(2) / sqrt(1.0*nl**3)

  END DO
     
   CONTAINS


  ! Evidence value (analytic).
  REAL(rk) FUNCTION ev_analytic()
    
    ev_analytic = log(fact(dim / 2)*((2*sigma**2)**(dim/2)))
    
  END FUNCTION ev_analytic
  
  
END PROGRAM test_sk_sys

