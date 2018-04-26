! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

! Include precision definitions
#include "defs.H90"

!! Module containing functions required for program.
MODULE test_sk_comp_mod

USE sampling
USE tree
USE rand

IMPLICIT NONE

INTEGER, PARAMETER :: ng = 5

REAL(rk), DIMENSION(ng) :: p1
REAL(rk), DIMENSION(ng) :: p2

REAL(rk), DIMENSION(ng) :: h
REAL(rk), DIMENSION(ng) :: w

! Set constants.
REAL(rk), PARAMETER :: pi = 3.1415926535
REAL(rk) :: min_inf = -1.0 * HUGE(pi)

INTEGER :: i

! Set model parameters.
INTEGER,  PARAMETER :: dim = 2


CONTAINS



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
  REAL(rk) FUNCTION like(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,n
    REAL(rk) :: r1, e1
    REAL(rk) :: r,e
    r = 0.0
    e = 0.0
    n = size(par)

    ! Find n-d radius.
    DO i=3,n
       r = r + (par(i)**2)
    END DO
    
    ! Return log likelihood.
    DO i=1,ng
       r1 = r + (par(1) - p1(i))**2 + (par(2) - p2(i))**2
       e = e + h(i)*exp(-1.0 * r1 / (2*w(i)**2))
    END DO

    like = log(e)
       
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
    
    INTEGER :: i

    ev_analytic = 0.0

    DO i=1,ng
       ev_analytic = ev_analytic + h(i) * fact(dim / 2)*((2*w(i)**2)**(dim/2))
    END DO

    ev_analytic = log(ev_analytic)

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
  
END MODULE test_sk_comp_mod


!! Calculate the example given in Skilling's Paper, modified to give
!! two clusters seperated by 2*hs.
PROGRAM test_sk_comp

  USE test_sk_comp_mod
  USE nested_cluster
  USE ellipsoid
  USE metropolis

  IMPLICIT NONE

  REAL(rk), DIMENSION(2) :: ev, e
  REAL(rk) :: ep2, ep, epv, e2
  TYPE(pset), POINTER :: cs

  INTEGER, PARAMETER :: nl = 25
  INTEGER :: is

  REAL(rk) :: en
  en = 1.0
  e2 = 0.0
  ev = 0.0
  e = 0.0

  p1 = (/0.45,-0.2,-0.35,-0.4,0.1/)
  p2 = (/0.1,0.15,0.2,-0.4,-0.15/)
  
  h = (/0.6,0.8,1.0,0.5,0.5/)
  w = (/0.05,0.03,0.01,0.01,0.02/)


  ! Set up sampling.
  CALL set_npoints(300)
  CALL set_nrepl(100000)
  CALL set_niter(5)
  CALL set_nevid(30)
  CALL set_volfrac(real(1e-4,rk))
  CALL set_cluster(5)
  CALL set_volred(real(0.3,rk))

  CALL set_clustsep(real(1.0,rk))


!  PRINT *,"=== Complicated Gaussians ==="
!  PRINT *,"Analytic:",ev_analytic()


  CALL init_ellipsoid(dim, en, like , prior)

  DO is=1,nl

     NULLIFY(cs)

     ! Sample
     CALL sample_points(dim, like, prior, &
          initial_sphere, replace_ellipsoid, cs)

     PRINT *,"Samp"
     e = evidence(cs)
    
     PRINT *,e
     
     PRINT *,"Simp"
     CALL evidence_simple(cs, ep)

     PRINT *,ep

     ev = ev + e
     e2 = e2 + e(1)**2

     epv = epv + ep
     ep2 = ep2 + ep**2
     
  END DO
  
  PRINT "(2f20.10)",(ev(1) / nl) - ev_analytic(),ev(2) / sqrt(1.0*nl**3)
  PRINT *,sqrt((e2 / nl) - (ev(1) / nl)**2), ev(2) / nl

  PRINT "(f20.10)",(epv / nl) - ev_analytic()
  PRINT "(f20.10)",sqrt((ep2 / nl) - (epv / nl)**2)

  PRINT *,"Efficiency"

!  ev = evidence(cs)
!  PRINT *,"Evidence",ev(1),"+/-",ev(2)
  
  ! Print sampling efficiency.
  CALL print_sampling_eff()
  PRINT *,"N_like:",r_lnum / nl
!  PRINT *,"Outputting point tree...."
  !CALL print_pset(cs, "pset_comp.dat")
!  PRINT *,"Done"

    
END PROGRAM test_sk_comp
