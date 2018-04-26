! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

! A generalised examples file, allowing multiple Gaussians in an n-d
! flat-prior space. Set the separations (ensuring all within unit
! sphere), widths, and heights, and it will tell you the analytic
! value, and perform a sampling run through the space.

! Include precision definitions
#include "defs.H90"

!! Module containing functions required for program.
MODULE test_sk_gen_mod

USE sampling
USE tree
USE rand

IMPLICIT NONE

! Set model parameters.
! Number of dimensions.
INTEGER,  PARAMETER :: dim = 2
! Number of Gaussians
INTEGER, PARAMETER :: ng = 5
! Positions
REAL(rk), DIMENSION(ng,dim) :: pos
! Heights of Gaussians
REAL(rk), DIMENSION(ng) :: h
! Widths of Gaussinas
REAL(rk), DIMENSION(ng) :: sig

! Set constants.
REAL(rk), PARAMETER :: pi = 3.1415926535
REAL(rk) :: min_inf = -1.0 * HUGE(pi)

INTEGER :: i

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


  ! General Gaussian Likelihood function.
  REAL(rk) FUNCTION like(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    INTEGER ::i,j,n
    REAL(rk) :: r1, e1
    REAL(rk) :: r,e
    r = 0.0
    e = 0.0

    ! Find radial distance from each Gaussian centre
    DO i=1,ng

       r1 = 0.0

       DO j=1,dim
          r1 = r1 + (par(j) - pos(i,j))**2
       END DO

       e = e + h(i)*exp(-1.0 * r1 / (2*sig(i)**2))

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

       ev_analytic = ev_analytic &
            + h(i) * fact(dim / 2)*((2*sig(i)**2)**(dim/2))

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
  
END MODULE test_sk_gen_mod


!! Calculate the example given in Skilling's Paper, modified to give
!! two clusters seperated by 2*hs.
PROGRAM test_sk_gen

  USE test_sk_gen_mod
  USE nested_cluster
  USE ellipsoid
  USE metropolis

  IMPLICIT NONE

  REAL(rk), DIMENSION(2) :: e
  REAL(rk) :: ep
  TYPE(pset), POINTER :: cs

  REAL(rk) :: en
  en = 1.0

  ! Set Gaussian positions.
  pos(:,1) = (/0.45,-0.2,-0.35,-0.4,0.1/)
  pos(:,2) = (/0.1,0.15,0.2,-0.4,-0.15/)
  
  ! Set heights.
  h = (/0.6,0.8,1.0,0.5,0.5/)
  ! Set widths.
  sig = (/0.05,0.03,0.01,0.01,0.02/)


  ! Set up sampling.
  ! Number of active points.
  CALL set_npoints(300)
  ! Replacement number upper bound.
  CALL set_nrepl(100000)
  ! Number of MH iterations in Prior volume.
  CALL set_niter(5)
  ! Number of chains of volume reductions for evidence estimation.
  CALL set_nevid(30)
  ! Volume fraction remaing stopping criterion.
  CALL set_volfrac(real(1e-4,rk))
  ! Check for clusters every n replacements
  CALL set_cluster(5)
  ! Cluster when volume reduced by factor of.
  CALL set_volred(real(0.3,rk))
  ! Cluster when fraction of line of centres is in ellipsoids.
  CALL set_clustsep(real(1.0,rk))


  PRINT *,"=== Generalised Gaussians ==="
  PRINT *,"Analytic:",ev_analytic()


  CALL init_ellipsoid(dim, en, like , prior)

  NULLIFY(cs)
  
  ! Sample
  CALL sample_points(dim, like, prior, &
       initial_sphere, replace_ellipsoid, cs)
  
  e = evidence(cs)
  PRINT *,"Evidence Sampled:",e
  
  CALL evidence_simple(cs, ep)  
  PRINT *,"Evidence Simple:",ep

  ! Print sampling efficiency.
  CALL print_sampling_eff()
  PRINT *,"N_like:",r_lnum
    
END PROGRAM test_sk_gen
