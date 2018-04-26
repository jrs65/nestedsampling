#define rk 8

! A minimal 2d Gaussian example.
MODULE min_mod

IMPLICIT NONE

  ! Set constants.
  REAL(rk), PARAMETER :: pi = 3.1415926535
  REAL(rk) :: min_inf = -1.0 * HUGE(pi)


CONTAINS

  ! Gaussian Likelihood function.
  REAL(rk) FUNCTION like(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    REAL(rk) :: r2
    
    r2 = sum(par**2)
    
    ! Return log likelihood.
    like = -1.0 * r2 / 0.002
    
  END FUNCTION like

  
  ! Uniform prior in the Unit circle,
  REAL(rk) FUNCTION prior(par)
    
    ! Parameters 
    REAL(rk), DIMENSION(:) :: par
    
    REAL(rk) :: r2

    r2 = sum(par**2)
    
    ! Check if in unit n-sphere.
    IF(r2 > 1.0) THEN
       prior = min_inf
    ELSE
       prior = -1.0 * log(pi)
    END IF
    
  END FUNCTION prior
  
    
END MODULE min_mod


PROGRAM min

  ! Include nested module.
  USE nested_cluster
  USE min_mod

  IMPLICIT NONE

  ! Array to store evidence values in.
  REAL(rk), DIMENSION(2) :: ev
  ! Point set to store samples points in.
  TYPE(pset), POINTER :: cs


  ! Initial widths for Metropolis Burn-in
  REAL(rk), DIMENSION(2,2), TARGET :: sc
  REAL(rk), DIMENSION(:,:), POINTER :: scp
  ! Initial parameters vector for Metropolis burn-in.
  REAL(rk), DIMENSION(2), TARGET :: ipars
  REAL(rk), DIMENSION(:), POINTER :: iparsp

  ! Set parameters for Metropolis burn-in.
  sc(1,:) = 0.1
  sc(2,:) = 1.0
  ipars = 0.1
  scp => sc
  iparsp => ipars

  ! Set initial widths and parameters.
  CALL set_initialwidths(scp)
  CALL set_initialpars(iparsp)
  ! Burn-in the Metropolis sampler within the Prior.
  ! Pass in, parameters dimension, prior function.
  CALL burnin_metropolis(2, prior)

  ! Set number of active points.
  CALL set_npoints(50)

  ! Initialise the ellipsoidal sampler.  
  ! Need to pass, par dimension and enlargement factor of ellipse.
  CALL init_ellipsoid(2, 1.0D0)

  ! Generate the set of sampled points (stored in cs). Need to pass in
  ! data array, data dim, par dim, likelihood function, prior
  ! function, initial sampling routine, replacement sampling routine,
  ! and a point set to store result in.
  CALL sample_points(2, like, prior, &
       initial_metropolis, replace_ellipsoid, cs)

  ! Calculate the evidence value, and its error.
  ev = evidence(cs)
  PRINT *,"Evidence",ev(1),"+/-",ev(2)

  ! Print out the set of points to the given file.
  CALL print_pset(cs, "min_2d.dat")
  
END PROGRAM min

