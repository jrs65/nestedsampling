PROGRAM test2d_el

! Include precision definitions
#include "defs.H90"

  USE nested_simple
  USE tree
  USE rand
  USE distr
  USE metropolis

  IMPLICIT NONE


  INTEGER, PARAMETER :: n = 200

  INTEGER :: i
  REAL(rk), DIMENSION(2,4), TARGET :: prior_pars
  REAL(rk), DIMENSION(:,:), POINTER :: prior_pars_p
  REAL(rk), DIMENSION(n) :: d1, d2
  REAL(rk), DIMENSION(2,n) :: data
  REAL(rk), DIMENSION(3,2) :: dgauss
  REAL(rk), DIMENSION(4), TARGET :: sc
  REAL(rk), DIMENSION(4) :: par
  REAL(rk), DIMENSION(4), TARGET :: ipar
  REAL(rk), DIMENSION(:), POINTER :: iparsp
  REAL(rk), DIMENSION(:), POINTER :: scp
  TYPE(pnode), DIMENSION(:), POINTER :: arr1, arr2
  
  ! Set up prior.
  prior_pars(1,1) = -10.0
  prior_pars(2,1) = 10.0
  prior_pars(1,2) = 0.0
  prior_pars(2,2) = 10.0
  prior_pars(1,3) = -10.0
  prior_pars(2,3) = 10.0
  prior_pars(1,4) = 0.0
  prior_pars(2,4) = 10.0

  prior_pars_p => prior_pars

  CALL set_prior_params(prior_pars_p)

  ipar = 1
  iparsp => ipar

  sc = 1
  scp => sc

  OPEN(1, FILE="data")

  CALL rng_gaussian_arr(d1, n)
  CALL rng_gaussian_arr(d2, n)

  DO i=1,n
     data(1,i) = d1(i)
     data(2,i) = d2(i)

     WRITE(1,*) data(:,i)
  END DO

  CALL datagaussian(data, 2, dgauss)
  
  CALL set_npoints(200)
  CALL set_nrepl(1000)
  CALL set_niter(50)
  CALL set_scale(scp)
  CALL set_initialpars(iparsp)

  CALL sample_points(dgauss, 2, 4, Lgaussian, Xuniform, \
                     initial_metropolis, replace_metropolis, arr1)

  OPEN(2, FILE="vals")
  CALL print_parray(arr1, 2)

  PRINT *,"Uniform"
  PRINT *,"Evidence",evidence(arr1)

  CALL dmean(dgauss)
  CALL mean_pars(arr1, 4, par)

  PRINT *,"Means",par

  PRINT *,""

END PROGRAM test2d_el
