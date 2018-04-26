PROGRAM test1d

! Include precision definitions
#include "defs.H90"

  USE nested_simple
  USE tree
  USE rand
  USE distr
  USE metropolis

  IMPLICIT NONE

  INTEGER, PARAMETER :: ndata = 500

  INTEGER :: i
  REAL(rk), DIMENSION(2,2), TARGET :: prior_pars
  REAL(rk), DIMENSION(:,:), POINTER :: prior_pars_p
  REAL(rk), DIMENSION(ndata) :: d1
  REAL(rk), DIMENSION(1,ndata) :: data
  REAL(rk), DIMENSION(3,1) :: dgauss
  REAL(rk), DIMENSION(2,2), TARGET :: sc
  REAL(rk), DIMENSION(2) :: par
  REAL(rk), DIMENSION(2), TARGET :: ipars
  REAL(rk), DIMENSION(:,:), POINTER :: scp
  REAL(rk), DIMENSION(:), POINTER :: iparsp
  REAL(rk), DIMENSION(2) :: ev
  TYPE(pnode), DIMENSION(:), POINTER :: arr1

  sc(1,1) = 0.1
  sc(2,1) = 20
  sc(1,2) = 0.1
  sc(2,2) = 10
  ipars = 1
  scp => sc
  iparsp => ipars

  ! Set up prior.
  prior_pars(1,1) = -10.0
  prior_pars(2,1) = 10.0
  prior_pars(1,2) = 0.0
  prior_pars(2,2) = 10.0

  prior_pars_p => prior_pars
  CALL set_prior_params(prior_pars_p)

  ! Set up, and write out data.
  OPEN(1, FILE="data")

  CALL rng_gaussian_arr(d1, ndata)

  DO i=1,ndata
     data(1,i) = d1(i)

     WRITE(1,*) data(:,i)
  END DO

  CALL datagaussian(data, 1, dgauss)
  

  ! Set up sampling.
  CALL set_npoints(100)
  CALL set_nrepl(1000)
  CALL set_niter(50)
  CALL set_initialwidths(scp)
  CALL set_initialpars(iparsp)

  ! Burnin.
  CALL burnin_metropolis(2, Xuniform)

  ! Sample!
  CALL sample_points(dgauss, 1, 2, Lgaussian, Xuniform, \
                     initial_metropolis, replace_metropolis, arr1)

  OPEN(2, FILE="vals")
  CALL print_parray(arr1, 2)

  PRINT *,"Uniform"
  ev = evidence(arr1)
  PRINT *,"Evidence",ev(1),"+/-",ev(2)

  CALL dmean(dgauss)
  CALL mean_pars(arr1, 2, par)

  PRINT *,"Means",par

  PRINT *,"Max",(arr1(size(arr1)) % pvec)

  PRINT *,""

  CALL print_sampling_eff()

END PROGRAM test1d
