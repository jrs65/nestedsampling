PROGRAM test_el

! Include precision definitions
#include "defs.H95"

  USE ellipsoidutil
  USE rand

  IMPLICIT NONE

  INTEGER :: i

  REAL(rk), DIMENSION(3,3), TARGET :: cov
  REAL(rk), DIMENSION(:,:), POINTER :: covp
  REAL(rk), DIMENSION(3) :: point

  cov = 0
  cov(1,1) = 17.0/2.0
  cov(1,2) = -15.0/2.0
  cov(2,1) = -15.0/2.0
  cov(2,2) = 17.0/2.0
  cov(3,3) = 1.0

  covp => cov

  CALL set_covar(covp)
  CALL print_cov()
  CALL init_ellipsoidutil(3)
  CALL transmatrix()
  CALL print_trans()

  OPEN(1, FILE="sphere")

  DO i=1,2000
     
     CALL rng_sphere_nd(point,3)
     WRITE(1,*) point

  END DO

  OPEN(2, FILE="ellipsoid")

  DO i=1,2000
     
     CALL rand_ellipsoid(point)
     WRITE(2,*) point

  END DO

END PROGRAM test_el
