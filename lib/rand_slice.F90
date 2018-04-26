!! A module for performing slice sampling, for drawing random numbers
!! form an arbritrary distribution. For details, see MacKay,
!! Information Theory, Inference and Learning Algorithms.
!! @author Richard Shaw
MODULE rand_slice

! Precision definitions.
#include "defs.H90"

USE rand

IMPLICIT NONE

!! Width interval to use.
REAL(rk), PRIVATE :: w = 0.025

CONTAINS



  !! Set the width for slice sampling.
  SUBROUTINE set_slice_width(wd)
    
    !! New width.
    REAL(rk), INTENT(IN) :: wd
    
    w = wd

  END SUBROUTINE set_slice_width
  


  !! Initialise a new slice sampling run.
  REAL(rk) FUNCTION init_slice(p)

    INTERFACE 
       !! Un-normalised probability distribution to draw from.
       REAL(rk) FUNCTION p(x)
         REAL(rk), INTENT(IN) :: x
       END FUNCTION p
    END INTERFACE

    INTEGER :: i

    init_slice = 0.5

    DO i=1,100

       init_slice = slice_sample1d(p, init_slice)

    END DO

  END FUNCTION init_slice



  !! 1-dimensional slice sampling routine. Pass in a distribution p,
  !! and current state x.
  REAL(rk) FUNCTION slice_sample1d(p, x)
    
    INTERFACE 
       !! Un-normalised probability distribution to draw from.
       REAL(rk) FUNCTION p(x)
         REAL(rk), INTENT(IN) :: x
       END FUNCTION p
    END INTERFACE
    
    !! Previous state to transition from.
    REAL(rk), INTENT(IN) :: x
    
    REAL(rk) :: xl, xr, xp, u, ps
    
    u = p(x)*rng_uniform()

    ! Create enclosing interval.
    xp  = rng_uniform()
    xl = x - xp * w
    xr = x + (1-xp) * w
    DO WHILE(p(xl) > u)
       xl = xl - w
    END DO

    DO WHILE(p(xr) > u)
       xr = xr + w
    END DO

    ! Loop until we've got out number.
    DO
       xp = xl + (xr - xl) * rng_uniform()
       ! If P*(x') > u' then exit, we've got out number!
       IF(p(xp) > u) THEN
          EXIT
       ! If not then change the interval (xl,xr)
       ELSE
          IF(xp < x) THEN
             xl = xp
          ELSE
             xr = xp
          END IF
       END IF
       
    END DO
          
    slice_sample1d = xp
  
  END FUNCTION slice_sample1d

END MODULE rand_slice
