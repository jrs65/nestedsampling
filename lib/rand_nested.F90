! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

!! Module containing some rand number generation routines for use
!! during nested sampling.
MODULE rand_nested

#include "defs.H90"

  USE rand
  USE rand_slice

  IMPLICIT NONE

  !! Slice sampler bootstrapped?
  LOGICAL, PRIVATE :: slice_init = .FALSE.

  !! Storage for generic volume reduction.
  INTEGER, PRIVATE :: ns, ms
  REAL(rk), PRIVATE :: xs
  CONTAINS



    !! Sample from distribution of t-factors. p(t) = N * t^(N-1)
    !! Gives cumulative distr F(t) = t^N.
    SUBROUTINE sample_volreduce(t_arr, N) 

      REAL(rk), INTENT(OUT), DIMENSION(:) :: t_arr
      INTEGER, INTENT(IN) :: N

      INTEGER :: num
      num = size(t_arr)

      CALL rng_uniform_arr(t_arr, num)

      t_arr = (t_arr**(1.0/N))

    END SUBROUTINE sample_volreduce



    !! Return the volume fraction of a partition of n points, into an
    !! m, and an (n-m) cluster. Returns fraction corresponding to m
    !! point cluster. Calculate using slice sampling, a MCMC
    !! technique.
    SUBROUTINE sample_volreduce_n(t_arr, n, m)

      REAL(rk), INTENT(OUT), DIMENSION(:) :: t_arr
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: m

      INTEGER, PARAMETER :: skip = 10

      REAL(rk) :: x
      INTEGER :: i,j

      ! Re-initialize sampler (if required).
      IF(ns /= n .OR. ms /= m .OR. .NOT. slice_init) THEN
         ns = n; ms = m
         x = init_slice(wvolreduce_n)
         slice_init = .TRUE.
      ELSE
         x = xs
      END IF

      DO i=1,size(t_arr)

         DO j=1,skip
            x = slice_sample1d(wvolreduce_n, x)
         END DO

         t_arr(i) = x

      END DO

      xs = x

    END SUBROUTINE sample_volreduce_n

      

    REAL(rk) PURE FUNCTION wvolreduce_n(x)

      REAL(rk), INTENT(IN) :: x

      IF(x > 1 .OR. x < 0) THEN
         wvolreduce_n = 0.0
      ELSE
         wvolreduce_n = ((1-x)**(ns-ms))*(x**(ms))
      END IF

    END FUNCTION wvolreduce_n

      

      

      
  END MODULE rand_nested
