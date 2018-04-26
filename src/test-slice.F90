PROGRAM test_slice

#include "defs.H90"
  
  USE rand_nested

  IMPLICIT NONE

  INTEGER, PARAMETER :: ns = 1500
  INTEGER, PARAMETER :: n = 50
  INTEGER, PARAMETER :: m = 25
  
  INTEGER :: i
  REAL(rk), DIMENSION(ns) :: arr

  CALL sample_volreduce_n(arr, n, m)

  OPEN(10, FILE="sl_test")

  DO i=1,ns
     WRITE(10,*) arr(i)
  END DO

END PROGRAM test_slice
