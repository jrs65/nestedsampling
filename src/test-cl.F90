PROGRAM test_cl


! Include precision definitions
#include "defs.H90"

USE rand
USE twomeans
USE tree
USE treedef

IMPLICIT NONE

INTEGER, PARAMETER :: n = 20
REAL(rk), PARAMETER :: hs = 4
TYPE(pnode), POINTER :: tr, t1, t2, t3, t4, t5, t6
INTEGER :: i,d,n1,n2
REAL(rk), DIMENSION(2) :: pv
REAL(rk) :: t
REAL(rk), DIMENSION(6*n) :: g_arr

pv = 0
d = 0 

NULLIFY(tr, t1, t2, t3, t4, t5, t6)

CALL init_cluster(2)

PRINT *,"=== 2 clusters ==="

! Add cluster 1 points.
CALL rng_gaussian_arr(g_arr, 4*n)
DO i=1,n
   t = i
   pv = g_arr(2*i-1:2*i)
   pv(1) = pv(1) + hs
   CALL insert_node(tr, t, t, pv)
END DO

! Add cluster 2 points.
DO i=n+1,2*n
   t = i+0.5
   pv = g_arr(2*i-1:2*i)
   pv(1) = pv(1) - hs
   CALL insert_node(tr, t, t, pv)
END DO

CALL print_ptree(tr, "c2_entire.dat")

IF(two_cluster(tr, t1, t2, n1, n2)) THEN
   PRINT *,"* Successfully clustered."
   CALL print_ptree(t1, "c2_t1.dat")
   CALL print_ptree(t2, "c2_t2.dat")
ELSE
   PRINT *,"* Clustering failed."   
END IF


PRINT *,"=== 3 way equilateral (two pass) ==="

NULLIFY(tr, t1, t2)

! Add cluster 1 points.
CALL rng_gaussian_arr(g_arr, 6*n)
DO i=1,n
   t = i
   pv = g_arr(2*i-1:2*i)
   pv(1) = pv(1) + hs
   CALL insert_node(tr, t, t, pv)
END DO

! Add cluster 2 points.
DO i=n+1,2*n
   t = i+0.5
   pv = g_arr(2*i:2*i+1)
   pv(1) = pv(1) - hs
   CALL insert_node(tr, t, t, pv)
END DO

! Add cluster 3 points.
DO i=2*n+1,3*n
   t = i+0.75
   pv = g_arr(2*i-1:2*i)
   pv(2) = pv(2) + hs*sqrt(3.0)
   CALL insert_node(tr, t, t, pv)
END DO
CALL print_ptree(tr, "c3_entire.dat")

IF(two_cluster(tr, t1, t2, n1, n2)) THEN
   PRINT *,"* Successfully clustered (step 1)."
   CALL print_ptree(t1, "c3_t1.dat")
   CALL print_ptree(t2, "c3_t2.dat")

   IF(two_cluster(t1, t3, t4, n1, n2)) THEN
      PRINT *,"* Successfully clustered (step 2i)."
      CALL print_ptree(t3, "c3_t3.dat")
      CALL print_ptree(t4, "c3_t4.dat")
   ELSE
      PRINT *,"* Clustering failed (step 2i)."   
   END IF

   IF(two_cluster(t2, t5, t6, n1, n2)) THEN
      PRINT *,"* Successfully clustered (step 2ii)."
      CALL print_ptree(t5, "c3_t5.dat")
      CALL print_ptree(t6, "c3_t6.dat")
   ELSE
      PRINT *,"* Clustering failed (step 2ii)."   
   END IF
ELSE
   PRINT *,"* Clustering failed (step 1)."   
END IF
   
END PROGRAM test_cl


