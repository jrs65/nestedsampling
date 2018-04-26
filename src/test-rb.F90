PROGRAM test_rb


! Include precision definitions
#include "defs.H90"

USE tree
USE treedef
USE treeutil

IMPLICIT NONE

INTEGER, PARAMETER :: n = 30

TYPE(pnode), POINTER :: node, n1
INTEGER :: i,d
REAL(rk), DIMENSION(2) :: pv
REAL(rk) :: t

pv = 0
d = 0 

NULLIFY(node)

PRINT *,"=== Short Test ==="

t = 1.0
CALL insert_node(node, t, t, pv)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="
t = 2.0
CALL insert_node(node, t, t, pv)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="
t = 3.0
CALL insert_node(node, t, t, pv)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="
CALL remove_smallest(node, n1)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="
t = 4.0
CALL insert_node(node, t, t, pv)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="
!CALL rotate_node(node % left, node)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="

NULLIFY(node)

PRINT *,"=== Long test ==="
DO i=1,n
   t = i

   CALL insert_node(node, t, t, pv)
   PRINT *,"=== Addition ===",i
   CALL check_tree(node)
   PRINT *,"=== End of Tree ==="
END DO

CALL walk_postorder(node, print_pnode)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="

NULLIFY(n1)

DO WHILE(next_node(n1, node))
   PRINT *,n1 % lval
END DO

DO i=1,(n/2)
   CALL remove_smallest(node, n1)
   CALL walk_postorder(node, print_pnode)
   CALL check_tree(node)
   PRINT *,"=== End of Tree ==="
END DO
CALL walk_postorder(node, print_pnode)
CALL check_tree(node)
PRINT *,"=== End of Tree ==="

NULLIFY(n1)

DO WHILE(next_node(n1, node))
   PRINT *,n1 % lval
END DO

END PROGRAM test_rb
