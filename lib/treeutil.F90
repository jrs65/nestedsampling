! Copyright (C) 2006 Richard Shaw <jrs65 at cam ac uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

!! Tree utility functions. Functions for tree traversal, printing and
!! checking Red-Black status, finding and removing smallest node.
!! @author Richard Shaw
MODULE treeutil

! Precision definitions.
#include "defs.H90"

USE tree

IMPLICIT NONE

INTEGER, PRIVATE :: nb = -1

PRIVATE :: wcheck_tree

CONTAINS



  ! Subroutine remove_smallest:
  !! Remove the smallest (by likelihood value) node
  !! and return it.
  SUBROUTINE remove_smallest(tree, node)
    
    ! Tree to search
    TYPE(pnode), POINTER :: tree
    ! Removed node
    TYPE(pnode), POINTER :: node

    CALL find_smallest(tree, node)
    CALL remove_node(node, tree)

  END SUBROUTINE remove_smallest



  ! Subroutine find_smallest:
  !! Find the smallest (by likelihood value) node
  !! and return it.
  SUBROUTINE find_smallest(tree, node)
    
    ! Tree to search
    TYPE(pnode), POINTER :: tree
    ! Removed node
    TYPE(pnode), POINTER :: node
    ! Temporary node pointers
    TYPE(pnode), POINTER :: current1, current2
    
    IF(.not. ASSOCIATED(tree)) THEN
       PRINT *,"Tree not associated"
    END IF

    current1 => tree
    ! Loop through to find bottom node.
    DO WHILE(ASSOCIATED(current1 % left))
       current2 => current1
       current1 => current1 % left
    END DO

    node => current1
       
  END SUBROUTINE find_smallest



  ! Subroutine find_largest:
  !! Find the largest (by likelihood value) node
  !! and return it.
  SUBROUTINE find_largest(tree, node)
    
    ! Tree to search
    TYPE(pnode), POINTER :: tree
    ! Removed node
    TYPE(pnode), POINTER :: node
    ! Temporary node pointers
    TYPE(pnode), POINTER :: current1, current2
    
    IF(.not. ASSOCIATED(tree)) THEN
       PRINT *,"Tree not associated"
    END IF

    current1 => tree
    ! Loop through to find bottom node.
    DO WHILE(ASSOCIATED(current1 % right))
       current2 => current1
       current1 => current1 % right
    END DO

    node => current1
       
  END SUBROUTINE find_largest



  ! Function next_node:
  !! Return true if node has a successor, and if so, updates the node
  !! pointer to it. Use for looping through a tree in order.
  LOGICAL FUNCTION next_node(node, tree)

    TYPE(pnode), POINTER :: node
    TYPE(pnode), POINTER :: tree

    TYPE(pnode), POINTER :: n1

    next_node = .TRUE.

    ! If node null, start at beginning.
    IF(.NOT. ASSOCIATED(node)) THEN
       CALL find_smallest(tree, node)
    ! If we have a right node then return smallest in right subtree.
    ELSE IF(ASSOCIATED(node % right)) THEN
       n1 => node % right
       CALL find_smallest(n1, node)
    ! If we are a left subtree, return parent.
    ELSE IF(node % lr == LEFT) THEN
       node => node % parent
    ! If we are a right subtree, move upwards to find successor.
    ELSE IF(node % lr == RIGHT) THEN
       n1 => node
       DO WHILE(n1 % lr == RIGHT)
          n1 => n1 % parent
       END DO
       IF(n1 % lr == LEFT) THEN
          node => n1 % parent
       ELSE
          node => NULL()
          next_node = .FALSE.
       END IF

    ! Else, we have no successor, return false.
    ELSE
       node => NULL()
       next_node = .FALSE.
    END IF

  END FUNCTION next_node



  ! Function prev_node:
  ! Return true if node has a successor, and if so, updates the node
  ! pointer to it.
  LOGICAL FUNCTION prev_node(node, tree)

    TYPE(pnode), POINTER :: node
    TYPE(pnode), POINTER :: tree

    TYPE(pnode), POINTER :: n1

    prev_node = .TRUE.

    ! If node null, start at beginning.
    IF(.NOT. ASSOCIATED(node)) THEN
       CALL find_largest(tree, node)
    ! If we have a right node then return smallest in right subtree.
    ELSE IF(ASSOCIATED(node % left)) THEN
       n1 => node % left
       CALL find_largest(n1, node)
    ! If we are a left subtree, return parent.
    ELSE IF(node % lr == RIGHT) THEN
       node => node % parent
    ! If we are a left subtree, move upwards to find successor.
!!$    ELSE IF(node % lr == LEFT) THEN
!!$       n1 => node
!!$       DO WHILE(n1 % lr == LEFT)
!!$          n1 => n1 % parent
!!$       END DO
!!$       IF(n1 % lr == RIGHT) THEN
!!$          !node => n1 % parent
!!$       ELSE
!!$          !node => NULL()
!!$          next_node = .FALSE.
!!$       END IF

    ! Else, we have no successor, return false.
    ELSE
       node => NULL()
       prev_node = .FALSE.
    END IF

  END FUNCTION prev_node



  ! Subroutine walk_inorder:
  !! Walk the tree (in order), and apply subroutine f to each node.
  RECURSIVE SUBROUTINE walk_inorder(tr, f)

    ! The tree to be processed.
    TYPE (pnode), POINTER :: tr
    
    INTERFACE

       ! The function f to be applied.
       SUBROUTINE f(node)
         USE treedef     
         TYPE (pnode), POINTER :: node
       END SUBROUTINE f

    END INTERFACE

    ! Check to see if tree associated.
    IF(.NOT. ASSOCIATED(tr)) THEN
       PRINT *,"Tree unassociated"
       RETURN
    END IF

    ! If left tree exists, walk it.
    IF(ASSOCIATED(tr % left)) THEN
       CALL walk_inorder(tr % left, f)
    END IF
    
    ! Apply f to current node.
    CALL f(tr)

    ! If right tree exists, walk it.
    IF(ASSOCIATED(tr % right)) THEN
       CALL walk_inorder(tr % right, f)
    END IF
    
  END SUBROUTINE walk_inorder




  ! Subroutine walk_postorder:
  !! Walk the tree, and apply subroutine f, in a post-order fashion. f
  !! takes one argument, the current node.
  RECURSIVE SUBROUTINE walk_postorder(tr, f)

    ! The tree to be processed.
    TYPE (pnode), POINTER :: tr
    
    INTERFACE

       ! The function f to be applied.
       SUBROUTINE f(node)
         USE treedef     
         TYPE (pnode), POINTER :: node
       END SUBROUTINE f

    END INTERFACE

    ! Check to see if tree associated.
    IF(.NOT. ASSOCIATED(tr)) THEN
       PRINT *,"Tree unassociated"
       RETURN
    END IF

    ! If left tree exists, walk it.
    IF(ASSOCIATED(tr % left)) THEN
       CALL walk_postorder(tr % left, f)
    END IF
    
    ! If right tree exists, walk it.
    IF(ASSOCIATED(tr % right)) THEN
       CALL walk_postorder(tr % right, f)
    END IF

    ! Apply f to current node.
    CALL f(tr)

    
  END SUBROUTINE walk_postorder



  ! Subroutine print_node:
  !! Print the point node to IO channel 10.
  SUBROUTINE print_node(tr)

    ! The tree to be processed.
    TYPE(pnode), POINTER :: tr
    
    ! Print parameters, prior probability, likelihood.
    CALL print_node_fh(tr, 10)

  END SUBROUTINE print_node



  ! Subroutine print_node:
  !! Print the point node to IO channel fh.
  SUBROUTINE print_node_fh(tr, fh)

    ! The tree to be processed.
    TYPE(pnode), POINTER :: tr
    ! FH
    INTEGER, INTENT(IN) :: fh
    
    ! Print parameters, prior probability, likelihood.
    WRITE(fh,"(30G18.6)") (tr % pvec), (tr % xval), (tr % axmval), (tr % w), (tr % lval)

  END SUBROUTINE print_node_fh



  ! Subroutine print_pnode:
  !! Print a point node to standard out.
  SUBROUTINE print_pnode(tr)

    ! The tree to be processed.
    TYPE(pnode), POINTER :: tr
    
    ! Print parameters, prior probability, likelihood.
    PRINT *,(tr % lval),tr % rb, tr % lr

  END SUBROUTINE print_pnode



  ! Subroutine print_ptree:
  !! Print the entire pnode tree into file.
  SUBROUTINE print_ptree(tr, file)
    ! Array of point nodes to use.
    TYPE(pnode), POINTER :: tr
    ! Filename to write into.
    CHARACTER(LEN=*), INTENT(IN) :: file


    OPEN(10, FILE=file)
    CALL walk_inorder(tr, print_node)
    CALL flush(10)
    CLOSE(10)

  END SUBROUTINE print_ptree



  ! Subroutine check_tree:
  !! Check the tree for violations of the Red-black properties, and
  !! give some detail on where they are.
  SUBROUTINE check_tree(node) 

    TYPE (pnode), POINTER :: node

    nb = -1

    CALL wcheck_tree(node, 0)

  END SUBROUTINE check_tree


  
  ! Subroutine wcheck_tree:
  !! A recursive private work routine for checking the consistency of
  !! the tree Red-Black properties.
  RECURSIVE SUBROUTINE wcheck_tree(node, cb)

    TYPE (pnode), POINTER :: node
    INTEGER, INTENT(IN) :: cb

    TYPE(pnode), POINTER :: n1
    INTEGER :: cb1

    cb1 = cb

    ! Check links are correct.
    IF(ASSOCIATED(node % parent)) THEN
       IF(node % lr == LEFT) THEN
          n1 => node % parent % left
       ELSE IF(node % lr == RIGHT) THEN
          n1 => node % parent % right
       END IF

       IF(n1 % lval /= node % lval) THEN
          PRINT *,"*** Error: broken links ***",node % lval, n1 % lval
       END IF

       ! Check for colour violations.
       IF(node % rb == RED .AND. node % parent % rb == RED) THEN
          PRINT *,"*** Error: Incorrect colours ***",node % lval
       END IF

    ELSE
       IF(node % rb == RED) THEN
          PRINT *,"*** Error: root node must be coloured black. ***",node % lval
       END IF
    END IF

    ! If current node black, increment black node counter.
    IF(node % rb == BLACK) THEN
       cb1 = cb1 + 1
    END IF
    
    IF((.NOT. ASSOCIATED(node % left)) .AND. &
         (.NOT. ASSOCIATED(node % right))) THEN

       ! Check for equal black nodes property.
       IF(nb /= -1) THEN
          IF(cb1 /= nb) THEN
             PRINT *,"*** Error: Not equal numbers of black nodes in path. ***",node % lval
          END IF
       ELSE
          nb = cb1
       END IF
       
    ELSE
       ! Recurse down subtrees
       IF(ASSOCIATED(node % left)) THEN
          CALL wcheck_tree(node % left, cb1)
       END IF

       IF(ASSOCIATED(node % right)) THEN
          CALL wcheck_tree(node % right, cb1)
       END IF
    END IF

  END SUBROUTINE wcheck_tree

END MODULE treeutil
