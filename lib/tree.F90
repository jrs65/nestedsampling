! General tree stuff. Functions for addition, removal, deletion
! etc. Now updated to use a Red-Black Balanced Binary tree.  For
! details see, any reference on Data Structures, or
! http://sage.mc.yu.edu/kbeen/teaching/algorithms/resources/red-black-tree.html
MODULE tree

! Include precision definitions
#include "defs.H90"

USE treedef

IMPLICIT NONE

PRIVATE :: fix_rb_add, fix_rb_del, rotate_node, sibling_node

CONTAINS


  ! Subroutine insertnode:
  ! Create and insert a node into the given tree, if tree unassociated,
  ! create tree (i.e return node).
  SUBROUTINE insert_node(tree, xval, lval, pvec)

    ! Tree to search
    TYPE(pnode), POINTER :: tree
    ! Prior value
    REAL(rk), INTENT(IN) :: xval
    ! Likelihood value
    REAL(rk), INTENT(IN) :: lval
    ! Parameters vector
    REAL(rk), INTENT(IN), DIMENSION(:) :: pvec

    TYPE(pnode), POINTER :: node

    ! Create pnode node
    ALLOCATE(node)
    node % lval = lval
    node % xval = xval
    ALLOCATE(node % pvec(size(pvec)))
    node % pvec = pvec

    ! Add the node to the tree.
    CALL add_node(node, tree)

  END SUBROUTINE insert_node


  
  ! Subroutine add_node:
  ! Add a node to the tree.
  SUBROUTINE add_node(node, tree)

    ! Node to add into tree.
    TYPE(pnode), POINTER :: node
    ! Tree root node.
    TYPE(pnode), POINTER :: tree
    
    TYPE(pnode), POINTER :: current
    node % rb = RED

    NULLIFY(node % parent, node % left, node % right)

    IF(ASSOCIATED(tree)) THEN
       current => tree

       ! Loop through to find insertion point.
       DO WHILE(ASSOCIATED(current))
          
          IF((node % lval) < (current % lval)) THEN
             IF(ASSOCIATED(current % left)) THEN
                current => current % left
             ELSE
                current % left => node
                node % parent => current
                node % lr = LEFT
                EXIT
             END IF
          ELSE
             IF(ASSOCIATED(current % right)) THEN
                current => current % right
             ELSE
                current % right => node
                node % parent => current
                node % lr = RIGHT
                EXIT
             END IF
          END IF
       END DO
    ELSE
       tree => node
       node % lr = ROOT
       NULLIFY(node % left, node % right, node % parent)
    END IF

    CALL fix_rb_add(node, tree)

  END SUBROUTINE add_node



  ! Subroutine remove_node:
  ! Remove node from tree,
  SUBROUTINE remove_node(node1, tree)

    ! Node to rotate on.
    TYPE(pnode), POINTER :: node1
    ! Tree root node.
    TYPE(pnode), POINTER :: tree

    ! Temp. node to use
    TYPE(pnode), POINTER :: t1, node, r, dummy

    INTEGER(1) :: rb, dm

    node => node1

    rb = node % rb

    ! If we have a left subtree
    IF(ASSOCIATED(node % left)) THEN

       ! and a right one, find the successor node
       IF(ASSOCIATED(node % right)) THEN
          t1 => node % right
          DO WHILE(ASSOCIATED(t1 % left))
             t1 => t1 % left
          END DO
          
          r => t1 % parent
          rb = t1 % rb

          IF(ASSOCIATED(t1 % right)) THEN
             t1 % right % parent => t1 % parent
             t1 % parent % left => t1 % right
          END IF
          
          t1 % left => node % left
          node % left % parent => t1
          
       ELSE
          t1 => node % left
          r => t1
       END IF

    ELSE
       
       IF(ASSOCIATED(node % right)) THEN
          t1 => node % right
          r => t1
       ELSE
          t1 => NULL()
          ! Errm, cunning hack. Need to fix on NULL, so rather than
          ! adding dummy black node, fix on node, and then cut final
          ! links (node % parent).
          r => node
       END IF

    
    END IF

    IF(ASSOCIATED(t1)) THEN
       IF(node % lr == LEFT) THEN
          node % parent % left => t1
          t1 % lr = LEFT
       ELSE IF(node % lr == RIGHT) THEN
          node % parent % right => t1
          t1 % lr = RIGHT
       ELSE
          tree => t1
          t1 % lr = ROOT
       END IF

       t1 % parent => node % parent
    ELSE
       IF(node % lr == LEFT) THEN
          node % parent % left => NULL()
       ELSE IF(node % lr == RIGHT) THEN
          node % parent % right => NULL()
       ELSE
          tree => NULL()
       END IF
    END IF

    IF(rb == BLACK) THEN
       CALL fix_rb_del(r, tree)
    END IF

    NULLIFY(node % parent, node % left, node % right)


  END SUBROUTINE remove_node



  ! Subroutine rotate_node:
  ! Node rotation. Operations on RB trees best expressed in terms of
  ! rotations around nodes.
  SUBROUTINE rotate_node(node1, tree)
    
    ! Node to rotate on.
    TYPE(pnode), POINTER :: node1
    ! Tree root node.
    TYPE(pnode), POINTER :: tree
    ! Temporary node pointer.
    TYPE(pnode), POINTER :: t1, node

    node => node1

    ! Only rotate if we have a left child.
    IF(.NOT. ASSOCIATED(node)) THEN
       PRINT *,"Cannot rotate on root node."
       RETURN
    END IF

    ! The new sub tree root.
    IF(ASSOCIATED(node % parent)) THEN
       t1 => node % parent
    ELSE
       nonfatal("Error: node does not have a valid parent")
       t1 => NULL()
    END IF
       
    ! If not the whole tree root.
    IF(ASSOCIATED(t1 % parent)) THEN

       ! Update link into subtree. Need to check if we are l/r.
       IF(t1 % lr == LEFT) THEN
          t1 % parent % left => node
       ELSE IF(t1 % lr == RIGHT) THEN
          t1 % parent % right => node
       END IF
       
       node % parent => t1 % parent
    ELSE
       NULLIFY(node % parent)
       tree => node
    END IF
       
    ! Check which sort of subtree we are, then rotate accordingly.
    IF(node % lr == LEFT) THEN

       ! Swap over subtree of new root, onto old root.
       IF(ASSOCIATED(node % right)) THEN
          t1 % left => node % right
          node % right % parent => t1
          t1 % left % lr = LEFT
       ELSE
          NULLIFY(t1 % left)
       END IF

       ! Update lr flags for nodes.
       node % lr = t1 % lr
       t1 % lr = RIGHT
       ! Update link between rotated nodes.
       node % right => t1
       t1 % parent => node
       
    ELSE IF(node % lr == RIGHT) THEN

       ! Swap over subtree of new root, onto old root.
       IF(ASSOCIATED(node % left)) THEN
          t1 % right => node % left
          node % left % parent => t1
          t1 % right % lr = RIGHT
       ELSE
          NULLIFY(t1 % right)
       END IF
       
       ! Update lr flags for nodes.
       node % lr = t1 % lr
       t1 % lr = LEFT
       ! Update link between rotated nodes.
       node % left => t1
       t1 % parent => node

    END IF

  END SUBROUTINE rotate_node
          


  ! Subroutine fix_rb_add:
  ! Private routine to Fix tree after insertion such that it is
  ! Red/Black again.
  RECURSIVE SUBROUTINE fix_rb_add(node, tree)

    ! Node to start at.
    TYPE(pnode), POINTER :: node
    ! Tree root node.
    TYPE(pnode), POINTER :: tree
    ! Temp node.
    TYPE(pnode), POINTER :: n1


    ! If node is black, what are we doing here?
    IF(node % rb == BLACK) THEN
       RETURN
    END IF

    ! If start node is root, change colour and exit.
    IF(node % lr == ROOT) THEN
       node % rb = BLACK
       RETURN
    END IF

    IF(.NOT. ASSOCIATED(node % parent)) THEN
       PRINT *,"Shouldn't get here"
       RETURN
    END IF


    ! Check if parent red. If not, nothing to do.
    IF(node % parent % rb == RED) THEN

       ! Check cases
       ! Case i: Node and parent have same sidededness.
       ! Fix with one rotation plus colour change. Recurse upwards.
       IF(node % parent % lr == node % lr) THEN
          
          n1 => node % parent

          ! Rotate, cunningly, node % parent remains the same
          CALL rotate_node(node % parent, tree)
          node % rb = BLACK

          ! Try to fix tree again.
          CALL fix_rb_add(n1, tree)

       ! Case ii: Node and parent have different sidededness.
       ! Fix with two rotations plus colour change. Recurse upwards.
       ELSE

          ! Fix colour first, as we lose which node to change w/
          ! rotations.
          node % parent % rb = BLACK

          n1 => node

          ! Rotate
          CALL rotate_node(node, tree)
          CALL rotate_node(node, tree)

          ! Start fixing another level up.
          CALL fix_rb_add(n1, tree)
          
       END IF

    END IF          
    
  END SUBROUTINE fix_rb_add



  ! Function sibling node:
  ! Return the sibling of a node, if it exists.
  FUNCTION sibling_node(node)

    ! Node to start at.
    TYPE(pnode), POINTER :: node
    TYPE(pnode), POINTER :: sibling_node
    TYPE(pnode), POINTER :: s1

    sibling_node => NULL()

    IF(node % lr == LEFT) THEN
       sibling_node => node % parent % right
    ELSE IF(node % lr == RIGHT) THEN
       sibling_node => node % parent % left
    END IF

  END FUNCTION sibling_node

    


  ! Subroutine fix_rb_del:
  ! Fix tree into Red/Black form after deletion.
  RECURSIVE SUBROUTINE fix_rb_del(node, tree)
    
    ! Node to start at.
    TYPE(pnode), POINTER :: node
    ! Tree root node.
    TYPE(pnode), POINTER :: tree

    ! Temp node
    TYPE(pnode), POINTER :: s1,s2

    ! Case i: node is red
    ! We fix just by swapping its colour.
    IF(node % rb == RED) THEN
       node % rb = BLACK
       RETURN
    END IF

    ! Case ii: node is the root node
    ! We've propagated to the root, tree is fixed.
    IF(node % lr == ROOT) THEN
       node % rb = BLACK
       RETURN
    END IF

    s1 => sibling_node(node)

    IF(ASSOCIATED(s1)) THEN

       IF(s1 % rb == BLACK) THEN
          
          ! Case iii: Sibling and both nephew black (remember null is black).
          ! Swap sibling colour to red.
          ! Rely heavily on .OR. not evaluating second argument if
          ! first is true.
          IF(((.NOT. ASSOCIATED(s1 % left)) .OR. s1 % left % rb == BLACK) .AND. &
               ((.NOT. ASSOCIATED(s1 % right)) .OR. s1 % right % rb == BLACK)) THEN

             ! Swap colour.
             s1 % rb = RED
             ! Recurse up to next level.
             CALL fix_rb_del(node % parent, tree)

          ELSE

             ! Case iv: One or more of nephews is red.
             ! Subcase i: Near nephew is red, colour black, rotate twice.
             ! Subcase ii: Far nephew is red, colour black, rotate once.
             IF(node % lr == LEFT) THEN              
                IF(ASSOCIATED(s1 % left) .AND. s1 % left % rb == RED) THEN
                   s2 => s1 % left
                   s2 % rb = node % parent % rb
                   node % parent % rb = BLACK
                   CALL rotate_node(s2, tree)
                   CALL rotate_node(s2, tree)                   
                ELSE
                   s2 => s1 % right
                   s2 % rb = BLACK
                   s1 % rb = node % parent % rb
                   node % parent % rb = BLACK
                   CALL rotate_node(s1, tree)
                END IF
             ELSE
                IF(ASSOCIATED(s1 % right) .AND. s1 % right % rb == RED) THEN
                   s2 => s1 % left
                   s2 % rb = BLACK
                   s1 % rb = node % parent % rb
                   node % parent % rb = BLACK
                   CALL rotate_node(s1, tree)
                ELSE
                   s2 => s1 % right
                   s2 % rb = node % parent % rb
                   node % parent % rb = BLACK
                   CALL rotate_node(s2, tree)
                   CALL rotate_node(s2, tree)
                END IF
             END IF
                
             
          END IF
             

       ! Case v: Sibling is red.    
       ! Rotate and colour parent red. Guarantees sibling and children
       ! now black.
       ELSE
          
          ! Change colour first.
          node % parent % rb = RED
          s1 % rb = BLACK
          ! Rotate on sibling node.
          CALL rotate_node(s1, tree)
          ! Recursively fix on modified tree.
          CALL fix_rb_del(node, tree)

       END IF

    END IF
    
  END SUBROUTINE fix_rb_del


  
  ! Subroutine del_node: 
  ! Deallocate an individual node. Takes no action on child nodes etc.
  SUBROUTINE del_node(node) 

    ! Node to deallocate.
    TYPE(pnode), POINTER :: node

    DEALLOCATE(node % pvec)
    DEALLOCATE(node)
    NULLIFY(node)

  END SUBROUTINE del_node


  
  ! Subroutine del_tree:
  ! Recursively deallocate a tree.
  RECURSIVE SUBROUTINE del_tree(tree)

    ! Root of tree to deallocate.
    TYPE(pnode), POINTER :: tree

    ! Deallocate left subtree.
    IF(ASSOCIATED(tree % left)) THEN
       CALL del_tree(tree % left)
    END IF

    ! Deallocate right subtree.
    IF(ASSOCIATED(tree % right)) THEN
       CALL del_tree(tree % right)
    END IF

    ! Deallocate me.
    CALL del_node(tree)

  END SUBROUTINE del_tree



END MODULE tree
