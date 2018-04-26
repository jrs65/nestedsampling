!! A module containing the pset type and related routines. It is
!! designed to hold the results of a sampling run.  
!! @author Richard Shaw
MODULE pset_mod

! Include precision definitions
#include "defs.H90"

USE tree
USE treedef
USE treeutil
USE nested
USE rand_nested

IMPLICIT NONE

! TYPE pset
!! A type for a cluster point set. Contains the point tree, as well as
!! metadata on how it was drawn (to allow evidence calculation).
TYPE, PUBLIC :: pset

   !! The root node of the point tree.
   TYPE(pnode), POINTER :: root
   !! Number of active points
   INTEGER :: npoints

   !! Remaining active points after termination.
   TYPE(pnode), POINTER :: rem

   !! Child cluster point sets
   TYPE(pset), POINTER :: c1, c2
   !! Child cluster (original) point numbers 
   INTEGER :: n1, n2

END TYPE pset


CONTAINS

  ! Subroutine evidence_pset:
  !! Take a set of an array of sampled points, and calculate
  !! the log-value of the evidence. Use the simple 1/N method.
  RECURSIVE SUBROUTINE evidence_pset(cs, x_arr, ev_arr)

    !! Array of point nodes to use.
    TYPE(pset), POINTER :: cs

    !! Array to output evidence x values in.
    REAL(rk), DIMENSION(:), INTENT(IN) :: x_arr

    !! Array to output evidence values in.
    REAL(rk), DIMENSION(:), INTENT(OUT) :: ev_arr

    ! Counter, size label
    INTEGER :: i,n,np
    ! Width working arrays, Evidence work array.
    REAL(rk), DIMENSION(:), ALLOCATABLE ::axm1, axm2, t, ev, ev1, ev2
    REAL(rk) :: average, scl,l1,f

    ! Temp node pointer.
    TYPE(pnode), POINTER :: node, tree

    tree => cs % root
    np = cs % npoints

    n = size(ev_arr)
    ALLOCATE(axm1(n), axm2(n), t(n), ev(n), ev1(n), ev2(n))

    ev = 0.0
    scl = 0.0

    axm1 = x_arr
    axm2 = x_arr

    !PRINT *,"X",sum(x_arr)/n

    IF(ASSOCIATED(tree)) THEN

       CALL find_largest(tree, node)
       
       ! Scale to reduce numerical inaccuracy on exponentiation,
       scl = -1.0 * (node % lval)
       
       NULLIFY(node)
       IF(.NOT. next_node(node, tree)) THEN
          fatal("Tree seems to be empty.")
       END IF
       
       ! Special case for first strip. Pin to (1,0).
       CALL sample_volreduce(t, np)
       axm1 = x_arr
       axm2 = axm1 * t
       ev = (axm1 - axm2) * 0.5 * EXP((node % lval) + scl)

       ! Save first likelihood
       l1 = node % lval
       
       ! Evaluate using the trapezium method. Points in replacement
       ! phase.
       DO WHILE(next_node(node, tree))
          
          ! Average scaled likelihood value.
          average = 0.5 * (EXP(l1 + scl) &
               + EXP((node % lval)+scl))
          ! Generate new set of t-factors.
          CALL sample_volreduce(t, np)
          ! X vals.
          axm1 = axm2
          axm2 = axm1 * t
          l1 = node % lval
          
          ev = ev + average * (axm1 - axm2)

       END DO

    END IF

    ! Check to see if we have sub clusters.
    ! If so, add in their contributions.
    IF(ASSOCIATED(cs % c1)) THEN
       ! TODO: Replace splitting with more sound method. DONE.
       CALL sample_volreduce_n(t, (cs % n1 + cs % n2), cs % n1)

       !PRINT *,"VR",sum(t)/n

       ! Recurse into sub clusters.
       CALL evidence_pset(cs % c1, axm2 * t, ev1)
       CALL evidence_pset(cs % c2, axm2 * (1-t), ev2)

       ev = EXP(ev1 + scl) + EXP(ev2 + scl) + ev

       !PRINT *,"Cl",sum(ev1)/n,sum(ev2)/n,sum(ev)/n

    ! In terminating cluster, need to include remaining active points.
    !ELSE IF(ASSOCIATED(cs % rem)) THEN

    ! Otherwise we are in a terminating cluster, need to pin to end.
    ELSE
       ! Special case for last strip. Pin to (1,0)
       average = EXP(l1 + scl)

       ev = ev + axm2 * average
    
    END IF

    ! Return the log evidence.
    ev_arr = log(ev) - scl

    !PRINT *,sum(ev_arr) / n
       
    DEALLOCATE(axm1, axm2, t, ev)
    
  END SUBROUTINE evidence_pset


  ! Function evidence:
  !! Give an estimate of the evidence (and uncertainty).  Evaluate
  !! using several sets {t}, and return means and standard deviation.
  FUNCTION evidence(cs)

    REAL(rk), DIMENSION(2) :: evidence

    ! Array of point nodes to use.
    TYPE(pset), POINTER :: cs

    REAL(rk), DIMENSION(nevid) :: ev, x

    REAL(rk) :: x1,x2
    INTEGER :: i

    x1 = 0
    x2 = 0

    x = 1.0
    CALL evidence_pset(cs, x, ev)

    DO i=1,nevid
       x1 = x1 + ev(i)
       x2 = x2 + ev(i)**2
    END DO

    evidence(1) = x1 / nevid
    evidence(2) = sqrt((x2 / nevid) - (x1 / nevid)**2)

  END FUNCTION evidence
    
  ! Subroutine average_pset:
  !! Walk through the entire pset, and apply a given routine to each
  !! node, passing its current weight along. Useful for performing
  !! averages.
  RECURSIVE SUBROUTINE average_pset(cs, x_arr, fx)

    !! Array of point nodes to use.
    TYPE(pset), POINTER :: cs

    !! Array to output evidence x values in.
    REAL(rk), DIMENSION(:), INTENT(IN) :: x_arr

    INTERFACE
       !! Function to apply to each node.
       SUBROUTINE fx(node, weight)

         USE treedef

         ! Node to operate on.
         TYPE(pnode), POINTER :: node, tree
         
         ! Weight of current point.
         REAL(rk), DIMENSION(:), INTENT(IN) :: weight

       END SUBROUTINE fx
    END INTERFACE
         
    ! Counter, size label
    INTEGER :: n, np
    ! Width working arrays, Evidence work array.
    REAL(rk), DIMENSION(:), ALLOCATABLE ::axm1, axm2, t, w
    REAL(rk) :: f

    ! Temp node pointer.
    TYPE(pnode), POINTER :: node, tree

    tree => cs % root
    np = cs % npoints

    n = size(x_arr)
    ALLOCATE(axm1(n), axm2(n), t(n), w(n))

    NULLIFY(node)

    axm1 = x_arr
    axm2 = x_arr

    ! Iterate over points in set, applying f to each.
    DO WHILE(next_node(node, tree))

       ! Generate new set of t-factors.
       CALL sample_volreduce(t, np)
       ! X vals.
       axm1 = axm2
       axm2 = axm1 * t
       w = (axm1 - axm2) * exp(node % lval)

       CALL fx(node, w)

    END DO

    ! Check to see if we have sub clusters.
    IF(ASSOCIATED(cs % c1)) THEN
       ! TODO: Replace splitting with more sound method. DONE.
       CALL sample_volreduce_n(t, (cs % n1 + cs % n2), cs % n1)
       ! Recurse into sub clusters.
       CALL average_pset(cs % c1, axm2 * t, fx)
       CALL average_pset(cs % c2, axm2 * (1-t), fx)

    END IF

    DEALLOCATE(axm1, axm2, t)
    
  END SUBROUTINE average_pset



  ! Subroutine print_pset:
  !! Print a pset into file, file.
  SUBROUTINE print_pset(cs, file)

    !! Array of point nodes to use.
    TYPE(pset), POINTER :: cs

   ! Filename to write into.
    CHARACTER(LEN=*), INTENT(IN) :: file

    OPEN(10, FILE=file)
    CALL print_pset_work(cs)
    CALL flush(10)
    CLOSE(10)

  CONTAINS

    ! Nested work subroutine.
    RECURSIVE SUBROUTINE print_pset_work(cset)
      !! Array of point nodes to use.
      TYPE(pset), POINTER :: cset
      
      ! Temp node pointer.
      TYPE(pnode), POINTER :: node, tree
      
      tree => cset % root

      IF(ASSOCIATED(tree)) THEN
         NULLIFY(node)
         
         WRITE(10,*) "### Cluster ###"
         
         ! Loop through points and print them out.
         DO WHILE(next_node(node, tree))
            CALL print_node(node)
         END DO
         
      END IF

      ! Check to see if we have sub clusters.
      ! If so, print them.
      IF(ASSOCIATED(cset % c1)) THEN

         WRITE(10,*) ""
         WRITE(10,*) ""
         
         CALL print_pset_work(cset % c1)
         CALL print_pset_work(cset % c2)
         
      ! Else print out remaining active points.
      ELSE IF(ASSOCIATED(cset % rem)) THEN

         WRITE(10,*) "### Remaining active points ###" 

         NULLIFY(node)
         DO WHILE(next_node(node, cset % rem)) 
            CALL print_node(node)
         END DO

         WRITE(10,*) ""
         WRITE(10,*) ""


      END IF
      
    END SUBROUTINE print_pset_work


  END SUBROUTINE print_pset



  RECURSIVE SUBROUTINE evidence_simple(cs, evo)

    !! Array of point nodes to use.
    TYPE(pset), POINTER :: cs

    !! Array to output evidence values in.
    REAL(rk), INTENT(OUT) :: evo

    REAL(rk) :: average, scl, ev1, ev2, ev

    INTEGER :: np

    ! Temp node pointer.
    TYPE(pnode), POINTER :: node, tree

    tree => cs % root
    np = cs % npoints

    ev = 0.0
    scl = 0.0

    IF(ASSOCIATED(tree)) THEN

       CALL find_largest(tree, node)
       
       ! Scale to reduce numerical inaccuracy on exponentiation,
       scl = -1.0 * (node % lval)
       
       NULLIFY(node)
       IF(.NOT. next_node(node, tree)) THEN
          fatal("Tree seems to be empty.")
       END IF
       

       DO WHILE(next_node(node, tree))
          
          ! Average scaled likelihood value.
          average = EXP((node % lval)+scl) * (node % w)
          
          ev = ev + average

       END DO

    END IF

    ! Check to see if we have sub clusters.
    ! If so, add in their contributions.
    IF(ASSOCIATED(cs % c1)) THEN

       ! Recurse into sub clusters.
       CALL evidence_simple(cs % c1, ev1)
       CALL evidence_simple(cs % c2, ev2)

       ev = EXP(ev1 + scl) + EXP(ev2 + scl) + ev

    END IF

    ! Return the log evidence.
    evo = log(ev) - scl

    !PRINT *,evo
       
  END SUBROUTINE evidence_simple


END MODULE pset_mod
