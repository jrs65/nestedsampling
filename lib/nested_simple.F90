! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

MODULE nested_simple

! Include precision definitions
#include "defs.H90"

USE nested

IMPLICIT NONE


CONTAINS



  ! Take a set of an array of sampled points, and calculate
  ! the log-value of the evidence. Use the simple 1/N method.
  REAL(rk) FUNCTION evidence_simple(array)

    ! Array of point nodes to use.
    TYPE(pnode), DIMENSION(:), POINTER :: array
    INTEGER :: i

    REAL(rk) :: ev, average, width, scl

    ! Scale to reduce numerical inaccuracy on exponentiation,
    scl = -1.0 * (array(size(array)) % lval)

    ! Special case for first strip. Pin to (1,0).
    width = 1 - exp(array(1) % axmval)
    ev = width * 0.5 * EXP((array(1) % lval) + scl)

    ! Evaluate using the trapezium method.
    DO i=1,(SIZE(array) - 1)

       ! Average scaled likelihood value.
       average = 0.5 * (EXP((array(i) % lval) + scl) + EXP((array(i+1) % lval)+scl))
       ! Width of trapezium.
       width = exp(array(i) % axmval) - exp(array(i+1) % axmval)

       ev = ev + average * width

    END DO
    
    ! Special case for last strip. Pin to (1,0)
    average = EXP((array(SIZE(array)) % lval) + scl)

    ev = ev + exp(array(SIZE(array)) % axmval) * average
    
    ! Return the log evidence.
    evidence_simple = log(ev) - scl
    
  END FUNCTION evidence_simple

  

    ! Take a set of an array of sampled points, and calculate
  ! the log-value of the evidence. Use the simple 1/N method.
  SUBROUTINE evidence_array(array, ev_arr)

    ! Array of point nodes to use.
    TYPE(pnode), DIMENSION(:), POINTER :: array

    ! Array to output evidence values in.
    REAL(rk), DIMENSION(:), INTENT(OUT) :: ev_arr

    ! Counter, size label
    INTEGER :: i,n
    ! Width working arrays, Evidence work array.
    REAL(rk), DIMENSION(:), ALLOCATABLE ::axm1, axm2, t, ev
    REAL(rk) :: average, scl

    ! Check number of points is correct.
    IF((nrepl + npoints) /= size(array)) THEN
       PRINT *,"Error: Incorrect number of points"
    END IF

    n = size(ev_arr)
    ALLOCATE(axm1(n), axm2(n), t(n), ev(n))

    ! Scale to reduce numerical inaccuracy on exponentiation,
    scl = -1.0 * (array(size(array)) % lval)

    ! Special case for first strip. Pin to (1,0).
    CALL sample_volreduce(t, npoints)
    axm1 = 1
    axm2 = t
    ev = (axm1 - axm2) * 0.5 * EXP((array(1) % lval) + scl)

    ! Evaluate using the trapezium method. Points in replacement
    ! phase.
    DO i=1,nrepl

       ! Average scaled likelihood value.
       average = 0.5 * (EXP((array(i) % lval) + scl) &
            + EXP((array(i+1) % lval)+scl))
       ! Generate new set of t-factors.
       CALL sample_volreduce(t, npoints)
       ! X vals.
       axm1 = axm2
       axm2 = axm1 * t

       ev = ev + average * (axm1 - axm2)

    END DO

    DO i=1,(npoints-1)

       ! Average scaled likelihood value.
       average = 0.5 * (EXP((array(nrepl+i) % lval) + scl) &
            + EXP((array(nrepl+i+1) % lval)+scl))
       ! Generate new set of t-factors.
       CALL sample_volreduce(t, (npoints-i+1))
       ! X vals.
       axm1 = axm2
       axm2 = axm1 * t

       ev = ev + average * (axm1 - axm2)

    END DO
    
    ! Special case for last strip. Pin to (1,0)
    average = EXP((array(SIZE(array)) % lval) + scl)

    ev = ev + axm2 * average
    
    ! Return the log evidence.
    ev_arr = log(ev) - scl

    PRINT *,ev_arr

    DEALLOCATE(axm1, axm2, t, ev)
    
  END SUBROUTINE evidence_array



  ! Give an estimate of the evidence (and uncertainty).  Evaluate
  ! using several sets {t}, and return means and standard deviation.
  FUNCTION evidence(array)

    REAL(rk), DIMENSION(2) :: evidence

    ! Array of point nodes to use.
    TYPE(pnode), DIMENSION(:), POINTER :: array

    REAL(rk), DIMENSION(nevid) :: ev

    REAL(rk) :: x1,x2
    INTEGER :: i

    x1 = 0
    x2 = 0

    CALL evidence_array(array, ev)

    DO i=1,nevid
       x1 = x1 + ev(i)
       x2 = x2 + ev(i)**2
    END DO

    evidence(1) = x1 / nevid
    evidence(2) = sqrt((x2 / nevid) - (x1 / nevid)**2)

  END FUNCTION evidence
    

  

  ! Print the array of points.
  SUBROUTINE print_parray(array, fd)
    
    ! Array of point nodes to use.
    TYPE(pnode), DIMENSION(:), POINTER :: array
    ! File handle to write to.
    INTEGER, INTENT(IN) :: fd
    INTEGER :: i

    DO i=1,(SIZE(array))
       WRITE(fd,*) (array(i) % pvec), (array(i) % xval), (array(i) % axmval), (array(i) % lval)
    END DO

  END SUBROUTINE print_parray
  



  ! Sample the Parameter space using a Nested Sampling Method.
  ! Returns an array of point nodes (pnode), sorted in ascending likelihood.
  SUBROUTINE sample_points(dpars, L, X, initial_points, replace_points, array)

    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars
    
    INTERFACE
       ! The Likelihood function
       REAL(rk) FUNCTION L(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION L
       
       ! The Prior function
       REAL(rk) FUNCTION X(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION X

       ! Select initial points from the prior volume.
       SUBROUTINE initial_points(tree, n, dpars, L, X)

         USE treedef
         
         ! The point tree to add to.
         TYPE(pnode), POINTER :: tree
         ! Number of points to add.
         INTEGER, INTENT(IN) :: n
         ! Dimensionality of parameter space.
         INTEGER, INTENT(IN) :: dpars
         
         INTERFACE
            ! The Likelihood function
            REAL(rk) FUNCTION L(ppars)
              REAL(rk), DIMENSION(:) :: ppars
            END FUNCTION L
            
            ! The Prior function
            REAL(rk) FUNCTION X(ppars)
              REAL(rk), DIMENSION(:) :: ppars
            END FUNCTION X

         END INTERFACE
         
       END SUBROUTINE initial_points

       ! Replace a point within the likelihood bound
       SUBROUTINE replace_points(tree, n_iter, lbound, axmval, dpars, L, X)

         USE treedef
         
         ! The point tree to add to.
         TYPE(pnode), POINTER :: tree
         ! Number of times to evolve.
         INTEGER, INTENT(IN) :: n_iter
         ! Likelihood lower bound.
         REAL(rk), INTENT(IN) :: lbound
         ! Accumulated Prior mass. Use to rescale the scale parameters.
         REAL(rk), INTENT(IN) :: axmval
         ! Dimensionality of parameter space.
         INTEGER, INTENT(IN) :: dpars
         
         INTERFACE
            ! The Likelihood function
            REAL(rk) FUNCTION L(ppars)
              REAL(rk), DIMENSION(:) :: ppars
            END FUNCTION L
            
            ! The Prior function
            REAL(rk) FUNCTION X(ppars)
              REAL(rk), DIMENSION(:) :: ppars
            END FUNCTION X
            
         END INTERFACE

       END SUBROUTINE replace_points

    END INTERFACE

    ! Array to output into.
    TYPE(pnode), DIMENSION(:), POINTER :: array

    ! Binary tree for initial population.
    TYPE(pnode), POINTER :: tree, node
    REAL(rk) :: t, axmval
    INTEGER :: i, p

    NULLIFY(tree)
    NULLIFY(node)


    ALLOCATE(array(npoints + nrepl))

    CALL initial_points(tree, npoints, dpars, L, X)
    
    ! Debug, print initial tree.
    !CALL print_ptree(tree, "initial")
    
    axmval = 0.0

    t = -1.0 / npoints

    p = nrepl / 100

    DO i=1,nrepl
       
       ! Remove the node with the smallest likelihood.
       CALL remove_smallest(tree, node)
       
       IF(MOD(i, p) == 0) THEN
          CALL progress_meter(i / p)
       END IF

       ! Set the Accumulated Prior mass value according to NS.
       ! Just use a crude estimate for setting in the pnode.
       axmval = axmval + t
       node % axmval = axmval
       ! Add the node to the array.
       array(i) = node

       ! Add a replacement point, using the given method.
       CALL replace_points(tree, niter, (node % lval), exp(axmval), dpars, L, X)

    END DO

    DO i=1,npoints
       IF(.NOT. ASSOCIATED(tree)) THEN
          PRINT *,"Problem"
          CALL exit(-1)
       END IF
       ! Set the AXM value according to 1/N prescription.
       axmval = axmval  - (1.0 / (npoints - i + 1))
       CALL remove_smallest(tree, node)
       node % axmval = axmval 
       array(i+nrepl) = node

       !PRINT *,"Removed"
       !CALL print_node(node)


    END DO

    !PRINT *,"Final"
    !CALL print_parray(array, 6)
    
  END SUBROUTINE sample_points

  ! Work out the mean of the parameter space.
  SUBROUTINE mean_pars(arr, dpars, par)
     
    TYPE(pnode), DIMENSION(:), POINTER :: arr
    INTEGER, INTENT(IN) :: dpars
     REAL(rk), INTENT(OUT), DIMENSION(:) :: par
    
    REAL(rk), DIMENSION(2) :: ev
    INTEGER :: i
    REAL(rk) :: width, scl

    ! Scale to reduce numerical inaccuracy.
    scl = -1.0 * (arr(size(arr)) % lval)
    PRINT *,"Scale",scl
    par = 0
    
    ! Sum over points, weighted by widths + likelihood
    DO i=2,size(arr)
       width = (arr(i-1) % axmval) - (arr(i) % axmval)
       
       par = par + ((arr(i) % pvec) * EXP((arr(i) % lval) + scl) * width)
       
    END DO
    
    ev = evidence(arr)

    par = par * (exp(-scl-ev(1)))
    
  END SUBROUTINE mean_pars

    
END MODULE nested_simple

      
