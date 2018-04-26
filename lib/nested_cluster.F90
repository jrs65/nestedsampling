! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

!! A module for performing Nested Sampling with recursive clustering.
!! @author Richard Shaw
!! @see nested
MODULE nested_cluster

! Include precision definitions
#include "defs.H90"

USE pset_mod
USE nested
USE twomeans
USE ellipsoid

IMPLICIT NONE

! Perform clustering (and recursive sampling)
INTEGER, PRIVATE :: cluster = 10

! Current number of branches.
INTEGER, PRIVATE :: n_branches = 0

! Maximum number of branchings.
INTEGER, PRIVATE :: m_branches = 10

! Output clusters to files
LOGICAL, PRIVATE :: cluster_output = .FALSE.

!! Filehandle for continuous output. If zero don't output.
INTEGER :: fh = 0

PRIVATE :: wnested_sample

CONTAINS


  !! Set continous IO filehandle.
  SUBROUTINE set_cont_output(n)
    INTEGER, INTENT(IN) :: n
    
    fh = n

  END SUBROUTINE set_cont_output


  !! Set clustering on or off. On, set non-zero n, where n is the
  !! clustering interval. Off, set = 0.
  SUBROUTINE set_cluster(n)

    INTEGER, INTENT(IN) :: n
    
    cluster = n

  END SUBROUTINE set_cluster


  !! Set output of clusters to files on/off
  SUBROUTINE set_cluster_output(tf)

    LOGICAL, INTENT(IN) :: tf
    
    cluster_output = tf

  END SUBROUTINE set_cluster_output

  
  !! Set the volume reduction factor.
  SUBROUTINE set_volred(v)

    REAL(rk), INTENT(IN) :: v
    
    vol_red = v

  END SUBROUTINE set_volred



  !! Set the volume reduction factor.
  SUBROUTINE set_clustsep(s)

    REAL(rk), INTENT(IN) :: s
    
    sep = s

  END SUBROUTINE set_clustsep



  !! Set the maximum number of branchings to perform.
  SUBROUTINE set_maxbranching(n)

    INTEGER, INTENT(IN) :: n
    
    m_branches = n

  END SUBROUTINE set_maxbranching


  ! Subroutine sample_points:
  !! Sample the Parameter space using a Nested Sampling Method.
  !! Returns an array of point nodes (pnode), sorted in ascending  likelihood.  
  SUBROUTINE sample_points(dpars, L, X, initial_points, replace_points, cs)

    !! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars
    
    INTERFACE
       !! The Likelihood function
       REAL(rk) FUNCTION L(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION L
       
       !! The Prior function
       REAL(rk) FUNCTION X(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION X

       !! Select initial points from the prior volume.
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

       !! Replace a point within the likelihood bound
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

    !! Point set to output into.
    TYPE(pset), POINTER :: cs

    ! Binary tree for initial population.
    TYPE(pnode), POINTER :: tree, node
    REAL(rk) :: t, axmval
    INTEGER :: i, p, n_rem

    NULLIFY(tree)
    NULLIFY(node)


    ALLOCATE(cs)

    CALL init_cluster(dpars)

    IF(.NOT. read_resume_initial(tree, npoints, dpars)) THEN
       CALL initial_points(tree, npoints, dpars, L, X)
       CALL write_resume_tree(tree)
    END IF

    ! Reset number of branches.
    n_branches = 0

    cs % root => tree
    cs % npoints = npoints
    NULLIFY(cs % c1, cs % c2)

    CALL wnested_sample(cs, npoints, REAL(1.0, rk), real(0.0, rk), dpars, L, X, replace_points)  
    
  END SUBROUTINE sample_points



  ! Subroutine wnested_sample:
  !! Private work routine for nested sampling. Using cluster method.
  RECURSIVE SUBROUTINE wnested_sample(cset, nc, xr, ec, dpars, L, X, replace_points)

    ! Current node tree.
    TYPE(pset), POINTER :: cset
    ! Current number of points in tree.
    INTEGER, INTENT(IN) :: nc
    ! accumulated X remaining.
    REAL(rk), INTENT(IN) :: xr
    ! current evidence (rough).
    REAL(rk), INTENT(IN) :: ec
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


    ! Binary tree for initial population.
    TYPE(pnode), POINTER :: tree, t1, t2, t_rem, node, node_resume

    ! Point sets for clustering.
    TYPE(pset), POINTER :: c1, c2

    ! Work integers.
    INTEGER :: i, p, j

    REAL(rk) :: tf, ax, f, e, ml

    CHARACTER(LEN=10) :: cfile

    ! Fractional decrease in X vol.
    tf = exp( -1.0 / cset % npoints)
    ax = xr
    e = ec

    IF(.NOT. ASSOCIATED(cset)) THEN
       fatal("Panic.")
    END IF

    ! Reset ellipsoid skip.
    CALL reset_ellskip()

    tree => cset % root

    NULLIFY(node)
    IF(ASSOCIATED(tree)) THEN
       CALL find_largest(tree, node)
       ml = node % lval
    END IF

    NULLIFY(t1, t2, node, t_rem)

    CALL find_smallest(tree, node)
    
    ! Progress
    ! p = ceiling((cset % npoints - nc) / 100.0)
    p = 5
!    PRINT *,"Topping up cluster."
!    PRINT *, "Make initial clustering attempt"
    ! Make initial clustering attempt.
    ! NOTE: ifort was being screwy and always evaluating two_cluster,
    ! have nested conditional to prevent this.
    IF(cluster /= 0 .AND. n_branches < m_branches) THEN
       IF(two_cluster(tree, t1, t2, cset % n1, cset % n2)) THEN
          
          NULLIFY(cset % root)
          PRINT *, "Call branch_cluster"
	  CALL branch_cluster(cset, t1, t2, ax, e)
          
          RETURN

       END IF
       
    END IF


    ! Top up the tree until it contains npoints points.
    DO i=1,(cset % npoints - nc)

       ! Check for a resume point.
       IF(resume .AND. read_resume_node(node_resume, dpars)) THEN
          CALL add_node(node_resume, tree)
       ELSE
          ! If not, add a replacement point, using the given method.
          CALL replace_points(tree, niter, (node % lval), 0.0D0, dpars, L, X)
       END IF

       !IF(i == 1 .AND. text_dialog("Ellipsoid print?")) THEN
       !   CALL ellipsoid_surface(100, "el.dat")
       !END IF

       ! Progress
       IF(MOD(i, p) == 0) THEN
 !         CALL progress_meter(i / p)
       END IF

    END DO

    ! Set up max likelihood
    IF(ASSOCIATED(node)) THEN
       CALL find_largest(tree, node)
       IF(node % lval > ml) THEN
          e = e * exp(ml - node % lval)
          ml = node % lval
       END IF
    ELSE
       CALL find_largest(tree, node)
       ml = node % lval
    END IF

!    PRINT *,"Starting replacements..."

    ! Progress
    !p = ceiling(nrepl / 1000.0)
    p = 5

    ! Start the replacement process
    DO i=1,nrepl
       
       ! Remove the node with the smallest likelihood.
       CALL remove_smallest(tree, node)

       ! Set the fractional volume change, and crudely estimate evidence.
       e = e + ax * (1 - tf) * exp(node % lval - ml)
       
       node % w = ax*(1 - tf)

       ax = ax * tf
       node % axmval = ax

       ! Output to continuous stream.
       IF(fh /= 0) THEN
          CALL print_node_fh(node, fh)
       END IF

       ! Add node into a further tree. 
       CALL add_node(node, t_rem)
       
       IF(resume .AND. read_resume_node(node_resume, dpars)) THEN
          CALL add_node(node_resume, tree)
       ELSE
          ! Add a replacement point, using the given method.
          CALL replace_points(tree, niter, (node % lval), 0.0D0, dpars, L, X)
       END IF
 
       ! Check stopping.
       IF(MOD(i, p) == 0) THEN
          CALL find_largest(tree, node)
          ! Check for new maximum likelihood.
          IF(node % lval > ml) THEN
             e = e * exp(ml - node % lval)
             ml = node % lval
          END IF
          
 !         CALL progress_meter1(ax / e, volfrac)

          ! Terminate run, add active points into pset.
          IF((ax / e) < volfrac) THEN
             PRINT *,"Found required volume fraction. Points:",i

             cset % rem => tree

             ! Assign estimated xvals
             NULLIFY(node)
             j = cset % npoints
             DO WHILE(next_node(node, tree))

                j = j - 1       
                node % axmval = j* ax / npoints
                node % w = ax / npoints
                
             END DO
             EXIT
          END IF
       END IF

       ! Try to cluster into two.
       ! NOTE: ifort was being screwy and always evaluating two_cluster,
       ! have nested conditional to prevent this.
       IF(cluster /= 0 .AND. n_branches < m_branches) THEN 
          IF(MOD(i,cluster) == 0) THEN 
             IF(two_cluster(tree, t1, t2, cset % n1, cset % n2)) THEN
       
                CALL branch_cluster(cset, t1, t2, ax, e)
                
                EXIT

             END IF

          END IF
          
       END IF

    END DO

    cset % root => t_rem


  CONTAINS 
    
    SUBROUTINE branch_cluster(cst, tree1, tree2, axc, ev)
      
      TYPE(pset), POINTER :: cst
      TYPE(pnode), POINTER :: tree1
      TYPE(pnode), POINTER :: tree2
      
      REAL(rk), INTENT(IN) :: axc
      REAL(rk), INTENT(IN) :: ev
      
      TYPE(pnode), POINTER :: node
      REAL(rk) :: sc1, sc2

      IF(cluster_output) THEN
         WRITE (cfile, '(A, I2.2, I1.1)') 'cluster', n_branches, 1
         CALL print_ptree(tree1, cfile)
         WRITE (cfile, '(A, I2.2, I1.1)') 'cluster', n_branches, 2
         CALL print_ptree(tree2, cfile)
      END IF
     
      ! Insert blank line to account for progress meter interruption.
      PRINT *,""
      PRINT *,"*** Clustered ***"
      
      ! Increment cluster counter.
      n_branches = n_branches + 1

      ! Create a cluster point set for cluster 1.
      ALLOCATE(cst % c1)
      NULLIFY(cst % c1 % c1, cst % c1 % c2)
      cst % c1 % root => tree1
      cst % c1 % npoints = cset % npoints
      
      ! Create a cluster point set for cluster 2.
      ALLOCATE(cst % c2)
      NULLIFY(cset % c2 % c1, cst % c2 % c2)
      cst % c2 % root => tree2
      cst % c2 % npoints = cst % npoints
      
      f = (1.0 * cst % n1) / (cst % n1 + cst % n2)
      
      PRINT '(a,f5.2,a)',"  * Cluster 1: ", 100*f," %"
      ! Output to continuous stream.
      IF(fh /= 0) THEN
         WRITE(fh,*)  "### Cluster 1 num:", cst % n1
      END IF

      ! Need to update ev scalings.
      CALL find_largest(tree1, node)
      sc1 =  exp(ml - node % lval)
      
      node => NULL()
      CALL find_largest(tree2, node)
      sc2 =  exp(ml - node % lval)


      ! Recursively sample cluster 1.
      CALL wnested_sample(cst % c1, cst % n1, axc * f, ev*sc1, &
           dpars, L, X, replace_points)
      
      PRINT '(a,f5.2,a)',"  * Cluster 2: ", 100*(1-f)," %"
      ! Output to continuous stream.
      IF(fh /= 0) THEN
         WRITE(fh,*)  "### Cluster 2 num:", cst % n2
      END IF

      ! Recursively sample cluster 2.
      CALL wnested_sample(cst % c2, cst % n2, axc * (1-f), ev*sc2, &
           dpars, L, X, replace_points)
      
    END SUBROUTINE branch_cluster
    
  END SUBROUTINE wnested_sample
  
  

END MODULE nested_cluster

      
