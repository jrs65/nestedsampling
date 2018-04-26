! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

!! Perform k-means (with k=2) to partition active points into two
!! clusters. May modify with scheme for more general k.
MODULE twomeans

#include "defs.H90"

USE rand
USE tree
USE treeutil

IMPLICIT NONE

! A doubly linked list for storing intermediate cluster states.
TYPE, PRIVATE :: plist
   
   ! Previous, next pointers for doubly linked list.
   TYPE(plist), POINTER :: prev, next
   
   ! Pointer to the node.
   TYPE(pnode), POINTER :: node
   
END TYPE plist

! List storing cluster 1.
TYPE(plist), POINTER, PRIVATE :: c1
! List storing cluster 2.
TYPE(plist), POINTER, PRIVATE :: c2

! Array storing cluster 1 centre.
REAL(rk), DIMENSION(:), PRIVATE, POINTER :: m1
! Array storing cluster 2 centre.
REAL(rk), DIMENSION(:), PRIVATE, POINTER :: m2

! Array storing cluster 1 work space.
REAL(rk), DIMENSION(:), PRIVATE, POINTER :: tm1
! Array storing cluster 2 work space.
REAL(rk), DIMENSION(:), PRIVATE, POINTER :: tm2

! Number of points in cluster 1.
INTEGER, PRIVATE :: n1 = 0
! Number of points in cluster 2.
INTEGER, PRIVATE :: n2 = 0


! Clustering iterations.
INTEGER, PRIVATE :: n_iter = 5

REAL(rk), PRIVATE :: rm1 = 0.0
REAL(rk) ::  vol_red = 0.2
REAL(rk) ::  sep = 0.4

REAL(rk), PRIVATE :: min_dist

INTEGER, PRIVATE :: nsurf_points = 200

! Array of surface points.
REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: surf_points
! Metric matrix.
REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: metric

! Set PRIVATE status on subroutines
PRIVATE :: cluster_work, add_plist, init_plist, cluster_volume

CONTAINS


  ! Set the number of clustering iterations to make.
  SUBROUTINE set_cluster_iterations(n)

    INTEGER, INTENT(IN) :: n

    n_iter = n

  END SUBROUTINE set_cluster_iterations



  ! Initialise the cluster work arrays etc.
  SUBROUTINE init_cluster(n) 
    
    INTEGER, INTENT(IN) :: n

    ALLOCATE(m1(n), m2(n), tm1(n), tm2(n))
    ALLOCATE(surf_points(nsurf_points,n))
    ALLOCATE(metric(n,n))
    NULLIFY(c1, c2)

  END SUBROUTINE init_cluster


  LOGICAL FUNCTION two_cluster(tree, tree1, tree2, num1, num2)

    ! The original tree to cluster on.
    TYPE(pnode), POINTER :: tree
    
    ! The tree for cluster 1.
    TYPE(pnode), POINTER :: tree1
    ! The tree for cluster 2.
    TYPE(pnode), POINTER :: tree2

    ! Number in tree 1.
    INTEGER, INTENT(OUT) :: num1
    ! Number in tree 2.
    INTEGER, INTENT(OUT) :: num2
  
    ! Volumes, and separations.
    REAL(rk) :: vol, vol1, vol2, min1, min2, dist

    TYPE(plist), POINTER :: t
    INTEGER :: i,j,k,n
    LOGICAL :: intersect

    intersect = .FALSE.

    NULLIFY(c1, c2)

    n = size(m1)

    n1 = 0
    CALL walk_inorder(tree, init_plist)
    
    vol = cluster_volume(c1, .FALSE.)

    m1 = c1 % node % pvec
    m2 = c2 % node % pvec
!    PRINT *,m1,m2
    NULLIFY(c2)

    DO i=1,n_iter
       
       CALL cluster_work()

    END DO

    vol1 = cluster_volume(c1, .TRUE.)
    min1 = min_dist
    vol2 = cluster_volume(c2, .FALSE.)
    min2 = min_dist

    COSMOS_OUT("*** Two means ***")
#ifdef _COSMOS_DBG
    PRINT "(A,G10.5)","Volumes:     ",(vol1 + vol2) / vol
    PRINT "(A,G10.5)","Separations: ",(min1+min2) / sqrt(distance(m1,m2))
#endif

    ! If we meet the seperation criteria, cluster.
    ! If not just leave the tree intact.
    IF((vol1 + vol2) < (vol_red * vol) .AND. n1 > n .AND. n2 > n &
         .AND. (min1+min2) < (sep*sqrt(distance(m1,m2)))) THEN

       ! Check for intersection.
       DO k=1,nsurf_points
          
          DO i=1,n
             DO j=1,n
                dist = dist + ((surf_points(k,i) - m2(i)) * metric(i,j) * (surf_points(k,j) - m2(j)))
                
             END DO
          END DO

          ! Set intersection to true
          IF(dist < 1.0) THEN
             intersect = .TRUE.
             
             COSMOS_OUT("*** Clusters intersect. ***")
             two_cluster = .FALSE.   
             COSMOS_OUT("Clustering failed.")

             EXIT
             
          END IF

       END DO

       IF(.NOT. intersect) THEN
          
          t => c1
          DO WHILE(ASSOCIATED(t)) 
             CALL add_node(t % node, tree1)
             t => t % next
          END DO
          
          t => c2
          DO WHILE(ASSOCIATED(t)) 
             CALL add_node(t % node, tree2)
             t => t % next
          END DO
          
          NULLIFY(tree)
          
          num1 = n1
          num2 = n2
          
          two_cluster = .TRUE.
          COSMOS_OUT("Clustering successful.")
       END IF
    ELSE
       
       two_cluster = .FALSE.   
       COSMOS_OUT("Clustering failed.")
    END IF
    
  END FUNCTION two_cluster



  SUBROUTINE init_plist(node)
    
    ! The node we're currently at.
    TYPE(pnode), POINTER :: node

    TYPE(plist), POINTER :: li
    
    ALLOCATE(li)

    li % node => node
    n1 = n1 + 1

    CALL add_plist(li, c2, c1)

  END SUBROUTINE init_plist


  ! Perform one 2-means iteration.
  SUBROUTINE cluster_work()

    TYPE(plist), POINTER :: p1, p2, pt1, pt2
    TYPE(plist), POINTER :: t,t2
    REAL(rk) :: d1, d2

    NULLIFY(p1, p2)
    ! Set temporary variables to zero.
    tm1 = 0; tm2 = 0; n1 = 0; n2 = 0;
    
    rm1 = 0

    ! Cycle through cluster 1, and reassign points.
    t => c1
    DO WHILE(ASSOCIATED(t))

       ! Calculate distances to centres.
       d1 = distance(m1, t % node % pvec)
       d2 = distance(m2, t % node % pvec)

       t2 => t % next

       ! Assign to each cluster.
       IF(d1 < d2) THEN
          CALL add_plist(t, pt1, p1)
          tm1 = tm1 + t % node % pvec
          n1 = n1 + 1
          IF(d1 > rm1) THEN 
             rm1 = d1
          END IF
       ELSE
          CALL add_plist(t, pt2, p2)
          tm2 = tm2 + t % node % pvec
          n2 = n2 + 1
       END IF

       t => t2

    END DO

    ! Cycle through cluster 1, and reassign points.
    t => c2
    DO WHILE(ASSOCIATED(t))

       ! Calculate distances to cluster centres.
       d1 = distance(m1, t % node % pvec)
       d2 = distance(m2, t % node % pvec)

       t2 => t % next

       ! Assign point to each cluster.
       IF(d1 < d2) THEN
          CALL add_plist(t, pt1, p1)          
          tm1 = tm1 + t % node % pvec
          n1 = n1 + 1
          ! Set the max radius of Cluster 1.
          IF(d1 > rm1) THEN 
             rm1 = d1
          END IF
       ELSE
          CALL add_plist(t, pt2, p2)
          tm2 = tm2 + t % node % pvec
          n2 = n2 + 1
       END IF

       t => t2
    END DO

    ! Re-calculate means based on new clusters.
    m1 = tm1 / n1
    m2 = tm2 / n2

    ! Re assign cluster pointers.
    c1 => p1
    c2 => p2

  END SUBROUTINE cluster_work



  SUBROUTINE add_plist(node, cp, list)

    TYPE(plist), POINTER :: node
    TYPE(plist), POINTER :: cp
    TYPE(plist), POINTER :: list

    IF(.NOT. ASSOCIATED(list)) THEN
       list => node
       cp => node
       NULLIFY(node % prev, node % next)
    ELSE IF(.NOT. ASSOCIATED(cp % next)) THEN
       cp % next => node
       node % prev => cp
       cp => node
       NULLIFY(node % next)
    ELSE
       fatal("Next pointer must be null.")
    END IF

  END SUBROUTINE add_plist

  
  ! Convenience function to calculate distance (squared) between two
  ! points in the parameter space.
  REAL(rk) FUNCTION distance(par1, par2)
    
    REAL(rk), DIMENSION(:), POINTER :: par1
    REAL(rk), DIMENSION(:), POINTER :: par2
    
    INTEGER :: i
    
    distance = 0.0
  
    ! Ensure that arrays are of same length.
    IF(size(par1) /= size(par2)) THEN
       fatal("Parameter array lengths must match.")
    END IF
    
    DO i = 1,size(par1)
       distance = distance + (par1(i) - par2(i))**2
    END DO
    
  END FUNCTION distance



  REAL(rk) FUNCTION cluster_volume(cl, flag) 

    TYPE(plist), POINTER :: cl

    ! Mode flag. If True, set array of points
    ! If false, set metric matrix.
    LOGICAL, INTENT(IN) :: flag

    TYPE(plist), POINTER :: cur

    REAL(rk), DIMENSION(:,:), ALLOCATABLE :: mat, met
    REAL(rk), DIMENSION(:), ALLOCATABLE :: means,work,evals,point,tv1
    REAL(rk) :: dist, max_dist, ts1, ts2

    INTEGER :: i,j,k,n,n_points,info

    n = size(cl % node % pvec)
    cur => cl
    ALLOCATE(mat(n,n), met(n,n), means(n), evals(n), work(n*70), point(n), tv1(n))
    n_points = 0
    means = 0
    mat = 0
    dist = 0
    max_dist = 0

    DO WHILE(ASSOCIATED(cur))
       
       n_points = n_points + 1
       means = means + (cur % node % pvec)
       
       DO i = 1,n
          DO j=1,n
             mat(i,j) = mat(i,j) + (cur % node % pvec(i))*(cur % node % pvec(j))
          END DO
       END DO

       cur => cur % next
    END DO

    ! Convert sums into moments.
    mat = mat / n_points
    means = means / n_points

    ! Generate covariances.
    DO i = 1,n
       DO j=1,n
          mat(i,j) = mat(i,j) - means(i) * means(j)
       END DO
    END DO

    CALL dsyev('V', 'U', n, mat, n, evals, work, n*70, info)

    DO i=1,n
       ! Set to zero for sum.
       met(i,:) = 0.0

       DO j=1,n
          DO k=1,n
             met(i,j) = met(i,j) + (mat(i,k) * mat(j,k) / evals(k))
          END DO
       END DO
    END DO

    cur => cl
    
    DO WHILE(ASSOCIATED(cur))

       dist = 0.0

       DO i=1,n
          DO j=1,n
             dist = dist + ((cur % node % pvec(i) - means(i)) * met(i,j) * (cur % node % pvec(j) - means(j)))
          END DO
       END DO
       
       ! If is greater than current max, update distance.
       IF(dist > max_dist) THEN
          max_dist = dist
       END IF
       
       cur => cur % next
    END DO

    !PRINT *,"Ev",evals,max_dist

    cluster_volume = max_dist**(n/2.0)

    DO i=1,n
       cluster_volume = cluster_volume * sqrt(evals(i))
    END DO

    DO i=1,n
       DO j=1,n
          mat(i,j) = mat(i,j) * sqrt(evals(j))
       END DO
    END DO

    ! Generate point array.
    IF(flag) THEN

       DO k=1,nsurf_points
          surf_points(k,:) = 0.0
       
          CALL rng_sphere_surf_nd(point, n)
          
          DO i=1,n
             DO j=1,n
                surf_points(k,i) = surf_points(k,i) + point(j) * mat(i,j)
             END DO
          END DO
          
          surf_points(k,:) = (sqrt(max_dist)*surf_points(k,:)) + means

       END DO

    ! Set metric.
    ELSE

       metric = met / max_dist

    END IF


    ! Zero sep variables.
    ts1 = 0.0
    ts2 = 0.0

    ! Calculate 'radius' of surface along line of centres.
    DO i=1,n
       
       ! Find norm of centre sep vector.
       ts1 = ts1 + (m1(i) - m2(i))**2

       tv1(i) = 0.0

       ! Transform into eigen basis, and scale.
       DO j=1,n
          tv1(i) = tv1(i) + mat(j,i) * (m1(j) - m2(j))
       END DO
       
       ! Find norm of scaled vector.
       ts2 = ts2 + tv1(i)**2

    END DO

    ! Scale by scaling factor, divide by norm of original vector.
    min_dist = sqrt(ts2 * max_dist / ts1)


    !PRINT *,"CV: ",cluster_volume
    !min_dist = sqrt(max_dist*evals(n))

  END FUNCTION cluster_volume

END MODULE twomeans
