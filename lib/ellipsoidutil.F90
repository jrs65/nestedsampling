! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

MODULE ellipsoidutil

! Include precision definitions
#include "defs.H90"

USE tree
USE treeutil
USE rand

IMPLICIT NONE



! Sigma enlargement factor.
REAL(rk), PRIVATE :: enlargement

! Internal store for the covariance matrix.
REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: cov_mat
! Internal store for transformation matrix.
REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: trans_mat
! Internal store for elliptical metric.
REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: metric_mat
! Internal store for the means.
REAL(rk), DIMENSION(:), POINTER, PRIVATE :: means
! Internal store for the eigenvalues (axes lengths).
!REAL(rk), DIMENSION(:), POINTER, PRIVATE :: evals
! Internal store for the number of points.
INTEGER, PRIVATE :: n_points
! The maximum square elliptical distance in the data set.
REAL(rk), PRIVATE :: max_dist
! The scaling required to include all points.
REAL(rk), PRIVATE :: scaling

! Storage for lapack routines.
REAL(rk), DIMENSION(:), POINTER, PRIVATE :: work
INTEGER, PRIVATE :: lwork

#ifdef _DEBUG_EL
INTEGER, PRIVATE :: num = 0
#endif

CONTAINS

  
  ! Set enlargement factor, relative to principal variances.
  SUBROUTINE set_enlargement(en)
    REAL(rk), INTENT(IN) :: en

    enlargement = en

  END SUBROUTINE set_enlargement

  
  
  ! Worker routine to walk the tree and do summations required for
  ! covariance matrix.
  SUBROUTINE cov_calc(node)

    TYPE(pnode), POINTER :: node

    INTEGER :: i,j

    n_points = n_points + 1
    means = means + (node % pvec)
    
    DO i = 1,size(node % pvec)
       DO j=1,size(node % pvec)

          cov_mat(i,j) = cov_mat(i,j) + (node % pvec(i))*(node % pvec(j))

       END DO

    END DO

  END SUBROUTINE cov_calc



  ! Worker routine to walk the tree and calculate the maximum
  ! elliptical distance.
  SUBROUTINE el_dist(node)

    TYPE(pnode), POINTER :: node

    INTEGER :: i,j,n
    REAL(rk) :: dist

    n = size(node % pvec)
    dist = 0.0

    ! Calculate inner product with elliptical metric.
    DO i=1,n
       DO j=1,n
          dist = dist + ((node % pvec(i) - means(i)) * metric_mat(i,j) * (node % pvec(j) - means(j)))
       END DO
    END DO

    ! If is greater than current max, update distance.
    IF(dist > max_dist) THEN
       max_dist = dist
    END IF

  END SUBROUTINE el_dist
       

  ! Calculate the covariance matrix of the nodes in tree
  ! and store the results in cov_mat
  SUBROUTINE covmatrix(tree)

    TYPE(pnode), POINTER :: tree

    INTEGER :: l,i,j
    
    l = size(tree % pvec)

    n_points = 0
    cov_mat = 0
    means = 0

    CALL walk_inorder(tree, cov_calc)

    ! Debug Section
    ! Second moments matrix.
    ! CALL print_cov()
    ! Means.
    ! PRINT *,means
    ! PRINT *,n_points

    ! Convert sums into moments.
    cov_mat = cov_mat / n_points
    means = means / n_points

    ! Generate covariances.
    DO i = 1,l
       DO j=1,l
          cov_mat(i,j) = cov_mat(i,j) - means(i) * means(j)

       END DO

    END DO

    !PRINT *,"Points:",n_points
    !CALL print_cov()
    !CALL flush()

  END SUBROUTINE covmatrix



  ! Calculate the point transformation matrix.
  SUBROUTINE transmatrix()

    INTEGER :: n, info
    REAL(rk), DIMENSION(:), ALLOCATABLE :: evals

    INTEGER :: i,j,k

    n = size(cov_mat(1,:))

    ALLOCATE(evals(n))
    trans_mat = cov_mat

    ! Find the eigenvalues and vectors of the covariance matrix.
    CALL dsyev('V', 'U', n, trans_mat, n, evals, work, lwork, info)

    IF(info /= 0) THEN
       PRINT *,"dsyev exited with code:",info
       CALL print_cov()
       CALL exit(-1)
    END IF

    ! Debug
#ifdef _DEBUG_EL

    IF(MOD(num, _DEBUG_EL) == 0) THEN
       !PRINT *,"Eccentricities",num
       !PRINT *,sqrt(1 - (evals(1:n-1)**2 / (evals(n)**2)))
       !PRINT *,"Ellipticities",num
       PRINT *,evals(1) / (evals(n))
    END IF

    num = num + 1
    
#endif


    ! Construct the metric matrix to allow computation of elliptical
    ! distance, for use in finding scaling of ellipse bound.
    DO i=1,n
       ! Set to zero for sum.
       metric_mat(i,:) = 0.0

       DO j=1,n
          DO k=1,n
             metric_mat(i,j) = metric_mat(i,j) &
                  + (trans_mat(i,k) * trans_mat(j,k) / evals(k))
          END DO
       END DO
    END DO
    
    ! Construct a transformation matrix. First rescale along principle
    ! coordinates (multiply by square of evals), then transform into
    ! real coordinate frame.
    DO i=1,n
       DO j=1,n
          trans_mat(i,j) = trans_mat(i,j) * sqrt(evals(j))
       END DO
    END DO

    DEALLOCATE(evals)

    ! Debug
    ! CALL print_trans()

  END SUBROUTINE transmatrix



  ! Calculate the scaling required to include all points in the active
  ! data set.
  SUBROUTINE calc_scaling(tree)

    TYPE(pnode), POINTER :: tree

    max_dist = 0.0

    CALL walk_inorder(tree, el_dist)

    scaling = sqrt(max_dist)

  END SUBROUTINE calc_scaling

  

  ! Init ellipsoid work space.
  SUBROUTINE init_ellipsoidutil(n)
    
    ! Dimension of parameter space.
    INTEGER, INTENT(IN) :: n

    ! Allocate a work array for lapack.
    ALLOCATE(work(n*70))
    lwork = n*70

    ! Allocate internal matrices and vectors.
    ALLOCATE(means(n))
    means = 0
    !ALLOCATE(evals(n))
    !evals = 0
    ALLOCATE(cov_mat(n,n))
    cov_mat = 0
    ALLOCATE(trans_mat(n,n))
    trans_mat = 0
    ALLOCATE(metric_mat(n,n))
    metric_mat = 0

    enlargement = 1.0
    scaling = 1.0

  END SUBROUTINE init_ellipsoidutil



  ! Generate a point randomly distributed in the ellipsoid.
  SUBROUTINE rand_ellipsoid(point)

    
    REAL(rk), DIMENSION(:), INTENT(OUT) :: point  

    INTEGER :: i,j,n

    REAL(rk), DIMENSION(:), ALLOCATABLE :: point2

    n = size(point)
    ALLOCATE(point2(n))
    point2 = 0

    CALL rng_sphere_nd(point, n)

    IF(n /= size(trans_mat(1,:))) THEN
       PRINT *,"Dimensions incorrect."
    END IF

    DO i=1,n
       DO j=1,n
          point2(i) = point2(i) + point(j) * trans_mat(i,j)
       END DO
    END DO

    point = (scaling*enlargement*point2) + means

    DEALLOCATE(point2)

  END SUBROUTINE rand_ellipsoid


  
  ! Test subr, set covar so don't need point tree.
  SUBROUTINE set_covar(cov)
    REAL(rk), DIMENSION(:,:), POINTER :: cov

    cov_mat => cov

  END SUBROUTINE set_covar



  ! Print the covariance matrix.
  SUBROUTINE print_cov()

    INTEGER :: i
    PRINT *,"Covariance matrix:"

    DO i=1,size(cov_mat(:,1))
       PRINT *,cov_mat(i,:)
    END DO

  END SUBROUTINE print_cov



  ! Print the transformation matrix.
  SUBROUTINE print_trans()

    INTEGER :: i

    PRINT *,"Transformation matrix:"
    
    DO i=1,size(trans_mat(:,1))
       PRINT *,trans_mat(i,:)
    END DO

  END SUBROUTINE print_trans


  ! Generate a point randomly distributed in the ellipsoid.
  SUBROUTINE ellipsoid_surface(num, file)

    INTEGER, INTENT(IN) :: num
    CHARACTER(LEN=*), INTENT(IN) :: file

    REAL(rk), DIMENSION(:), ALLOCATABLE :: point  

    INTEGER :: i,j,n,k

    REAL(rk), DIMENSION(:), ALLOCATABLE :: point2

    OPEN(8,FILE=file)

    n = size(trans_mat(1,:))
    ALLOCATE(point(n))
    ALLOCATE(point2(n))

    DO k=1,num
       point2 = 0
       
       CALL rng_sphere_surf_nd(point, n)

       DO i=1,n
          DO j=1,n
             point2(i) = point2(i) + point(j) * trans_mat(i,j)
          END DO
       END DO
       
       point = (enlargement*scaling*point2) + means

       WRITE(8,*) point,"0"

    END DO

    CLOSE(8)

    DEALLOCATE(point, point2)

  END SUBROUTINE ellipsoid_surface



  ! Return the determinant of the covariance matrix (product of
  ! eigenvalues).
  REAL(rk) FUNCTION det_covmatrix() 

    INTEGER :: i

    det_covmatrix = 1
    !DO i=1,size(evals)
       !det_covmatrix = det_covmatrix * evals(i)
    !END DO

  END FUNCTION det_covmatrix

END MODULE ellipsoidutil


