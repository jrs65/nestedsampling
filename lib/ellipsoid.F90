! A sampler using an ellipsoid bounding method for
! replacements. Initial sampler is the simple initial_metropolis. For
! replacement, a principle ellipsoid is found from the covariance
! matrix and scaled to include all points. We then perform MH, within
! this hard bound. See Mukherjee et al. 2006 (arXiv:astro-ph/0508461)
MODULE ellipsoid

! Include precision definitions
#include "defs.H90"

USE sampling
USE metropolis, ONLY: initial_metropolis, burnin_metropolis
USE tree
USE treeutil
USE rand
USE ellipsoidutil
!USE nested

IMPLICIT NONE

! Enlargement method
INTEGER, PRIVATE :: use_scaling = 1

! Number of replacements so far.
INTEGER, PRIVATE :: n = 0

! Evaluate ellipse every n_el iterations.
INTEGER, PRIVATE :: n_el = 10



! Initialise the ellipsoid sampling routine.
!
! Burn-in Metropolis for initial selection. Set-up the work
! arrays. Set an enlargement factor.
INTERFACE init_ellipsoid

   ! New version
   MODULE PROCEDURE init_ellipsoid_short
   ! For compatibility.
   MODULE PROCEDURE init_ellipsoid_long

END INTERFACE



CONTAINS



  ! Set the enlargement method. 0 is no scaling (other than
  ! enlargement factor), any other val is max_dist (include all
  ! points).
  SUBROUTINE set_scalingmethod(val)

    INTEGER, INTENT(IN) :: val

    use_scaling = val

  END SUBROUTINE set_scalingmethod


 
  ! Evaluate ellipse only every 10 iterations.
  SUBROUTINE set_ellskip(n)

    INTEGER, INTENT(IN) :: n

    n_el = n

  END SUBROUTINE set_ellskip


  
  ! Force re-evalutation of ellipse.
  SUBROUTINE reset_ellskip()

    n = 0

  END SUBROUTINE reset_ellskip



  ! Initialise the ellipsoid sampling routine.
  !
  ! Burn-in Metropolis for initial selection. Set-up the work
  ! arrays. Set an enlargement factor.
  SUBROUTINE init_ellipsoid_long(dpars, en, L, X)

    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars
    ! Enlargement factor
    REAL(rk) :: en

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

    CALL init_ellipsoidutil(dpars)

    CALL set_enlargement(en)

  END SUBROUTINE init_ellipsoid_long


  SUBROUTINE init_ellipsoid_short(dpars, en)

    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars
    ! Enlargement factor
    REAL(rk) :: en

    CALL init_ellipsoidutil(dpars)

    CALL set_enlargement(en)

  END SUBROUTINE init_ellipsoid_short





  ! Select replacement point by the Metropolis method with an
  ! elliptical bound defined by the covariance ellipsoid. Within the
  ! ellipsoid we evolve n_iter times from a randomly selected point.
  !
  ! Effectively this is MH, with a fixed proposal distribution uniform
  ! within an ellipse.
  SUBROUTINE replace_ellipsoid(tree, n_iter, lbound, axmval, dpars, L, X)
    
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

    ! Vectors for evolving parameters
    REAL(rk), DIMENSION(dpars) :: rvec, par1, par2
    REAL(rk) :: like, x1, x2, r, sc
    INTEGER :: i

    TYPE(pnode), POINTER :: node

    ! For debugging
#ifdef _DEBUG_EL
    CHARACTER(LEN=10) :: file
#endif

    i = 0

    IF(MOD(n,n_el) == 0) THEN
       ! Generate the covariance matrix of the active points.
       CALL covmatrix(tree)
       ! Generate the transformation matrix to transform sphere->ellipse
       CALL transmatrix()
       
       ! Calculate the scaling to include all points, if use_scaling set.
       IF(use_scaling /= 0) THEN
          CALL calc_scaling(tree)
       END IF
    END IF
       

    ! Initial select any random point in ellipse.
    CALL rand_ellipsoid(par1)
    x1 = X(par1)
       
    DO WHILE(i < n_iter)
       
       CALL rand_ellipsoid(par2)

       ! Increment transition counter
       r_tnum = r_tnum + 1
       
       x2 = X(par2)
    
       ! Standard Metropolis selection.
       IF (x2 > x1) THEN
          par1 = par2
          x1 = x2
          
          i = i + 1
          
       ELSE
          r = rng_uniform()
          
          IF(r < EXP(x2 -x1)) THEN
             par1 = par2
             x1 = x2
             
             i = i + 1
             
          END IF
       END IF
        
    END DO

    like = L(par1)
    ! Increment likelihood counter
    r_lnum = r_lnum + 1
!    WRITE(*,"(A,G20.10)")   "Bound",lbound
!    WRITE(*,"(A,20G20.10)") "Trial",par1,like

    ! Keep within likelihood contour.
    ! If not repeat process.
    ! Double negative to get around NaNs
    DO WHILE(like < lbound .OR. isnan(like))

       CALL rand_ellipsoid(par2)

       ! Increment transition counter
       r_tnum = r_tnum + 1
       
       x2 = X(par2)
    
       ! Standard Metropolis selection.
       IF (x2 > x1) THEN
          par1 = par2
          x1 = x2
          
          i = i + 1

          like = L(par1)
          ! Increment likelihood counter
          r_lnum = r_lnum + 1

!          WRITE(*,"(A,20G20.10)") "Trial",par1,like
          
       ELSE
          r = rng_uniform()
          
          IF(r < EXP(x2 -x1)) THEN
             par1 = par2
             x1 = x2
             
             i = i + 1

             like = L(par1)
             ! Increment likelihood counter
             r_lnum = r_lnum + 1

!             WRITE(*,"(A,20G20.10)") "Trial",par1,like

          END IF
       END IF

    END DO
    
!    PRINT *,"Accepted"

    ! Debug
#ifdef _DEBUG_EL
    IF(MOD(n, _DEBUG_EL) == 0) THEN
       WRITE (file, '(A, I5.5)') 'inter', n
!       CALL print_ptree(tree, file)
       WRITE (file, '(A, I5.5)') 'bell', n
       CALL ellipsoid_surface(100, file)
     END IF
#endif
     
     n = n+1
     
     ALLOCATE(node)
     ALLOCATE(node % pvec(dpars))
     
     node % pvec = par1
     node % lval = like
     node % xval = x1
     CALL add_node(node, tree)
     CALL write_resume_node(node)
!    CALL check_tree(tree)
    ! Increment replacement counter.
    r_num = r_num + 1
    ! Increment accepted transitions counter.
    r_anum = r_anum + n_iter

  END SUBROUTINE replace_ellipsoid



  ! Return the volume of the ellipsoid (less the n-d angular factor).
  REAL(rk) FUNCTION ellipsoid_volume(tree)

    TYPE(pnode), POINTER :: tree

    ! Generate the covariance matrix of the active points.
    CALL covmatrix(tree)
    ! Generate the transformation matrix to transform sphere->ellipse
    CALL transmatrix()

    ! Return the volume.
    ellipsoid_volume = det_covmatrix()


  END FUNCTION ellipsoid_volume
    
END MODULE ellipsoid
