! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

! A simple sampler based on Metropolis-Hastings. Contains, an initial
! sampler (including a burn-in method, and a less useful bounded
! replacement sampler.
MODULE metropolis

! Include precision definitions
#include "defs.H90"

USE rand
USE tree
USE treeutil
USE sampling

IMPLICIT NONE


CONTAINS



  ! The Likelihood function
  REAL(rk) FUNCTION dummy_like(ppars)
    REAL(rk), DIMENSION(:) :: ppars

    dummy_like = 0.0

  END FUNCTION dummy_like

  
  ! Burn in. Adjust width parameters to maintain rejection rate.
  ! Make 5 iterations of 200 runs to adjust widths.
  SUBROUTINE burnin_metropolis(dpars, X)

    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: dpars

    INTERFACE
    
       ! The Prior function
       REAL(rk) FUNCTION X(ppars)
         REAL(rk), DIMENSION(:) :: ppars
       END FUNCTION X

    END INTERFACE

    ! The point tree to add to.
    TYPE(pnode), POINTER :: tree

    INTEGER :: i


    REAL(rk), DIMENSION(dpars) :: sc1, sc2, sc
    REAL(rk) :: ar1, ar2, ar

    sc1 = initial_widths(1,:)
    sc2 = initial_widths(2,:)

    NULLIFY(tree)
    IF(ASSOCIATED(scalelength)) THEN
       DEALLOCATE(scalelength)
    END IF

    ALLOCATE(scalelength(dpars))    
    ! Assume larger proposal gives greater rejection.
    scalelength = sc1
    CALL initial_metropolis(tree, 200, dpars, dummy_like, X)
    ar1 = 1.0*i_anum / i_tnum

    i_tnum = 0
    i_anum = 0

    DEALLOCATE(tree)
    NULLIFY(tree)

    scalelength = sc2
    CALL initial_metropolis(tree, 200, dpars, dummy_like, X)
    ar2 = 1.0*i_anum / i_tnum
    
    i_tnum = 0
    i_anum = 0

    DEALLOCATE(tree)
    NULLIFY(tree)

    DO i=1,5
       sc = (sc1 + sc2) / 2.0
       scalelength = sc
       NULLIFY(tree)

       CALL initial_metropolis(tree, 200, dpars, dummy_like, X)
       ! If ar > 0.5 move to larger proposal size.
       ar = 1.0*i_anum / i_tnum
       IF(ar > 0.5) THEN
          sc1 = sc
          ar1 = ar
       ELSE
          sc2 = sc
          ar2 = ar
       END IF

       DEALLOCATE(tree)

       i_tnum = 0
       i_anum = 0

    END DO

    scalelength = (sc1 + sc2) / 2.0
    NULLIFY(tree)

    CALL initial_metropolis(tree, 1000, dpars, dummy_like, X)

    DEALLOCATE(tree)

    i_tnum = 0
    i_anum = 0

  END SUBROUTINE burnin_metropolis



  ! Select initial points by the Metropolis-Hastings method
  ! Add the initial n points.
  SUBROUTINE initial_metropolis(tree, n, dpars, L, X)
    
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

    ! Vectors for evolving parameters
    REAL(rk), DIMENSION(dpars) :: rvec, par1, par2
    REAL(rk) :: like, x1, x2, r
    INTEGER :: i

    i = 0

    par1 = initial_par
    x1 = X(par1)
    ! PRINT *,scalelength
    ! Loop to add inital points. Use standard metropolis selection.
    DO WHILE(i < n)
       CALL rng_gaussian_arr(rvec, dpars)
       par2 = par1 + (rvec * scalelength)

       x2 = X(par2)

       ! Increment transition counter
       i_tnum = i_tnum + 1
       
       IF (x2 > x1) THEN
          par1 = par2
          x1 = x2
          like = L(par1)

          ! Insert point into tree.
          !PRINT *,i,x1,like
          CALL insert_node(tree, x1, like, par1)
 !         CALL check_tree(tree)

          i = i + 1
          
          
       ELSE
          r = rng_uniform()
          IF(r < EXP(x2 -x1)) THEN
             par1 = par2
             x1 = x2
             like = L(par1)

             ! Insert point into tree.
             !PRINT *,i,x1,like
             CALL insert_node(tree, x1, like, par1)
!             CALL check_tree(tree)

             i = i + 1
          END IF
       END IF
    END DO

    ! Increment accepted transition counter.
    i_anum = i_anum + n

  END SUBROUTINE initial_metropolis


             
  ! Select replacement point by the Metropolis-Hastings method.
  ! Evolve n_iter times from a randomly selected point.
  SUBROUTINE replace_metropolis(tree, n_iter, lbound, axmval, dpars, L, X)
    
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

    i = 0

    sc = axmval ** (1.0/dpars)

    IF(ASSOCIATED(tree)) THEN
       par1 = tree % pvec
    ELSE
       par1 = initial_par
    END IF

    x1 = X(par1)
       
    DO WHILE(i < n_iter)
       CALL rng_gaussian_arr(rvec, dpars)
       
       par2 = par1 + (rvec * sc* scalelength)

       ! Increment transition counter
       r_tnum = r_tnum + 1
       
       x2 = X(par2)
    
       like = L(par2)
       ! Increment likelihood counter
       r_lnum = r_lnum + 1
       
       ! Keep within likelihood contour.
       ! If not repeat process.
       IF (like > lbound) THEN
       
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
       END IF
       
    END DO

    CALL insert_node(tree, x1, like, par1)

    ! Increment replacement counter.
    r_num = r_num + 1
    ! Increment accepted transitions counter.
    r_anum = r_anum + n_iter

  END SUBROUTINE replace_metropolis             
      

END MODULE
