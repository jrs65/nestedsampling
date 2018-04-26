! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

MODULE treedef

! Include precision definitions
#include "defs.H90"
  
  IMPLICIT none

  ! Define RED/BLACK parameters
  INTEGER(KIND=1), PARAMETER :: RED = 0
  INTEGER(KIND=1), PARAMETER :: BLACK = 1
  
  ! Define LEFT RIGHT parameters.
  INTEGER(KIND=1), PARAMETER :: LEFT = -1
  INTEGER(KIND=1), PARAMETER :: ROOT = 0
  INTEGER(KIND=1), PARAMETER :: RIGHT = 1
  

  ! pnode
  ! Point tree type, store sample points in binary 
  ! tree structure for fast searching, insertion etc.
  TYPE, PUBLIC :: pnode
     
     ! left,right subtrees
     TYPE(pnode), POINTER :: left, right, parent

     ! Red/Black state. 0 if red. 1 if black.
     INTEGER(KIND=1) :: rb

     ! Left/Right state. -1 if left. 1 if right. 0 if root.
     INTEGER(KIND=1) :: lr     

     ! Likelihood value
     REAL(rk) :: lval
     
     ! Prior value
     REAL(rk) :: xval
     
     ! Accumulated Prior Mass
     REAL(rk) :: axmval

     ! Weight
     REAL(rk) :: w
     
     ! Vector of point coordinates
     REAL(rk), DIMENSION(:), POINTER :: pvec
     
  END TYPE pnode
  
END MODULE treedef
