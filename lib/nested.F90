! Copyright (C) 2006 Richard Shaw <jrs65 at cam.ac.uk>
! Licensed under GNU GPL v2. See COPYING or 
! http://www.gnu.org/copyleft/gpl.html

MODULE nested

! Include precision definitions
#include "defs.H90"

USE treedef
USE tree
USE treeutil

USE rand
USE rand_nested

IMPLICIT NONE


!! Number of points within likelihood contour at any one time.
INTEGER :: npoints = 500
!! Number of remove-replace cycles.
INTEGER :: nrepl = 1000
!! Number of iterations to create replacement.
INTEGER :: niter = 10
!! Number of evidence estimates to calculate.
INTEGER :: nevid = 10


!! Fraction of volume remaining for auto stop.
REAL(rk) :: volfrac = 10e-3
  
CONTAINS


  !! Set the number of points
  SUBROUTINE set_npoints(np)
    INTEGER :: np
    
    npoints = np

  END SUBROUTINE set_npoints


  !! Set the number of replacements
  SUBROUTINE set_nrepl(nr)
    INTEGER :: nr
    
    nrepl = nr

  END SUBROUTINE set_nrepl


  !! Set the number of iterations/replacement
  SUBROUTINE set_niter(ni)
    INTEGER :: ni
    
    niter = ni

  END SUBROUTINE set_niter


  !! Set the number of evidence calculations to make.
  SUBROUTINE set_nevid(ne)
    INTEGER :: ne
    
    nevid = ne

  END SUBROUTINE set_nevid




  !! A progress meter drawing routine. Prints the percent completion
  !! of a task. Need to call for every percent update you wish to
  !! make.
  SUBROUTINE progress_meter(percent) 

    !! Percent completed.
    INTEGER, INTENT(IN) :: percent

    WRITE(*,'(a,i0,2a)',advance='no') '    Progress:   ',percent,' %',char(13) 
    
    IF(percent == 100) THEN
       PRINT *,'    Done    : '
    END IF

  END SUBROUTINE progress_meter



  !! A progress meter drawing routine. Prints the percent completion
  !! of a task. Need to call for every percent update you wish to
  !! make.
  SUBROUTINE progress_meter1(af, tf) 

    REAL(rk), INTENT(IN) :: af
    REAL(rk), INTENT(IN) :: tf

#ifndef _COSMOS_OUT
    WRITE(*,'(a,f10.3,a,f10.3,a)',advance='no') '    Target:  ',tf,'     Current:  ',af,char(13)
#else
    PRINT *,'    Target:  ',tf,'     Current:  ',af
#endif
    IF(af <= tf) THEN
       PRINT *,' '
    END IF

  END SUBROUTINE progress_meter1


  
  !! Set the fraction of remaining volume to stop stop sampling at.
  SUBROUTINE set_volfrac(vf) 

    REAL(rk), INTENT(IN) :: vf

    volfrac = vf

  END SUBROUTINE set_volfrac


  LOGICAL FUNCTION  text_dialog(txt)

    ! Text to display.
    CHARACTER(LEN=*), INTENT(IN) :: txt
    CHARACTER(LEN=1) :: input

    PRINT *,""
    PRINT *,txt,"[y/N]"
    READ(*,*) input

    IF(input == "y" .OR. input == "Y") THEN
       text_dialog = .TRUE.
    ELSE
       text_dialog = .FALSE.
    END IF

  END FUNCTION text_dialog

    
END MODULE nested

      
