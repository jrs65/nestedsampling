MODULE sampling


! Include precision definitions
#include "defs.H90"

USE tree
USE treeutil

IMPLICIT NONE


! Scalelengths vector
REAL(rk), DIMENSION(:), POINTER :: scalelength
! Initial Parameters
REAL(rk), DIMENSION(:), POINTER :: initial_par
! Initial widths (1 min width, 2 max width)
REAL(rk), DIMENSION(:,:), POINTER :: initial_widths

! Filehandle to read an write resuming data to.
INTEGER :: resume_fh
! To resume, or not to resume!
LOGICAL :: resume = .FALSE.
! Continue reading from resume file, true, then read. False, then write.
LOGICAL :: resume_read = .TRUE.


! Efficiency monitoring
!
! Initial
! Total and accepted transitions respectively.
INTEGER :: i_tnum = 0
INTEGER :: i_anum = 0

! Replacement
! Total and accepted transitions respectively.
INTEGER :: r_tnum = 0
INTEGER :: r_anum = 0

! Total number points replaced
INTEGER :: r_num = 0
! Number of likelihood evaluations
INTEGER :: r_lnum = 0


CONTAINS



  ! Set the scale parameter vector.
  ! Scales the jumps in each direction.
  SUBROUTINE set_scale(scale) 
    REAL(rk), DIMENSION(:), POINTER :: scale

    scalelength => scale

  END SUBROUTINE set_scale



  ! Set the scale parameter vector.
  ! Scales the jumps in each direction.
  SUBROUTINE set_initialpars(pars) 
    REAL(rk), DIMENSION(:), POINTER :: pars

    initial_par => pars

  END SUBROUTINE set_initialpars



  ! Set the min/max width vector.
  ! Used to burn in for acceptable acceptance rate (~0.5)
  SUBROUTINE set_initialwidths(widths) 
    REAL(rk), DIMENSION(:,:), POINTER :: widths

    initial_widths => widths

  END SUBROUTINE set_initialwidths



  ! Reset the sampling efficiency data.
  SUBROUTINE reset_sampling_eff() 
    
    i_tnum = 0
    i_anum = 0
    r_tnum = 0
    r_anum = 0
    r_num = 0
    
  END SUBROUTINE reset_sampling_eff

  
  ! Print out data on sampling efficiency.
  SUBROUTINE print_sampling_eff()

    PRINT *,"Sampling efficiency"
    
    PRINT *,"Initial"
    PRINT *,i_tnum,"Transitions"
    PRINT *,(1.0*i_anum)/i_tnum,"Accepted transitions rate"

    PRINT *,"Replacements"
    PRINT *,r_tnum,"Transitions"
    PRINT *,(1.0*r_anum)/r_tnum,"Accepted transitions"
    PRINT *,r_num,"Replacements"
    PRINT *,(1.0*r_lnum)/r_num,"L per repl"

  END SUBROUTINE print_sampling_eff



  ! Set the filehandle for resuming, must be rw.
  SUBROUTINE set_resume_fh(fh)

    ! Filehandle to set.
    INTEGER, INTENT(IN) :: fh

    resume = .TRUE.
    resume_fh = fh

  END SUBROUTINE set_resume_fh


  ! Read in node from the resume file.
  LOGICAL FUNCTION read_resume_node(node, n)

    ! Node to return.
    TYPE(pnode), POINTER :: node
    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: n

    INTEGER :: is

    IF(.NOT. resume_read .OR. .NOT. resume) THEN
       read_resume_node = .FALSE.
       RETURN
    END IF

    ALLOCATE(node)
    ALLOCATE(node % pvec(n))
    
    ! Read in node.
    READ (resume_fh,*,IOSTAT=is) node % pvec, node % xval, node % lval
    
    ! Check status.
    ! No more resume data.
    IF(is < 0) THEN
       read_resume_node = .FALSE.
       resume_read = .FALSE.
    ! Read error.
    ELSE IF(is > 0) THEN
       PRINT *,"Error reading node."
       resume_read = .FALSE.
       read_resume_node = .FALSE.
    ! Read finished correctly.
    ELSE
       read_resume_node = .TRUE.
    END IF

  END FUNCTION read_resume_node


  ! Write node data to resume file.
  SUBROUTINE write_resume_node(node)

    TYPE(pnode), POINTER :: node

    IF(resume_read .OR. .NOT. resume) THEN
       RETURN
    END IF

    WRITE (resume_fh,*) node % pvec, node % xval, node % lval

  END SUBROUTINE write_resume_node

  
  ! Initial point resume. Read in initial n points. If succesful
  ! return true, if less than n points, return number remaining n_rem,
  ! and return false.
  LOGICAL FUNCTION read_resume_initial(tree, np, n)

    ! Tree of nodes to return.
    TYPE(pnode), POINTER :: tree

    ! Number of points to read.
    INTEGER, INTENT(IN) :: np
    ! Dimensionality of parameter space.
    INTEGER, INTENT(IN) :: n

    INTEGER :: i

    TYPE(pnode), POINTER :: node

    read_resume_initial = .TRUE.

    IF(.NOT. resume_read .OR. .NOT. resume) THEN
       read_resume_initial = .FALSE.
       RETURN
    END IF

    DO i=1,np

       IF(read_resume_node(node, n)) THEN
          CALL add_node(node, tree)
       ELSE
          read_resume_initial = .FALSE.
          PRINT *,"Error: could not read all initial points."
          EXIT
       END IF

    END DO

  END FUNCTION read_resume_initial


  SUBROUTINE write_resume_tree(tree)

    ! Tree of nodes to write.
    TYPE(pnode), POINTER :: tree

    TYPE(pnode), POINTER :: node

    IF(resume_read .OR. .NOT. resume) THEN
       RETURN
    END IF

    NULLIFY(node)

    DO WHILE(next_node(node, tree)) 

       CALL write_resume_node(node)

    END DO

  END SUBROUTINE write_resume_tree

END MODULE sampling
