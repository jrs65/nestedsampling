!! A module wrapping around the generator provided by rand_mt. Provides
!! functions for single, and arrays of uniform random numbers. Plus
!! functions for generating gaussian random variables, and points
!! uniformly within, and on the surface of an n-sphere.
MODULE rand

! Include precision definitions
#include "defs.H90"

  USE mt95

  IMPLICIT NONE

  INTEGER, PRIVATE :: init = 0

  INTEGER, PRIVATE :: gflag
  REAL(rk), PRIVATE :: rng_g

  CONTAINS



    !! Initialise Random number stuff
    SUBROUTINE init_rng()
      
      REAL :: rs
      INTEGER, DIMENSION(8) :: is
      INTEGER :: seed
      
      CALL date_and_time(values = is)

      seed = is(7) + is(8)
      PRINT *,"Seed",seed

      CALL genrand_init(seed)

      init = 1

    END SUBROUTINE init_rng



    !! Generate a uniform random number in the range [0,1]
    REAL(rk) FUNCTION rng_uniform()

      IF(init == 0) THEN
         CALL init_rng()
      END IF

      CALL genrand_real1(rng_uniform)

    END FUNCTION rng_uniform



    !! Generate an array of uniform random numbers.
    SUBROUTINE rng_uniform_arr(arr, len)

      REAL(rk), DIMENSION(:), INTENT(OUT) :: arr
      INTEGER, INTENT(IN) :: len

      IF(init == 0) THEN
         CALL init_rng()
      END IF

      CALL genrand_real1(arr)

    END SUBROUTINE rng_uniform_arr
    


    !! Generate a gaussian distributed random number.
    REAL(rk) FUNCTION rng_gaussian()

      REAL(rk) :: x1,x2,r,f
      REAL(rk), DIMENSION(2) :: t

      IF(gflag == 0) THEN
         
         IF(init == 0) THEN
            CALL init_rng()
         END IF
         
         DO
            CALL genrand_real1(t)

            x1 = 2.0 * t(1) - 1.0
            x2 = 2.0 * t(2) - 1.0
            
            r = x1 * x1 + x2 * x2
            
            IF(r < 1.0 .AND. r /= 0) EXIT
         END DO
         
         f = SQRT(-2.0 * LOG(r) / r)
         
         rng_g = x1 * f
         gflag = 1
         rng_gaussian = x2 * f
         
      ELSE

         rng_gaussian = rng_g
         gflag = 0

      END IF

    END FUNCTION rng_gaussian

      

    !! Generate a gaussian distributed random number.
    SUBROUTINE rng_gaussian_arr(arr, len)

      REAL(rk), INTENT(OUT), DIMENSION(:) :: arr
      INTEGER, INTENT(IN) :: len

      INTEGER :: l1,i
      REAL(rk) :: x1,x2,r,f
      REAL(rk), DIMENSION(2) :: t

      l1  = (len / 2)

      IF(init == 0) THEN
         CALL init_rng()
      END IF
      

      DO i=1,l1
         
         DO   
            CALL genrand_real1(t)

            x1 = 2.0 * t(1) - 1.0
            x2 = 2.0 * t(2) - 1.0
            
            r = x1 * x1 + x2 * x2
               
            IF(r < 1.0 .AND. r /= 0) EXIT
         END DO
         
         f = SQRT(-2.0 * LOG(r) / r)
         
         arr(2*i-1) = x1 * f
         arr(2*i) = x2 * f
         
      END DO

      IF(2*l1 < len) THEN
         IF(gflag == 1) THEN
            arr(len) = rng_g
            gflag = 0
         ELSE
            arr(len) = rng_gaussian()
         END IF
      END IF
         
      
    END SUBROUTINE rng_gaussian_arr

      

    !! Select a point randomnly distributed within the 
    !! boundary of an n-d unit sphere.
    !! 
    !! First select a point uniformly distributed in the surface
    !! of an n-sphere. Then scale for a point in the interior.
    SUBROUTINE rng_sphere_nd(point, n)

      REAL(rk), DIMENSION(:), INTENT(OUT) :: point
      
      INTEGER, INTENT(IN) :: n

      INTEGER :: i
      REAL(rk) :: t

      t = 0

      CALL rng_gaussian_arr(point, n)

      DO i=1,n
         t = t + point(i)**2
      END DO

      ! Map point to surface.
      t = sqrt(t)
      ! Map to interior point. Power ensures uniform in volume.
      t = (rng_uniform()**(1.0/n)) / t

      point = point * t

    END SUBROUTINE rng_sphere_nd



    !! Select a point randomnly distributed on the surface of an n-d
    !! unit sphere.
    SUBROUTINE rng_sphere_surf_nd(point, n)

      REAL(rk), DIMENSION(:), INTENT(OUT) :: point
      
      INTEGER, INTENT(IN) :: n

      INTEGER :: i
      REAL(rk) :: t

      t = 0

      CALL rng_gaussian_arr(point, n)

      DO i=1,n
         t = t + point(i)**2
      END DO

      ! Map point to surface.
      t = sqrt(t)
      ! Map to interior point. Power ensures uniform in volume.

      point = point / t

    END SUBROUTINE rng_sphere_surf_nd
      
END MODULE rand

