MODULE distr

#include "defs.H90"

  IMPLICIT NONE

  ! A set of parameters for calculating the prior.
  REAL(rk), DIMENSION(:,:), POINTER, PRIVATE :: prior_params

  ! Pi
  REAL(rk), PARAMETER :: pi = 3.1415926535
  ! Smallest representable number. For excluded volumes.
  REAL(rk) :: min_inf = -1.0 * HUGE(pi)

  REAL(rk), DIMENSION(:), PRIVATE, POINTER :: cent
  REAL(rk), PRIVATE :: width = 1

  PRIVATE :: fact

  CONTAINS


    ! Top hat
    SUBROUTINE set_tophat(w, c)
      
      REAL(rk), INTENT(IN) :: w
      REAL(rk), DIMENSION(:), INTENT(IN) :: c

      ALLOCATE(cent(size(c)))

      cent = c
      width = w

    END SUBROUTINE set_tophat



    ! Log Likelihood function for 2D Gaussian Data
    ! Use sufficient statistics x1 = mean, x2 = var 
    ! Pack in data arr, (n, x1, x2) repeat for each dim.
    REAL(rk) FUNCTION Lgaussian(data , par)
      ! The array of data
      REAL(rk), DIMENSION(:,:) :: data
      ! Parameters 
      REAL(rk), DIMENSION(:) :: par

      ! Array shape
      INTEGER, DIMENSION(2) :: sh
      INTEGER :: i,s
      REAL(rk) :: m1,s1,t,x1,x2,n

      sh = (SHAPE(data))
      s = sh(2)

      ! Check array dimensions.
      IF(size(par) /= 2*s .OR. sh(1) /= 3)  THEN
         PRINT *,"Arrays incorrect dimensions."
      END IF

      t = 0.0

      ! Iterate over dimensions.
      DO i=1,s

         m1 = par(2*i-1)
         s1 = par(2*i)
         n  = data(1,i)
         x1 = data(2,i)
         x2 = data(3,i)

         t = t - n * 0.5 *( LOG(2*pi*s1*s1) + (x2 - 2*x1*m1 + m1**2) / (s1**2))

      END DO

      Lgaussian = t

    END FUNCTION Lgaussian



    ! Take m n-dimensional points and generate the data array
    ! to be used by Lgaussian.
    SUBROUTINE datagaussian(data, ddata, out)
      ! Array of data points to process.
      REAL(rk), DIMENSION(:,:), INTENT(IN) :: data
      ! Dimensionality of data.
      INTEGER, INTENT(IN) :: ddata
      ! Array to output formatted statistics in.
      REAL(rk), DIMENSION(:,:), INTENT(OUT) :: out

      INTEGER :: i,j,n 
      REAL(rk) :: x1=0.0, x2=0.0

      ! Number of data points.
      n = size(data(1,:))

      ! Calculate statistics.
      DO i=1,ddata
         out(1,i) = n

         DO j=1,size(data(1,:))
            x1 = x1 + data(i,j)
            x2 = x2 + data(i,j)**2
         END DO

         ! Write out array values.
         out(2,i) = x1 / n
         out(3,i) = x2 / n

      END DO

    END SUBROUTINE datagaussian



    ! Print the data statistics from the gaussian array.
    SUBROUTINE dmean(data)
      
      REAL(rk), DIMENSION(:,:) :: data
      INTEGER :: i
      
      DO i=1,size(data(1,:))
         PRINT *,"Dim",i
         PRINT *,"Mean",data(2,i)
         PRINT *," ","Sigma",(data(1,i)*data(3,i) - (data(2,i))**2) / (data(1,i) - 1)
      END DO
      
    END SUBROUTINE dmean
    


    ! Set the prior parameters.
    SUBROUTINE set_prior_params(pars)

      REAL(rk), DIMENSION(:,:), POINTER :: pars

      prior_params => pars

    END SUBROUTINE set_prior_params



    ! An n-dimensional gaussian prior. Mean and std-dev
    ! from prior_params (1,2) respectively.
    REAL(rk) FUNCTION Xgaussian(par)
      ! Parameter array to check.
      REAL(rk), DIMENSION(:) :: par

      REAL(rk) :: p, t=0.0
      INTEGER :: i

      ! Check to see if dimensionailty consistent.
      IF(size(prior_params(1,:)) /= size(par) \
         .OR. size(prior_params(:,1)) /= 2) THEN
         
         PRINT *,"Incorrect array dimensions."
      END IF

      DO i=1,size(prior_params(1,:))
         p = par(i)

         ! Exponent
         t = t - 0.5 * ((p - prior_params(1,i)) / prior_params(2,i))**2
         ! Normalisation
         t = t - 0.5 * LOG(2*pi*prior_params(2,i)**2)

      END DO

      Xgaussian = t

    END FUNCTION Xgaussian



    ! A prior uniform in n-dimensional space between boundarys 
    ! set by set_prior_params.
    REAL(rk) FUNCTION Xuniform(par)

      REAL(rk), DIMENSION(:) :: par

      REAL(rk) :: p
      REAL(rk) :: v
      INTEGER :: i
      v = 1.0
      
      ! Check to see if dimensionailty consistent.
      IF(size(prior_params(1,:)) /= size(par) \
         .OR. size(prior_params(:,1)) /= 2) THEN
         
         PRINT *,"Incorrect array dimensions."
      END IF

      DO i=1,size(prior_params(1,:))
         p = par(i)

         ! Check to see if point in accepted volume.
         IF(p < prior_params(1,i) .OR. p > prior_params(2,i)) THEN
            Xuniform = min_inf
            RETURN
         END IF

         ! Calculate volume of prior space.
         v = v * (prior_params(2,i) - prior_params(1,i))

      END DO

      Xuniform = -1.0 * LOG(v)

    END FUNCTION Xuniform



    ! Uniform prior in the Unit sphere,
    REAL(rk) FUNCTION Xtophat(par)
      
      ! Parameters 
      REAL(rk), DIMENSION(:) :: par
      
      INTEGER ::i,n
      REAL(rk) :: r2
      r2 = 0.0
      n = size(par)
      
      ! Find n-d radius.
      DO i=1,n
         r2 = r2 + ((par(i) - cent(i))**2)
      END DO
      
      ! Check if in unit n-sphere.
      IF(r2 > 1.0) THEN
         Xtophat = min_inf
      ELSE
         Xtophat = log(fact(n/2) / (pi**(n/2) * width**n))
      END IF
      
    END FUNCTION Xtophat


    ! Factorial function.
    INTEGER FUNCTION fact(n) 
      
      INTEGER :: i,n
      fact = 1
      IF(n < 1) THEN 
         PRINT *,"n must be a positive integer."
      END IF
      
      DO i=1,n
         fact = fact * i
      END DO
      
    END FUNCTION fact
    

  END MODULE distr
