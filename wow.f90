PROGRAM wow

implicit none

 REAL(8), PARAMETER :: Tmin = 0.000833333D0 ! min signal duration =72 seconds
 REAL(8), PARAMETER :: Tmax = 0.36787944117144233D0 ! max signal duration =e^-1 days
 REAL(8), PARAMETER :: Lmin = 0.01831563888873418D0 ! min repeat rate =e^-4 days
 REAL(8), PARAMETER :: Lmax = 7.38905609893065D0 ! max repeat rate =e^2
 INTEGER, PARAMETER :: kmax = 100

 REAL(8) :: T_draw, L_draw
 INTEGER :: obscounts, k1, k2, sumcounts, n, ngaveup
 REAL(8), PARAMETER :: precision2_threshold = 10.0D0
 INTEGER, PARAMETER :: maxcounts = 1000 ! maximum # of counts for iterationcount subroutine
 REAL(8) :: precision2
 
 ! Open up jammi.log
 open(2,FILE='output.log',FORM='FORMATTED',ACCESS='APPEND')
 
 DO k1=1,150
  DO k2=1,160
	 write(*,*) ' '
	 write(*,*) '<===> Realization k1,k2 = ',k1,k2,' <===>'
	 T_draw = DEXP( DLOG(Tmin) + (k1-1)*(DLOG(Tmax)-DLOG(Tmin))/DBLE(kmax-1) )
	 !L_draw = loguniform_draw(Lmin,1.0D0/T_draw)
	 L_draw = DEXP( DLOG(Lmin) + (k2-1)*(DLOG(Lmax)-DLOG(Lmin))/DBLE(kmax-1) )
	 !write(*,*) '[T,L,1/L] = ',T_draw,L_draw,1.0D0/L_draw
	 
	 sumcounts = 0 ! sum of trials needed
	 n = 0 ! counter for number of loops run thus far
	 ngaveup = 0
	 precision2 = 0.0
	 
	 ! start a loop which continues until we've measured p to sufficient precision
	 IF( DLOG(L_draw) .LE. (-3.057394436628309 - 1.2724645048999665*DLOG(T_draw)) ) THEN
		 ! the above if statements checks for cases which are found from preliminary runs to
		 ! a) always fail b) take an inordinate amount of time to compute
	   DO WHILE( precision2 .LT. precision2_threshold )
		   n = n + 1
       call iterationcount(maxcounts,T_draw,L_draw,obscounts)
		   IF( obscounts .EQ. maxcounts ) THEN ! => the iterationsubroutine gave up
			   ngaveup = ngaveup + 1
		   END IF
		   sumcounts = sumcounts + obscounts
		   precision2 = precision_squared(n,sumcounts)
		   !write(*,*) 'Loop #',n,' p = ',p_mode(n,sumcounts),' precision^2 = ',precision2
	   END DO
	 ELSE
		 ngaveup = 11
		 n = 11
		 sumcounts = 11000
		 precision2 = precision_squared(n,sumcounts)
	 END IF
	 write(*,*) 'T_draw, L_draw = ',T_draw,L_draw,n,sumcounts
	 write(*,*) 'p = ',p_mode(n,sumcounts),"+/-",p_mode(n,sumcounts)/precision2
	 write(*,*) 'failures = ',ngaveup
   write(2,*) T_draw,L_draw,p_mode(n,sumcounts),precision2,ngaveup!,DBLE(ngaveup/n)
  END DO
 END DO
 
 ! =======================================================
 CONTAINS
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE iterationcount(max_counts,Tdraw,Ldraw,final_counts)

 implicit none
 
 ! Inputs
 REAL(8), INTENT(IN) :: Tdraw, Ldraw
 INTEGER, INTENT(IN) :: max_counts ! limit to counter
 
 ! Intermediates
 INTEGER :: counts, oktostop
 INTEGER :: max_len, nsignals
 REAL(8), DIMENSION(:), ALLOCATABLE :: t_signals
 INTEGER :: i, d, hits
 LOGICAL :: successA, successB
 
 ! Output
 INTEGER, INTENT(OUT) :: final_counts
 
 ! Fixed parameters
 REAL(8), PARAMETER :: baseline = 2673.0D0 ! baseline ~7 years
 REAL(8), PARAMETER :: tgap = 0.001875D0 ! time between horns, =72+90 seconds
 REAL(8), PARAMETER :: W = 0.000833333333333333D0 ! obs window
 INTEGER, PARAMETER :: ndays = 90 ! 90 good observation days
 
 ! Observation days
 REAL(8), DIMENSION(ndays), PARAMETER :: days = (/ 0.5, &
 1.5, 2.5, 4.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, &
 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 52.5, 78.5, 240.5, &
 241.5, 247.5, 254.5, 255.5, 256.5, 257.5, 258.5, 261.5, 262.5, &
 380.5, 381.5, 382.5, 383.5, 384.5, 385.5, 386.5, 387.5, 1968.5, &
 1969.5, 1970.5, 1973.5, 1979.5, 1981.5, 1988.25, 1988.75, 1989.5, &
 1990.5, 1991.5, 1994.5, 1995.5, 1999.25, 1999.75, 2001.5, 2002.5, &
 2003.5, 2004.5, 2009.5, 2011.5, 2012.5, 2013.5, 2014.5, 2015.5, &
 2016.5, 2017.5, 2018.5, 2019.5, 2020.5, 2021.5, 2022.5, 2023.5, &
 2024.5, 2025.5, 2026.5, 2032.5, 2033.5, 2034.5, 2037.5, 2038.5, &
 2039.5, 2045.5, 2046.5, 2048.5, 2051.5, 2053.5, 2054.5, 2082.5, &
 2095.5, 2672.5 /)
 
 ! initial conditions
 counts = 0
 oktostop = 0
 max_len = INT(DBLE(baseline*L_draw) + 10.0D0*DSQRT(baseline*L_draw) )
 ALLOCATE (t_signals(max_len))
 
 ! dowhile loop
 DO WHILE( oktostop .EQ. 0 ) ! stop condition is exactly 1 hit found, or counts >= max_counts
	 counts = counts + 1 ! increase the counter
	 
	 ! generate a poisson process
	 call poisson_process(L_draw,baseline,max_len,t_signals,nsignals)
	 !write(*,*) 'Over B0, we observed ',actual_len,', expected ',INT(B0*L_draw),' and capped at ',max_len
	 
	 ! test how many signals fall into an observing window
	 hits = 0
	 DO d=1,ndays ! go through each day
		 ! observation A
		 successA = .FALSE.
		 DO i=1,nsignals
		   IF( (days(d)-0.5*W).GT.(t_signals(i)-0.5*T_draw) ) THEN
				 IF( (days(d)+0.5*W).LT.(t_signals(i)+0.5*T_draw) ) THEN
				   successA = .TRUE.
				 END IF
			 END IF
		 END DO
		 ! was observation A a success?
		 IF( successA ) THEN
			 hits = hits + 1
		 END IF
		 ! observation B
		 successB = .FALSE.
		 DO i=1,nsignals
		   IF( (days(d)+tgap-0.5*W).GT.(t_signals(i)-0.5*T_draw) ) THEN
				 IF( (days(d)+tgap+0.5*W).LT.(t_signals(i)+0.5*T_draw) ) THEN
				   successB = .TRUE.
				 END IF
			 END IF
		 END DO
		 ! was observation B a success?
		 IF( successB ) THEN
			 hits = hits + 1
		 END IF
	 END DO
	 !write(*,*) '#',z,' finds ',hits,' hit days'
	 IF( hits .EQ. 1) THEN
		 write(*,*) 'We did it! hits = ',hits,' for iteration ',counts
	 END IF
	 IF( hits .EQ. 1 .OR. counts .GE. max_counts ) THEN
		 oktostop = 1 ! stop condition satisfied
	 END IF
 END DO
 final_counts = counts
 !write(*,*) 'Probability of this realization = 1/',z

 END SUBROUTINE iterationcount
 ! =======================================================
 
 ! =======================================================
 FUNCTION uniform_draw(xmin,xmax)
	
 REAL(8) :: uniform_draw
 REAL(8) :: xmin, xmax
 
 uniform_draw = xmin + (xmax-xmin)*random_uniform()
 
 RETURN
 
 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION loguniform_draw(xmin,xmax)
	
 REAL(8) :: loguniform_draw
 REAL(8) :: xmin, xmax
 
 loguniform_draw = DEXP( DLOG(xmin) + (DLOG(xmax)-DLOG(xmin))*random_uniform() )
 
 RETURN
 
 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION exponential_draw(lambda)
	
 REAL(8) :: exponential_draw
 REAL(8) :: lambda
 
 exponential_draw = -(1.0D0/lambda)*DLOG( 1.0D0 - random_uniform() )
 
 RETURN
 
 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION random_uniform()

 REAL(8) :: random_uniform
 REAL(8) :: seeda

 call random_seed()
 call random_number(seeda)
 random_uniform=seeda

 RETURN

 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION precision_squared(m,s)

 REAL(8) :: precision_squared
 INTEGER :: m, s

 precision_squared = m**2*(2.0D0+m+s)**2*(3.0D0+m+s)
 precision_squared = precision_squared/( (m+s)**2*(1.0D0+m)*(1.0D0+s) )

 RETURN

 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION p_mode(m,s)

 REAL(8) :: p_mode
 INTEGER :: m, s

 p_mode = DBLE(m)/DBLE(m+s)

 RETURN

 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE poisson_process(lambda,tmax,lentimes,times,obslen)

 implicit none
 
 REAL(8), INTENT(IN) :: lambda, tmax
 INTEGER, INTENT(IN) :: lentimes
 REAL(8), DIMENSION(lentimes), INTENT(OUT) :: times
 REAL(8) :: tlast
 INTEGER :: j
 INTEGER, INTENT(OUT):: obslen

 tlast = 0.0D0
 DO WHILE( tlast .LT. tmax )
	 j = 0
	 tlast = 0.0D0
   DO WHILE( tlast .LT. tmax .AND. j .LT. lentimes )
	   j = j + 1
	   times(j) = tlast + exponential_draw(lambda)
	   tlast = times(j)
   END DO
   obslen = j
 END DO
 
 END SUBROUTINE poisson_process
 ! =======================================================

END PROGRAM wow
