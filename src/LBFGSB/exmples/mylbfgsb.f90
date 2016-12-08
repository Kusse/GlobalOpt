
MODULE mylbfgsb
!    subroutine LBFGS_BC(func, dfunc, x, lb, ub, fx, gx, FEV, GEV)
!        use data_type
!        use config, only: LBFGS_Memory,LBFGSB_Factr, gradient_tolerance, &
!        function_tolerance, max_gradient_evaluations, max_iterations, iPrint
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP=KIND(1.0D0)
CONTAINS
    SUBROUTINE lbfgsbc(IND1, IND2, gtol)
!        use MYDATA_TYPE, only: DP, I4B
        USE COMMONS, ONLY: NDIM, LOWERBOUND, UPPERBOUND, LBFGS_MEMORY, IPRINT
        USE OBJECTIVE, ONLY: OBJ_FUNC
        USE GA_POPULATION, ONLY: MYGA_POP_ENERGY, MYGA_POP_COORDS, MYGA_POP_COORDS
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IND1, IND2
        REAL(KIND=DP), INTENT(IN) :: Gtol
        INTEGER :: NFG
        REAL(KIND=dP)     :: f
        REAL(KIND=dp), DIMENSION(NDIM) :: x
        REAL(KIND=dp), DIMENSION(NDIM) :: g
!        integer, parameter :: maxfun = NDIM
        INTEGER, DIMENSION(NDIM) :: NBOUNDS
        REAL(KIND=dp), ALLOCATABLE          :: wa(:), iwa(:)
        REAL(KIND=dp)          :: dsave(29)
        CHARACTER(LEN=60)      :: task, csave
        LOGICAL                :: lsave(4)
        INTEGER                :: isave(44)
        INTEGER                :: I, J, n, m
        REAL(KIND=DP), PARAMETER :: FACTR = 1.D+1
        REAL(KIND=DP), PARAMETER :: PGTOL = 1.D-12
!            real(kind=DP), parameter :: LBFGSB_LA_FACTR = 1.d+12
        n = NDIM
        m = LBFGS_MEMORY
        ALLOCATE(iwa(3*n), wa(2*m*n + 5*n + 11*m*m + 8*m))
        NBOUNDS = 2 ! Bounded both above and below
        NFG = 0
        DO I=IND1, IND2
            task = 'START';
            X = MYGA_POP_COORDS(:, I)
            DO WHILE(task(1:2).EQ.'FG'.OR.task.EQ.'NEW_X'.OR. &
                task.EQ.'START')
                CALL setulb(NDIM, M, X, LOWERBOUND, UPPERBOUND, NBOUNDS, &
                    f, g, factr, pgtol, wa, iwa, task, iPrint, csave, lsave, isave, dsave)
                IF (task(1:2) .EQ. 'FG') THEN
                    CALL OBJ_FUNC(X, F, G, NDIM)
                    NFG = NFG + 1
                END IF
                IF (task(1:5) .EQ. 'NEW_X') THEN
                    IF (dsave(13)  .LE. 1.0E-10*(1.0_DP + abs(f))) THEN
                        task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
                        EXIT
                    ENDIF
                    IF (isave(34) .GE. 3*NDIM) THEN
                        task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
                        EXIT
                    ENDIF
                    IF (task(1:4) .EQ. 'STOP') THEN
                        WRITE (6,*) task
                        WRITE (6,*) 'Final X='
                        WRITE (6,'((1x,1p, 6(1x,d11.4)))') (x(J),J = 1,n)
                        EXIT
                    ENDIF
                    IF(MAXVAL(ABS(G)) < GTOL) EXIT
                ENDIF
            ENDDO
            MYGA_POP_COORDS(:, I) = X
            MYGA_POP_ENERGY(I) = F
        ENDDO
        DEALLOCATE (iwa, wa)
    END  SUBROUTINE lbfgsbc

END MODULE mylbfgsb















!                                                                                      
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
!  or “3-clause license”)                                                              
!  Please read attached file License.txt                                               
!                                        
!
!                      DRIVER1 in Fortran 90
!     --------------------------------------------------------------
!
!        L-BFGS-B is a code for solving large nonlinear optimization
!             problems with simple bounds on the variables.
!
!        The code can also be used for unconstrained problems and is
!        as efficient for these problems as the earlier limited memory
!                          code L-BFGS.
!
!        This is the simplest driver in the package. It uses all the
!                    default settings of the code.
!
!
!     References:
!
!        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!        memory algorithm for bound constrained optimization'',
!        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!        Subroutines for Large Scale Bound Constrained Optimization''
!        Tech. Report, NAM-11, EECS Department, Northwestern University,
!        1994.
!
!
!          (Postscript files of these papers are available via anonymous
!           ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                              *  *  *
!
!         March 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales
!
!     --------------------------------------------------------------
!             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
!     --------------------------------------------------------------
!
!     n is an INTEGER variable that must be set by the user to the
!       number of variables.  It is not altered by the routine.
!
!     m is an INTEGER variable that must be set by the user to the
!       number of corrections used in the limited memory matrix.
!       It is not altered by the routine.  Values of m < 3  are
!       not recommended, and large values of m can result in excessive
!       computing time. The range  3 <= m <= 20 is recommended. 
!
!     x is a DOUBLE PRECISION array of length n.  On initial entry
!       it must be set by the user to the values of the initial
!       estimate of the solution vector.  Upon successful exit, it
!       contains the values of the variables at the best point
!       found (usually an approximate solution).
!
!     l is a DOUBLE PRECISION array of length n that must be set by
!       the user to the values of the lower bounds on the variables. If
!       the i-th variable has no lower bound, l(i) need not be defined.
!
!     u is a DOUBLE PRECISION array of length n that must be set by
!       the user to the values of the upper bounds on the variables. If
!       the i-th variable has no upper bound, u(i) need not be defined.
!
!     nbd is an INTEGER array of dimension n that must be set by the
!       user to the type of bounds imposed on the variables:
!       nbd(i)=0 if x(i) is unbounded,
!              1 if x(i) has only a lower bound,
!              2 if x(i) has both lower and upper bounds, 
!              3 if x(i) has only an upper bound.
!
!     f is a DOUBLE PRECISION variable.  If the routine setulb returns
!       with task(1:2)= 'FG', then f must be set by the user to
!       contain the value of the function at the point x.
!
!     g is a DOUBLE PRECISION array of length n.  If the routine setulb
!       returns with taskb(1:2)= 'FG', then g must be set by the user to
!       contain the components of the gradient at the point x.
!
!     factr is a DOUBLE PRECISION variable that must be set by the user.
!       It is a tolerance in the termination test for the algorithm.
!       The iteration will stop when
!
!        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!       where epsmch is the machine precision which is automatically
!       generated by the code. Typical values for factr on a computer
!       with 15 digits of accuracy in double precision are:
!       factr=1.d+12 for low accuracy;
!             1.d+7  for moderate accuracy; 
!             1.d+1  for extremely high accuracy.
!       The user can suppress this termination test by setting factr=0.
!
!     pgtol is a double precision variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       The user can suppress this termination test by setting pgtol=0.
!
!     wa is a DOUBLE PRECISION  array of length 
!       (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace.
!       This array must not be altered by the user.
!
!     iwa is an INTEGER  array of length 3nmax used as
!       workspace. This array must not be altered by the user.
!
!     task is a CHARACTER string of length 60.
!       On first entry, it must be set to 'START'.
!       On a return with task(1:2)='FG', the user must evaluate the
!         function f and gradient g at the returned value of x.
!       On a return with task(1:5)='NEW_X', an iteration of the
!         algorithm has concluded, and f and g contain f(x) and g(x)
!         respectively.  The user can decide whether to continue or stop
!         the iteration. 
!       When
!         task(1:4)='CONV', the termination test in L-BFGS-B has been 
!           satisfied;
!         task(1:4)='ABNO', the routine has terminated abnormally
!           without being able to satisfy the termination conditions,
!           x contains the best approximation found,
!           f and g contain f(x) and g(x) respectively;
!         task(1:5)='ERROR', the routine has detected an error in the
!           input parameters;
!       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
!         contains additional information that the user can print.
!       This array should not be altered unless the user wants to
!          stop the run for some reason.  See driver2 or driver3
!          for a detailed explanation on how to stop the run 
!          by assigning task(1:4)='STOP' in the driver.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     csave  is a CHARACTER working array of length 60.
!
!     lsave is a LOGICAL working array of dimension 4.
!       On exit with task = 'NEW_X', the following information is
!         available:
!       lsave(1) = .true.  the initial x did not satisfy the bounds;
!       lsave(2) = .true.  the problem contains bounds;
!       lsave(3) = .true.  each variable has upper and lower bounds.
!
!     isave is an INTEGER working array of dimension 44.
!       On exit with task = 'NEW_X', it contains information that
!       the user may want to access:
!         isave(30) = the current iteration number;
!         isave(34) = the total number of function and gradient
!                         evaluations;
!         isave(36) = the number of function value or gradient
!                                  evaluations in the current iteration;
!         isave(38) = the number of free variables in the current
!                         iteration;
!         isave(39) = the number of active constraints at the current
!                         iteration;
!
!         see the subroutine setulb.f for a description of other 
!         information contained in isave
!
!     dsave is a DOUBLE PRECISION working array of dimension 29.
!       On exit with task = 'NEW_X', it contains information that
!         the user may want to access:
!         dsave(2)  = the value of f at the previous iteration;
!         dsave(5)  = the machine precision epsmch generated by the code;
!         dsave(13) = the infinity norm of the projected gradient;
!
!         see the subroutine setulb.f for a description of other 
!         information contained in dsave
!
!     --------------------------------------------------------------
!           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
!     --------------------------------------------------------------
!    
!     program driver
!DULE LBFGS_B
!    This simple driver demonstrates how to call the L-BFGS-B code to
!      solve a sample problem (the extended Rosenbrock function
!      subject to bounds on the variables). The dimension n of this
!      problem is variable.
!    implicit none
!    Declare variables and parameters needed by the code.
!      Note thar we wish to have output at every iteration.
!         iprint=1
!      We also specify the tolerances in the stopping criteria.
!         factr  = 1.0d+7, pgtol  = 1.0d-5
!      A description of all these variables is given at the beginning
!      of  the driver
!  CONTAINS
!  SUBROUTINE GA_LBFGSB_RUN(IND, factr, pgtol)
!      USE MYDATA_TYPE, ONLY: DP, I4B
!      USE COMMONS, ONLY: NDIM, LOWERBOUND, UPPERBOUND, LBFGSB_MEMORY, IPRINT
!      USE OBJECTIVEFUNCTION, ONLY: DFUNC, FUNC
!      USE GA_POPULATION, ONLY: POP_ENERGY, POP_FORCES, POP_COORDS
!      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: IND
!      REAL(KIND=DP), INTENT(IN) :: FACTR, pgtol
!     integer,  parameter    :: n = 25, m = 5, iprint = 1
!     integer,  parameter    :: dp = kind(1.0d0)
!     real(dp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-5
!    character(len=60)      :: task, csave
!    logical                :: lsave(4)
!    integer                :: isave(44)
!    real(dp)               :: f
!    real(dp)               :: dsave(29)
!    integer,  allocatable  :: iwa(:), wa(:)   ! nbd(:),
!     real(dp), allocatable  :: x(:), l(:), u(:), g(:), wa(:)
!    REAL(KIND=DP), DIMENSION(NDIM) :: X, G, NBD
!    Declare a few additional variables for this sample problem
!    real(dp)               :: t1, t2
!    integer                :: i, m, n
!    Allocate dynamic arrays
!    N = NDIM
!    M = LBFGSB_MEMORY
!     LB = LOWER
!     allocate ( nbd(NDIM), x(NDIM), l(n), u(n), g(n) )
!    allocate ( iwa(3*NDIM))
!    allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
!     do 10 i=1, n, 2
!        nbd(i) = 2
!        l(i)   = 1.0d0
!        u(i)   = 1.0d2
! 10  continue
!    Next set bounds AND DEFINE STARTING POINT
!  DO I=1, NDIM
!      NBD(I) = 2
!      X(I) = POP_COORDS(I, IND)
!      G(I) = 1.0D10
!  ENDDO
!rint*, "nbd = ", nbd
!top
!     do 12 i=2, n, 2
!        nbd(i) =  2
!        l(i)   = -1.0d2
!        u(i)   =  1.0d2
! 12   continue
!    We now define the starting point.
!     do 14 i=1, n
!        x(i) = 3.0d0
! 14  continue
!    write (6,16)
!16  format(/,5x, 'Solving sample problem.', &
!           /,5x, ' (f = 0.0 at the optimal solution.)',/)
!    We start the iteration by initializing task.
!    task = 'START'
!    The beginning of the loop
!    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
!             task.eq.'START')
!    This is the call to the L-BFGS-B code.
!       call setulb (NDIM, lbfgsb_memory, x, LOWERBOUND, UPPERBOUND, nbd, &
!          f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
!       if (task(1:2) .eq. 'FG') then
!          call func(x, f)
!           f=.25d0*( x(1)-1.d0 )**2
!           do 20 i=2, n
!              f = f + ( x(i)-x(i-1 )**2 )**2
!20         continue
!           f = 4.d0*f
!    Compute gradient g for the sample problem.
!          call dfunc(x, g)
!          POP_FORCES(:, IND) = G
!          POP_ENERGY(IND) = F
!          POP_COORDS(:, IND) = X
!           t1   = x(2) - x(1)**2
!           g(1) = 2.d0*(x(1) - 1.d0) - 1.6d1*x(1)*t1
!           do 22 i=2, n-1
!              t2   = t1
!              t1   = x(i+1) - x(i)**2
!              g(i) = 8.d0*t2 - 1.6d1*x(i)*t1
!22         continue
!           g(n) = 8.d0*t1
!       end if
!       if (maxval(abs(g)) < 1.0e-3) return
!    end do
!END SUBROUTINE GA_LBFGSB_RUN
!    end of loop do while
!D MODULE LBFGS_B
!     end program driver
