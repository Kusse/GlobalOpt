MODULE commons
    IMPLICIT NONE
    INTEGER, PARAMETER :: NATOMS = 309 
    INTEGER, PARAMETER :: NDIM = 3*NATOMS
    INTEGER, SAVE :: MYUNIT = 16
    INTEGER, PARAMETER :: NTYPEA = NATOMS       !NATOMS
    INTEGER, PARAMETER :: GATOM  = 23           ! ATOM TYPE
    INTEGER, SAVE :: NPCALL=0                   ! The number of potential calls
    INTEGER, SAVE :: LBFGS_FUNCALLS=0
    INTEGER, SAVE :: TOTAL_FUNCALLS=0
    INTEGER, SAVE :: MCSTEPS(3) = 1
    LOGICAL, SAVE :: HIT
    LOGICAL, SAVE :: BLJCLUSTER = .FALSE.
    LOGICAL, SAVE :: BGUPTAT = .FALSE.
    LOGICAL, SAVE :: BLNT = .FALSE.
    DOUBLE PRECISION, DIMENSION(NATOMS),SAVE :: ATMASS = 1.0D0
    DOUBLE PRECISION, DIMENSION(3),SAVE :: STEP = (/0.1d0, 0.1d0, 0.1d0 /)
    DOUBLE PRECISION, DIMENSION(3,NATOMS), SAVE :: COORDS
    DOUBLE PRECISION, DIMENSION(3,NATOMS), SAVE :: COORDSO
    DOUBLE PRECISION,SAVE :: TSTART
    DOUBLE PRECISION,SAVE :: BQMAX = 1.0E-3
    DOUBLE PRECISION,SAVE :: CQMAX = 1.0E-5
    DOUBLE PRECISION, SAVE :: LBFGS_MEMORY = 10
!    DOUBLE PRECISION, SAVE :: SPHERE_RADIUS = 0.4D0*(NATOMS**(1.0D0/3.0D0))
    DOUBLE PRECISION, DIMENSION(NDIM), SAVE ::  LOWERBOUND = -0.80D0*(NATOMS**(1.0D0/3.0D0))  !-2.0D0*SPHERE_RADIUS
    DOUBLE PRECISION, DIMENSION(NDIM), SAVE ::  UPPERBOUND =  0.80D0*(NATOMS**(1.0D0/3.0D0))  ! 2.0D0*SPHERE_RADIUS
    INTEGER, SAVE :: IPRINT = -10
    CHARACTER(LEN=1), ALLOCATABLE,SAVE :: BEADLETTER(:)
    INTEGER, DIMENSION(:), ALLOCATABLE,SAVE :: RANDOMSEED
!    CHARACTER(LEN=256),SAVE :: CMD_LINE

!CONTAINS
!    SUBROUTINE DEFAULT_VALUES()
!        IMPLICIT NONE
!        NDIM = 3*NATOMS
!        NTYPEA = NATOMS
!        GATOM = 23
!        ALLOCATE(COORDS(3, NATOMS), COORDSO(3, NATOMS), LOWERBOUND(NDIM), UPPERBOUND(NDIM))
!        LOWERBOUND = -2.0D0*(NATOMS**(1.0D0/3.0D0))
!        UPPERBOUND = 2.0D0*(NATOMS**(1.0D0/3.0D0))
!    END SUBROUTINE DEFAULT_VALUES

END MODULE commons