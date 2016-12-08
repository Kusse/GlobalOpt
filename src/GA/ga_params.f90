!
! Parameters for genetic algorithm searches
!
MODULE GA_PARAMS
   IMPLICIT NONE
   INTEGER, PARAMETER   :: MYGA_NSTRUC = 20    ! Population size
   INTEGER, PARAMETER   :: MYGA_NOFF = 1       ! Number of offspring
   INTEGER, SAVE        :: MYGA_NMUT = 0       ! Number of mutants. Re-sets every generation
   INTEGER, PARAMETER   :: MYGA_GENS = 20000   ! Number of generations
   INTEGER, SAVE        :: MAXFUNCALLS
   INTEGER, SAVE        :: MAXRUNS = 30
   INTEGER, PARAMETER   :: MYGA_TOURN_SIZE = 3
   LOGICAL, SAVE    :: MYGA_L_EPOCH=.FALSE.
   INTEGER,SAVE     :: MYGA_EPOCH_SAVE = 0
   LOGICAL,SAVE     :: MYGA_EPOCH_DUMP = .FALSE.
   INTEGER,SAVE     :: MYGA_COUNT_EPOCH = 1
   DOUBLE PRECISION,PARAMETER :: MYGA_MUT_RATE=0.3 ! Mutation rate
   INTEGER,PARAMETER :: MYGA_CROSS = 1 ! n-point crossover
   DOUBLE PRECISION,SAVE :: MYGA_BQMAX, MYGA_CQMAX !Loose & tight convergence thresholds
   DOUBLE PRECISION,SAVE :: MYGA_DUPLICATE_ETHRESH=0.0 !Energy threshold for duplicate predator
!  DOUBLE PRECISION,SAVE :: MYGA_DUPLICATE_GTHRESH=10*4.d0*DATAN(1.D0)/180d0 !Geometry threshold for duplicate predator
   DOUBLE PRECISION,PARAMETER :: MYGA_DUPLICATE_GTHRESH=0.174532925D0
   DOUBLE PRECISION,SAVE :: MYGA_EPOCH_THRESH=0.01D0
   DOUBLE PRECISION,SAVE :: MYGA_LAST_ENERGY=1D10
   LOGICAL,SAVE :: MYGA_L_ROUL=.FALSE.
   LOGICAL,SAVE :: MYGA_L_CHAIN=.FALSE.!Generate random structures as a continuous chain
   LOGICAL,SAVE :: MYGA_L_SPHERE=.FALSE.!
   DOUBLE PRECISION,SAVE :: MYGA_SPHERE_RADIUS= 3.0D0   !0.8*(55**(1.0D0/3.0D0)) !
   INTEGER,SAVE :: MYGA_COUNT_MIN = 0
   INTEGER,SAVE :: CURR_GEN=0
   DOUBLE PRECISION,SAVE :: MYGA_BH_INIT !Initial # of basin-hopping steps for each individual
   DOUBLE PRECISION,SAVE :: MYGA_BH_INCR=0 !Increment for basin-hopping steps after each generation
   DOUBLE PRECISION,SAVE :: MYGA_BH_STEPS !Current # of basin-hopping steps for each individual
   LOGICAL,SAVE :: MYGA_DUMP_POP=.FALSE. !Dump genomes of whole population to files after each generation
   LOGICAL,SAVE :: MYGA_TWIN=.FALSE. !Allow twinning moves (mating with two copies of one parent)
END MODULE GA_PARAMS
!
