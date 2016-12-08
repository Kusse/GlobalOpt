!mo361
!Modules for genetic algorithms

! Module to hold whole population of genetic algorithm
!
MODULE GA_POPULATION
   USE GA_PARAMS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE :: MYGA_POP_FOUND(:)
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_ENERGY(:)
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_GENOME(:,:) !3*natoms,popsize
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_COORDS(:,:) !3*natoms,popsize
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_FITNESS(:) !Fitness for roulette selection
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_MAXFORCE(:)
   INTEGER, ALLOCATABLE :: MYGA_POP_TYPE(:,:) !Atom types for multi-component clusters.
   INTEGER, ALLOCATABLE :: MYGA_POP_MOVE(:,:) !Attempted move record 3,popsize
   INTEGER, ALLOCATABLE :: MYGA_MOVE_HISTORY(:,:) !Move history 4,natoms (1=Attempted crossovers, 2=accepted crossovers, 3=attempted mutations, 4=accepted mutations)
END MODULE GA_POPULATION
