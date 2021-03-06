!mo361
!Core genetic algorithm routines for GMIN
!
!Main routine for GA

SUBROUTINE MYGA_RUN()
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY : MYUNIT,MCSTEPS,TSTART,HIT, LBFGS_FUNCALLS
   USE MYLBFGSB, ONLY: LBFGSBC
   IMPLICIT NONE
   INTEGER,SAVE:: I,J, RUN=0
   INTEGER CURR_STRUC, MY_UNIT1, MY_UNIT2
   INTEGER:: REMOVED=0 !Count of structures removed by predator in current generation
   INTEGER GOOD_MUT,GOOD_CROSS !Count of new structures getting into next generation
   DOUBLE PRECISION DPRAND, AVG, STD
   DOUBLE PRECISION TIME
   DOUBLE PRECISION MYGA_MEAN_ENERGY
   CHARACTER(LEN=50) :: FILENAME
   RUN = RUN+1
   LBFGS_FUNCALLS = 0
   PRINT*, "SIMULATION STARTED."
   WRITE(FILENAME, "(A,I0,A)") "Emin_", RUN,".txt"
   OPEN(NEWUNIT=MY_UNIT1, FILE=ADJUSTL(TRIM(FILENAME)), ACTION="WRITE", STATUS="REPLACE")
   WRITE(FILENAME, "(A,I0,A)") "Xmin_",RUN,".txt"
   OPEN(NEWUNIT=MY_UNIT2, FILE=ADJUSTL(TRIM(FILENAME)), ACTION="WRITE", STATUS="REPLACE")
   CALL MYGA_SETUP()
   !Minimise initial population
   PRINT*, "INITIAL OPTIMIZATION."
   CALL LBFGSBC(1, MYGA_NSTRUC, 1.0D-2)

!   CALL MYGA_OPTIMISE_POP(1,MYGA_NSTRUC,0,MYGA_BQMAX)
   MYGA_COUNT_MIN=MYGA_COUNT_MIN+MYGA_NSTRUC
   !Sort initial population
   PRINT*, "SORTING INITIAL POPULATION."
   CALL MYGA_SORT(MYGA_NSTRUC)
   IF (MYGA_DUMP_POP) CALL MYGA_DUMP_IC()
   ! Main GA loop starts here
!   PRINT*, "MAIN LOOP.", MYGA_POP_ENERGY(1), MYGA_POP_ENERGY(MYGA_NSTRUC)
!   STOP
   WRITE(MY_UNIT1, "(7(A10,2X),/)") "GENERATION", "FUNCALLS", "MIN_ENERGY", &
	 "MAX_ENERGY", "AVG_ENERGY", "STD_ENERGY", "MAXFORCE"
   WRITE(MY_UNIT2, "(4(A10,2X))") "GENERATION", "FUNCALLS", "ENERGY", "COORDS"
   WRITE(MY_UNIT2, *) 0, LBFGS_FUNCALLS, MYGA_POP_ENERGY(1), MYGA_POP_COORDS(:,1)
   DO CURR_GEN = 1,MYGA_GENS
!	DO I=1, MYGA_NSTRUC
	  AVG = MYGA_MEAN_ENERGY()
	  STD = SQRT(SUM((MYGA_POP_ENERGY(1:MYGA_NSTRUC) - AVG)**2)/DBLE(MYGA_NSTRUC))
          WRITE(MY_UNIT1, *) CURR_GEN-1, LBFGS_FUNCALLS, MYGA_POP_ENERGY(1), &
		MYGA_POP_ENERGY(MYGA_NSTRUC), AVG, STD, MYGA_POP_MAXFORCE(1)
!	DO I=1, MYGA_NSTRUC
!     	ENDDO
        print*, "gen = ", CURR_GEN, LBFGS_FUNCALLS, myga_pop_energy(1)
!      IF (MYGA_MEAN_ENERGY() < -52.0D0) THEN
!        CALL ANALYZE(1, MYGA_NSTRUC, CURR_GEN)
!      ENDIF
      !Check for stagnation and start new epoch if needed
      IF (MYGA_L_EPOCH.AND.(MYGA_LAST_ENERGY-MYGA_MEAN_ENERGY()).LT.MYGA_EPOCH_THRESH) THEN
         CALL MYGA_NEW_EPOCH()
         MYGA_COUNT_MIN=MYGA_COUNT_MIN+MYGA_NSTRUC-MYGA_EPOCH_SAVE
      ENDIF
      MYGA_LAST_ENERGY=MYGA_MEAN_ENERGY()
      ! Mark all structures as old
      DO CURR_STRUC=1,MYGA_NSTRUC
         MYGA_POP_MOVE(1,CURR_STRUC)=0
      ENDDO
      IF (MYGA_L_ROUL)THEN !Generate wheel for roulette selection
         CALL MYGA_FIT_TANH()
      ENDIF
      !Make offspring
      DO CURR_STRUC= MYGA_NSTRUC+1,MYGA_NSTRUC+MYGA_NOFF
         CALL MYGA_SELECT(I,J) !Select parents
         CALL MYGA_MATE(I,J,CURR_STRUC)
         MYGA_POP_FOUND(CURR_STRUC)=CURR_GEN
      ENDDO
      ! Optimise all offpring using loose convergence criterion
      CALL LBFGSBC(MYGA_NSTRUC+1, MYGA_NSTRUC+MYGA_NOFF, 1.0D-2)

!      CALL MYGA_OPTIMISE_POP(MYGA_NSTRUC+1,MYGA_NSTRUC+MYGA_NOFF,MCSTEPS(1),MYGA_BQMAX)
      !Make mutants
      MYGA_NMUT=0
      DO I=1,MYGA_NSTRUC+MYGA_NOFF
         IF (DPRAND().LT.MYGA_MUT_RATE) THEN
            MYGA_NMUT=MYGA_NMUT+1
            CALL MYGA_MUTATE(I,MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT)
         ENDIF
      ENDDO
      MYGA_COUNT_MIN=MYGA_COUNT_MIN+(MYGA_NOFF+MYGA_NMUT)*(MCSTEPS(1)+1)
      ! Optimise all mutants using loose convergence criterion
      CALL LBFGSBC(MYGA_NSTRUC+MYGA_NOFF+1,MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT, 1.0D-2)

!      CALL MYGA_OPTIMISE_POP(MYGA_NSTRUC+MYGA_NOFF+1,MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT,MCSTEPS(1),MYGA_BQMAX)
      !All structures generated
      !Sort population
      CALL MYGA_SORT(MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT)
      IF (MYGA_DUPLICATE_ETHRESH.GT.0) THEN
         CALL MYGA_REMOVE_DUPLICATES(MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT,REMOVED)
      ENDIF
      ! Optimise best structure with tight convergence criterion
      CALL LBFGSBC(1,1,1.0D-3)
!      CALL MYGA_OPTIMISE_POP(1,1,0,MYGA_CQMAX)
      IF (MYGA_POP_ENERGY(1).LT.MYGA_POP_ENERGY(0)) THEN
         CALL COPY_STRUCTURE(1,0)
      ENDIF
      ! Increment number of basin-hopping steps for BHGA hybrid
      MYGA_BH_STEPS=MYGA_BH_STEPS+MYGA_BH_INCR
      MCSTEPS(1)=INT(MYGA_BH_STEPS)
      ! Update move history
      GOOD_MUT=0
      GOOD_CROSS=0
      DO CURR_STRUC = 1,MYGA_NSTRUC
         IF(MYGA_POP_MOVE(1,CURR_STRUC).EQ.1) THEN !Offspring structure in new population
            GOOD_CROSS=GOOD_CROSS+1
            DO I=1,MYGA_CROSS
               MYGA_MOVE_HISTORY(2,MYGA_POP_MOVE(I+1,CURR_STRUC))=MYGA_MOVE_HISTORY(2,MYGA_POP_MOVE(I+1,CURR_STRUC))+1
            ENDDO
         ELSE IF(MYGA_POP_MOVE(1,CURR_STRUC).EQ.2) THEN !Mutant structure in new population
            GOOD_MUT=GOOD_MUT+1
            MYGA_MOVE_HISTORY(4,MYGA_POP_MOVE(2,CURR_STRUC))=MYGA_MOVE_HISTORY(4,MYGA_POP_MOVE(2,CURR_STRUC))+1
         ENDIF
      ENDDO
      ! Print out summary of last generation
      CALL MYCPU_TIME(TIME)
      WRITE(MYUNIT,'(A,I5,4(A,F15.6))') &
      "GA> Generation=",CURR_GEN,&
      " Best=",MYGA_POP_ENERGY(1),&
      " Mean=",MYGA_MEAN_ENERGY(),&
      " Worst=",MYGA_POP_ENERGY(MYGA_NSTRUC),&
      " Time=",TIME-TSTART
      WRITE(MYUNIT,'(A,I5,4(A,F15.6))') &
      "GA> Generation=",CURR_GEN,&
      " Duplicates=",DBLE(REMOVED)/DBLE(MYGA_NSTRUC+MYGA_NOFF),&
      " Offspring in new pop.= ",DBLE(GOOD_CROSS)/DBLE(MYGA_NOFF),&
      " Mutants in new pop.= ",DBLE(GOOD_MUT)/DBLE(MYGA_NMUT)
      IF (MYGA_DUMP_POP) CALL MYGA_DUMP_IC()
      !WRITE(*,*) "Generation",CURR_GEN
      !CALL MYGA_PRINT_POP()
      !Terminate search if we've found the target structure
      IF (HIT) THEN
         WRITE(MYUNIT,'(A,I5)')"Target hit at generation",CURR_GEN
         GOTO 10
      ENDIF
      DO I=1, MYGA_NSTRUC
         IF(MYGA_POP_ENERGY(I) < -52) THEN
           WRITE(MY_UNIT2, *) CURR_GEN-1, &
           LBFGS_FUNCALLS, MYGA_POP_ENERGY(I), MYGA_POP_COORDS(:,I)
         ELSE
           EXIT
         ENDIF
      ENDDO
   ENDDO
   CURR_GEN=CURR_GEN-1 !Correct generation count when search runs out of steps
10  CALL lbfgsbc(1,MYGA_NOFF,1.0D-8)
   CALL lbfgsbc(1,MYGA_NOFF,1.0D-8)
!   DO I=1, MYGA_NSTRUC
!	IF(MYGA_POP_ENERGY(I) < -52) THEN
!          WRITE(MY_UNIT2, *) CURR_GEN-1, &
!	  LBFGS_FUNCALLS, MYGA_POP_ENERGY(I), MYGA_POP_COORDS(:,I)
!	ELSE
!	  EXIT
!	ENDIF
!  ENDDO
CLOSE(MY_UNIT1)
CLOSE(MY_UNIT2)
! GA finished
!10  CALL lbfgsbc(1,MYGA_NOFF,1.0D-8)
!CALL MYGA_OPTIMISE_POP(1,MYGA_NOFF,0,MYGA_CQMAX)
   CALL MYGA_FINALIO()
   PRINT*, "FMIN = ", MYGA_POP_ENERGY(1)
   RETURN
END
!
!Set up population for GA.
!
SUBROUTINE MYGA_SETUP()
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS
   !, ONLY : NATOMS,MYUNIT,BQMAX,CQMAX,BLNT,MCSTEPS,BLJCLUSTER,BGUPTAT,STEP
   IMPLICIT NONE
   INTEGER I
   INTEGER POPSIZE
!   MYGA_SPHERE_RADIUS  = 1.0D0*(NATOMS**(1.0D0/3.0D0)) ! BY KUSSE
   ! Make arrays large enough to store population + offspring + mutants + best structures
   POPSIZE=MYGA_NSTRUC*3+MYGA_NOFF*2
   !Population arrays start at -2 to hold temporary structures.
!   PRINT*, "POPSIZE = ", POPSIZE, NATOMS
   IF(.NOT. ALLOCATED(MYGA_POP_ENERGY)) ALLOCATE(MYGA_POP_ENERGY(-2:POPSIZE))
   MYGA_POP_ENERGY=1D10
   IF(.NOT. ALLOCATED(MYGA_POP_FITNESS)) ALLOCATE(MYGA_POP_FITNESS(-2:POPSIZE))
   IF(.NOT. ALLOCATED(MYGA_POP_MAXFORCE)) ALLOCATE(MYGA_POP_MAXFORCE(-2:POPSIZE))
   IF(.NOT. ALLOCATED(MYGA_POP_FOUND))	ALLOCATE(MYGA_POP_FOUND(-2:POPSIZE))
   MYGA_POP_FOUND=0
   ! Attempted moves for structures in current population
   ! Size=MYGA_CROSS+2 to prevent problems when doing uniform crossover (MYGA_CROSS=0)
   IF(.NOT. ALLOCATED(MYGA_POP_MOVE)) ALLOCATE(MYGA_POP_MOVE(MYGA_CROSS+2,-2:POPSIZE))
   MYGA_POP_MOVE=0
   IF(.NOT. ALLOCATED(MYGA_POP_GENOME)) ALLOCATE(MYGA_POP_GENOME(3*NATOMS,-2:POPSIZE))
   IF(.NOT. ALLOCATED(MYGA_POP_COORDS)) ALLOCATE(MYGA_POP_COORDS(3*NATOMS,-2:POPSIZE))
   ! Attempted & accepted moves for each gene in genome.
   ! Index starts at 0 to prevent over-writing when doing uniform crossover (MYGA_CROSS=0)
   IF(.NOT. ALLOCATED(MYGA_MOVE_HISTORY)) ALLOCATE(MYGA_MOVE_HISTORY(4,0:NATOMS))
   MYGA_MOVE_HISTORY=0
   IF (BLJCLUSTER.OR.BGUPTAT) THEN
      IF(.NOT.ALLOCATED(MYGA_POP_TYPE)) ALLOCATE(MYGA_POP_TYPE(NATOMS,-2:POPSIZE))
   ENDIF
   ! Set up random structure generator
   CALL RANDOM_SEED(SIZE=I)
   IF(.NOT. ALLOCATED(RANDOMSEED)) ALLOCATE(RANDOMSEED(I))
   CALL RANDOM_SEED(GET=RANDOMSEED)
   !First check if default options are overrridden
   IF (MYGA_L_SPHERE.OR.MYGA_L_CHAIN) THEN
      !Do nothing
   ELSE
      IF (BLNT) THEN
         MYGA_L_CHAIN=.TRUE.
      ELSE
         MYGA_L_SPHERE=.TRUE.
      ENDIF
   ENDIF
!   PRINT*, "Fill population with random structures"
   !Fill population with random structures
   DO I=1,MYGA_NSTRUC
      CALL MYGA_RANDOM_STRUC(I)
   ENDDO
!   PRINT*, "Copy convergence thresholds"
   ! Copy convergence thresholds
   MYGA_BQMAX=BQMAX
   MYGA_CQMAX=CQMAX
!   PRINT*, "Print search parameters"
   !Print search parameters
!   WRITE(MYUNIT,*)       "CMD_ARGS    = ", CMD_LINE
   WRITE(MYUNIT,*)       "Random seed = ", RANDOMSEED
   WRITE(MYUNIT,'(A,I5)')"GA> Population=",MYGA_NSTRUC
   WRITE(MYUNIT,'(A,I5)')"GA> Offspring=",MYGA_NOFF
   WRITE(MYUNIT,'(A,F10.3)')"GA> Mutation rate=",MYGA_MUT_RATE
   IF (STEP(1).EQ.0.AND.MYGA_MUT_RATE.GT.0) THEN
      WRITE(MYUNIT,'(A)')"GA> WARNING: Attempting mutation moves of size 0. Increase STEP."
   ENDIF
   IF (MYGA_L_ROUL) THEN
      WRITE(MYUNIT,'(A)')"GA> Roulette selection"
   ELSE
      WRITE(MYUNIT,'(A,I5)')"GA> Tournament selection, size=",MYGA_TOURN_SIZE
   ENDIF
   IF (MYGA_L_EPOCH) THEN
      WRITE(MYUNIT,'(A,D10.3,A,I5)')"GA> Epoch convergence threshold=",MYGA_EPOCH_THRESH," survival rate=",MYGA_EPOCH_SAVE
   ENDIF
   IF (MYGA_DUPLICATE_ETHRESH.NE.0) THEN
      WRITE(MYUNIT,'(A,D10.3)')"GA> Duplicate predator energy threshold=",MYGA_DUPLICATE_ETHRESH
   ENDIF
   IF (MYGA_CROSS.EQ.0) THEN
      WRITE(MYUNIT,'(A,I1,A)')"GA> Using uniform crossover"
   ELSE
      WRITE(MYUNIT,'(A,I1,A)')"GA> Using ",MYGA_CROSS,"-point crossover"
   ENDIF
   ! Set up initial basin-hopping GA steps
   MYGA_BH_INIT=MCSTEPS(1)
   MYGA_BH_STEPS=MCSTEPS(1)
   IF(MYGA_BH_INIT.GT.0.OR.MYGA_BH_INCR.GT.0) THEN
      WRITE(MYUNIT,'(A,F10.3)')"GA> Initial basin-hopping steps=",MYGA_BH_INIT
      WRITE(MYUNIT,'(A,F10.3)')"GA> Basin-hopping step increment=",MYGA_BH_INCR
   ENDIF
   RETURN
END
!
! Call appropriate mating routine for this system.
!
SUBROUTINE MYGA_MATE(I,J,CURR_STRUC)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY : BLNT
   IMPLICIT NONE
   INTEGER I,J,CURR_STRUC
!   IF (BLNT) THEN
!      IF (MYGA_CROSS.GT.0) THEN
!         CALL MYGA_BLN_MATE(I,J,CURR_STRUC)
!      ELSE
!         CALL MYGA_BLN_MATE_UNIFORM(I,J,CURR_STRUC)
!      ENDIF
!   ELSE
      CALL MYGA_CLUSTER_MATE(I,J,CURR_STRUC)
!   ENDIF
   RETURN
END
!
!Call appropriate mutation routine for this system
!Make a copy of structure at index INSTRUC and place mutant at OUTSTRUC
!
SUBROUTINE MYGA_MUTATE(INSTRUC,OUTSTRUC)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS , ONLY : BLNT
   IMPLICIT NONE
   INTEGER INSTRUC,OUTSTRUC
!   IF (BLNT) THEN
!      CALL MYGA_BLN_MUTATE_TOR(INSTRUC,OUTSTRUC)
!   ELSE ! Mutate by displacing atoms
      !CALL MYGA_CLUSTER_RANDOM(OUTSTRUC)
      CALL MYGA_CLUSTER_MUTATE_MOVE(INSTRUC,OUTSTRUC)
!   ENDIF
   RETURN
END
!
! Generate random structure at index I.
! Call the appropriate random structure generator for this system.
!
SUBROUTINE MYGA_RANDOM_STRUC(I)
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
   INTEGER I
!   IF (MYGA_L_CHAIN) THEN
!      CALL MYGA_BLN_RANDOM(I)
!   ELSE
      CALL MYGA_CLUSTER_RANDOM(I)
!   ENDIF
   MYGA_POP_ENERGY(I)=1D10
   RETURN
END
!
!Minimize a set of structures running from index I1 to I2
!Make an MPI version of this one day
!
SUBROUTINE MYGA_OPTIMISE_POP(I1,I2,STEPS,THRESH)
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
   INTEGER STEPS !Basin-hopping steps
   INTEGER I,I1,I2
   DOUBLE PRECISION THRESH !Convergence threshold
   PRINT*, "IN MYGA_OPTIMISE_POP BEFORE MYGA_OPT"
   DO I=I1,I2
      CALL MYGA_OPT(I,STEPS,THRESH)
   ENDDO
   PRINT*, "IN MYGA_OPTIMISE_POP AFTER MYGA_OPT"

   RETURN
END
!
!Optimize one structure
!This is the GA/GMIN interface
!
SUBROUTINE MYGA_OPT(STRUC,STEPS,THRESH)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS , ONLY : NATOMS,COORDS,BQMAX,BLNT
   USE QMODULE, ONLY : QMIN,QMINP
   USE MYLBFGSB, ONLY: LBFGSBC
   IMPLICIT NONE
   INTEGER STRUC !Index of structure in population arrays
   INTEGER I
   INTEGER STEPS !Basin-hopping steps
   DOUBLE PRECISION SCREENC(3*NATOMS)
   DOUBLE PRECISION THRESH !Convergence threshold
   ! First, make sure no coords or energies are left behind from old optimisations
   PRINT*, "   ! First, make sure no coords or energies are left behind from old optimisations"
   QMIN=1d10
   QMINP=0d0
   SCREENC=0d0
   ! Set convergence threshold
   BQMAX=THRESH
   ! Copy coordinates to GMIN array
   PRINT*, "Copy coordinates to GMIN array = "
   DO I=1,3*NATOMS
      COORDS(I,1)=MYGA_POP_COORDS(I,STRUC)
   ENDDO
   ! Optimise
   PRINT*, "IN MYGA_OPT BEFOR LBFGS = "
   CALL lbfgsbc(STRUC, STRUC, THRESH)
   PRINT*, "IN MYGA_OPT AFTER LBFGS = "
!   CALL MC(STEPS,1d0,SCREENC)
   ! Copy optimised coordinates & energy back to population array
   MYGA_POP_ENERGY(STRUC)=QMIN(1)
   DO I=1,3*NATOMS
      MYGA_POP_COORDS(I,STRUC)=(QMINP(1,I))
   ENDDO
   ! Only do this for BLN
!   IF (BLNT) THEN
!      CALL CART2IC(STRUC)
!   ENDIF
   RETURN
END
!
!Copy structure from index in to index out
!
SUBROUTINE COPY_STRUCTURE(IN,OUT)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY : NATOMS,BLJCLUSTER,BGUPTAT
   IMPLICIT NONE
   INTEGER IN, OUT,I, K
   MYGA_POP_ENERGY(OUT)=MYGA_POP_ENERGY(IN)    ! Copy energy
   MYGA_POP_FOUND(OUT)=MYGA_POP_FOUND(IN)
   DO I=1,MYGA_CROSS+2
      MYGA_POP_MOVE(I,OUT)=MYGA_POP_MOVE(I,IN) ! Copy move details
   ENDDO
   DO K=1,3*NATOMS
   ! Copy coordinates
      MYGA_POP_COORDS(K,OUT)=MYGA_POP_COORDS(K,IN)
   !Copy genome
      MYGA_POP_GENOME(K,OUT)=MYGA_POP_GENOME(K,IN)
   END DO
   IF (BLJCLUSTER.OR.BGUPTAT) THEN !Copy atom types for binary clusters
      DO K=1,NATOMS
         MYGA_POP_TYPE(K,OUT)=MYGA_POP_TYPE(K,IN)
      ENDDO
   ENDIF
   RETURN
END
!
!Sort population in energy order.
!Sort from index=1 to index=sortsize
!
SUBROUTINE MYGA_SORT(SORTSIZE)
USE GA_PARAMS
USE GA_POPULATION
IMPLICIT NONE
! INTERNAL INTEGER VALUES
INTEGER IR   !
INTEGER I    !
INTEGER J    !
INTEGER L    !
INTEGER SORTSIZE
! INTERNAL REAL VALUES
! (ORIGINAL COMMENT)
!     SORTS THE ENERGY ARRAY IN ORDER OF INCREAING ENERGY
!     AND MAKES THE CORRESPONDING RE-ARRANGMENT OF OTHER ARRAYS
!     USES THE NUMERICAL RECIPES HEAP SORT
!     START OF HEAP SORT
! POPULATION INDEX -1 IS TEMPRARY HOLDER FOR SWAPPING
L = (SORTSIZE)/2+1    ! MIDDLE POINT
IR = SORTSIZE    ! DUPLICATED
10   CONTINUE
IF(L.GT.1)THEN      ! IF ARRAY SIZE IS GREATER THAN 1
   L=L-1        ! TAKE ONE BELOW MIDDLE POINT
   CALL COPY_STRUCTURE(L,-1)
ELSE        ! ELSE IF ARRAY SIZE IS SMALLER THAN 1

   CALL COPY_STRUCTURE(IR,-1)
   CALL COPY_STRUCTURE(1,IR)
   IR=IR-1    ! SUBTRACT 1 FROM IR
   IF(IR.EQ.1)THEN  ! CHECK IF WE ARE AT LAST VALUE
      CALL COPY_STRUCTURE(-1,1)
      ! Check if best structure in current epoch is better than best in all epochs
      RETURN    ! RETURN FUNCTION
   ENDIF
ENDIF
I=L    ! SET I TO MIDDLE-1
J=L+L    ! SET J TO MIDDLE
20 IF(J.LE.IR)THEN      ! IF MIDDLE <= LENGTH
   IF(J.LT.IR)THEN    ! IF MIDDLE < LENGTH
      IF(MYGA_POP_ENERGY(J).LT.MYGA_POP_ENERGY(J+1))J=J+1    ! SWAP J IF ENERGY OF NEXT IS GREATER
   ENDIF
   IF(MYGA_POP_ENERGY(-1).LT.MYGA_POP_ENERGY(J))THEN
      ! IF TEMP ENERGY IS GREATER THAN THAT AT MIDDLE(+1)
      CALL COPY_STRUCTURE(J,I)
      I=J    ! MIDDLE(+1)
      J=J+J  ! DOUBLE J
   ELSE
                ! OTHERWISE J = LENGTH + 1
         J=IR+1
   ENDIF
        ! RETURN TO START OF LOOP
   GO TO 20
ENDIF
CALL COPY_STRUCTURE(-1,I)
! RETURN TO START OF SUBROUTINE
GO TO 10
END
!
!Look through population and remove duplicate structures
!LENGTH=maximum index to check for duplicates
!REMOVED=number of duplicated removed
!
SUBROUTINE MYGA_REMOVE_DUPLICATES(LENGTH,REMOVED)
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
! INTERNAL INTEGER VALUES
   INTEGER LENGTH
   INTEGER I
   INTEGER, ALLOCATABLE ::  MARK(:) ! TEMPORARY ARRAY TO STORE DELETED POSITIONS
   INTEGER COUNT  ! TEMPORARY COUNTER
   INTEGER REMOVED ! COUNT OF REMOVED DUPLICATES
   LOGICAL ARE_DUPLICATES !Function
   LOGICAL UNIQUE_STRUCTURE !Function
   ALLOCATE(MARK(LENGTH))
   REMOVED=0
!     REMOVES DUPLICATE STRINGS FROM THE POPULATION
   DO I=1,LENGTH
! FOR ALL POSITIONS SET REFERENCE TO 0
      MARK(I) = 0
   END DO
   COUNT = 1  ! COUNTER SET TO START
   MARK(1) = 1 ! FIRST POINT IN ARRAY SET TO START
   DO I=2,LENGTH
! FOR ENTIRE LENGTH OF ARRAY BEYOND FIRST ATOM
! (ORIGINAL COMMENT)
!     IDENTIFY UNIQUE STRINGS
      IF (UNIQUE_STRUCTURE(I)) THEN
! IF ENERGIES ARE ESSENTIALLY THE SAME AS PREVIOUS POINT
! NOTE: THIS IS USED AFTER ARRAYS HAVE BEEN SORTED
         COUNT = COUNT + 1  ! INCREASE COUNT OF IDENTICAL STRINGS BY 1
         MARK(COUNT) = I  ! SET NEXT DUPLICATE POSITION TO COUNTER POINT
      ELSE
         REMOVED=REMOVED+1
      END IF
   END DO
! (ORIGINAL COMMENT)
!     REMOVE NON-UNIQUE STRINGS
   DO I = MARK(1),COUNT
! FOR POINT NEED TO BE REMOVED TO END OF COUNTER
! LOOP THROUGH ALL POSSIBLE OPTIONS AND MOVE THEM FORWARD BY ONE
      CALL COPY_STRUCTURE(MARK(I),I)
   END DO
   DO I=COUNT + 1, MYGA_NSTRUC + MYGA_NOFF + MYGA_NMUT
      MYGA_POP_ENERGY(I) = 1D10
   ENDDO
! MIGHT NEED SOMETHING HERE TO CATCH CASES WHERE LOTS OF STRUCTURES ARE REMOVED AND POPULATION BECOMES TOO SMALL.
   RETURN  ! END OF SUBROUTINE
END
!
!Are structures I and I the same?
!For BLN proteins, this checks all the torsion angles
!For other systems, this just checks the energies of I and J
!
FUNCTION ARE_DUPLICATES(I,J)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY : NATOMS,BLNT
   IMPLICIT NONE
   INTEGER I,J,K
   INTEGER IATOM
   LOGICAL ARE_DUPLICATES
   DOUBLE PRECISION :: TOR
!   DOUBLE PRECISION, PARAMETER :: PI = 4.d0*DATAN(1.D0)
   DOUBLE PRECISION, PARAMETER :: PI = 3.141592654
   ! CHECK FOR CLOSE ENERGIES
   IF(DABS(MYGA_POP_ENERGY(I)-MYGA_POP_ENERGY(J)).LT.MYGA_DUPLICATE_ETHRESH) THEN
   ! Geometry check
      IF (BLNT) THEN
         DO IATOM =4,NATOMS
            TOR=DABS(MYGA_POP_GENOME(3*IATOM,I)-MYGA_POP_GENOME(3*IATOM,J))
            IF (TOR.GT.PI) THEN
               TOR =2*PI-TOR
            ENDIF
            ! Check for close torsion angles
            IF (TOR.GT.MYGA_DUPLICATE_GTHRESH) THEN
               ARE_DUPLICATES=.FALSE.
               !WRITE(*,*)"DEBUG:",I,J,K,TOR/DEG2RAD,MYGA_ENERGIES(I),MYGA_ENERGIES(J)
               RETURN
            ENDIF
         ENDDO
       ENDIF
       !Passed geometry tests
       ARE_DUPLICATES=.TRUE.
   ELSE
      ARE_DUPLICATES=.FALSE.
   ENDIF
   RETURN
END
!
!Is structure i unique compared to all more stable members of population?
!Only works properly after population has been sorted
!
FUNCTION UNIQUE_STRUCTURE(I)
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
   INTEGER I,J,K
   LOGICAL UNIQUE_STRUCTURE,ARE_DUPLICATES
   ! LOOP THROUGH STRUCTURES, STARTING WITH MOST STABLE
   DO J=1,I-1
      ! CHECK FOR CLOSE ENERGY
      IF(DABS(MYGA_POP_ENERGY(I)-MYGA_POP_ENERGY(J)).LT.MYGA_DUPLICATE_ETHRESH) THEN
         IF (ARE_DUPLICATES(I,J)) THEN
            UNIQUE_STRUCTURE=.FALSE. ! DUPLICATE FOUND
            ! Check for earliest instance of this structure (matters most for loose convergence of minimisation)
            IF (MYGA_POP_FOUND(I).LT.MYGA_POP_FOUND(J)) THEN !New structure is slight improvement on old structure
                MYGA_POP_FOUND(J)=MYGA_POP_FOUND(I) ! Keep earliest hit of this structure
                MYGA_POP_MOVE(1,J)=0 ! Mark as old structure
                DO K=1,MYGA_CROSS
                   MYGA_POP_MOVE(1+K,J)=MYGA_POP_MOVE(1+K,I)
                ENDDO
            ENDIF
            RETURN
         ENDIF
      ENDIF
   ENDDO
   UNIQUE_STRUCTURE=.TRUE. !RUN OUT OF STRUCTURES, NO DUPLICATES FOUND
   RETURN
END
!
!Calculate the mean energy of all members of the population
!
FUNCTION MYGA_MEAN_ENERGY()
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
   INTEGER I
   DOUBLE PRECISION MYGA_MEAN_ENERGY
   MYGA_MEAN_ENERGY=0D0
   DO I=1,MYGA_NSTRUC
      MYGA_MEAN_ENERGY=MYGA_MEAN_ENERGY+MYGA_POP_ENERGY(I)
   ENDDO
   MYGA_MEAN_ENERGY=MYGA_MEAN_ENERGY/DBLE(MYGA_NSTRUC)
   RETURN
END
!
!Start a new epoch
!
SUBROUTINE MYGA_NEW_EPOCH()
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY : MYUNIT,MCSTEPS
   IMPLICIT NONE
   INTEGER I
   INTEGER REMOVED
   CHARACTER*50 FILENAME
   IF (MYGA_EPOCH_DUMP) THEN
      WRITE(FILENAME,*) MYGA_COUNT_EPOCH
      FILENAME="epoch."//TRIM(ADJUSTL(FILENAME))
      CALL MYGA_WRITE_LOWEST(filename)
   ENDIF
   WRITE(MYUNIT,'(A)') "GA> Population converged. Starting new epoch."
   !Sort population+saved structures
   CALL MYGA_SORT(MYGA_NSTRUC*3+MYGA_NOFF*2)
   CALL MYGA_REMOVE_DUPLICATES(MYGA_NSTRUC*3+MYGA_NOFF*2,REMOVED)
   ! Save best structures
   DO I=1,MYGA_NSTRUC
      CALL COPY_STRUCTURE(I,I+2*(MYGA_NSTRUC+MYGA_NOFF))
   ENDDO
   !Make new population of random structures
   DO I=MYGA_EPOCH_SAVE+1,MYGA_NSTRUC
      CALL MYGA_RANDOM_STRUC(I)
   ENDDO
   !Optimise new population
   CALL MYGA_OPTIMISE_POP(MYGA_EPOCH_SAVE+1,MYGA_NSTRUC,0,MYGA_BQMAX)
   !Reset empty structures
   DO I=MYGA_NSTRUC+1,(MYGA_NSTRUC+MYGA_NOFF)*2
      MYGA_POP_ENERGY(I)=1D10
   ENDDO
   CALL MYGA_SORT(MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT)
   MYGA_COUNT_EPOCH=MYGA_COUNT_EPOCH+1
! Reset basin-hopping length
   MYGA_BH_STEPS=MYGA_BH_INIT
   MCSTEPS(1)=INT(MYGA_BH_STEPS)
   RETURN
END
!
!Write best structures to file "lowest"
!
SUBROUTINE MYGA_FINALIO()
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY :NATOMS,MYUNIT, NPCALL
   IMPLICIT NONE
   INTEGER I,J1,J2,J3,JJ
   DOUBLE PRECISION TIME
   INTEGER REMOVED
   !Sort population+saved structures
   CALL MYGA_SORT(MYGA_NSTRUC*3+MYGA_NOFF*2)
   CALL MYGA_REMOVE_DUPLICATES(MYGA_NSTRUC*3+MYGA_NOFF*2,REMOVED)
   ! Save best structures
   !DO I=1,MYGA_NSTRUC
   !   CALL COPY_STRUCTURE(I,I+2*(MYGA_NSTRUC+MYGA_NOFF))
   !ENDDO
   CALL MYGA_WRITE_LOWEST("lowest")
   ! Write GA move summary
   WRITE(MYUNIT,'(A)')"GA> Genetic operations summary:"
   DO I=1,NATOMS
      WRITE(MYUNIT,'(5(A,I8))')"GA> Position=",I,&
      " Crossovers=",MYGA_MOVE_HISTORY(1,I),&
      " Accepted=",MYGA_MOVE_HISTORY(2,I),&
      " Mutations=",MYGA_MOVE_HISTORY(3,I),&
      " Accepted=",MYGA_MOVE_HISTORY(4,I)
   ENDDO
   CALL MYCPU_TIME(TIME)
   WRITE(MYUNIT,'(4(A,I12),A,F15.3)')&
   "GA> Run finished: Evaluations=",NPCALL,&
   " Minimisations=",MYGA_COUNT_MIN,&
   " Generations=",CURR_GEN,&
   " Epochs=",MYGA_COUNT_EPOCH,&
   " Time=",TIME
   RETURN
END
!
! Write all structures in current population to file
! At the end of a search, this is used to write the "lowest" file
!
SUBROUTINE MYGA_WRITE_LOWEST(FILENAME)
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS, ONLY :NATOMS,BEADLETTER,BLNT,BLJCLUSTER,BGUPTAT
   IMPLICIT NONE
   CHARACTER*(*) FILENAME
   INTEGER I,J1,J2,J3,JJ
   DOUBLE PRECISION TIME
   CHARACTER*2 S2
   INTEGER REMOVED
   OPEN(25,FILE=TRIM(ADJUSTL(FILENAME)))
   DO J1=1,MYGA_NSTRUC
      JJ=J1
      WRITE(25,'(I12)')NATOMS
      WRITE(25,10)J1,MYGA_POP_ENERGY(JJ),MYGA_POP_FOUND(JJ)
10    FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at generation ',I8)
      IF (BLNT) THEN
         DO J2=1,NATOMS
            WRITE(25,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(MYGA_POP_COORDS((3*(J2-1)+J3),JJ),J3=1,3)
         ENDDO
      ELSEIF (BLJCLUSTER.OR.BGUPTAT) THEN
         DO J2=1,NATOMS
            IF (MYGA_POP_TYPE(J2,JJ).EQ.1) THEN
               S2="LA"
            ELSE
               S2="LB"
            ENDIF
            WRITE(25,'(A,3F20.10)') S2,(MYGA_POP_COORDS((3*(J2-1)+J3),JJ),J3=1,3)
         ENDDO
      ELSE
         DO J2=1,NATOMS
            WRITE(25,'(A,3F20.10)') 'X ',(MYGA_POP_COORDS((3*(J2-1)+J3),JJ),J3=1,3)
         ENDDO
      ENDIF
   ENDDO
   CLOSE(25)
   RETURN
END
!
!Print details of whole population for debugging
!
SUBROUTINE MYGA_PRINT_POP
   USE GA_PARAMS
   USE GA_POPULATION
   IMPLICIT NONE
   INTEGER I,J
   WRITE(*,*)"ENERGY SUMMARY"
   DO I=1,MYGA_NSTRUC+MYGA_NOFF+MYGA_NMUT
      WRITE(*,*)I,MYGA_POP_ENERGY(I),MYGA_POP_FOUND(I),(MYGA_POP_MOVE(J,I),J=1,MYGA_CROSS+1)
   ENDDO
   RETURN
END
!
! Dump genome of whole population to file
! Works for BLN proteins
!
SUBROUTINE MYGA_DUMP_IC()
   USE GA_PARAMS
   USE GA_POPULATION
   USE COMMONS , ONLY : NATOMS
   IMPLICIT NONE
   INTEGER I,J,CURR_STRUC
   CHARACTER*10 FILENAME
   WRITE(FILENAME,'(I10)')CURR_GEN
   OPEN (50,FILE="POP."//ADJUSTL(FILENAME))
   DO CURR_STRUC=1,MYGA_NSTRUC
      DO J=4,NATOMS
         WRITE(50,*)J,MYGA_POP_GENOME(3*J,CURR_STRUC),MYGA_POP_GENOME(3*J,1)
      ENDDO
      WRITE(50,*)
   ENDDO
   CLOSE(50)
   RETURN
END
!
!SUBROUTINE ANALYZE(I1, I2, IGEN)
!    USE COMPARATOR, ONLY: CHECK_UNIQUENESS
!    USE GA_POPULATION, ONLY: MYGA_POP_COORDS, MYGA_POP_ENERGY
!    USE COMMONS, ONLY: NATOMS
!    USE GA_PARAMS, ONLY: MYGA_DUPLICATE_ETHRESH
!    IMPLICIT NONE
!    INTEGER, INTENT(IN) :: I1, I2, IGEN
!    INTEGER, SAVE :: I, J1, MYUNIT1, MYUNIT2, NUM_UNIQUE_MINIMA = 0, DUPLICATE_INDX = 0
!    REAL(KIND=8), DIMENSION(3*NATOMS, 10000), SAVE :: UNIQUE_COORDS != 1.0D10
!    REAL(KIND=8), DIMENSION(10000),SAVE :: UNIQUE_ENERGY != 1.0D10
!    REAL(KIND=8), DIMENSION(10000, 4), SAVE :: MINIMA_STAT != 1.0D0
!    LOGICAL, SAVE :: ISUNIQUE, INIT = .TRUE.
!    IF (INIT) THEN
!        OPEN(NEWUNIT=MYUNIT1, FILE="minstat.txt",ACTION="WRITE", STATUS="REPLACE")
!        OPEN(NEWUNIT=MYUNIT2, FILE="stat.txt", action = "write", status="replace")
!        INIT = .FALSE.
!    ENDIF
!    IF(NUM_UNIQUE_MINIMA == 0) THEN
!        NUM_UNIQUE_MINIMA = 1
!        UNIQUE_COORDS(:,1) = MYGA_POP_COORDS(:, I1)
!        UNIQUE_ENERGY(1)   = MYGA_POP_ENERGY(I1)
!        J1 = I1+1
!    ELSE
!        J1 = I1
!    ENDIF
!    DO I=J1, I2
!        IF(isnan(MYGA_POP_ENERGY(I))) THEN
!            PRINT*, "error in analyze subroutine"
!            STOP
!        ENDIF
!        CALL CHECK_UNIQUENESS(MYGA_POP_ENERGY(I), UNIQUE_ENERGY(1:NUM_UNIQUE_MINIMA), &
!            MYGA_POP_COORDS(:,I), UNIQUE_COORDS(:,1:NUM_UNIQUE_MINIMA), 3*NATOMS, &
!            NUM_UNIQUE_MINIMA, DUPLICATE_INDX, ISUNIQUE, MYGA_DUPLICATE_ETHRESH+1.0E-5)
!        IF (ISUNIQUE .AND. DUPLICATE_INDX == 0) THEN
!            NUM_UNIQUE_MINIMA = NUM_UNIQUE_MINIMA + 1
!            UNIQUE_ENERGY(NUM_UNIQUE_MINIMA) = MYGA_POP_ENERGY(I)
!            UNIQUE_COORDS(:,NUM_UNIQUE_MINIMA) = MYGA_POP_COORDS(:,I)
!            MINIMA_STAT(NUM_UNIQUE_MINIMA, 1) = IGEN
!            MINIMA_STAT(NUM_UNIQUE_MINIMA, 2) = IGEN
!            MINIMA_STAT(NUM_UNIQUE_MINIMA, 3) = 1
!            MINIMA_STAT(NUM_UNIQUE_MINIMA, 4) = MYGA_POP_ENERGY(I)
!            WRITE(MYUNIT1, *) IGEN, NUM_UNIQUE_MINIMA, INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 1)),&
!              INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 2)),  INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 3)),&
!              (MINIMA_STAT(NUM_UNIQUE_MINIMA, 4))
!            WRITE(MYUNIT2, *)
!            WRITE(MYUNIT2, *), "ENERGY                    = ", MYGA_POP_ENERGY(I)
!            WRITE(MYUNIT2, *), "FIRST FOUND AT GENERATION = ", INT(IGEN)
!            WRITE(MYUNIT2, *), "NUMBER OF TIMES VISITED   = ", INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 3))
!        ELSEIF(DUPLICATE_INDX /= 0) THEN
!            MINIMA_STAT(NUM_UNIQUE_MINIMA, 2) = IGEN
!            MINIMA_STAT(DUPLICATE_INDX, 3) = MINIMA_STAT(DUPLICATE_INDX, 3) + 1
!!            PRINT*, "I, DUPLICATE_INDX,  ISUNIQUE = ",I, DUPLICATE_INDX, ISUNIQUE, MYGA_POP_ENERGY(I)
!            WRITE(MYUNIT2, *)
!            WRITE(MYUNIT2, *), "CURRENT ENERGY            = ", MYGA_POP_ENERGY(I)
!            WRITE(MYUNIT2, *), "FIRST FOUND AT GENERATION = ", INT(MINIMA_STAT(DUPLICATE_INDX, 1))
!            WRITE(MYUNIT2, *), "LAST FOUND AT GENERATION  = ", INT(MINIMA_STAT(DUPLICATE_INDX, 2))
!            WRITE(MYUNIT2, *), "NUMBER OF TIMES VISITED   = ", INT(MINIMA_STAT(DUPLICATE_INDX, 3))
!            WRITE(MYUNIT2, *), "LAST FOUND ENERGY         = ", UNIQUE_ENERGY(DUPLICATE_INDX)
!            WRITE(MYUNIT1, *) IGEN, NUM_UNIQUE_MINIMA, INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 1)),&
!              INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 2)),  INT(MINIMA_STAT(NUM_UNIQUE_MINIMA, 3)),&
!              (MINIMA_STAT(NUM_UNIQUE_MINIMA, 4))
!
!        ELSE
!            PRINT*, "ERROR : DUPLICATE_INDX = ", DUPLICATE_INDX
!            PRINT*, "ISUNIQUE               = ", ISUNIQUE
!            STOP
!        ENDIF
!    ENDDO
!    PRINT*, "NUMBER OF UNIQUE MINIMA SO FAR = ", NUM_UNIQUE_MINIMA
!    PRINT*,
!END SUBROUTINE ANALYZE
