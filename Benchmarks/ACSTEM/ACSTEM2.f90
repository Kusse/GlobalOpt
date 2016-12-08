MODULE ACSTEM2
    IMPLICIT NONE
    INTEGER, PARAMETER, PRIVATE :: DP     = KIND(1.0D0)
    INTEGER, PARAMETER, PRIVATE :: NATOMS = 309                     ! NUMBER OF ATOMS
    INTEGER, PARAMETER, PUBLIC  :: NDIMS  = 3*NATOMS + 2            ! NUMBER OF DIMENSIONS
    INTEGER, PARAMETER, PRIVATE :: NROWS  = 50                      ! NUMBER OF ROWS, IN PIXELS, IN THE IMAGE
    INTEGER, PARAMETER, PRIVATE :: NCOLS  = 50                      ! NUMBER OF COLUMNS, IN PIXELS, IN THE IMAGE
    REAL(KIND=DP), PARAMETER, PRIVATE  :: IMAGE_WIDTH = 24.0D0      ! IN NANOMETERS
    REAL(KIND=DP), PARAMETER, PRIVATE  :: IMAGE_HEIGHT= 24.0D0      ! IN NANOMETERS

    REAL(KIND=DP), PARAMETER, PRIVATE  :: LBX = -IMAGE_WIDTH/2.0D0
    REAL(KIND=DP), PARAMETER, PRIVATE  :: UBX =  IMAGE_WIDTH/2.0D0
    REAL(KIND=DP), PARAMETER, PRIVATE  :: LBY = -IMAGE_HEIGHT/2.0D0
    REAL(KIND=DP), PARAMETER, PRIVATE  :: UBY =  IMAGE_HEIGHT/2.0D0
    REAL(KIND=DP), DIMENSION(NDIMS, NDIMS), SAVE :: HESSIAN          ! HESSIAN MATRIX OF THE OBJECTIVE FUNCTION
    REAL(KIND=DP), DIMENSION(NROWS, NCOLS), SAVE :: MODEL_IMAGE
    REAL(KIND=DP), DIMENSION(NROWS, NCOLS), PRIVATE , SAVE :: EXPT_IMAGE

    ! POTENTIAL PARAMETERS
    INTEGER, PARAMETER, PRIVATE :: ATOMTYPE = 8                 ! GOLD
    INTEGER, PARAMETER, PRIVATE :: ATOMZ    = 79                ! ATOMIC NUMBER
    REAL(KIND=DP), PARAMETER, PRIVATE       :: RHO = 1.0D0
    CHARACTER(LEN=5), PARAMETER, PRIVATE    :: POTTYPE = "MORSE"

    ! GAUSSIAN SETTINGS
    REAL(KIND=DP), DIMENSION(NDIMS), PROTECTED, SAVE, PUBLIC :: STEMLB = MIN(LBX, LBY)       ! LOWER LIMIT ON THE VALUE OF COORDINATES
    REAL(KIND=DP), DIMENSION(NDIMS), PROTECTED, SAVE, PUBLIC :: STEMUB = MAX(UBX, UBY)       ! UPPER LIMIT ON THE VALUE OF THE COORDINATES
    REAL(KIND=DP), PARAMETER, PRIVATE       :: GAUSS_WIDTH_DEFAULT = 1.0D0        ! IN NANOMETERS
    REAL(KIND=DP), PARAMETER, PRIVATE       :: GAUSS_HEIGHT_DEFAULT = 1.0D0     ! IN NANOMETERS
    CHARACTER(LEN=100), PARAMETER, PRIVATE  :: EXPT_IMFILE = "img.txt"

    LOGICAL, SAVE :: INIT=.TRUE.
    CONTAINS
    SUBROUTINE SETBOUNDS()
        IMPLICIT NONE
        INTEGER :: I
        PRINT*, "SETING UP"
        !STOP
        DO I=1, NATOMS
            STEMLB(3*I - 2) = -IMAGE_WIDTH/2.0D0
            STEMLB(3*I - 1) = -IMAGE_HEIGHT/2.0D0
            STEMLB(3*I)     = -MAX(IMAGE_WIDTH/2.0D0, IMAGE_HEIGHT/2.0D0)
            STEMUB(3*I - 2) = IMAGE_WIDTH/2.0D0
            STEMUB(3*I - 1) = IMAGE_HEIGHT/2.0D0
            STEMUB(3*I)     = MAX(IMAGE_WIDTH/2.0D0, IMAGE_HEIGHT/2.0D0)
        ENDDO
        ! BOUNDS ON GAUSSAIN WIDTH
        STEMLB(NDIMS - 1)  = 0.0D0
        STEMUB(NDIMS - 1)  = 2.0D0
        ! BOUNDS ON GAUSSIAN HEIGHT
        STEMLB(NDIMS)      = 0.0D0
        STEMUB(NDIMS)      = 0.0D0
        IF(INIT) THEN
            !EXPT_IMAGE = 0.0D0
            CALL LOADIMAGE(EXPT_IMFILE, NROWS, NCOLS, EXPT_IMAGE)
            INIT = .FALSE.
            PRINT*, "EXPERIMENTAL IMAGE LOADED."
            !STOP
        ENDIF
    END SUBROUTINE SETBOUNDS

    SUBROUTINE MYOBJFUNC(ND, XSTEM, ENERGY, GRAD, GTEST, HTEST)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ND
        LOGICAL, INTENT(IN) :: GTEST, HTEST
        REAL(KIND=8), INTENT(IN),  DIMENSION(ND):: XSTEM
        REAL(KIND=8), INTENT(OUT), DIMENSION(ND):: GRAD
        REAL(KIND=8), INTENT(OUT) :: ENERGY

        REAL(KIND=8), DIMENSION(ND-2)   :: GPOT
        REAL(KIND=8), DIMENSION(ND)     :: GSTEM
        REAL(KIND=8) :: ESTEM, EPOT, W
        IF(INIT) THEN
            !EXPT_IMAGE = 0.0D0
            CALL LOADIMAGE(EXPT_IMFILE, NROWS, NCOLS, EXPT_IMAGE)
            INIT = .FALSE.
            PRINT*, "EXPERIMENTAL IMAGE LOADED."
            !STOP
        ENDIF
        CALL MY_ACSTEM_MODEL(ND, XSTEM, ESTEM, GSTEM, GTEST, HTEST)
        IF (ESTEM < 0.1) THEN
            CALL POTENTIAL(NATOMS, XSTEM(1:ND-2), EPOT, GPOT, GTEST, HTEST)
        ELSE
            epot = 0.0d0
        ENDIF
        IF (ESTEM < 0.01) THEN
            PRINT*, "SIMULATION CONVERGED!"
            PRINT*, "ERROR  = ", ESTEM
            PRINT*, "ENERGY = ", EPOT
            PRINT*, "FORCE  = ", NORM2(GPOT)
        ENDIF
        W = 0.0D0   !E1/(E2 + 1.0E-20)
        PRINT*, "ESTEM = ", ESTEM
        PRINT*, "EPOT = ", EPOT
        IF (ESTEM <= 0.1) THEN
            ENERGY = ESTEM
        ELSE
            ENERGY = ESTEM + W*EPOT
        ENDIF
        GRAD = GSTEM
        GRAD(1:ND-2) = GRAD(1:ND-2) + GPOT
    END SUBROUTINE MYOBJFUNC

    SUBROUTINE MY_ACSTEM_MODEL(ND, XSTEM, ESTEM, GSTEM, GTEST, HTEST)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ND
        REAL(KIND=DP), INTENT(IN), DIMENSION(ND)  :: XSTEM
        REAL(KIND=DP), INTENT(OUT), DIMENSION(ND) :: GSTEM
        REAL(KIND=DP), INTENT(OUT)  :: ESTEM
        LOGICAL, INTENT(IN)         :: GTEST, HTEST

        REAL(KIND=DP), DIMENSION(NATOMS) :: IMIJ, DX, DY
        REAL(KIND=DP) :: XSTEP, YSTEP, HEIGHT, ALPHA, SIGMA, BETA
        REAL(KIND=DP) :: XJ, YI, ERRIJ, TMP, DHDB
        INTEGER :: I, J, K

        XSTEP = (UBX - LBX)/DBLE(NCOLS)
        YSTEP = (UBY - LBY)/DBLE(NROWS)
        ESTEM = 0.0D0
        GSTEM = 0.0D0
        SIGMA = XSTEM(ND-1)
        ALPHA = 1.0D0/(2.0D0*SIGMA*SIGMA)
        BETA = XSTEM(ND)
        HEIGHT = (DBLE(ATOMZ))**BETA
        DHDB = HEIGHT*LOG(DBLE(ATOMZ))
        IF (GTEST) THEN
            DO J=1, NCOLS
                XJ = LBX + XSTEP*DBLE(J)
                DO I=1, NROWS
                    YI = LBY + YSTEP*DBLE(I)
                    DO K= 1, NATOMS
                        DX(K) = XJ - XSTEM(3*K-2)
                        DY(K) = YI - XSTEM(3*K-1)
                        IMIJ(K) = HEIGHT*EXP(-ALPHA*(DX(K)**2 + DY(K)**2))
                    ENDDO
                    MODEL_IMAGE(I,J) = SUM(IMIJ)
                    ERRIJ = EXPT_IMAGE(I,J) - MODEL_IMAGE(I,J)
                    ESTEM = ESTEM + ERRIJ*ERRIJ
                    TMP = -4.0D0*ALPHA*ERRIJ
                    DO K=1, NATOMS
                        GSTEM(3*K - 2) = GSTEM(3*K - 2) + TMP*DX(K)*IMIJ(K)
                        GSTEM(3*K - 1) = GSTEM(3*K - 1) + TMP*DY(K)*IMIJ(K)
                        GSTEM(ND-1)    = GSTEM(ND-1) - 2.0D0*ERRIJ*(DX(K)**2 + DY(K)**2)*IMIJ(K)*ALPHA/SIGMA
                    ENDDO
                    GSTEM(ND)   = GSTEM(ND) - 2.0*ERRIJ*MODEL_IMAGE(I,J)*DHDB/HEIGHT
                ENDDO
            ENDDO
        ELSE IF (HTEST) THEN
            HESSIAN = 0.0D0
            PRINT*, "The second derivative function for stem model image is not implemented."
        ELSE
            DO J=1, NCOLS
                XJ = LBX + XSTEP*DBLE(J)
                DO I=1, NROWS
                    YI = LBY + YSTEP*DBLE(I)
                    DO K= 1, NATOMS
                        DX(K) = XJ - XSTEM(3*K-2)
                        DY(K) = YI - XSTEM(3*K-1)
                        IMIJ(K) = HEIGHT*EXP(-ALPHA*(DX(K)**2 + DY(K)**2))
                    ENDDO
                    MODEL_IMAGE(I,J) = SUM(IMIJ)
                    ERRIJ = EXPT_IMAGE(I,J) - MODEL_IMAGE(I,J)
                    ESTEM = ESTEM + ERRIJ*ERRIJ
                ENDDO
            ENDDO
        ENDIF
    END SUBROUTINE MY_ACSTEM_MODEL

    SUBROUTINE MY_IMAGE(N, COORDS, NR, NC, SIGMA, BETA, IMAGE)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NR, NC
        REAL(KIND=DP), INTENT(IN) :: SIGMA, BETA
        REAL(KIND=DP), DIMENSION(N,3), INTENT(IN) :: COORDS
        REAL(KIND=DP), DIMENSION(NR,NC), INTENT(OUT) :: IMAGE
        REAL(KIND=DP) :: DX, DY, XJ, YI, HEIGHT, ALPHA, IM(N), XSTEP, YSTEP
        INTEGER :: I, J, K
        XSTEP = (UBX - LBX)/NC
        YSTEP = (UBY - LBY)/NR
        HEIGHT = (DBLE(ATOMZ))**BETA
        ALPHA = 1.0D0/(2.0D0*SIGMA*SIGMA)
        DO J=1, NC
            XJ = LBX + DBLE(J)*XSTEP
            DO I=1, NR
                YI = LBY + DBLE(I)*YSTEP
                DO K=1, N
                    DX = XJ - COORDS(K,1)
                    DY = YI - COORDS(K,2)
                    IM(K) = HEIGHT*EXP(-ALPHA*(DX*DX + DY*DY))
                ENDDO
                IMAGE(I,J) = SUM(IM)
            ENDDO
        ENDDO
    END SUBROUTINE MY_IMAGE

    SUBROUTINE LOADIMAGE(FILENAME, NROWS, NCOLS, IMAGE)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NROWS, NCOLS
        REAL(8), INTENT(OUT), DIMENSION(NROWS, NCOLS) :: IMAGE
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER :: I, J, IUNIT, ERR
        OPEN(NEWUNIT=IUNIT, FILE=ADJUSTL(TRIM(FILENAME)), &
                STATUS="OLD", ACTION="READ", IOSTAT=ERR)
        IF (ERR /= 0) THEN
            PRINT*, "ERROR! Problem opening file " , &
                ADJUSTL(TRIM(FILENAME)), ". Please check if the file exists"
            CLOSE(IUNIT)
            STOP
        ELSE
            DO I=1, NROWS
                PRINT*, "I= ", I
                READ(IUNIT,*) IMAGE(I, :)
            ENDDO

            CLOSE(IUNIT)
        ENDIF
    END SUBROUTINE LOADIMAGE

    SUBROUTINE WRITEIMAGE(FILENAME, NROWS, NCOLS, IMAGE)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NROWS, NCOLS
        REAL(8), INTENT(IN), DIMENSION(NROWS, NCOLS) :: IMAGE
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER :: I, J, IUNIT, ERR
        OPEN(NEWUNIT=IUNIT, FILE=ADJUSTL(TRIM(FILENAME)), STATUS="UNKNOWN", ACTION="WRITE", IOSTAT=ERR)
        IF (ERR /= 0) THEN
            PRINT*, "ERROR! Problem opening file for writing " , &
                ADJUSTL(TRIM(FILENAME)), ". Please check if the file already exists"
            CLOSE(IUNIT)
            STOP
        ELSE
            DO I=1, NROWS
                WRITE(IUNIT,*) IMAGE(I, :)
            ENDDO
            CLOSE(IUNIT)
        ENDIF
    END SUBROUTINE WRITEIMAGE

    SUBROUTINE READCOORDINATES(FILENAME, N, FILETYPE, COORDS)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        CHARACTER(LEN=*), INTENT(IN) :: FILETYPE
        REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS, 3) :: COORDS
        CHARACTER(LEN=4) :: FTYPE
        CHARACTER(LEN=2), DIMENSION(N) :: ATOMIC_SYMBOL
        INTEGER :: I, J, IUNIT, ERR
        OPEN(NEWUNIT=IUNIT, FILE=ADJUSTL(TRIM(FILENAME)), STATUS="OLD", ACTION="READ", IOSTAT=ERR)
        IF (ERR /= 0) THEN
            PRINT*, "ERROR! Problem opening file." , &
                ADJUSTL(TRIM(FILENAME)), ". Please check if the file exists", ERR
            CLOSE(IUNIT)
            STOP
        ENDIF
        FTYPE = ADJUSTL(TRIM(FILETYPE))
        CALL TO_UPPER(FTYPE)
        SELECT CASE (FTYPE)
            CASE ("TXT")
                DO I=1, N
                    READ(IUNIT,*) COORDS(I, :)
                ENDDO
                CLOSE(IUNIT)
            CASE ("XYZ")
                READ(IUNIT, *) J
                IF(J /= N) WRITE(*, *) "WARNING! number of atoms is not equal to number of data rows"
                READ(IUNIT, *)
                DO I=1, N
                    PRINT*, I
                    READ(IUNIT, *) ATOMIC_SYMBOL(I), COORDS(I, :)
                ENDDO
            CASE DEFAULT
                PRINT*, "File type not recognized."
        END SELECT
    END SUBROUTINE READCOORDINATES

    SUBROUTINE To_upper(str)
        CHARACTER(*), INTENT(IN OUT) :: str
        INTEGER :: i
        DO i = 1, len(str)
            SELECT CASE(str(i:i))
                CASE("a":"z")
                    str(i:i) = achar(iachar(str(i:i))-32)
            END SELECT
        END DO
    END SUBROUTINE To_upper
    SUBROUTINE To_lower(str)
        CHARACTER(*), INTENT(IN OUT) :: str
        INTEGER :: i
        DO i = 1, len(str)
            SELECT CASE(str(i:i))
                CASE("A":"Z")
                    str(i:i) = achar(iachar(str(i:i))+32)
            END SELECT
        END DO
    END SUBROUTINE To_Lower

    SUBROUTINE GETENERGY(PIXLECOORDS, ATOMCOORDS, ATOMZ, BETA, SIGMA, ENERGY)
        IMPLICIT NONE
        REAL(KIND=DP), INTENT(IN), DIMENSION(:) :: PIXLECOORDS
        REAL(KIND=DP), INTENT(IN), DIMENSION(:) :: ATOMCOORDS
        REAL(KIND=DP), INTENT(IN) :: BETA, SIGMA
        REAL(KIND=DP), INTENT(IN) :: ATOMZ
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
        REAL(KIND=DP) :: DX, DY
        DX = PIXLECOORDS(1) - ATOMCOORDS(1)
        DY = PIXLECOORDS(2) - ATOMCOORDS(2)
        ENERGY = (ATOMZ**BETA)*EXP(-(DX*DX + DY*DY)/(SIGMA*SIGMA))
    END SUBROUTINE GETENERGY

    SUBROUTINE GETGRADIENT(PIXLECOORDS, ATOMCOORDS, ATOMZ, BETA, SIGMA, ENERGY, GRAD)
        IMPLICIT NONE
        REAL(KIND=DP), INTENT(IN), DIMENSION(:) :: PIXLECOORDS
        REAL(KIND=DP), INTENT(IN), DIMENSION(:) :: ATOMCOORDS
        REAL(KIND=DP), INTENT(IN) :: ATOMZ, BETA, SIGMA
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
        REAL(KIND=DP), INTENT(OUT), DIMENSION(:) :: GRAD
        REAL(KIND=DP) :: DFDX, DFDY, DFDZ, DFDB, DFDS, SIGMA2, SIGMA3
        CALL GETENERGY(PIXLECOORDS, ATOMCOORDS, ATOMZ, BETA, SIGMA, ENERGY)
        SIGMA2 = SIGMA*SIGMA
        SIGMA3= SIGMA*SIGMA*SIGMA
        DFDX = 2.0D0*((PIXLECOORDS(1) - ATOMCOORDS(1))/SIGMA2) * ENERGY
        DFDY = 2.0D0*((PIXLECOORDS(2) - ATOMCOORDS(2))/SIGMA2) * ENERGY
        DFDZ = 0.0D0
        DFDS = (1.0D0/SIGMA3) * ENERGY
        DFDB = LOG(BETA) * ENERGY
        GRAD(1) = DFDX
        GRAD(2) = DFDY
        GRAD(3) = DFDS
        GRAD(4) = DFDB
    END SUBROUTINE GETGRADIENT

    !
    !    SUBROUTINE GETHESSIAN()
    !
    !    END SUBROUTINE GETHESSIAN

    SUBROUTINE POTENTIAL(NATOMS, COORDS, ENERGY, GRAD, GTEST, HTEST)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NATOMS
        REAL(KIND=DP), INTENT(IN), DIMENSION(3*NATOMS) :: COORDS
        REAL(KIND=DP), INTENT(OUT), DIMENSION(3*NATOMS) :: GRAD
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
!        CHARACTER(LEN=5), INTENT(IN) :: POTTYPE
        CHARACTER(LEN=5) :: POT
        LOGICAL, INTENT(IN) :: GTEST, HTEST

        INTERFACE
            SUBROUTINE GUPTA(NATOMS,COORDS, ENERGY, GRAD, GTEST, ATOMTYPE)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: NATOMS, ATOMTYPE
                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: COORDS
                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: GRAD
                REAL(KIND=8), INTENT(OUT) :: ENERGY
                LOGICAL, INTENT(IN) :: GTEST
            END SUBROUTINE GUPTA

            SUBROUTINE MORSEGH(NATOMS, COORDS, ENERGY, GRAD, RHO, GTEST, HTEST)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: NATOMS
                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: COORDS
                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: GRAD
                REAL(KIND=8), INTENT(IN) :: RHO
                REAL(KIND=8), INTENT(OUT) :: ENERGY
                LOGICAL, INTENT(IN) :: GTEST, HTEST
            END SUBROUTINE MORSEGH

            SUBROUTINE LJGH(NATOMS, COORDS, ENERGY, GRAD, GTEST, HTEST)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: NATOMS
                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: COORDS
                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: GRAD
                REAL(KIND=8), INTENT(OUT) :: ENERGY
                LOGICAL, INTENT(IN) :: GTEST, HTEST
            END SUBROUTINE LJGH
        END INTERFACE
        POT = POTTYPE
        CALL TO_UPPER(POT)
        POT = ADJUSTL(TRIM(POT))
        SELECT CASE(POT)
            CASE("GUPTA")
                CALL GUPTA(NATOMS,COORDS, ENERGY, GRAD, GTEST, ATOMTYPE)
            CASE("MORSE")
                CALL MORSEGH(NATOMS, COORDS, ENERGY, GRAD, RHO, GTEST, HTEST)
            CASE("LJ")
                CALL LJGH(NATOMS, COORDS, ENERGY, GRAD, GTEST, HTEST)
            CASE DEFAULT
                PRINT*, "POTENTIAL NOT IMPLEMENTED"
        END SELECT
    END SUBROUTINE POTENTIAL
END MODULE ACSTEM2



!    SUBROUTINE SETBOUND()
!        IMPLICIT NONE
!        INTEGER :: I
!        DO I=1, NATOMS
!            LB(3*I - 2) = -IMAGE_WIDTH/2.0D0
!            LB(3*I - 1) = -IMAGE_HEIGHT/2.0D0
!            LB(3*I)     = -MAX(IMAGE_WIDTH, IMAGE_HEIGHT)
!            UB(3*I - 2) = IMAGE_WIDTH/2.0D0
!            UB(3*I - 1) = IMAGE_HEIGHT/2.0D0
!            UB(3*I)     = MAX(IMAGE_WIDTH, IMAGE_HEIGHT)
!        ENDDO
!        ! BOUNDS ON GAUSSAIN WIDTH
!        LB(ND - 1)  = 0.2D0*GAUSS_WIDTH_DEFAULT
!        UB(ND - 1)  = GAUSS_WIDTH_DEFAULT + 0.2D0*GAUSS_WIDTH_DEFAULT
!        ! BOUNDS ON GAUSSIAN HEIGHT
!        LB(ND)      = 1.0D0
!        UB(ND)      = 3.0D0
!    END SUBROUTINE SETBOUND


!        INTERFACE
!            SUBROUTINE GUPTA(NATOMS,X2, E2, G2, GTEST, ATOMTYPE)
!                IMPLICIT NONE
!                INTEGER, INTENT(IN) :: NATOMS, ATOMTYPE
!                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: X2
!                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: G2
!                REAL(KIND=8), INTENT(OUT) :: E2
!                LOGICAL, INTENT(IN) :: GTEST
!            END SUBROUTINE GUPTA
!
!            SUBROUTINE MORSEGH(NATOMS, X2, E2, G2, RHO, GTEST, HTEST)
!                IMPLICIT NONE
!                INTEGER, INTENT(IN) :: NATOMS
!                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: X2
!                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: G2
!                REAL(KIND=8), INTENT(IN) :: RHO
!                REAL(KIND=8), INTENT(OUT) :: E2
!                LOGICAL, INTENT(IN) :: GTEST, HTEST
!            END SUBROUTINE MORSEGH
!
!            SUBROUTINE LJGH(NATOMS, X2, E2, G2, GTEST, HTEST)
!                IMPLICIT NONE
!                INTEGER, INTENT(IN) :: NATOMS
!                REAL(KIND=8), INTENT(IN), DIMENSION(NATOMS) :: X2
!                REAL(KIND=8), INTENT(OUT), DIMENSION(NATOMS) :: G2
!                REAL(KIND=8), INTENT(OUT) :: E2
!                LOGICAL, INTENT(IN) :: GTEST, HTEST
!            END SUBROUTINE LJGH
!        END INTERFACE

!        CALL GUPTA(NATOMS,X2,E2,G2,GTEST,ATOMTYPE)


!    REAL(KIND=DP), PARAMETER, PRIVATE       :: HEIGHT_DEFAULT  = ATOMIC_NUMBER**2


!    REAL(KIND=DP), PARAMETER :: SIGMA   = 0.25
!    REAL(KIND=DP), PARAMETER :: SIGMAX  = SIGMA
!    REAL(KIND=DP), PARAMETER :: SIGMAY  = SIGMA
!    REAL(KIND=DP), PARAMETER :: SIGMAX2 = SIGMAX*SIGMAX
!    REAL(KIND=DP), PARAMETER :: SIGMAY2 = SIGMAY*SIGMAY
!    REAL(KIND=DP), PARAMETER :: SIGMAX4 = SIGMAX2*SIGMAX2
!    REAL(KIND=DP), PARAMETER :: SIGMAY4 = SIGMAY2*SIGMAY2
!    REAL(KIND=DP), PARAMETER :: SIGMAX2SIGMAY2 = SIGMAX2*SIGMAY2

