MODULE OBJECTIVE

    IMPLICIT NONE
    INTEGER, PARAMETER, PRIVATE :: dp       = KIND(1.0D0)
    INTEGER, PARAMETER, PUBLIC  :: numatoms = 55        !309   !55
    INTEGER, PARAMETER, PUBLIC  :: numdims  = 3*numatoms + 2
    INTEGER, PARAMETER, PUBLIC  :: func_num = 1
    INTEGER, PARAMETER, PUBLIC  :: NUM_OBJFUNC = 2
    REAL(KIND=DP), DIMENSION(NUM_OBJFUNC), SAVE :: WEIGHTS
    REAL(KIND=DP), PARAMETER :: LAMBDA = 1.0E-4
    CHARACTER(LEN=10),  PUBLIC, PARAMETER :: func_name = "ACSTEM"
    REAL(KIND=dp), PARAMETER, PUBLIC :: VTR = 0.0D0 ! VALUE TO REACH
    REAL(KIND=dp), DIMENSION(numdims), PROTECTED, SAVE, PUBLIC :: lowerbound
    REAL(KIND=dp), DIMENSION(numdims), PROTECTED, SAVE, PUBLIC :: upperbound
    REAL(KIND=dp), DIMENSION(numatoms), PUBLIC, SAVE :: interatomic_distance
!    CHARACTER(LEN=100), PARAMETER, PRIVATE :: experimental_image_file = "img.txt" 
    CHARACTER(LEN=100), PARAMETER, PRIVATE :: experimental_image_file = "Au_20000000X_0020.txt"
    CHARACTER(LEN=255), SAVE, PRIVATE :: MYRESULTSPATH=""

    INTEGER, PARAMETER, PUBLIC :: numrows  = 128        !50
    INTEGER, PARAMETER, PUBLIC :: numcols  = 128        !50
    REAL(KIND=dp), PARAMETER, PRIVATE  :: image_width = 20.0D0         ! IN ANGSTROM
    REAL(KIND=dp), PARAMETER, PRIVATE  :: image_height= 20.0D0         ! IN ANGSTROM

    REAL(KIND=dp), PARAMETER, PRIVATE  :: lowerbound_x = 0.0D0
    REAL(KIND=dp), PARAMETER, PRIVATE  :: upperbound_x = image_width
    REAL(KIND=dp), PARAMETER, PRIVATE  :: lowerbound_y = 0.0D0
    REAL(KIND=dp), PARAMETER, PRIVATE  :: upperbound_y = image_height
    REAL(KIND=dp), PARAMETER, PRIVATE  :: lowerbound_z = 0.0d0  !MIN(lowerbound_x, lowerbound_y)
    REAL(KIND=dp), PARAMETER, PRIVATE  :: upperbound_z = MAX(upperbound_x, upperbound_y)
    REAL(KIND=dp), PARAMETER, PRIVATE  :: lowerbound_gausswidth = 1.0D0
    REAL(KIND=dp), PARAMETER, PRIVATE  :: upperbound_gausswidth = 2.0D0
    REAL(KIND=dp), PARAMETER, PRIVATE  :: lowerbound_gaussheight = 0.0D0
    REAL(KIND=dp), PARAMETER, PRIVATE  :: upperbound_gaussheight = 2.0D0
    REAL(KIND=dp), DIMENSION(numdims, numdims), SAVE :: hessian
    REAL(KIND=dp), DIMENSION(numrows, numcols), PUBLIC, SAVE :: simulated_image
    REAL(KIND=dp), DIMENSION(numrows, numcols), PRIVATE , SAVE :: experimental_image

    REAL(KIND=DP), DIMENSION(NUM_OBJFUNC), SAVE, PRIVATE :: MY_ENERGIES
    INTEGER, SAVE, PRIVATE :: NUM_ACSTEM_CALLS, THIS_RUN=0, THIS_GENERATION=0, ISTORE
    ! POTENTIAL PARAMETERS
    INTEGER, PARAMETER, PRIVATE             :: atomtype         = 23    ! GARZON
    INTEGER, PARAMETER, PRIVATE             :: atomic_number    = 79
    REAL(KIND=dp), PARAMETER, PRIVATE       :: morse_rho        = 1.0D0
    CHARACTER(LEN=5), PARAMETER, PRIVATE    :: potential_type   = "GUPTA"

    LOGICAL, SAVE, PRIVATE  ::  init=.TRUE.
    LOGICAL, SAVE, PRIVATE  ::  debug_acstem=.TRUE.

    CONTAINS
    SUBROUTINE SETUP_OBJECTIVE()
        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: CURRENT_RUN
!        CHARACTER(LEN=*), INTENT(IN) :: WHERE_TO_SAVE
        INTEGER :: i, myunit, seedsize
        PRINT*, "Setting the upper and lower bounds of the search space ...."
!        THIS_RUN = THIS_RUN + 1
!        THIS_RUN = CURRENT_RUN
        DO i=1, numatoms
            lowerbound(3*i - 2) = lowerbound_x
            lowerbound(3*i - 1) = lowerbound_y
            lowerbound(3*i)     = lowerbound_z
            upperbound(3*i - 2) = upperbound_x
            upperbound(3*i - 1) = upperbound_y
            upperbound(3*i)     = upperbound_z
        ENDDO
        ! BOUNDS ON GAUSSAIN WIDTH
        lowerbound(numdims - 1)  = lowerbound_gausswidth
        upperbound(numdims - 1)  = upperbound_gausswidth
        ! BOUNDS ON GAUSSIAN HEIGHT
        lowerbound(numdims)      = lowerbound_gaussheight
        upperbound(numdims)      = upperbound_gaussheight
        PRINT*, "Search bounds set successfully."
        PRINT*, "Loading experimental image .... "
        IF(init) THEN
            ! LOAD EXPERIMENTAL IMAGE
            CALL LOADIMAGE(experimental_image_file, numrows, numcols, experimental_image)
            experimental_image = experimental_image/MAXVAL(experimental_image)
            init = .FALSE.
            PRINT*, "Experimental image loaded successfully."
        ENDIF
    END SUBROUTINE SETUP_OBJECTIVE

    SUBROUTINE OBJFUNC(nd, COORDS, ENERGY, GRAD, gtest, htest)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ND
        LOGICAL, INTENT(IN) :: gtest, htest
        REAL(KIND=8), INTENT(IN),  DIMENSION(nd):: COORDS
        REAL(KIND=8), INTENT(OUT)  :: ENERGY
        REAL(KIND=8), INTENT(OUT), DIMENSION(nd):: GRAD
        REAL(KIND=8) :: STEM_ENERGY, POT_ENERGY
        REAL(KIND=8), DIMENSION(ND) :: GRAD_STEM, GRAD_POT
        INTEGER, SAVE:: N3, INIT=.TRUE.
        IF(INIT) THEN
            CALL SETUP_OBJECTIVE()
            INIT = .FALSE.
        ENDIF
        GRAD_POT = 0.0D0
        N3 = 3*NUMATOMS
        CALL MY_ACSTEM_MODEL(nd, COORDS, STEM_ENERGY, GRAD_STEM, gtest, htest)
        CALL POTENTIAL(numatoms, COORDS(1:3*NUMATOMS), POT_ENERGY, GRAD_POT(1:N3), gtest, htest)
        NUM_ACSTEM_CALLS = NUM_ACSTEM_CALLS + 1
        WEIGHTS(2) = EXP(-LAMBDA*DBLE(NUM_ACSTEM_CALLS))
        ENERGY = STEM_ENERGY + WEIGHTS(2)*POT_ENERGY !+ WEIGHT(3)*NOISE
        GRAD = GRAD_STEM + WEIGHTS(2)*GRAD_POT
    END SUBROUTINE OBJFUNC

    SUBROUTINE MY_ACSTEM_MODEL(nd, xstem, estem, gstem, gtest, htest)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nd
        REAL(KIND=dp), INTENT(IN), DIMENSION(nd)  :: xstem
        REAL(KIND=dp), INTENT(OUT), DIMENSION(nd) :: gstem
        REAL(KIND=dp), INTENT(OUT)  :: estem
        LOGICAL, INTENT(IN)         :: gtest, htest

        REAL(KIND=dp), DIMENSION(numatoms) :: imij, dx, dy
        REAL(KIND=dp) :: xstep, ystep, height, alpha, sigma, beta
        REAL(KIND=dp) :: xi, yj, errij, tmp, dhdb
        INTEGER :: i, j, k

        xstep = (upperbound_x - lowerbound_x)/DBLE(numrows)
        ystep = (upperbound_y - lowerbound_y)/DBLE(numcols)
        estem = 0.0D0
        gstem = 0.0D0
        sigma = XSTEM(nd-1)
        alpha = 1.0D0/(2.0D0*sigma*sigma)
        beta = XSTEM(nd)
        height = (DBLE(atomic_number))**beta
        dhdb = height*LOG(DBLE(atomic_number))
        IF (gtest) THEN
            DO j=1, numcols
                yj = lowerbound_y + ystep*DBLE(j)
                DO i=1, numrows
                    xi = lowerbound_x + xstep*DBLE(i)
                    DO k= 1, numatoms
                        dx(k) = xi - XSTEM(3*k-2)
                        dy(k) = yj - XSTEM(3*k-1)
                        imij(k) = height*EXP(-alpha*(DX(k)**2 + DY(k)**2))
                    ENDDO
                    simulated_image(i,j) = SUM(imij)
                    errij = EXPERIMENTAL_IMAGE(i,j) - SIMULATED_IMAGE(i,j)
                    estem = estem + errij*errij
                    tmp = -4.0D0*alpha*errij
                    DO k=1, numatoms
                        gstem(3*k - 2) = GSTEM(3*k - 2) + tmp*DX(k)*IMIJ(k)
                        gstem(3*k - 1) = GSTEM(3*k - 1) + tmp*DY(k)*IMIJ(k)
                        gstem(nd-1)    = GSTEM(nd-1) - 2.0D0*errij*(DX(k)**2 + DY(k)**2)*IMIJ(k)*alpha/sigma
                    ENDDO
                    gstem(nd)   = GSTEM(nd) - 2.0*errij*SIMULATED_IMAGE(i,j)*dhdb/height
                ENDDO
            ENDDO
        ELSE IF (htest) THEN
            hessian = 0.0D0
            PRINT*, "The second derivative function for stem model image is not implemented."
        ELSE
            DO j=1, numcols
                yj = lowerbound_y + ystep*DBLE(j)
                DO i=1, numrows
                    xi = lowerbound_x + xstep*DBLE(i)
                    DO k= 1, numatoms
                        dx(k) = xi - XSTEM(3*k-2)
                        dy(k) = yj - XSTEM(3*k-1)
                        imij(k) = height*EXP(-alpha*(DX(k)**2 + DY(k)**2))
                    ENDDO
                    simulated_image(i,j) = SUM(imij)
                    errij = EXPERIMENTAL_IMAGE(i,j) - SIMULATED_IMAGE(i,j)
                    estem = estem + errij*errij
                ENDDO
            ENDDO
        ENDIF
    END SUBROUTINE MY_ACSTEM_MODEL

    SUBROUTINE MY_IMAGE(n, coords, nr, nc, sigma, beta, image)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, nr, nc
        REAL(KIND=dp), INTENT(IN) :: sigma, beta
        REAL(KIND=dp), DIMENSION(n,3), INTENT(IN) :: coords
        REAL(KIND=dp), DIMENSION(nr,nc), INTENT(OUT) :: image
        REAL(KIND=dp) :: dx, dy, xi, yj, height, alpha, im(n), xstep, ystep
        INTEGER :: i, j, k
        xstep = (upperbound_x - lowerbound_x)/nc
        ystep = (upperbound_y - lowerbound_y)/nr
        height = (DBLE(atomic_number))**beta
        alpha = 1.0D0/(2.0D0*sigma*sigma)
        DO j=1, nc
            yj = lowerbound_y + DBLE(j)*ystep
            DO i=1, nr
                xi = lowerbound_x + DBLE(i)*xstep
                DO k=1, n
                    dx = xi - COORDS(k,1)
                    dy = yj - COORDS(k,2)
                    im(k) = height*EXP(-alpha*(dx*dx + dy*dy))
                ENDDO
                image(i,j) = SUM(im)
            ENDDO
        ENDDO
    END SUBROUTINE MY_IMAGE

    SUBROUTINE LOADIMAGE(filename, nrows, ncols, image)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nrows, ncols
        REAL(8), INTENT(OUT), DIMENSION(nrows, ncols) :: image
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: i, j, iunit, err
        OPEN(NEWUNIT=iunit, FILE=adjustl(TRIM(filename)), STATUS="OLD", ACTION="READ", IOSTAT=ERR)
        IF (err /= 0) THEN
            PRINT*, "ERROR! Problem opening file " , ADJUSTL(TRIM(filename)), ". Please check if the file exists"
            CLOSE(iunit)
            STOP
        ELSE
            DO i=1, nrows
                !PRINT*, "I= ", I
                READ(iunit,*) IMAGE(i, :)
            ENDDO

            CLOSE(iunit)
        ENDIF
    END SUBROUTINE LOADIMAGE

    SUBROUTINE WRITEIMAGE(filename, nrows, ncols, image)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nrows, ncols
        REAL(8), INTENT(IN), DIMENSION(nrows, ncols) :: image
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: i, j, iunit, err
        OPEN(NEWUNIT=iunit, FILE=adjustl(TRIM(filename)), STATUS="UNKNOWN", ACTION="WRITE", IOSTAT=ERR)
        IF (err /= 0) THEN
            PRINT*, "ERROR! Problem opening file for writing " , &
                ADJUSTL(TRIM(filename)), ". Please check if the file already exists"
            CLOSE(iunit)
            STOP
        ELSE
            DO i=1, nrows
                WRITE(iunit,*) IMAGE(i, :)
            ENDDO
            CLOSE(iunit)
        ENDIF
    END SUBROUTINE WRITEIMAGE

    SUBROUTINE READCOORDINATES(filename, natoms, filetype, coords)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: natoms
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=*), INTENT(IN) :: filetype
        REAL(KIND=8), INTENT(OUT), DIMENSION(natoms, 3) :: coords
        CHARACTER(LEN=4) :: ftype
        CHARACTER(LEN=2), DIMENSION(natoms) :: atomic_symbol
        INTEGER :: i, j, iunit, err
        OPEN(NEWUNIT=iunit, FILE=adjustl(TRIM(filename)), STATUS="OLD", ACTION="READ", IOSTAT=ERR)
        IF (err /= 0) THEN
            PRINT*, "ERROR! Problem opening file." , &
                ADJUSTL(TRIM(filename)), ". Please check if the file exists", err
            CLOSE(iunit)
            STOP
        ENDIF
        ftype = ADJUSTL(TRIM(filetype))
        CALL TO_UPPER(ftype)
        SELECT CASE (ftype)
            CASE ("TXT")
                DO i=1, natoms
                    READ(iunit,*) COORDS(i, :)
                ENDDO
                CLOSE(iunit)
            CASE ("XYZ")
                READ(iunit, *) J
                IF(j /= natoms) WRITE(*, *) "WARNING! number of atoms is not equal to number of data rows"
                READ(iunit, *)
                DO i=1, natoms
                    PRINT*, i
                    READ(iunit, *) ATOMIC_SYMBOL(i), COORDS(i, :)
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

    SUBROUTINE GETENERGY(pixlecoords, atomcoords, atomz, beta, sigma, energy)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: pixlecoords
        REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: atomcoords
        REAL(KIND=dp), INTENT(IN) :: beta, sigma
        REAL(KIND=dp), INTENT(IN) :: atomz
        REAL(KIND=dp), INTENT(OUT) :: energy
        REAL(KIND=dp) :: dx, dy
        dx = PIXLECOORDS(1) - ATOMCOORDS(1)
        dy = PIXLECOORDS(2) - ATOMCOORDS(2)
        energy = (atomz**beta)*EXP(-(dx*dx + dy*dy)/(sigma*sigma))
    END SUBROUTINE GETENERGY

    SUBROUTINE GETGRADIENT(pixlecoords, atomcoords, atomz, beta, sigma, energy, grad)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: pixlecoords
        REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: atomcoords
        REAL(KIND=dp), INTENT(IN) :: atomz, beta, sigma
        REAL(KIND=dp), INTENT(OUT) :: energy
        REAL(KIND=dp), INTENT(OUT), DIMENSION(:) :: grad
        REAL(KIND=dp) :: dfdx, dfdy, dfdz, dfdb, dfds, sigma2, sigma3
        CALL GETENERGY(pixlecoords, atomcoords, atomz, beta, sigma, energy)
        sigma2 = sigma*sigma
        sigma3= sigma*sigma*sigma
        dfdx = 2.0D0*((PIXLECOORDS(1) - ATOMCOORDS(1))/sigma2) * energy
        dfdy = 2.0D0*((PIXLECOORDS(2) - ATOMCOORDS(2))/sigma2) * energy
        dfdz = 0.0D0
        dfds = (1.0D0/sigma3) * energy
        dfdb = LOG(beta) * energy
        grad(1) = dfdx
        grad(2) = dfdy
        grad(3) = dfds
        grad(4) = dfdb
    END SUBROUTINE GETGRADIENT

    SUBROUTINE POTENTIAL(natoms, coords, energy, grad, gtest, htest)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: natoms
        REAL(KIND=dp), INTENT(IN), DIMENSION(3*natoms) :: coords
        REAL(KIND=dp), INTENT(OUT), DIMENSION(3*natoms) :: grad
        REAL(KIND=dp), INTENT(OUT) :: energy
!        CHARACTER(LEN=5), INTENT(IN) :: POTTYPE
        CHARACTER(LEN=5) :: pot
        LOGICAL, INTENT(IN) :: gtest, htest

        INTERFACE
            SUBROUTINE GUPTA(natoms,coords, energy, grad, gtest, atomtype)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: natoms, atomtype
                REAL(KIND=8), INTENT(IN), DIMENSION(natoms) :: coords
                REAL(KIND=8), INTENT(OUT), DIMENSION(natoms) :: grad
                REAL(KIND=8), INTENT(OUT) :: energy
                LOGICAL, INTENT(IN) :: gtest
            END SUBROUTINE GUPTA

            SUBROUTINE MORSEGH(natoms, coords, energy, grad, morse_rho, gtest, htest)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: natoms
                REAL(KIND=8), INTENT(IN), DIMENSION(natoms) :: coords
                REAL(KIND=8), INTENT(OUT), DIMENSION(natoms) :: grad
                REAL(KIND=8), INTENT(IN) :: morse_rho
                REAL(KIND=8), INTENT(OUT) :: energy
                LOGICAL, INTENT(IN) :: gtest, htest
            END SUBROUTINE MORSEGH

            SUBROUTINE LJGH(natoms, coords, energy, grad, gtest, htest)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: natoms
                REAL(KIND=8), INTENT(IN), DIMENSION(natoms) :: coords
                REAL(KIND=8), INTENT(OUT), DIMENSION(natoms) :: grad
                REAL(KIND=8), INTENT(OUT) :: energy
                LOGICAL, INTENT(IN) :: gtest, htest
            END SUBROUTINE LJGH
        END INTERFACE
        pot = potential_type
        CALL TO_UPPER(pot)
        pot = ADJUSTL(TRIM(pot))
        SELECT CASE(pot)
            CASE("GUPTA")
                CALL GUPTA(natoms,coords, energy, grad, gtest, atomtype)
            CASE("MORSE")
                CALL MORSEGH(natoms, coords, energy, grad, morse_rho, gtest, htest)
            CASE("LJ")
                CALL LJGH(natoms, coords, energy, grad, gtest, htest)
            CASE DEFAULT
                PRINT*, "POTENTIAL NOT IMPLEMENTED"
        END SELECT
    END SUBROUTINE POTENTIAL


      SUBROUTINE CHECK_AND_SAVE_OBJECTIVE_SETUP(MYUNIT, CURRENT_RUN)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: MYUNIT, CURRENT_RUN
        LOGICAL, SAVE :: setup_error = .FALSE.
        IF (numdims /= 3*numatoms+2) THEN
            PRINT*, "WARNING! Number of atoms and number of dimensions not matching."
            PRINT*, "NUMDIMS = ", numdims
            PRINT*, "NUMATOMS = ", numatoms
        ENDIF
        IF (numrows == 0) THEN
            PRINT*, "ERRO!  IMAGE CAN'T HAVE ZERO ROWS"
            setup_error = .TRUE.
        ELSEIF (numcols == 0) THEN
            PRINT*, "ERRO!  IMAGE CAN'T HAVE ZERO ROWS"
            setup_error = .TRUE.
        ENDIF
        IF (image_width == 0 .OR. image_width < 0) THEN
            PRINT*, "IMAGE_WIDTH CAN'T BE ZERO OR NEGATIVE."
            setup_error = .TRUE.
        ELSEIF(image_height == 0 .OR. image_height < 0) THEN
            PRINT*, "IMAGE_HEIGHT CAN'T BE ZERO OR NAGATIVE"
            setup_error = .TRUE.
        ENDIF
        IF (ANY(LOWERBOUND(1:3*numatoms:3) /= lowerbound_x)) THEN
            PRINT*, "ERROR! LOWERBOUND NOT SET PROPERLY FOR X COORDINATES"
            setup_error = .TRUE.
        ELSEIF(ANY(LOWERBOUND(2:3*numatoms:3) /= lowerbound_y)) THEN
            PRINT*, "ERROR! LOWERBOUND NOT SET PROPERLY FOR Y COORDINATES"
            setup_error = .TRUE.
        ELSEIF (ANY(LOWERBOUND(3:3*numatoms:3) /= lowerbound_z)) THEN
            PRINT*, "WARNING! LOWERBOUND NOT SET PROPERLY FOR Z COORDINATES"
            setup_error = .TRUE.
        ELSEIF (ANY(UPPERBOUND(1:3*numatoms:3) /= upperbound_x)) THEN
            PRINT*, "ERROR! UPPERBOUND NOT SET PROPERLY FOR X COORDINATES"
            setup_error = .TRUE.
        ELSEIF(ANY(UPPERBOUND(2:3*numatoms:3) /= upperbound_y)) THEN
            PRINT*, "ERROR! UPPERBOUND NOT SET PROPERLY FOR Y COORDINATES"
            setup_error = .TRUE.
        ELSEIF (ANY(UPPERBOUND(3:3*numatoms:3) /= upperbound_z)) THEN
            PRINT*, "WARNING! UPPERBOUND NOT SET PROPERLY FOR Z COORDINATES"
            setup_error = .TRUE.
        ENDIF
        WRITE(myunit, "(/)")
        WRITE(myunit, "(2x, A)") repeat("=", 72)
        WRITE(myunit, "(2x, 3A)") "=",repeat(" ", 70), "="
        WRITE(myunit, "(2x, 5A)") "=",repeat(" ", 19), "SETUP OF THE OBJECTIVE FUNCTION", repeat(" ", 20), "="
        WRITE(myunit, "(2x, 3A)") "=",repeat(" ", 70), "="
        WRITE(myunit, "(2x, A)") repeat("=", 72)
            WRITE(myunit, "(2x, A, I8)")    "Number of atoms                                   :   ", numatoms
            WRITE(myunit, "(2x, A, I8)")    "Number of additional tunable parameters           :   ", 2
            WRITE(myunit, "(2x, A, I8)")    "Total Number of Variables = 3*natoms + tunables   :   ", numdims
            WRITE(MYUNIT, "(2x, A, I8)")    "Number of objective functions, in case of MOO     :   ",  2
            WRITE(myunit, "(2x, A, f8.3)")  "Image width in Angstrom                           :   ", image_width
            WRITE(myunit, "(2x, A, f8.3)")  "Image height in Angstrom                          :   ", image_width
            WRITE(myunit, "(2x, A, I8)")    "Number of rows in the image                       :   ", numrows
            WRITE(myunit, "(2x, A, I8)")    "Number of columns in the image                    :   ", numcols
            WRITE(myunit, "(2x, A, f8.3)")  "Lower bound on X-coordinates                      :   ", lowerbound_x
            WRITE(myunit, "(2x, A, f8.3)")  "Lower bound on Y-coordinates                      :   ", lowerbound_y
            WRITE(myunit, "(2x, A, f8.3)")  "Lower bound on Z-coordinates                      :   ", lowerbound_z
            WRITE(myunit, "(2x, A, f8.3)")  "Upper bound on X-coordinates                      :   ", upperbound_x
            WRITE(myunit, "(2x, A, f8.3)")  "Upper bound on Y-coordinates                      :   ", upperbound_y
            WRITE(myunit, "(2x, A, f8.3)")  "Upper bound on Z-coordinates                      :   ", upperbound_z
            WRITE(myunit, "(2x, A, f8.3)")  "Lower bound on Gaussian width                     :   ", lowerbound_gausswidth
            WRITE(myunit, "(2x, A, f8.3)")  "Lower bound on Gaussian height power              :   ", lowerbound_gaussheight
            WRITE(myunit, "(2x, A, f8.3)")  "Upper bound on Gaussian width                     :   ", upperbound_gausswidth
            WRITE(myunit, "(2x, A, f8.3)")  "Upper bound on Gaussian height power              :   ", upperbound_gaussheight
            WRITE(myunit, "(2x, A, a)")     "Type of interactive potential used                :   ", potential_type
            WRITE(myunit, "(2x, A, I8)")    "INDEX of atoms in the nanocluster                 :   ", atomtype
            WRITE(myunit, "(2x, A, I8)")    "Atomic numbers of the atom                        :   ", atomic_number
            WRITE(myunit, "(2x, A, f8.3)")  "Decay constant of the weight function (lambda)    :   ", lambda
            WRITE(myunit, "(2x, A, A)")     "Experimental image used in the simulation         :   ", experimental_image_file
            WRITE(myunit, "(/)")
        IF(setup_error) THEN
            PRINT*, "PLEASE CHECK THE SETUP FILE IN THE CURRENT DIERECTORY."
            STOP
        ENDIF
    END SUBROUTINE  CHECK_AND_SAVE_OBJECTIVE_SETUP

END MODULE OBJECTIVE





























!    SUBROUTINE SAVE_ACSTEM_RESULTS(WHERE_TO_SAVE, CURRENT_RUN, CURRENT_GENERATION)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: CURRENT_RUN, CURRENT_GENERATION
!        CHARACTER(LEN=*), INTENT(IN) :: WHERE_TO_SAVE
!        INTEGER, SAVE :: I, MYUNIT, IRUN
!        LOGICAL :: MY_INIT=.TRUE., NEW_RUN=.FALSE.
!        CHARACTER(LEN=255), SAVE :: MY_FILE_NAME
!        IF(MY_INIT) THEN
!            IRUN = CURRENT_RUN
!            NEW_RUN = .TRUE.
!            MY_INIT = .FALSE.
!        ELSEIF(CURRENT_RUN == IRUN) THEN
!            NEW_RUN = .FALSE.
!        ELSEIF(CURRENT_RUN > IRUN) THEN
!            NEW_RUN = .TRUE.
!            IRUN = CURRENT_RUN
!        ELSE
!            PRINT*,"ERROR IN SAVING ACSTEM RESULTS"
!            STOP
!        ENDIF
!
!        IF(ADJUSTL(TRIM(FUNC_NAME)) == "") THEN
!            WRITE(MY_FILE_NAME, "(a,3(i0,a))") "ObjValues_F",FUNC_NUM,"_",NUMDIMS,"D_R", CURRENT_RUN, ".txt"
!        ELSE
!            WRITE(MY_FILE_NAME, "(3a,2(i0,a))") "ObjValues_",&
!                ADJUSTL(TRIM(FUNC_NAME)),"_",NUMDIMS,"D_R", CURRENT_RUN, "_out.txt"
!        ENDIF
!        MY_FILE_NAME = ADJUSTL(TRIM(WHERE_TO_SAVE))//"/"//ADJUSTL(TRIM(MY_FILE_NAME))
!        IF(NEW_RUN) THEN
!            OPEN(NEWUNIT=MYUNIT, FILE=ADJUSTL(TRIM(MY_FILE_NAME)), STATUS="REPLACE", ACTION="WRITE")
!            !!!WRITE(MYUNIT, *) "STMCALLS", "Σ(I0-I)^2", "Σ(V(R))", "WEIGHT", "Σ(I0-I)^2 + WEIGHT * Σ(V(R))"
!            WRITE(MYUNIT, *) "First column is the EA generation."
!            WRITE(MYUNIT, *) "Second column is the number of calls to ACSTEM."
!            WRITE(MYUNIT, *) "Third column is the value of f1 (sum of squared error between true and simulated images)."
!            WRITE(MYUNIT, *) "Forth column is the value of f2 (potential energy of the resulting cluster). "
!            WRITE(MYUNIT, *) "Fifth column is the weight by which f2 is multiplied."
!            WRITE(MYUNIT, *) "Sixth column is the value of the objective function (f1 + weight*f2)."
!            WRITE(MYUNIT, *)
!        ELSE
!            OPEN(NEWUNIT=MYUNIT, FILE=ADJUSTL(TRIM(MY_FILE_NAME)), ACCESS="APPEND", ACTION="WRITE")
!        ENDIF
!        WRITE(MYUNIT, *) CURRENT_GENERATION, NUM_ACSTEM_CALLS, MY_ENERGIES(1), MY_ENERGIES(2), WEIGHTS(2), SUM(MY_ENERGIES)
!!        WRITE(MYUNIT, *) "ACSTEM CALLS                                              =   ",MY_ENERGIES(1)
!!        WRITE(MYUNIT, *) "SUM OF SQUARED ERROR BETWEEN TARGET AND SIMULATED IMAGES  =   ",MY_ENERGIES(2)
!!        WRITE(MYUNIT, *) "POTENTIAL ENERGY OF THE NANOCLUSTER                       =   ",MY_ENERGIES(3)
!!        WRITE(MYUNIT, *) "WEIGHT BY WHICH POTENTIAL ENERGY WAS MULTIPLIED           =   ",MY_ENERGIES(4)
!!        WRITE(MYUNIT, *) "VALUE OF THE OBJECTIVE FUNCTION AFTER WEIGHT              =   ",MY_ENERGIES(5)
!        CLOSE(MYUNIT)
!        WRITE(MY_FILE_NAME, "(3(a,i0),a)")"Simulated_image_",NUMDIMS,"D_R", &
!            CURRENT_RUN,"_G", CURRENT_GENERATION, ".txt"
!        MY_FILE_NAME = ADJUSTL(TRIM(WHERE_TO_SAVE))//"/"//ADJUSTL(TRIM(MY_FILE_NAME))
!        OPEN(NEWUNIT=MYUNIT, FILE=ADJUSTL(TRIM(MY_FILE_NAME)), STATUS="REPLACE", ACTION="WRITE")
!        DO I=1, NUMROWS
!            WRITE(MYUNIT,*) SIMULATED_IMAGE(I,:)
!        ENDDO
!        CLOSE(MYUNIT)
!    END SUBROUTINE SAVE_ACSTEM_RESULTS
