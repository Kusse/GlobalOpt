PROGRAM main
    IMPLICIT NONE
    INTEGER :: I
    DO I=1, 10
     CALL MYGA_RUN()
    ENDDO
!CONTAINS
!    SUBROUTINE READ_CMD_ARGS()
!        USE COMMONS
!        USE GA_PARAMS
!        IMPLICIT NONE
!        integer :: nargs, i, j, STAT, CLEN
!        CHARACTER(LEN=50) :: buffer, fn
!        CALL DEFAULT_VALUES()
!        !=================================================================
!        ! Open the input file, read in the values and assign them to the
!        ! appropraite varaibles
!        !=================================================================
!
!        NARGS = COMMAND_ARGUMENT_COUNT()
!        CALL get_command(CMD_LINE, CLEN, STAT)
!        DO I=1, NARGS
!            CALL GET_COMMAND_ARGUMENT(i, buffer)
!            buffer = TRIM(buffer)
!            IF(I==1) THEN
!                READ(buffer,    *) NATOMS
!            ELSEIF(I==2) THEN
!                READ(buffer,    *) NDIM
!            ELSEIF(I==3) THEN
!                READ(buffer,    *) GATOM
!            ELSEIF(i==4) THEN
!                READ(buffer,    *) MYGA_NSTRUC
!            ELSEIF(i==5) then
!                READ(buffer,    *) MYGA_NOFF
!            ELSEIF(i==6) then
!                READ(buffer,    *) MYGA_CROSS
!            ELSEIF(i==7) then
!                READ(buffer,   FMT="(f15.8)") MYGA_MUT_RATE
!            ELSEIF(i==8) then
!                READ(buffer,    FMT="(f15.8)") STEP
!            ELSEIF(i==9) then
!                READ(buffer,   *) MAXFUNCALLS
!            ELSEIF(i==10) then
!                READ(buffer, *) MAXRUNS
!            ELSE
!                PRINT*, "Last command line argument is not recognized."
!                stop
!            ENDIF
!        ENDDO
!    END SUBROUTINE READ_CMD_ARGS
END PROGRAM main
