MODULE cna
    implicit none
!    INTEGER :: I, MYUNIT
!    real(8), dimension(13,3) :: coords
!    INTEGER, dimension(13, 32) :: NEIGBOR_LIST
!    OPEN(NEWUNIT=MYUNIT, FILE="Au13.xyz", action="read")
!    read(myunit, *)
!    read(myunit, *)
!    do i=1, 13
!        read(myunit, *) coords(i,:)
!        print*, i," : ", coords(i, :)
!    enddo
!    call MYCNA(13, COORDS, 1.40d0, NEIGBOR_LIST)
!    do i=1, 13
!        print*, I,":    ", NEIGBOR_LIST(i, :)
!    enddo

    contains

    subroutine MYCNA(NATOMS, COORDS, CUTOFF, NEIGHBOR_LIST) !CNA_INDICES)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NATOMS
        REAL(KIND=8), DIMENSION(NATOMS, 3), INTENT(IN) :: COORDS
!        REAL(KIND=8), DIMENSION(NATOMS, 3), INTENT(OUT) :: CNA_INDICES
        REAL(KIND=8), INTENT(IN) :: CUTOFF
        REAL(KIND=8) :: RC2, DELTA, DIST
        INTEGER, PARAMETER :: MAX_NNB = 32
        INTEGER :: I, J, K, M, N, NNB(NATOMS), NCOMMON, NBOND, NLCHAIN
        INTEGER, DIMENSION(NATOMS, MAX_NNB) :: NEIGHBOR_LIST, CMN
        REAL(8), DIMENSION(NATOMS, NATOMS) :: DISTMATRIX
        RC2 = CUTOFF*CUTOFF
        ! CULCULATE THE FULL NEIGHBOUR LIST OF EACH ATOM
        NEIGHBOR_LIST = 0
        DO I=1, NATOMS
            NNB(I) = 0
            DO J=1, NATOMS
                IF(I /= J) THEN
                    DIST = NORM2(COORDS(I, :) - COORDS(J, :))
                    DISTMATRIX(I, J) = DIST
                    IF(DIST < CUTOFF) THEN
!                        PRINT*, "DITST ",I,J," : ", DIST
                        NNB(I) = NNB(I) + 1
                        IF(NNB(I) < MAX_NNB) THEN
                            NEIGHBOR_LIST(I, NNB) = J
!                            NEIGHBOR_LIST(J, NNB) = I
                        ELSE
                            PRINT*, "ERROR. ATOM HAS MORE THAN ", NNB, "NEAREST NEIGHBORS."
                            STOP
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
!            PRINT*, DISTMATRIX(I, :)
!        ENDDO
!        PRINT*, "NNB = ", NNB
        DO M=1, NNB(I)
            J = NEIGHBOR_LIST(I, M)
            PRINT*, "J = ", J, " NNB(J) ", NNB(J)
            NCOMMON = 0
            DO N=1, NNB(J)
                K = NEIGHBOR_LIST(J, N)
                PRINT*, I, J, "K = ", K
                IF(J == K) THEN
                    NCOMMON = NCOMMON + 1
                    CMN(I, NCOMMON) = NEIGHBOR_LIST(M, N)
                ENDIF
            ENDDO
            PRINT*, I, " CMN = ", NCOMMON
!            STOP
        ENDDO
!        PRINT*, I, " CMN = ", NCOMMON   !,   " VALS = ", CMN(I, NCOMMON)
    ENDDO


    end subroutine mycna

end MODULE cna
