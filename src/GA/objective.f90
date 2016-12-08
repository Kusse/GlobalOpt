MODULE OBJECTIVE
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE OBJ_FUNC(X, F, G, ND)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ND
        REAL(KIND=8), INTENT(IN), DIMENSION(ND) :: X
        REAL(KIND=8), INTENT(OUT), DIMENSION(ND) :: G
        REAL(KIND=8), INTENT(OUT) :: F
        INTEGER, PARAMETER :: ATOMTYPE = 23
        INTEGER :: NUM_ATOMS
        NUM_ATOMS = ND/3
        CALL GUPTA(ND/3, X, G, F, .TRUE., ATOMTYPE)
    END SUBROUTINE OBJ_FUNC
END MODULE OBJECTIVE