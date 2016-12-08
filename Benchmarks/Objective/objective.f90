MODULE OBJECTIVE
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE OBJ_FUNC(X, F, G, ND)
	USE MORSE, ONLY: MORSEGH
	USE LJ, ONLY: LJGH

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ND
        REAL(KIND=8), INTENT(IN), DIMENSION(ND) :: X
        REAL(KIND=8), INTENT(OUT), DIMENSION(ND) :: G
        REAL(KIND=8), INTENT(OUT) :: F
	LOGICAL, PARAMETER :: GTEST=.TRUE., HTEST=.FALSE.
        INTEGER, PARAMETER :: ATOMTYPE = 23
	REAL(KIND=8), PARAMETER :: RHO = 1.0D0
	CHARACTER(10), PARAMETER :: POTENTIAL_TYPE = "GUPTA"
        INTEGER :: NTOMS
        NTOMS = ND/3
	IF (POTENTIAL_TYPE == "GUPTA") THEN
        	CALL GUPTA(NTOMS, X, F, G, GTEST, ATOMTYPE)
   	ELSEIF(POTENTIAL_TYPE == "MORSE") THEN
		CALL MORSEGH(NTOMS, X, F, G, GTEST, HTEST, RHO)
	ELSEIF(POTENTIAL_TYPE == "LJ") THEN
		CALL LJGH(NTOMS, X, F, G,  GTEST, HTEST)
	ELSE
		PRINT*, "POTEITNAL FUNCTION NOT AVAILABLE."
		STOP
	ENDIF
 	END SUBROUTINE OBJ_FUNC
END MODULE OBJECTIVE
