program main
    implicit none
    character(len=100), parameter :: xyz_infile = "Au55chiral.xyz"
    character(len=100), parameter :: r_infile="Au55chiral.r"
    INTEGER, PARAMETER :: NATOMS = 55
    INTEGER, PARAMETER :: NATOM_TYPES = 1
    real(8), parameter :: Bohr2Angstrom=0.529177249D0
    CHARACTER(LEN=2) :: SYMBOL="Au"
    CHARACTER(LEN=50) :: RCRD
    REAL(KIND=8) :: X, Y, Z
    REAL(KIND=8), DIMENSION(3*NATOMS) :: COORDS 
    REAL(KIND=8), DIMENSION(3) :: MAXBOX=(/8.0, 8.0, 8.0/)
    REAL(KIND=8), DIMENSION(3) :: MINBOX=(/-8.0, -8.0, -8.0 /)
    integer :: indx, i, j, MYUNIT1, MYUNIT2
   OPEN(NEWUNIT=MYUNIT1, FILE=ADJUSTL(TRIM(r_infile)), ACTION="WRITE", STATUS="REPLACE")
   WRITE(MYUNIT1, *) SYMBOL
   WRITE(MYUNIT1, *) NATOMS, DBLE(NATOM_TYPES), DBLE(90)
   WRITE(MYUNIT1, *) MAXBOX
   WRITE(MYUNIT1, *) MINBOX
   WRITE(MYUNIT1, *) 0.2796684300000000E-02, 13 
!   WRITE(MYUNIT1, *)
   OPEN(NEWUNIT=MYUNIT2 , FILE=ADJUSTL(TRIM(xyz_infile)), ACTION="READ", STATUS="OLD")
   READ(MYUNIT2, *) RCRD
   READ(MYUNIT2, *)
   print*, "record = ", RCRD
   DO I=1, NATOMS
	READ(MYUNIT2, *) RCRD, J, COORDS(3*I-2:3*I), J 
	WRITE(MYUNIT1, "(3F22.16)") Bohr2Angstrom*COORDS(3*I-2:3*I)
        WRITE(MYUNIT1, "(3F22.16)") 0.0, 0.0, 0.0
	WRITE(MYUNIT1, "(I1)") 1 
    ENDDO
   CLOSE(MYUNIT1)
   CLOSE(MYUNIT2)


end program main
