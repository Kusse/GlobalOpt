/*****************************************************************************
 *$$$$$$$$$$$$$$$$$$$$$$ Common Neighbor Analysis Code $$$$$$$$$$$$$$$$$$$$$$*
 *****************************************************************************
 *  AUTHOR:                                                                  *
 *     Srinivasan G. Srivilliputhur,  T-11;   sgsrini@lanl.gov               *
 *     Srinivasan G. Srivilliputhur,  University of Washington, Seattle      *
 *                                                                           *
 *                  *** NOT FOR GENERAL RELEASE ***                          *
 *                                                                           *
 * DISCLAIMER : This is prerelease software.  We reserve the right to make   *
 * change(s). Do NOT USE or DISTRIBUTE it without our written consent.       *
 *****************************************************************************
 * Details of CNA is in:                                                     *
 *          1. Clarke and Jonsson, PRE 47, 3975 (1991)                       *
 *          2. Faken and Jonsson,  Comp. Mater. Sci., (1994)                 *
 *****************************************************************************
 * $Author: sgsrini $ : Srinivasan G. Srivilliputhur (T-11, LANL)
 * $Date: 2000/01/28 16:03:15 $
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNmacros.h,v $
 * $Revision: 1.3 $
 * $Log: CNmacros.h,v $
 * Revision 1.3  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
 *
 * Revision 1.2  2000/01/27 04:17:09  sgsrini
 * Changed CN index table to ushort
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *              $$$$ macros.h: Vector and Matrix Operations                  *
 *****************************************************************************
 * THIS FILE CONTAINS CODE FRAGMENTS OF SOME GENERAL PURPOSE MACROS USED     *
 * FOR MATRIX AND ROUTINE ARITHMETIC. EX: MULTIPLICATION, INVERSE ETC. THE   *
 * 3 X 3 matrices are saved in 1-D arrays. Relation between the 1-D array    *
 * indices and the 3 X 3 matrix indices are as Follows:                      *
 *          '0' = '11'; '1' = '12'; '2' = '13';                              *
 *          '3' = '21'; '4' = '22'; '5' = '23';                              *
 *          '6' = '31'; '7' = '32'; '8' = '33'.                              *
 *****************************************************************************
 *  A. Macros for Some Routine Arithmetic operations- Absolute Values etc.   *
 *****************************************************************************/
/* Absolute Value of a Number */
#define ABS(a)          ( (a) < 0.0 ? -(a) : (a) ) 
/*---------------------------------------------------------------------------*/
/* Square of a number */
#define SQR(x)          ((x) * (x))
/*---------------------------------------------------------------------------*/
/* Smaller of the Two Numbers x and y:  */
#define MIN(x,y)        ( (x) < (y) ? (x) : (y) )
/*---------------------------------------------------------------------------*/
/* Larger of the Two Numbers x and y:  */
#define MAX(x,y)        ( (x) > (y) ? (x) : (y) )
/*---------------------------------------------------------------------------*/
/* Cube of a number */
#define CUBE(x)         ((x) * (x) * (x))
/*---------------------------------------------------------------------------*/
/* Larger of the 3 Numbers x and y */
/* #define MAX3(x,y,z)     MAX(x, MAX(y,z)) */
/*---------------------------------------------------------------------------*/
/* Smaller of the 3 Numbers x and y */
/* #define MIN3(x,y,z)     MIN(((x) < (y) ? (x) : (y)),z) */
/*---------------------------------------------------------------------------*/
/* Square of Distance/Distance between two points */
/* #define DIST_SQ(x,y)   (SQR(x[0]-y[0])+SQR(x[1]-y[1])+SQR(x[2]-y[2])) */
/* #define DIST(x,y)      sqrt( (double) (DIST_SQ( x, y )) ) */
/*---------------------------------------------------------------------------*/
                  /*********************************************
		   B. Macros for Matrix Arithmetic operations
		   *********************************************/
/* Add two 3x3 matrix "B" and "C";  i.e. A = B + C */
#define ADD_MATRICES( A, B, C )                                           \
{ A[0] = B[0] + C[0];    A[1] = B[1] + C[1];   A[2] = B[2] + C[2];        \
  A[3] = B[3] + C[3];    A[4] = B[4] + C[4];   A[5] = B[5] + C[5];        \
  A[6] = B[6] + C[6];    A[7] = B[7] + C[7];   A[8] = B[8] + C[8];        \
}
/*---------------------------------------------------------------------------*/
/* Assign 3x3 matrix "B" to "A"; i.e. A = B */
#define ASSIGN_MATRICES( A, B )       {                                   \
  A[0] = B[0];            A[1] = B[1];          A[2] = B[2];              \
  A[3] = B[3];            A[4] = B[4];          A[5] = B[5];              \
  A[6] = B[6];            A[7] = B[7];          A[8] = B[8];              \
}
/*---------------------------------------------------------------------------*/
/* Invert 3x3 matrix "A"; result in "B"; RETURNS determinant(A). Matrices
   A and B must be different; i.e. A-inverse cannot be assigned to A itself */
#define INVERT_MATRIX( B, A, tmpVol )  {                                  \
  register int4  i;                                                       \
  tmpVol = 0.0;                              /* initialize to zero */     \
  /* trap error if same variables are used for A and B */                 \
  if( B == A ) Error("INVERT()/macros.h: Don't use same A and B");        \
  B[0] = A[4] * A[8] - A[5] * A[7];   B[1] = A[2] * A[7] - A[1] * A[8];   \
  B[2] = A[1] * A[5] - A[2] * A[4];   B[3] = A[5] * A[6] - A[3] * A[8];   \
  B[4] = A[0] * A[8] - A[2] * A[6];   B[5] = A[2] * A[3] - A[0] * A[5];   \
  B[6] = A[3] * A[7] - A[4] * A[6];   B[7] = A[1] * A[6] - A[0] * A[7];   \
  B[8] = A[0] * A[4] - A[1] * A[3];                                       \
  tmpVol = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];                             \
  for(i = 0; i < 9; i++) B[i] /= tmpVol;                                  \
}
/*---------------------------------------------------------------------------*/
/* Multiply two 3x3 matrices "B" and "C", put result in "A".    i.e. A=BC */
#define MULTIPLY_MATRIX_MATRIX( A, B, C )  {                              \
  /* trap error if same variables are used for A and B */                 \
  if( (B==A) || (A==C) )                                                  \
    Error("MULTIPLY_MATRIX_MATRIX()/macros.h: Don't use same A and B/C"); \
  A[0] = B[0]*C[0] + B[1]*C[3] + B[2]*C[6];                               \
  A[1] = B[0]*C[1] + B[1]*C[4] + B[2]*C[7];                               \
  A[2] = B[0]*C[2] + B[1]*C[5] + B[2]*C[8];                               \
  A[3] = B[3]*C[0] + B[4]*C[3] + B[5]*C[6];                               \
  A[4] = B[3]*C[1] + B[4]*C[4] + B[5]*C[7];                               \
  A[5] = B[3]*C[2] + B[4]*C[5] + B[5]*C[8];                               \
  A[6] = B[6]*C[0] + B[7]*C[3] + B[8]*C[6];                               \
  A[7] = B[6]*C[1] + B[7]*C[4] + B[8]*C[7];                               \
  A[8] = B[6]*C[2] + B[7]*C[5] + B[8]*C[8];                               \
}
/*---------------------------------------------------------------------------*/
/* Divide 3x3 matrix "A" by Scalar "m"; store result in 3X3 "C" 
   June 16, 1999: Use Multiply Instead of Divide */
#define SCALAR_DIVIDE_MATRIX( A, SCALAR, B )  {                           \
  Real invSCAL = ((Real) 1.0/(Real) SCALAR);                              \
  A[0] = B[0]*invSCAL;    A[1] = B[1]*invSCAL;    A[2] = B[2]*invSCAL;    \
  A[3] = B[3]*invSCAL;    A[4] = B[4]*invSCAL;    A[5] = B[5]*invSCAL;    \
  A[6] = B[6]*invSCAL;    A[7] = B[7]*invSCAL;    A[8] = B[8]*invSCAL;    \
}
/*---------------------------------------------------------------------------*/
/* Multiply 3x3 matrix "A" by Scalar "m"; store result in 3X3 "C" */
#define SCALAR_MULTIPLY_MATRIX( A, SCALAR, B )  {                         \
  A[0] = SCALAR * B[0];  A[1] = SCALAR * B[1];  A[2] = SCALAR * B[2];     \
  A[3] = SCALAR * B[3];  A[4] = SCALAR * B[4];  A[5] = SCALAR * B[5];     \
  A[6] = SCALAR * B[6];  A[7] = SCALAR * B[7];  A[8] = SCALAR * B[8];     \
}
/*---------------------------------------------------------------------------*/
/* Subtracts 3x3 matrix "B" and "C";  i.e. A = B - C */
#define SUBTRACT_MATRICES( A, B, C )          {                           \
  A[0] = B[0] - C[0];    A[1] = B[1] - C[1];   A[2] = B[2] - C[2];        \
  A[3] = B[3] - C[3];    A[4] = B[4] - C[4];   A[5] = B[5] - C[5];        \
  A[6] = B[6] - C[6];    A[7] = B[7] - C[7];   A[8] = B[8] - C[8];        \
}
/*---------------------------------------------------------------------------*/
/* Transpose 3x3 matrix "A"; store result in "B" */
#define TRANSPOSE_MATRIX( B, A )                    {                     \
  if( ((char*) B) == ((char*) A) ) {                                      \
   Error("TRANSPOSE(): matrices MUST be different\n");                    \
   }                                                                      \
  B[0] = A[0];   B[1] = A[3];   B[2] = A[6];                              \
  B[3] = A[1];   B[4] = A[4];   B[5] = A[7];                              \
  B[6] = A[2];   B[7] = A[5];   B[8] = A[8];                              \
  }
/*---------------------------------------------------------------------------*/
#define ZERO_THE_MATRIX( matrix )                                         \
{ matrix[0] = matrix[1] = matrix[2] = matrix[3] = 0.0;                    \
  matrix[4] = matrix[5] = matrix[6] = matrix[7] = matrix[8] = 0.0;        \
}
/*---------------------------------------------------------------------------*/
                /*****************************************************
		 C. Macros for Vector Arithmetic operations (3D)
		 *****************************************************/
/* Assign Vector "B" to "A"; i.e. A = B */
#define ASSIGN_3D_VECTORS( A, B)  { A.x = B.x;  A.y = B.y;  A.z = B.z; }
/*---------------------------------------------------------------------------*/
/* Add two 3x1 vectors A and B; result in 3X1 C */
#define ADD_3D_VECTORS( C, B, A)    {                                     \
  (C.x) = (B.x) + (A.x);                                                  \
  (C.y) = (B.y) + (A.y);                                                  \
  (C.z) = (B.z) + (A.z);                                                  \
  }
/*---------------------------------------------------------------------------*/
/* Post-Multiply 3x3 matrix "B" with 3X1 Vector "C"; Vec(A) != Vec(C) */
#define MATRIX_3D_VECTOR_MULTIPLY( A, B, C)     {                         \
  /* trap error if same variables are used for A and B */                 \
  char *p1, *p2;                                                          \
  p1 = (char *) &A; p2 = (char *) &C;                                     \
  if(p1==p2) Error("MULTIPLY_VECTOR()/macros.h: Don't use same A and C"); \
  (A.x) = B[0] * (C.x) + B[1] * (C.y) + B[2] * (C.z);                     \
  (A.y) = B[3] * (C.x) + B[4] * (C.y) + B[5] * (C.z);                     \
  (A.z) = B[6] * (C.x) + B[7] * (C.y) + B[8] * (C.z);                     \
  }
/*---------------------------------------------------------------------------*/
/* Multiply 3x1 vector "A" by Scalar "m" */
#define SCALAR_3D_VECTOR_MULTIPLY( VECTOR2, SCALAR, VECTOR1)   {          \
  VECTOR2.x = SCALAR * (VECTOR1.x);                                       \
  VECTOR2.y = SCALAR * (VECTOR1.y);                                       \
  VECTOR2.z = SCALAR * (VECTOR1.z);                                       \
  }
/*---------------------------------------------------------------------------*/
/* Subtract 3x1 vector "A" from 3X1 Vector "B" */
#define SUBTRACT_3D_VECTORS( C, B, A )                  {                 \
  (C.x) = (B.x) - (A.x);                                                  \
  (C.y) = (B.y) - (A.y);                                                  \
  (C.z) = (B.z) - (A.z);                                                  \
  }
/*---------------------------------------------------------------------------*/
/* Square of a 'vector'  */
#define SUMSQ(vec)    ((vec.x)*(vec.x)+(vec.y)*(vec.y)+(vec.z)*(vec.z))  
/*---------------------------------------------------------------------------*/
/* DOT or SCALAR PRODUCT of 3x1 vectors "A" and "B"; C = (B * A) = (A * B)
    NOTE: This is different from "SCALAR_MULTIPLY()"  */
#define VECTOR_DOT_PRODUCT_3D( C, B, A )             {                    \
  C = ((B.x)*(A.x)) + ((B.y)*(A.y)) + ((B.z)*(A.z));                      \
  }
/*---------------------------------------------------------------------------*/
/* Cross or Vector product of 2 vectors: C = A x B */
#define VECTOR_CROSS_PRODUCT_3D( C, A, B )                {                \
  /* trap error if same variables are used for C and A OR C and B */      \
  char *p1, *p2, *p3;                                                     \
  p1 = (char *) &A; p2 = (char *) &B; p3 = (char *) &C;                   \
  if( (p1==p3) || (p2==p3) )                                              \
    Error("VEC_CROSS_PROD()/macros.h: Don't use same A/C or B/C");        \
  C.x = ( ((A.y)*(B.z)) - ((A.z)*(B.y)) );                                \
  C.y = ( ((A.z)*(B.x)) - ((A.x)*(B.z)) );                                \
  C.z = ( ((A.x)*(B.y)) - ((A.y)*(B.x)) );                                \
}
/*---------------------------------------------------------------------------*/
#define ZERO_3D_VECTOR_( A ) {  A.x = A.y = A.z = 0.0;   } 
/*---------------------------------------------------------------------------*/
