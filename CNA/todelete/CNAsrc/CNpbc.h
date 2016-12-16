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
 * $Date: 2000/01/27 04:17:09 $
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNpbc.h,v $
 * $Revision: 1.2 $
 * $Log: CNpbc.h,v $
 * Revision 1.2  2000/01/27 04:17:09  sgsrini
 * Changed CN index table to ushort
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *               $$$$ CNpbc.h: Link-List based neighbor lists                *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
#define OUT_OF_BOUNDS   -99
#define INSIDE_THE_BOX  0
#define CROSS_BOUNDARY  1
/*---------------------------------------------------------------------------*/
/* 1a. Apply PERIODIC BOUNDARY CONDITIONS in 1-Dimension */
#define APPLY_PBC_1D(x, ax, a2x, minax, Flag)                            \
{                                                                        \
  if( x >= 0.0 && x < ax )  {                                            \
    if( Flag != OUT_OF_BOUNDS )    Flag = INSIDE_THE_BOX;                \
  }                                                                      \
  else if( x < 0.0 && x >= minax ) {                                     \
    x    += ax;                                                          \
    if( Flag != OUT_OF_BOUNDS )    Flag  = CROSS_BOUNDARY;               \
  }                                                                      \
  else if( x >= ax && x < a2x )   {                                      \
    x    -= ax;                                                          \
    if( Flag != OUT_OF_BOUNDS )    Flag  = CROSS_BOUNDARY;               \
  }                                                                      \
  /* If 'x' is very large, ERROR */                                      \
  else if( x >= a2x  )   {                                               \
    Flag = OUT_OF_BOUNDS;                                                \
    /* June 6, 2001: Set to largest float                                \
    x = 0.99;   */                                                       \
    /* Flag = INSIDE_THE_BOX;*/                                          \
  }                                                                      \
  else if( x < minax )  {                                                \
    Flag = OUT_OF_BOUNDS;                                                \
    /* June 6, 2001: Set to smallest float                             \
    x = 0.0001;    */                                                    \
  }                                                                      \
 }
/*---------------------------------------------------------------------------*/
/* 1b. Apply PERIODIC BOUNDARY CONDITIONS in 3-Dimension */
#define APPLY_PBC(s_Vector, edgeA, edge2A, edgeMinA, Flag)               \
{    APPLY_PBC_1D(s_Vector.x, edgeA[0], edge2A[0], edgeMinA[0], Flag)    \
     APPLY_PBC_1D(s_Vector.y, edgeA[1], edge2A[1], edgeMinA[1], Flag)    \
     APPLY_PBC_1D(s_Vector.z, edgeA[2], edge2A[2], edgeMinA[2], Flag)    \
}
/*---------------------------------------------------------------------------*/
