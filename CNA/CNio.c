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
 * $Date: 2002/08/13 20:36:38 $
 * $Source: /home/srini/LJ_EAM/CN2000/RCS/CNio.c,v $
 * $Revision: 1.1 $
 * $Log: CNio.c,v $
 * Revision 1.1  2002/08/13 20:36:38  sgsrini
 * Initial revision
 *
 * Revision 1.6  2001/06/06 22:34:03  sgsrini
 * Error in computing box edge length for Baskes's file
 *
 * Revision 1.5  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
 *
 * Revision 1.4  2000/01/27 21:45:31  sgsrini
 * *** empty log message ***
 *
 * Revision 1.3  2000/01/27 21:39:21  sgsrini
 * Added functions to compute and write CNN pair statistics
 *
 * Revision 1.2  2000/01/27 04:16:44  sgsrini
 * Changed CN index table to ushort
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *              $$$$ CNio.c: Input-Output module for CNA code                *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/*---------------------------------------------------------------------------*/
#include "CNmain.h"
#include "CNmacros.h"
#include "CNpairStat.h"
#include "CNpbc.h"
/*---------------------------------------------------------------------------*/
/* For Mike Baskes's file format */
/* #define FREE_SFCE_X 1
   #define FREE_SFCE_Y 1
   #define FREE_SFCE_Z 1 */

#define FREE_SFCE_X 0
#define FREE_SFCE_Y 0
#define FREE_SFCE_Z 0
/*---------------------------------------------------------------------------*/
/****** extern global function declarations *****************************/
extern void Error( char*, ... );
extern void LCInitializeSlots(int4, int4, int4, uint4, int4 );
extern int LCInsertAtomIntoSlot(float*,  float*, Vector_t*, 
				float,   uint4, int4, char*, int4 );
extern LCType *LCgetLinkCellP( void );
/*---------------------------------------------------------------------------*/
/****** extern global variable declarations *****************************/
extern Atom_t *AtomP;
/*---------------------------------------------------------------------------*/
#define ARY_SZ            9             /* Buffer Array Sizes     */
/* #define N_TYPE            6 */
/*---------------------------------------------------------------------------*/
static float        SnapTime = 0.0;
static float        BoxVectors[9];
static float        BoxVecDeri[9];
static char         Buf_CH[ARY_SZ];
static double       Buf_DO[ARY_SZ];
static float        Buf_FL[ARY_SZ];
static unsigned int Buf_UI[ARY_SZ];
static double       AtomMass[N_TYPE];

/* 06/28/03: AtomicNumber[] was declared along with AtomMass[] as a double */
static unsigned int        AtomicNumber[N_TYPE];
/*---------------------------------------------------------------------------*/
/* Mike's dunamo code has origin at box center => while writing out file
   we need to move origin to box center for dynamo format files */
int                 Ntype                     = 1;
static int          BaskesSurfaceConfigFlag   = NO;
static int          MoveOriginToBoxCenterFlag = NO;
static float        InitCooAtom01[3];
static float        MaxBoxEdges[3];
static float        MinBoxEdges[3];
/*---------------------------------------------------------------------------*/
struct FHead  {
  double  TimeOfSnap;             /* Total Timeoffset included in Time */
  uint4 Dimension;               /* 2 or 3                            */
  uint4 TotNumAtoms;             /* Total Number of Atoms             */
  uint4 TotBluAtoms;             /* Total Num Blue  Atoms             */
  uint4 TotRedAtoms;             /* Total Num Red   Atoms             */
  uint4 TotGrnAtoms;             /* Total Num Green Atoms             */
  uint4 CoordFlag;               /* Precision: 0/Single/Double        */
  uint4 VelFlag;                 /* Precision: 0/Single/Double        */
  uint4 PEFlag;                  /* Precision: 0/Single/Double        */
  uint4 StressFlag;              /* Precision: 0/Single/Double        */
  uint4 NstressComp;             /* Number of Stress Components       */
  uint4 TypeFlag;                /* TypeInfo  Absent(0)/Present(1)    */
  uint4 ColorFlag;               /* ColorInfo Absent(0)/Present(1)    */
  uint4 MoveFlag;                /* MoveInfo  Absent(0)/Present(1)    */

  uint4 InitCoordFlag;           /* Coord of Start system: 0/Single(1)*/
  uint4 SpinFlag;                /* WSpin Vector Coord: 0/1/2         */
  uint4 padding;                 /*Make header multiple sizeof(double)*/
};
/*---------------------------------------------------------------------------*/
static struct FHead FileHeader;
static struct FHead *FileHeaderP    = &FileHeader;
/*---------------------------------------------------------------------------*/
#define FREAD_STRUCT( filHdrP )                       {                   \
  if( fread( filHdrP, sizeof(FileHeader), 1, ifp ) != 1 )                 \
    Error("Sa. ReadStruct(): Cannot read header structure");              \
}
/*---------------------------------------------------------------------------*/
#define FREAD_UINT( nElem )                          {                   \
  if( nElem >  ARY_SZ ) Error("Ua. ReadUint(): UINT Array Overflow");    \
  if( fread( Buf_UI, sizeof(unsigned int), nElem, ifp ) != nElem )       \
    Error("Ub. ReadUint(): File read error");                            \
}
/*---------------------------------------------------------------------------*/
#define FREAD_SINGLE_PRECISION( nElem )               {                   \
  if( nElem >  ARY_SZ ) Error("Fa. ReadFloat(): Float Array Overflow");   \
  if( fread( Buf_FL, sizeof(float), nElem, ifp ) != nElem )               \
    Error("Fb. ReadFloat(): File read error");                            \
}
/*---------------------------------------------------------------------------*/
#define FREAD_SINGLE_PRECISION_SPASM( nElem )        {                   \
  if( nElem >  ARY_SZ ) Error("Fa. ReadFloat(): Float Array Overflow");   \
  if( fread( Buf_FL, sizeof(float), nElem, ifp ) != nElem ) break;        \
}
/*---------------------------------------------------------------------------*/

#define FREAD_DOUBLE_PRECISION( nElem )               {                   \
  if( nElem >  ARY_SZ ) Error("Da. ReadDouble(): Double Array Overflow"); \
  if( fread( Buf_DO, sizeof(double), nElem, ifp ) != nElem )              \
    Error("Db. ReadDouble(): File read error");                           \
}
/*---------------------------------------------------------------------------*/
#define FREAD_CHARACTERS( nElem )                     {                   \
  if( nElem >  ARY_SZ ) Error("Ca. ReadChar(): Char Array Overflow");     \
  if( fread( Buf_CH, sizeof(char), nElem, ifp ) != nElem )                \
    Error("Cb. ReadChar(): File read error");                            \
}
/*---------------------------------------------------------------------------*/
#define FWRITE_STRUCT( filHdrP )                      {                   \
  if( fwrite( filHdrP, sizeof(FileHeader), 1, ofp ) != 1 )    {           \
    if( feof(ofp) ) Error("Sa. WriteStruct(): Premature end of file");    \
    else           Error("Sb. WriteStruct(): File write error");          \
  }                                                                       \
}
/*---------------------------------------------------------------------------*/
#define FWRITE_UINT( nElem )                         {                   \
  if( nElem >  ARY_SZ ) Error("Ua. WriteUint(): UINT Array Overflow");   \
  if( fwrite( Buf_UI, sizeof(unsigned int), nElem, ofp ) != nElem )      \
    Error("Ub. WriteUint(): File write error");                          \
}
/*---------------------------------------------------------------------------*/
#define FWRITE_SINGLE_PRECISION( nElem )              {                   \
  if( nElem >  ARY_SZ ) Error("Fa. WriteFloat(): Float Array Overflow");  \
  if( fwrite( Buf_FL, sizeof(float), nElem, ofp ) != nElem )              \
    Error("Fb. WriteFloat(): File write error");                          \
}
/*---------------------------------------------------------------------------*/
#define FWRITE_DOUBLE_PRECISION( nElem )              {                   \
  if( nElem >  ARY_SZ ) Error("Da. WriteDouble(): Double Array Overflow");\
  if( fwrite( Buf_DO, sizeof(double), nElem, ofp ) != nElem )             \
    Error("Db. WriteDouble(): File write error");                         \
}
/*---------------------------------------------------------------------------*/
#define FWRITE_CHARACTERS( nElem )                    {                   \
  if( nElem >  ARY_SZ ) Error("Ca. WriteChar(): Char Array Overflow");    \
  if( fwrite( Buf_CH, sizeof(char), nElem, ofp ) != nElem )               \
    Error("Cb. WriteChar(): File write error");                           \
}
/*---------------------------------------------------------------------------*/
void DbgPrintReadInAtoms( uint4 nAtoms )
{
  int4 ia = 0;

  for( ia = 0; ia < nAtoms; ia++ )  {
    fprintf( stderr,
	    "$$$ DbgAtoms(): ID(%7u)/X(%f %f %f)\n",
	    AtomP[ia].id, (float) AtomP[ia].coo.x, (float) AtomP[ia].coo.y, 
	    (float) AtomP[ia].coo.z );

    fflush( stderr );
  }
}
/*---------------------------------------------------------------------------*/
/* 0a. TIME FUNCTIONS ---*/
static void GetStartTime( FILE *fp )
{
  static int    firstTime = 1;  /* first time function is called */
  static time_t tNow;           /* Time Now                      */

  if( firstTime ) {
    firstTime = 0;

    /* Get time */
    if( time( &tNow ) == -1 )
      Error("$$$GetStartTime(): Time not available?");
  } /* if( firstTime ) */

  /* Convert to Ascii Format and print to stdout ot to file fp*/
  if( fp == NULL ) {
    fprintf( stdout, "#\tMD Simulation Started At:\t%s\n", ctime( &tNow ) );
    fflush( stdout );
  }
  else {
    fprintf( fp, "\t****************************\n\tCN-Anal Started At:\t%s", 
            ctime( &tNow ) );
    fprintf( fp, "\t****************************\n" );
    fflush( fp );
  }
} /* void GetStartTime() */
/*---------------------------------------------------------------------------*/
/* 0b. Call at the end of a run (outside the while loop) */
static void GetFinishTime( FILE *fp )
{
  static time_t tNow;   /* Time Now */

  /* Get time */
  if( time( &tNow ) == -1 )
    Error("$$$GetFinishTime(): Time not available?");

  /* Convert to Ascii Format and print to stdout ot to file fp*/
  if( fp == NULL ) {
    /* Convert to Ascii Format and print to stdout */
    fprintf( stdout, "#\tCN-Anal Finished At:\t%s\n", ctime( &tNow ) );
    fflush( stdout );
  }
  else {
    fprintf( fp, "\t****************************\n\tCN-Anal Finished At:\t%s", 
            ctime( &tNow ) );
    fprintf( fp, "\t****************************\n\n" );
    fprintf( fp, "\t\n############################\n\n" );
    fflush( fp );
  }
} /* void GetFinishTime( FILE *fp ) */
/*---------------------------------------------------------------------------*/
/* 0c. Print run-control params, coordination data, pairs-type etc. info */
void WritePairsInfo(CNstat_t *cnStatTblP, char *ifile, float *boxP, float rcut,
		    uint4 nAtoms, int4 verbose )
{ 
  LCType     *lcP        = NULL;
  FILE       *fp         = NULL;
  uint4      ia          = 0;
  uint4      l           = 0;
  uint4      len         = 0;
  uint8      totNumPairs = 0;
  int4       itmp, i, j, k;
  int4       nType       = 1;
  uint4      totNbcc     = 0,
             totNfcc     = 0,
             totNhcp     = 0,
             totNico     = 0,
             totNdisIco  = 0;

  /* Append Mode */
  if (!( fp = fopen("CNanalysis.log", "a+")))
    Error("\nWriteCoordAndPairsInfo():\tCan't open 'CNanalysis.log'\n");

  /* Get Start Time for CN analysis */
  GetStartTime( fp );

  if( Ntype == 2 ) nType = 2;
  else             nType = 1;

  fprintf(fp, "\n*************** SYSTEM PARAMETERS ***************\n");
  /* 1. CutOff distances for forming pairs */
  fprintf(fp, "CutOff Distance for Neighbor Pairs\t:\t%f\n", (float) rcut );

  /* 2. Name of Atom Configuration file */
  fprintf(fp, "Input Atom config file\t\t\t:\t%s\n", ifile );

  /* 3a. Box Vector "A" */
  fprintf(fp, "Box Vector (A)\t\t\t\t:\t(%f, %f, %f)\n", 
	  (float) boxP[0], (float) boxP[3], (float) boxP[6] );

  /* 3b. Box Vector "B" */
  fprintf(fp, "Box Vector (B)\t\t\t\t:\t(%f, %f, %f)\n", 
	  (float) boxP[1], (float) boxP[4], (float) boxP[7] );

  /* 3c. Box Vector "B" */
  fprintf(fp, "Box Vector (C)\t\t\t\t:\t(%f, %f, %f)\n", 
	  (float) boxP[2], (float) boxP[5], (float) boxP[8] );

  /* 4. Get slot etc. info */
  lcP = (LCType*) LCgetLinkCellP();
  fprintf(fp, "Slots\t\t\t\t\t:\t(%d, %d, %d)\n", 
	  lcP->NslotsX, lcP->NslotsY, lcP->NslotsZ );

  fprintf(fp, "Min/Max Pair Separation\t\t\t:\t(%f, %f)\n", 
	  (float) lcP->MinRij, (float) lcP->MaxRij );

  fprintf(fp, "Min/Max Number of Nbors\t\t\t:\t(%d, %d)\n", 
	  lcP->MinNbors, lcP->MaxNbors );

  fprintf(fp, "*****************************************************\n\n");

  fprintf(fp, "************ # BCC, FCC, HCP, ICO Atoms *************\n");

  fprintf( fp,"# Atoms In the System\t\t\t:\t%u\n", nAtoms );

  /* %% prints percent symbol */
  for( i = 0; i < nType; i++ ) {
    totNbcc     += cnStatTblP->Nbcc[i];
    totNfcc     += cnStatTblP->Nfcc[i];
    totNhcp     += cnStatTblP->Nhcp[i];
    totNico     += cnStatTblP->Nico[i];
    totNdisIco  += cnStatTblP->Nicodist[i];

    if( cnStatTblP->Nbcc[i] > 0 )
      fprintf(fp,
	      "# Type-%d BCC Atoms In the Box\t\t:\t%5u (%.3f  %%) @Time= %.8e\n", 
	      i+1, cnStatTblP->Nbcc[i],
	      ((((float) cnStatTblP->Nbcc[i]) * 100.0)/((float)  nAtoms)),
	      SnapTime );
    
    /*   %% prints percent symbol */
    if( cnStatTblP->Nfcc[i] > 0 )
      fprintf(fp,
	      "# Type-%d FCC Atoms In the Box\t\t:\t%5u (%.3f  %%)  @Time= %.8e\n", 
	      i+1, cnStatTblP->Nfcc[i], 
	      ((((float) cnStatTblP->Nfcc[i]) * 100.0)/((float)  nAtoms)),
	      SnapTime);
    
    /* %% prints percent symbol */
    if( cnStatTblP->Nhcp[i] > 0 )
      fprintf(fp,
	      "# Type-%d HCP Atoms In the Box\t\t:\t%5u (%.3f  %%)  @Time= %.8e\n", 
	      i+1, cnStatTblP->Nhcp[i],
	      ((((float) cnStatTblP->Nhcp[i]) * 100.0)/((float)  nAtoms)),
	      SnapTime);
    
    /* %% prints percent symbol */
    if( cnStatTblP->Nico[i] > 0 )
      fprintf(fp,
	      "# Type-%d ICOSAHEDRA Coord Atom\t\t:\t%5u (%.3f  %%)  @Time= %.8e\n", 
	      i+1, cnStatTblP->Nico[i],
	      ((((float) cnStatTblP->Nico[i]) * 100.0)/((float)  nAtoms)),
	      SnapTime);
    
    if( cnStatTblP->Nicodist[i] > 0 )
      fprintf(fp,
	      "# Type-%d DISTORTED-ICO Coord Atom\t:\t%5u (%.3f  %%)  @Time= %.8e\n",
	      i+1, cnStatTblP->Nicodist[i],
	      ((((float) cnStatTblP->Nicodist[i]) * 100.0)/((float) nAtoms)),
	      SnapTime);
  }
  
  fprintf(fp, "*****************************************************\n");

#ifdef SRINI_COMMENT_NOV_2004
  if( cnStatTblP->NatmInFullICO > 0 )
    fprintf(fp,"# Atoms in a FULL ICOSAHEDRA\t\t:\t%u (%.3f  %%)  @ Time= %.8e\n",
	    cnStatTblP->NatmInFullICO,
	    ((((float) cnStatTblP->NatmInFullICO) * 100.0)/((float)  nAtoms)),
	    SnapTime);
  
  
  if( cnStatTblP->NatmInDistICO > 0 )
    fprintf(fp,"# Atoms in a DISTORTED-ICO\t\t:\t%u (%.3f  %%)  @ Time= %.8e\n",
	    cnStatTblP->NatmInDistICO,
	    ((((float) cnStatTblP->NatmInDistICO) * 100.0)/((float) nAtoms)),
	    SnapTime);
  fprintf(fp, "*****************************************************\n\n");
#endif
  
    fprintf(fp, "***************** COORDINATION DATA *****************\n");

    /* NOV 12, 2004: This is commented out because it's wrongly computed */
    /* fprintf( fp,"Average Number of Neighbors\t\t:\t%.4f\n",
       (float) cnStatTblP->avneb ); */

  /* Input file has NO pairPE info */
  fprintf( fp,"CoordNum 'i'\t# of Atoms with CN 'i'\t%% Atoms With CN 'i'\n" );
  fprintf( fp,"------------\t----------------------\t-------------------\n" );
  for( i = 0; i < NNBR; i++ )  {
    if( cnStatTblP->NcoordnumbM[i] <= 0 ) continue;
    fprintf( fp,"\t%d\t\t%u\t\t%.4f\n", i, cnStatTblP->NcoordnumbM[i], 
	    (float) cnStatTblP->PercoordnAA[i] );
  } /* for( i ) */

  if( totNbcc > 0 )
    fprintf( fp,"\tBCC\t\t%u\t\t%.4f\n", totNbcc,
	    ((((float) totNbcc) * 100.0)/((float)  nAtoms)) );

  if( totNfcc > 0 )
    fprintf( fp,"\tFCC\t\t%u\t\t%.4f\n", 
             totNfcc, 
             ((((float) totNfcc) * 100.0)/((float)  nAtoms)) );

  if( totNhcp > 0 )
    fprintf( fp,"\tHCP\t\t%u\t\t%.4f\n", 
             totNhcp,
             ((((float) totNhcp) * 100.0)/((float)  nAtoms)) );

  if( totNico > 0 )
    fprintf( fp,"\tICO\t\t%u\t\t%.4f\n", 
             totNico,
             ((((float) totNico) * 100.0)/((float)  nAtoms)) );

  if( totNdisIco > 0 )
    fprintf( fp,"\tDIC\t\t%u\t\t%.4f\n", 
             totNdisIco,
             ((((float) totNdisIco) * 100.0)/((float)  nAtoms)) );
  fprintf(fp, "*****************************************************\n\n");

  fprintf(fp, "************** PAIRS WITH DIFFERENT INDICES **************\n");

  /* %% prints percent symbol */
  fprintf( fp,"Pair Index (J, K, L)\t# of (J,K,L) Pairs\t%% of Total Pairs\n");
  fprintf( fp,"--------------------\t------------------\t-----------------\n");

  /* August 19, 1997: First count total num pairs */
  totNumPairs = 0;
  for( i = 0; i < NNBR; i++ )
    for( j = 0; j < NNBR; j++ )
      for( k = 0; k < NNBR; k++ )  {
        if( cnStatTblP->NumDiffPairs[i][j][k] <= 0 ) continue;
        totNumPairs += ((uint4) cnStatTblP->NumDiffPairs[i][j][k]);
      }  /* for(i:j:k) */


  /* Now compute percent of various pairs */
  for( i = 0; i < NNBR; i++ )
    for( j = 0; j < NNBR; j++ )
      for( k = 0; k < NNBR; k++ )  {
	long double invPairs     = 1.0/((long double) totNumPairs);
	long double percentPair  = 0.0;

        if( cnStatTblP->NumDiffPairs[i][j][k] <= 0 ) continue;

	percentPair = ( ((long double) (cnStatTblP->NumDiffPairs[i][j][k])) * 
		       invPairs ) * 100.00;

        /* percentPair  = (((float)  (cnStatTblP->NumDiffPairs[i][j][k]))/
	   ((float)  totNumPairs) );
	   percentPair *= ((float)  100.0); */

        fprintf( fp,"\t(%d, %d, %d)\t\t%u\t\t\t\t%5.2f\n", i, j, k, 
                 cnStatTblP->NumDiffPairs[i][j][k], (float) percentPair );
      }  /* for(i:j:k) */

  fprintf(fp, "*****************************************************\n\n");


  l = 0;
  /* Number of ICOSAHEDRA Atoms  */
  for( ia = 0; ia < nAtoms; ia++ )   {
    register Atom_t *ap = AtomP + ia;


    if( ap->tags[2] == 'I' )  {
      
      /* Print header: ONCE */
      if( !l )
        fprintf(fp,"************%u-ATOMS IN FULL ICOSAHEDRA***********\nID = ",
                 cnStatTblP->NatmInFullICO );

      l++;
      fprintf( fp, "%u\t", (uint4) ap->id );
      if(!(l%6) ) fprintf( fp, "\n" );  /* break into columns of 6 */
    }
  }  /* for( ia = 0; ia < nAtoms; ia++ ) */

  if( l )
    fprintf(fp, "\n*****************************************************\n\n");

  l = 0;

#ifdef REMOVE_THIS
  /* Number of DISTORTED ICOSAHEDRA Atoms  */
  for( ia = 0; ia < nAtoms; ia++ )   {
    register Atom_t *ap = AtomP + ia;

    if( ap->tags[2] == 'D' ) {

      /* Print header: ONCE */
      if( !l ) {
        fprintf(fp,"*********%u-ATOMS IN DISTORTED ICOSAHEDRA*********\nID = ",
		cnStatTblP->NatmInDistICO );
      }

      l++;
      fprintf( fp, "%u\t", (uint4) ap->id );
      if(!(l%6) ) fprintf( fp, "\n" );  /* break into columns of 6 */
    }
  }  /* for( ia = 0; ia < nAtoms; ia++ ) */

  if( l )
    fprintf(fp, "\n*****************************************************\n\n");
#endif

  /* Get Finish Time for CN analysis */
  GetFinishTime( fp );

  fflush(fp);
  fclose(fp);
}  /* void WriteCoordAndPairsInfo() */
/*---------------------------------------------------------------------------*/
/***************************************************************************
 * (1a) Compute UnitVector (u) normal to Plane containing Vectors 'a' and  *
 *       'b'. Here, a and b are simulation box vectors.                    *
 ***************************************************************************/
static void Normal( float *a, float *b, float *u )
{
  float p = 0.0;
  u[0]  = a[1]*b[2] - a[2]*b[1];
  u[1]  = a[2]*b[0] - a[0]*b[2];
  u[2]  = a[0]*b[1] - a[1]*b[0];
  p     = sqrt( (double) (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) );
  u[0] /= p, u[1] /= p, u[2] /= p;
} /* Normal() */
/*---------------------------------------------------------------------------*/
/***************************************************************************
 * (1b) Compute the number of slots along X/Y/Z directions                 *
 ***************************************************************************
 *          rcut - Value of cut-off (= slot size)                          *
 ***************************************************************************/
static void NumSlotsXYZ(float *boxP, float rcut, 
			int4 *sxP, int4 *syP, int4 *szP, int4 verbose )
{
  float a[3], b[3], c[3];
  float v[3], dp;
  int4  sX, sY, sZ;


  /* COLUMNS of boxDatap->hpresent[] are A[], B[] and C[], the  Edge vectors.
     THE ROWS WERE ASSIGNED (ERRONEOUSLY) BEFORE AND HELL BROKE LOOSE */
  a[0] = boxP[0];                                         /* Box A-Vector */
  a[1] = boxP[3];
  a[2] = boxP[6];

  b[0] = boxP[1];                                         /* Box B-Vector */
  b[1] = boxP[4];
  b[2] = boxP[7];

  c[0] = boxP[2];                                         /* Box C-Vector */
  c[1] = boxP[5];
  c[2] = boxP[8];

  /* Num. Cells Along normal to b and c */
  Normal( b, c, v );
  dp = ( a[0]*v[0] + a[1]*v[1] + a[2]*v[2] ) / rcut;
  sX = (int4) (dp < 0.0 ? -dp: dp);

  /* Num. Cells Along normal to c and a */
  Normal( c, a, v );
  dp = ( b[0]*v[0] + b[1]*v[1] + b[2]*v[2] ) / rcut;
  sY = (int4) (dp < 0.0 ? -dp: dp);

  /* Num. Cells Along normal to a and b */
  Normal( a, b, v );
  dp = ( c[0]*v[0] + c[1]*v[1] + c[2]*v[2] ) / rcut;
  sZ = (int4) (dp < 0.0 ? -dp: dp);

#ifdef DYNAMIC_MEM_ALLOC
  if((sX < 2) || (sX >= MAX_CELLS_X) ||
     (sY < 2) || (sY >= MAX_CELLS_Y) ||
     (sZ < 2) || (sZ >= MAX_CELLS_Z) )
    Error("$$$ NumSlotsXYZ(): Invalid number of cells (%d %d %d)\n", 
	  sX, sY, sZ );
#endif


  if((sX < 2) || (sY < 2) || (sZ < 2)  )
    Error("$$$ NumSlotsXYZ(): Invalid number of cells (%d %d %d)\n", 
	  sX, sY, sZ );



  /* Initialize BoxData.slotsX/Y/Z fields */
  (*sxP) = sX;
  (*syP) = sY;
  (*szP) = sZ;

  if( verbose )  {
    printf("$$$ NumSlotsXYZ(): SlotsX/Y/Z\t\t\t=\t(%3d %3d %3d)\n",
	   sX, sY, sZ);
    printf("$$$ NumSlotsXYZ(): RcutMax\t\t\t=\t%7.4f\n", (float) rcut);
  }

}  /* static void NumSlotsXYZ(float *boxP, float rcut, 
                              int4 *sxP, int4 *syP, int4 *szP, int4 verbose) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read LASS Checkpoint/Snap file. These files have following fields:   *
 *         ** Single/DOUBLE precision: (X,Y,Z,Vx,Vy,Vz)                     *
 *         ** uint: AtomicNumber                                            *
 *         ** Char:  (ColorTag, MovementTag)                                *
 *         ** UINT: AtomID                                                  *
 *                                                                          *
 * (B) Check-Point File is used to restart a simulation.                    *
 ****************************************************************************/
int4 ReadConLASS(char *ifnameP, float *boxP, float rcut, int4 verbose )
{
  float      boxVec[9];
  float      boxInv[9];
  float      boxVol            =  1.0;
  float      snaptime          = -1.0;

  char        inpBuf[BUF_SIZE];
  FILE        *ifp              = NULL;    /* Turab Input file ptr */ 

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0;
  int4        numAtomsProcessed = 0;

  int4        totNumAtoms  = 0;
  int4        numBluAtoms  = 0;
  int4        numRedAtoms  = 0;
  int4        numGrnAtoms  = 0;

  int4        nSlotsX      = 0,
              nSlotsY      = 0, 
              nSlotsZ      = 0;

  /* Read Input file (LASS .CK/.SN  Format) */
  if( !(ifp = fopen(ifnameP, "rb") ) )
    Error("1. ReadConLASS(): Unable to Open InFile: %s\n", ifnameP );


  /* Read Input *.CK/*.SN File Header */
  FREAD_STRUCT( FileHeaderP );

  /* Read Header and save for debug testing */
  snaptime    = (float) FileHeader.TimeOfSnap;
  totNumAtoms = FileHeader.TotNumAtoms;
  numBluAtoms = FileHeader.TotBluAtoms;
  numRedAtoms = FileHeader.TotRedAtoms;
  numGrnAtoms = FileHeader.TotGrnAtoms;

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConLASS(): Time(%lf)/Dim(%u)/TotAtoms(%u)/Blu(%u)/Red(%u)/Grn(%u)\n",
	    FileHeader.TimeOfSnap,  FileHeader.Dimension, 
	    FileHeader.TotNumAtoms, FileHeader.TotBluAtoms,
	    FileHeader.TotRedAtoms, FileHeader.TotGrnAtoms);
    
    fprintf(stderr, 
	    "$$$ ReadConLASS(): CooFlag(%u)/Vflag(%u)/PEflag(%u)/StrFlag(%u)/NstrComp(%u)\n",
	    FileHeader.CoordFlag, FileHeader.VelFlag, FileHeader.PEFlag,
	    FileHeader.StressFlag, FileHeader.NstressComp);

    fprintf(stderr, 
	    "$$$ ReadConLASS(): TypeFlag(%u)/Colorflag(%u)/Moveflag(%u)/InitCoordFlag(%u)\n",
	    FileHeader.TypeFlag, FileHeader.ColorFlag, FileHeader.MoveFlag,
	    FileHeader.InitCoordFlag );

    fprintf(stderr, "$$$ ReadConLASS(): SpinFlag(%u)/padding(%u)\n",
	    FileHeader.SpinFlag, FileHeader.padding );
    fflush( stderr );
  }

  if( FileHeader.StressFlag )
    Error("$$$ ReadConLASS(): CNA code cannot handle stress in SN/CK file\n");

  /* Read Box Edge Vectors */
  FREAD_DOUBLE_PRECISION( 9 ); 

  for( i = 0; i < 9; i++ )  BoxVectors[i] = (float) Buf_DO[i];

  /* TRUE box matrix = Transpose of the matrix read in */
  TRANSPOSE_MATRIX( boxP, BoxVectors );

  ASSIGN_MATRICES( boxVec, boxP )

  /* Compute Box matrix Inverse */
  INVERT_MATRIX( boxInv, boxVec, boxVol );

  /* Compute TRUE number of slots */
  NumSlotsXYZ(boxVec, rcut, &nSlotsX, &nSlotsY, &nSlotsZ, verbose );

  /* Initialize true number of slots */
  LCInitializeSlots(nSlotsX, nSlotsY, nSlotsZ, totNumAtoms, verbose );

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConLASS(): Hvec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    Buf_DO[0],Buf_DO[1],Buf_DO[2],Buf_DO[3],Buf_DO[4],Buf_DO[5],
	    Buf_DO[6],Buf_DO[7],Buf_DO[8]);
    fflush( stderr );
  }

  /* Read Box Edge Vectors Derivatives */
  FREAD_DOUBLE_PRECISION( 9 ); 
  for( i = 0; i < 9; i++ )  BoxVecDeri[i] = (float) Buf_DO[i];

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConLASS(): Hdot\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    Buf_DO[0],Buf_DO[1],Buf_DO[2],Buf_DO[3],Buf_DO[4],Buf_DO[5],
	    Buf_DO[6],Buf_DO[7],Buf_DO[8]);
    fflush( stderr );
  }
  
  
  /* Loop, and Read Parameters */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    float    pe;
    Vector_t xyz;                       /* FLOAT coord of grid points   */
    Vector_t xyz0;                      /* Initial coordinates          */
    Vector_t vel;

    uint4    atomicNum;
    uint4    id;                        /* Atom ID                      */
    uint4    arrIn;                     /* Index of atom's position 
					   in AtomP[]                   */
    char     tags[4] = {'M', 'B'};      /* [0] = movetag, [1] = color;  */
    
    numAtomsProcessed++;

    if( numAtomsProcessed%500000 == 0 )
      fprintf( stderr, "$$$ReadConSPaSM(): Read %u atoms\n",
	      numAtomsProcessed );

    /* A1. Read X/Y/Z in DOUBLE PRECISION */
    switch( FileHeader.CoordFlag )   {
    case 1:
      FREAD_SINGLE_PRECISION( 3 );
      xyz.x = (float) Buf_FL[0];
      xyz.y = (float) Buf_FL[1];
      xyz.z = (float) Buf_FL[2];
      break;
    case 2:
      FREAD_DOUBLE_PRECISION( 3 );
      xyz.x = (float) Buf_DO[0];
      xyz.y = (float) Buf_DO[1];
      xyz.z = (float) Buf_DO[2];
      break;
    default: Error("4a. ReadConLASS(): Invalid Precision for Coord");
    }
    
    /* A2. Read initial CARTESIAN coord X0/Y0/Z0 in SINGLE PRECISION */
    switch( FileHeader.InitCoordFlag )   {
    case 0:
      ;                                /* DO NOTHING */
      break;
    case 1:                           /* Single Precision */
      FREAD_SINGLE_PRECISION( 3 );
      xyz0.x = (float) Buf_FL[0];
      xyz0.y = (float) Buf_FL[1];
      xyz0.z = (float) Buf_FL[2];
      break;
    default: Error("4b. ReadConLASS(): NOT Single Precision Coord");
    }
    
    /* A3. Read Sx/Sy/Sz  */
    switch( FileHeader.SpinFlag )   {
    case 0:
      ;                                /* DO NOTHING */
      break;
    default: Error("4c. ReadConLASS(): Invalid SpinFlag");
    }
    
    /* B1. Read VX/Y/Z in DOUBLE PRECISION */
    switch( FileHeader.VelFlag )   {
    case 0:
      ;         /* Do nothing */
      break;
    case 2:
      FREAD_DOUBLE_PRECISION( 3 );
      vel.x = (float) Buf_DO[0];
      vel.y = (float) Buf_DO[1];
      vel.z = (float) Buf_DO[2];
      break;
    default: Error("5a. ReadConLASS(): NOT Double Precision Velocity");
    }

    /* B2. Read PotentialEnergy in DOUBLE PRECISION */
    switch( FileHeader.PEFlag )   {
    case 0:
      ;         /* Do nothing */
      break;
    case 1:
      FREAD_SINGLE_PRECISION( 1 );
      pe = (float) Buf_FL[0];
      break;
    case 2:
      FREAD_DOUBLE_PRECISION( 1 );
      pe = (float) Buf_DO[0];
      break;
    default: Error("5b. ReadConLASS(): NOT Single/Double Precision PE");
    }

    /* C. Read UINT AtomicNumber Flag */
    switch( FileHeader.TypeFlag )   {
    case 1:
      FREAD_UINT( 1 );
      atomicNum = Buf_UI[0];
      break;
    default: Error("6. ReadConLASS(): Need to Read AtomicNumber");
    }
    
    /* D. Read CHAR Color  Tags */
    switch( FileHeader.ColorFlag )   {
    case 1:
      FREAD_CHARACTERS( 1 );
      /* see we use [1] instead of [0] */
      tags[1] = Buf_CH[0];
      break;
    default: Error("7b. ReadConLASS(): Need to Read Color Tag");
    }
    
    /* E. Read CHAR Movement Tags */
    switch( FileHeader.MoveFlag )   {
    case 0:
      break;                 /* No flag => skip */
    case 1:
      FREAD_CHARACTERS( 1 );
      /* see we use [0]instead of [1]  */
      tags[0] = Buf_CH[0];

      break;
    default: Error("8. ReadConLASS(): Need to Read Movement Flag");
    }
    
    /* F. Read Atom ID : UINT */
    FREAD_UINT( 1 );
    id = Buf_UI[0];

    /****** DEBUG ********
    fprintf( stderr, "ReadLassCon(): Id=%u\n", id );
    fflush( stderr );
    ****** DEBUG ********/
    
    /* Insert into link-cell slots: RETURN value is index of atom's
       position in AtomP[] (ignored for now) */
    arrIn = LCInsertAtomIntoSlot(boxVec, boxInv, &xyz, pe, id, atomicNum,
				 tags, verbose );

  }   /* Loop Over Atoms.               */
  
  fclose(ifp);


  /* if( verbose ) DbgPrintReadInAtoms( totNumAtoms ); */

  printf("$$$ ReadConLASS():TotAtoms/Nblu/red/green=(%d %d %d %d)\n",
	 totNumAtoms, numBluAtoms, numRedAtoms, numGrnAtoms );

  return( numAtomsProcessed );
} /* int ReadConLASS(char *ifnameP, float *boxP, float rcut, int4 verbose ) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Checkpoint file Function. These files have following fields:         *
 *         ** DOUBLE precision: (X,Y,Z)                                     *
 *         ** uint: AtomicNumber                                            *
 *         ** Char:  (ColorTag, MovementTag)                                *
 *         ** UINT: AtomID                                                  *
 *                                                                          *
 * (B) Check-Point File is used to restart a simulation.                    *
 ****************************************************************************/
void WriteConLASS( char *outFil, int verbose )
{
  float        snaptime          = -1.0;
  char         inpBuf[1024];
  FILE         *ofp              = NULL;    /* LASS  CHKPT file ptr */
  FILE         *tfp              = NULL;    /* Temp Text file ptr */
  int4         ai                = 0;
  int4         i                 = 0;
  int4         j                 = 0;
  int4         movement_tag      = 1;
  int4         numAtomsProcessed = 0;
  int4         totNumAtoms       = 0;
  int4         numBluAtoms       = 0;
  int4         numRedAtoms       = 0;
  int4         numGrnAtoms       = 0;

  if( verbose )   {
    if( FileHeader.TotNumAtoms < 100000 ) {
      if( !(tfp = fopen("DbgCNAfile.txt", "w") ) )    {
	Error("3. WriteConLASS(): Unable to Open DbgCNAfile.txt file\n");
      }
    }
  }

  /* Write output file (CK Format) after rotating Blue Atoms */
  if( !(ofp = fopen(outFil, "wb") ) )    {
    fprintf( stderr, "WriteCkPt(): Unable to Open Output File: %s\n", outFil );
    fflush( stderr );
    Error("3. WriteConLASS(): Unable to open file");
  }

  /********* Always write SINGLE precision Coord ***********/
  FileHeader.CoordFlag      = 1;

  /* Velocity, Initial Atom Coord, Spins, Stress are never printed out here */
  FileHeader.InitCoordFlag  = 0;         
  FileHeader.SpinFlag       = 0;
  FileHeader.VelFlag        = 0;   
  FileHeader.StressFlag     = 0;
  FileHeader.NstressComp    = 0;

  /* For small systems PE, atomicNum, color and move tags may be written.
     For NOT small systems, these flags are NOT written */
#ifndef SMALL_SYSTEM
  FileHeader.PEFlag         = 0;    
#endif

  /* Write the File Header */
  FWRITE_STRUCT( FileHeaderP );
    
  /* Read Header and save for debug testing */
  snaptime    = (float) FileHeader.TimeOfSnap;
  totNumAtoms = FileHeader.TotNumAtoms;
  numBluAtoms = FileHeader.TotBluAtoms;
  numRedAtoms = FileHeader.TotRedAtoms;
  numGrnAtoms = FileHeader.TotGrnAtoms;

  if( verbose && (tfp != NULL) )   {
    fprintf(tfp, "Time(%lf)/Dim(%u)/TotAtoms(%u)/Blu(%u)/Red(%u)/Grn(%u)\n",
	    FileHeader.TimeOfSnap,  FileHeader.Dimension, 
	    FileHeader.TotNumAtoms, FileHeader.TotBluAtoms,
	    FileHeader.TotRedAtoms, FileHeader.TotGrnAtoms);

    fprintf(tfp, "CooFlag(%u)/Vflag(%u)/PEflag(%u)/StrFlag(%u)/NstrComp(%u)\n",
	    FileHeader.CoordFlag, FileHeader.VelFlag, FileHeader.PEFlag,
	    FileHeader.StressFlag, FileHeader.NstressComp);
    
    fprintf(tfp, "TypeFlag(%u)/Colorflag(%u)/Moveflag(%u)/InitCoordFlag(%u)\n",
	    FileHeader.TypeFlag, FileHeader.ColorFlag, FileHeader.MoveFlag,
	    FileHeader.InitCoordFlag );
    
    fprintf(tfp, "SpinFlag(%u)/padding(%u)\n",
	    FileHeader.SpinFlag, FileHeader.padding );
    fflush( tfp );
  }

  for( i = 0; i < 9; i++ )  Buf_DO[i] = (double) BoxVectors[i];

  /* Write Box Edge Vectors */
  FWRITE_DOUBLE_PRECISION( 9 );

  if( verbose && (tfp != NULL) )   {
    fprintf(tfp,
	    "Hvec\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    Buf_DO[0],Buf_DO[1],Buf_DO[2],Buf_DO[3],Buf_DO[4],Buf_DO[5],
	    Buf_DO[6],Buf_DO[7],Buf_DO[8]);
  }

  /* Write Box Edge Vectors Derivatives */
  for( i = 0; i < 9; i++ )  Buf_DO[i] = (double) BoxVecDeri[i];
  FWRITE_DOUBLE_PRECISION( 9 );

  if( verbose && (tfp != NULL) )   {
    fprintf(tfp, 
	    "Hdot\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    Buf_DO[0],Buf_DO[1],Buf_DO[2],Buf_DO[3],Buf_DO[4],Buf_DO[5],
	    Buf_DO[6],Buf_DO[7],Buf_DO[8]);
  }
  
  /* Loop and Print out Parameters */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    numAtomsProcessed++;

    if( verbose && (tfp != NULL) )   {
      int4 il;
      
      int4 pairs = (int4) AtomP[ai].pairs;

      (void) fprintf(tfp, "%7d %2d", AtomP[ai].id, (int4) AtomP[ai].nNbors);
      (void) fprintf(tfp, " %8.4f %8.4f %8.4f  ", (float) AtomP[ai].coo.x, 
		    (float) AtomP[ai].coo.y, (float) AtomP[ai].coo.z );
      for (il = 0; il < pairs; il++)       {         
	(void) fprintf(tfp, "%1d",   (int4) AtomP[ai].index[il][0]);
	(void) fprintf(tfp, "%1d",   (int4) AtomP[ai].index[il][1]);
	(void) fprintf(tfp, "%1d  ", (int4) AtomP[ai].index[il][2]);
      }
      (void) fprintf(tfp, "\n");
      fflush( tfp );
    }

    /* A1. Print X/Y/Z in SINGLE PRECISION */
    switch( FileHeader.CoordFlag )   {
    case 1:
      Buf_FL[0] = (float) AtomP[ai].coo.x;
      Buf_FL[1] = (float) AtomP[ai].coo.y;
      Buf_FL[2] = (float) AtomP[ai].coo.z;
      FWRITE_SINGLE_PRECISION( 3 );
      break;
    default: Error("A1. WriteConLASS(): NOT SINGLE Precision Coord");
    }
    
    /* A2. Print initial CARTESIAN coord X0/Y0/Z0 in SINGLE PRECISION */
    switch( FileHeader.InitCoordFlag )   {
    case 0:
      ;                                /* DO NOTHING */
      break;
    default: Error("A2. WriteConLASS(): Init Coord Should not occur");
    }
    
    /* A3. Print Sx/Sy/Sz  */
    switch( FileHeader.SpinFlag )   {
    case 0:
      ;                                /* DO NOTHING */
      break;
    default: Error("A3. WriteConLASS(): Spins Should NOT occur");
    }
    
    /* B1. Print VX/Y/Z in DOUBLE PRECISION */
    switch( FileHeader.VelFlag )   {
    case 0:
      ;         /* Do nothing */
      break;
    default: Error("B1. WriteConLASS(): Velocity Should NOT occur");
    }

    /* B2. Read PotentialEnergy in DOUBLE PRECISION */
    switch( FileHeader.PEFlag )   {
    case 0:
      ;         /* Do nothing */
      break;
    case 1:
#ifdef SMALL_SYSTEM
      Buf_FL[0] = (float) AtomP[ai].potE;
      FWRITE_SINGLE_PRECISION( 1 );
#else
      Error("B2_1. WriteConLASS(): PE NOT saved for large systems\n");
#endif
      break;
    default: Error("B2. ReadConLASS(): Invalid PEFlag (must be 0 or 1)");
    }


    /* C. Print UINT AtomicNumber Flag */
    switch( FileHeader.TypeFlag )   {
    case 0:
      break;                      /* Do nothing */
    case 1:
      /* Buf_UI[0] = (uint4) AtomP[ai].atomicNum; */
      Buf_UI[0] = (uint4) AtomP[ai].cnTag;
      FWRITE_UINT( 1 );
      break;
    default: Error("C. WriteConLASS(): Need to Print AtomicNumber/CNtag");
    }
    
    /* D. Print CHAR Color  Tags */
    switch( FileHeader.ColorFlag )   {
    case 0:
      break;                      /* Do nothing */
    case 1:
      /* note it is: tags[1] */
      Buf_CH[0] = AtomP[ai].tags[1];
      FWRITE_CHARACTERS( 1 );
      break;
    default: Error("D. WriteConLASS(): Need to Print Color Tag");
    }
    
    /* E. Print CHAR Movement Tags */
    switch( FileHeader.MoveFlag )   {
    case 0:
      break;                 /* No flag => skip */
    case 1:
      /* note it is: tags[0] */
      Buf_CH[0] = AtomP[ai].tags[0];
      FWRITE_CHARACTERS( 1 );
      break;
    default: Error("E. WriteConLASS(): Invalid Movement Flag");
    }
    
    /* L. Print Atom ID : UINT */
    Buf_UI[0] = AtomP[ai].id;
    FWRITE_UINT( 1 );
  }   /* Loop Over Atoms.               */
  
  fclose(ofp);

  if( verbose && (tfp != NULL) ) fclose(tfp);

  if( numAtomsProcessed != totNumAtoms )   {
    char buf[512];
    
    sprintf(buf, "rm -f %s\0", outFil );
    system( buf );           /* Remove file */
    
    Error("******* Error: (B) Write(): numAtomsProcessed != totNumAtoms");
  }

  if( verbose )
    fprintf(stderr, 
	    "$$$ WriteConLASS():TotAtoms/Nblu/red/green=(%d %d %d %d)\n",
	    totNumAtoms, numBluAtoms, numRedAtoms, numGrnAtoms );

} /*void WriteConLASS() */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read SPaSM file. These files have following fields:                  *
 *         ** SpaSM-1 Format: int id + 6 doubles (X,Y,Z,Vx,Vy,Vz)           *
 *         ** SpaSM-2 Format: int id + 6 doubles (X,Y,Z,X0,Y0,Z0Vx,Vy,Vz)   *
 ****************************************************************************/
void ReadConSPaSM(char *ifnameP, float *boxP, float rcut, 
		  uint4 nAtoms,  int4  fileFormat, char *coordToCarveDomainP,
		  float minCoo, float maxCoo, int4 verbose )
{
  float      boxVec[9];
  float      boxInv[9];
  float      boxVol            =  1.0;

  char        inpBuf[BUF_SIZE];
  FILE        *ifp              = NULL;    /* Turab Input file ptr */ 

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0;
  int4        numAtomsProcessed = 0;
  uint4       numInsertedInList = 0;

  int4        totNumAtoms       = nAtoms;

  int4        nSlotsX           = 0,
              nSlotsY           = 0, 
              nSlotsZ           = 0;

  int4        removeFlag         = 0;                  /* Don't carve */

  /* Set file header for writing out CK/SN files */
  FileHeader.TimeOfSnap     = 0.0;
  FileHeader.Dimension      = 3;
  FileHeader.TotNumAtoms    = nAtoms;
  FileHeader.TotBluAtoms    = nAtoms;
  FileHeader.TotRedAtoms    = 0;
  FileHeader.TotGrnAtoms    = 0;
  FileHeader.CoordFlag      = 1;        /* Single precision coord */
  FileHeader.VelFlag        = 0;   
  FileHeader.PEFlag         = 0;    
  FileHeader.StressFlag     = 0;
  FileHeader.NstressComp    = 0;
  FileHeader.TypeFlag       = 1;        /* Print CN tag           */
  FileHeader.ColorFlag      = 0;
  FileHeader.MoveFlag       = 0;
  FileHeader.InitCoordFlag  = 0;         
  FileHeader.SpinFlag       = 0;
  FileHeader.padding        = 0;


  /* Read Input file (LASS .CK/.SN  Format) */
  if( !(ifp = fopen(ifnameP, "rb") ) )
    Error("1. ReadConSPaSM(): Unable to Open InFile: %s\n", ifnameP );

  /* Box Vectors */
  for( i = 0; i < 9; i++ )  {

    BoxVecDeri[i] = 0.0;
    BoxVectors[i] = (float) boxP[i];
    
    /* if(i == 1 || i == 2 || i == 3 ||
       i == 5 || i == 6 || i == 7 )
       if( ((float) boxP[i]) != 0.0 )
       Error("ReadConSPaSM(): Not a orthogonal box\n"); */
  }


  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConSPaSM(): BoxVec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    boxP[0],boxP[1],boxP[2],boxP[3],boxP[4],boxP[5],
	    boxP[6],boxP[7],boxP[8]);
    fflush( stderr );
  }



  /* TRUE box matrix = Transpose of the matrix read in */
  TRANSPOSE_MATRIX( boxVec, BoxVectors );

  /* Compute Box matrix Inverse */
  INVERT_MATRIX( boxInv, boxVec, boxVol );

  /* Compute TRUE number of slots */
  NumSlotsXYZ(boxVec, rcut, &nSlotsX, &nSlotsY, &nSlotsZ, verbose );

  /* Initialize true number of slots */
  LCInitializeSlots(nSlotsX, nSlotsY, nSlotsZ, totNumAtoms, verbose );

  /* Loop, and Read Parameters */
  for(; !feof(ifp); )   {
    /* for( ai = 0; ai < totNumAtoms; ai++ )   { */
    float    pe         = 0.0;
    Vector_t xyz;                       /* double coordinates of grid points */
    Vector_t xyz0;                      /* double coordinates of grid points */
    Vector_t vel;                       /* double coordinates of grid points */
    uint4    atomicNum  = 0;
    uint4    id         = 0;            /* Atom ID                           */
    uint4    arrIn      = 0;            /* Index of atom position in AtomP[] */
    char     tags[4]    = {'M', 'B'};   /* [0] = movetag, [1] = color;       */
    
    if( numAtomsProcessed%500000 == 0 )
      fprintf( stderr, "$$$ReadConSPaSM(): Read %u atoms\n",
	      numAtomsProcessed );



#ifdef DBG
    /* Read Atom ID */
    FREAD_UINT( 1 );
    id = Buf_UI[0];
#endif
    
    id = numAtomsProcessed+1;

    if( id > nAtoms ) break;


    switch( fileFormat )   {
    case SPASM_1_FORMAT:
      FREAD_SINGLE_PRECISION_SPASM( 5 );
      FileHeader.CoordFlag      = 1;
      xyz.x = (float) Buf_FL[0];        /* get X/Y/Z                         */
      xyz.y = (float) Buf_FL[1];
      xyz.z = (float) Buf_FL[2];
      vel.x = (float) Buf_FL[3];        /* get Vx                            */
      pe    = (float) Buf_FL[4];
      /* vel.z = (float) Buf_DO[5]; */

      switch( coordToCarveDomainP[0] )   {
      case 'N':                                 /* Don't Carve */
	removeFlag = 0;
	break;

      case 'X':                                 /* Use X-coord */
      case 'x':

	removeFlag = 1;

	if( Buf_FL[0] > minCoo && Buf_FL[0] <= maxCoo ) removeFlag = 0;
	break;

      case 'Y':                                 /* Use X-coord */
      case 'y':

	removeFlag = 1;

	if( Buf_FL[1] > minCoo && Buf_FL[1] <= maxCoo ) removeFlag = 0;
	break;

      case 'Z':                                 /* Use X-coord */
      case 'z':

	removeFlag = 1;

	if( Buf_FL[2] > minCoo && Buf_FL[2] <= maxCoo ) removeFlag = 0;
	break;
	
      default: Error("$$$ReadConSPaSM(): Unknown carve region tag\n");
      }
      break;

    case SPASM_2_FORMAT:
      FREAD_UINT( 2 );                  /* Tags are read and ignored         */
      FREAD_DOUBLE_PRECISION( 9 );
      xyz0.x  = (float) Buf_DO[0];      /* get X/Y/Z                         */
      xyz0.y  = (float) Buf_DO[1];
      xyz0.z  = (float) Buf_DO[2];
      xyz.x   = (float) Buf_DO[3];      /* get X0/Y0/Z0                      */
      xyz.y   = (float) Buf_DO[4];
      xyz.z   = (float) Buf_DO[5];
      vel.x   = (float) Buf_DO[6];      /* get Vx/Vy/Vz                      */
      vel.y   = (float) Buf_DO[7];
      vel.z   = (float) Buf_DO[8];
      break;
    break;
    default:
      Error("$$$ReadConSPaSM(): Unknown file format\n");
    }

    /* Insert into link-cell slots: RETURN value is index of atom's
       position in AtomP[] (ignored for now) */
    if( !removeFlag )  {
      numInsertedInList++;
      arrIn = LCInsertAtomIntoSlot(boxVec, boxInv, &xyz, pe, id, atomicNum,
				   tags, verbose );
    }

    numAtomsProcessed++;
  }   /* Loop Over Atoms.               */
  
  fclose(ifp);


  /* if( verbose ) DbgPrintReadInAtoms( totNumAtoms ); */

  printf("$$$ ReadConSPaSM():TotAtoms/Nblu=(%d %d)\n",
	 totNumAtoms, FileHeader.TotBluAtoms );

  if( numInsertedInList != nAtoms )
    Error("$$$ ReadConSPaSM(): numInsertedInList(%u) != nAtoms(%u)\n",
	  numInsertedInList, nAtoms );

} /* int ReadConSPaSM(char *ifnameP, float *boxP, float rcut, int4 verbose ) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read Mike Baskes's Snap file (Text Format). This is a surface        *
 *     configuration and we need to massage the data to reduce the number   *
 *     of linkcells(to decrease memory allocation). It has following fields:*
 *         ** First Line is Header (String)                                 *
 *         ** 2nd line: (NumAtoms, NumTypes, Float, Float)                  *
 *         ** 3rd line: (Max values of Box A/B/C edges)                     *
 *         ** 4th line: (Min values of Box A/B/C edges)                     *
 *         ** 5th/6th Lines: Read and IGNORE                                *
 *         ** From 7th Line Onwards:                                        *
 *               Line-1: X/Y/Z, Line-2: Vx/Vy/Vz, Line-3: AtomType(int)     *
 ****************************************************************************/
static int4 ReadConBaskesSurface(char *ifnameP, float rcut, int4 verbose )
{
  Vector_t    *cooDP  = NULL;
  Vector_t    *velDP  = NULL;
  int4        *typeIP = NULL;

  Vector_t    maxCoo  = {-999999., -999999., -999999.};
  float       edgeA[3], edge2A[3], edgeMinA[3];
  float       maxXYZ[3];
  float       minXYZ[3];

  float       boxVec[9];
  float       newBoxCenter[3];

  float       timeF             = -1.0;
  float       angF              = -1.0;
  float       massF[N_TYPE];
  int4        type[N_TYPE];

  char        inpBuf[1024];
  char        comment[1024];
  FILE        *ifp              = NULL; 
  FILE        *ofp              = NULL; 

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0;
  int4        numAtomsProcessed = 0;

  uint4       totNumAtoms       = 0;
  uint4       nType             = 1;

  /* Read Input file (LASS .CK/.SN  Format) */
  if( !(ifp = fopen(ifnameP, "r") ) )
    Error("1a. ReadConBaskesSfce(): Unable to Open InFile: %s\n", ifnameP );

  /* Scan First SIX Lines */
#ifdef DBG
  /* ** First Line is Header (String);                                *
   * ** 2nd line: (NumAtoms, NumTypes, Float, Float)                  *
   * ** 3rd line: (Max values of Box A/B/C edges)                     *
   * ** 4th line: (Min values of Box A/B/C edges)                     *
   * ** 5th/6th Lines: Read and IGNORE                                */
#endif

  /* Loop over 1st 4 lines in file */
  for( i = 1; i <= 4; i++ )    {

    /* Line-1-6: Get Header and Discard */
    fgets(inpBuf, sizeof(inpBuf), ifp );
  
    switch( i )    {
    case 1:                  /* 1st Line is comment                          */
      sscanf( inpBuf, "%s", comment );
      break;

    case 2:                  /* 2nd Line: (NumAtoms, NumTypes, Float, Float) */
      /* Line-2:  NumAtoms/Type/Float/Float*/
      sscanf( inpBuf, "%u %u %f %f", &totNumAtoms, &nType, &SnapTime, &angF );
      break;
    case 3:                  /* 3rd line: (Max values of Box A/B/C edges)    */
      sscanf( inpBuf, "%f %f %f", &maxXYZ[0], &maxXYZ[1], &maxXYZ[2] );      
      break;
    case 4:                  /* 4th line: (Min values of Box A/B/C edges)    */
      sscanf( inpBuf, "%f %f %f", &minXYZ[0], &minXYZ[1], &minXYZ[2] );
      break;
    default: Error("$$$ReadConBaskesSfce(): Problem Reading File Header\n");
    }
  }  /* for(i) */


  /* This is dynamo file => box origin is at box center */
  BaskesSurfaceConfigFlag   = YES;
  MoveOriginToBoxCenterFlag = YES;
  for( i = 0; i < 3; i++ ) {
    MaxBoxEdges[i] = maxXYZ[i];
    MinBoxEdges[i] = minXYZ[i];
  }

  /* Loop through Num Atom Types and scan mass/type info */
  if( nType >=  N_TYPE )
    Error("ReadConBaskesSfce(): Too many atomtypes (%d)\n", nType);
  for( i = 0; i < nType; i++ )  {
     fgets(inpBuf, sizeof(inpBuf), ifp );
     sscanf( inpBuf, "%f %d", &massF[i], &type[i]  );
  }

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConBaskesSfce(): Time(%lf)/TotAtoms(%u)\n",
	    SnapTime,  totNumAtoms );
    
    fprintf(stderr, "$$$ ReadConBaskesSfce(): MinBoxXYZ(%f %f %f)\n",
	    minXYZ[0], minXYZ[1], minXYZ[2] );

    fprintf(stderr, "$$$ ReadConBaskesSfce(): MaxBoxXYZ(%f %f %f)\n",
	    maxXYZ[0], maxXYZ[1], maxXYZ[2] );

    fflush( stderr );
  }


  /* Initialize Box Vectors */
  for( i = 0; i < 9; i++ )    {
    float lenB;

    boxVec[i] = 0.0;

    switch( i )    {
    case 0: 
      lenB      = (float) (minXYZ[0] - maxXYZ[0]);
      boxVec[0] = (float) (ABS(lenB));

      edgeA[0]     = boxVec[0];
      edge2A[0]    = 2.0*edgeA[0];
      edgeMinA[0]  = - edgeA[0];

      /* BoxVectors[0] = (float) (ABS(minXYZ[0]) + ABS(maxXYZ[0])); */
      break;

    case 4: 
      lenB      = (float) (minXYZ[1] - maxXYZ[1]);
      boxVec[4] = (float) (ABS(lenB));

      edgeA[1]     = boxVec[4];
      edge2A[1]    = 2.0*edgeA[1];
      edgeMinA[1]  = - edgeA[1];

      /* BoxVectors[4] = (float) (ABS(minXYZ[1]) + ABS(maxXYZ[1])); */
      break;

    case 8: 
      lenB      = (float) (minXYZ[2] - maxXYZ[2]);
      boxVec[8] = (float) (ABS(lenB));

      edgeA[2]     = boxVec[8];
      edge2A[2]    = 2.0*edgeA[2];
      edgeMinA[2]  = - edgeA[2];

      /* BoxVectors[8] = (float) (ABS(minXYZ[2]) + ABS(maxXYZ[2])); */
      break;
    }
  }


  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConBaskesSfce(): Hvec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    (float) boxVec[0], (float) boxVec[1], (float) boxVec[2],
	    (float) boxVec[3], (float) boxVec[4], (float) boxVec[5],
	    (float) boxVec[6], (float) boxVec[7], (float) boxVec[8]);
    fflush( stderr );
  }

  /* Allocate memory for atoms */
  cooDP = (Vector_t*) calloc( totNumAtoms, sizeof(Vector_t) );
  if( !cooDP )
      Error("ConBaskesSfce(): Cannot Allocate %u MBytes Memory for atom Coo\n",
	    ((totNumAtoms * sizeof(Vector_t))/1000000 + 1) );

  velDP = (Vector_t*) calloc( totNumAtoms, sizeof(Vector_t) );
  if( !velDP )
      Error("ConBaskesSfce(): Cannot Allocate %u MBytes Memory for atom Vel\n",
	    ((totNumAtoms * sizeof(Vector_t))/1000000 + 1) );

  typeIP = (int4*) calloc( totNumAtoms, sizeof(int4) );
  if( !typeIP )
      Error("ConBaskesSfce(): Cannot Allocate %u MBytes Memory for atomType\n",
	    ((totNumAtoms * sizeof(int4))/1000000 + 1) );


  /* Loop, and Read Parameters */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    Vector_t xyz    = {0., 0., 0.};     /* FLOAT coord of grid points   */
    double   coo[3] = {0., 0., 0.};     /* Tmp     coordinates          */
    Vector_t vel    = {0., 0., 0.};

    int4     pbcFlag[3] = {INSIDE_THE_BOX,INSIDE_THE_BOX,INSIDE_THE_BOX};
    uint4    arrIn;                     /* Index of atom's position 
					   in AtomP[]                   */
    
    numAtomsProcessed++;

    if( numAtomsProcessed%5000 == 0 )
      fprintf( stderr, "$$$ReadConBaskesSfce(): Read %u atoms\n",
	      numAtomsProcessed );

    /* A. Read X/Y/Z  */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%lf %lf %lf",  &coo[0], &coo[1], &coo[2] );

    /* Assuming the origin is at box center */
    /* This is dynamo file => box origin is at box center */
    if( ai == 0 ) {
      InitCooAtom01[0] = (float) coo[0];
      InitCooAtom01[1] = (float) coo[1];
      InitCooAtom01[2] = (float) coo[2];
    }

    /* Coord values must be between 0 and Box edge => Move origin to box 
       corner */
    cooDP[ai].x = ((float) -minXYZ[0]) + ((float) coo[0]);
    cooDP[ai].y = ((float) -minXYZ[1]) + ((float) coo[1]);
    cooDP[ai].z = ((float) -minXYZ[2]) + ((float) coo[2]);

    /* Atoms close to the edge could have escaped across the face (even though
       this is a free-surface system => translate all coordinates by +50 
       angstrom and then apply PBC to bring all atoms close to each other */
    cooDP[ai].x += 50.0;
    cooDP[ai].y += 50.0;
    cooDP[ai].z += 50.0;

    /* Apply Periodic Boundary Condition */
    APPLY_PBC_1D(cooDP[ai].x, edgeA[0], edge2A[0], edgeMinA[0], pbcFlag[0]);
    APPLY_PBC_1D(cooDP[ai].y, edgeA[1], edge2A[1], edgeMinA[1], pbcFlag[1]);
    APPLY_PBC_1D(cooDP[ai].z, edgeA[2], edge2A[2], edgeMinA[2], pbcFlag[2]);


    /* Get maximum  coordinate values */
    maxCoo.x = MAX( maxCoo.x, cooDP[ai].x );
    maxCoo.y = MAX( maxCoo.y, cooDP[ai].y );
    maxCoo.z = MAX( maxCoo.z, cooDP[ai].z );

    /* B. Read VX/Y/Z  */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%f %f %f", (float*) &velDP[ai].x, 
	   (float*) &velDP[ai].y, (float*) &velDP[ai].z );

    /* C. Read UINT Atom Type */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%u", (uint4*) &typeIP[ai] );

  }   /* Loop Over Atoms.               */
  
  fclose(ifp);

  
  /* Maximum value for box edge= 50+ maxCoo and still retain free-surface.
     Origin of new coordinate system is at this box center */
  newBoxCenter[0] = 0.5*(maxCoo.x + 50.0);
  newBoxCenter[1] = 0.5*(maxCoo.y + 50.0);
  newBoxCenter[2] = 0.5*(maxCoo.z + 50.0);

  /* Write NEW ouput file (BASKES'S  Format) */
  if( !(ofp = fopen("BaskesTmp.atm", "w") ) )
    Error("1b. ReadConBaskesSfce(): Unable to Open OutFile: BaskesTmp.atm\n" );


   /* Line-01: Comment line */
   fprintf(ofp, "%s\n", comment );

   /* Line-02: (NumAtoms, NumTypes, time, box-angle) */
   fprintf(ofp, "%10d%10d%15.8lE%10.5lf\n", totNumAtoms, nType,
           (double) SnapTime, (double) angF  );

   /* Line-03: (Max values of Box A/B/C edges) */
   fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", (double) newBoxCenter[0],
           (double) newBoxCenter[1], (double) newBoxCenter[2] );

   /* Line-04: (Min values of Box A/B/C edges) */
   fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", (double) -newBoxCenter[0],
           (double) -newBoxCenter[1], (double) -newBoxCenter[2] );


   /* Line-05: (Atomic mass (amu*1.0365e-4) and atomic number */
  /* Loop through Num Atom Types and scan mass/type info */
  for( i = 0; i < nType; i++ )  {
     fprintf( ofp, "%25.16lE%10d\n", (double) massF[i], type[i]  );
  }

  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
     double x, y, z;

     x = (double)  (cooDP[ai].x - newBoxCenter[0]);
     y = (double)  (cooDP[ai].y - newBoxCenter[1]);
     z = (double)  (cooDP[ai].z - newBoxCenter[2]);


     /* A. Write X/Y/Z  */
     fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", x, y, z );
     
     /* B. Write VX/Y/Z  */
     fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", (double) velDP[ai].x, 
	     (double) velDP[ai].y, (double) velDP[ai].z );
     
     /* C. Write Atom Type */
     fprintf(ofp, "%10d\n", typeIP[ai] );
  }

  fclose(ofp);

  /* Rename this file as the new atom file */
  /* sprintf(inpBuf, "mv baskes.tmp %s\0", ifnameP);
     system(inpBuf); */
  /* Use BaskesTmp.atm as atom file */
  sprintf( ifnameP, "BaskesTmp.atm\0" );

  /* Free Memory */
  free(cooDP );
  free(velDP );
  free(typeIP);

  return( numAtomsProcessed );
} /* int ReadConBaskesSfce(char *ifnameP, int4 verbose)          */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read Mike Baskes's Snap file (Text Format). It has following fields: *
 *         ** First Line is Header (String)                                 *
 *         ** 2nd line: (NumAtoms, NumTypes, Float, Float)                  *
 *         ** 3rd line: (Max values of Box A/B/C edges)                     *
 *         ** 4th line: (Min values of Box A/B/C edges)                     *
 *         ** 5th/6th Lines: Read and IGNORE                                *
 *         ** From 7th Line Onwards:                                        *
 *               Line-1: X/Y/Z, Line-2: Vx/Vy/Vz, Line-3: AtomType(int)     *
 ****************************************************************************/
int4 ReadConBaskes(char *ifnameP, float *boxP, float rcut, int4 verbose )
{
  float       maxXYZ[3];
  float       minXYZ[3];
  float       boxVec[9];
  float       boxInv[9];
  float       boxVol            =  1.0;
  float       tmpF              = -1.0;
  uint4       tmpI              = 0;
  const int4  zero              = 0;

  char        inpBuf[1024];
  FILE        *ifp              = NULL;    /* Turab Input file ptr */ 

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0;
  int4        numAtomsProcessed = 0;

  uint4       totNumAtoms       = 0;
  uint4       nType             = 1;


  int4        nSlotsX           = 0,
              nSlotsY           = 0, 
              nSlotsZ           = 0;

  /* For mike Baskes's Free-Surface system we need to reduce the number 
     of link-cells to lowere memory consumption. We first read in  the 
     original Baskes con file and adjust for minimum image criterion and
     shrink the box.  The new coordinate values are rewritten into the same 
     input atom filename  */
#if FREE_SFCE_X && FREE_SFCE_Y && FREE_SFCE_Z
  ReadConBaskesSurface(ifnameP, rcut, verbose );
#endif

  for(i = 0; i < 9; i++ ) {
    boxVec[i] = boxInv[i] = 0.0;
  }

  /* Read Input file (LASS .CK/.SN  Format) */
  if( !(ifp = fopen(ifnameP, "r") ) )
    Error("1. ReadConBaskes(): Unable to Open InFile: %s\n", ifnameP );

  /* Set file header for writing out CK/SN files */
  memset((void*) &FileHeader, zero, (size_t)sizeof(struct FHead) );

  FileHeader.Dimension      = 3;
  FileHeader.CoordFlag      = 1;        /* Single precision coord */
  FileHeader.TypeFlag       = 1;        /* Print CN tag           */


  /* Scan First SIX Lines */

#ifdef DBG
  /* ** First Line is Header (String);                                *
   * ** 2nd line: (NumAtoms, NumTypes, Float, Float)                  *
   * ** 3rd line: (Max values of Box A/B/C edges)                     *
   * ** 4th line: (Min values of Box A/B/C edges)                     *
   * ** 5th/6th Lines: (atommass, atomic number)                      */
#endif

  /* Loop over 1st 4 lines in file */
  for( i = 1; i <= 4; i++ )    {

    /* Line-1-6: Get Header and Discard */
    fgets(inpBuf, sizeof(inpBuf), ifp );
  
    switch( i )    {
    case 1:                  /* 1st Line is Header Ignored                   */
      break;

    case 2:                  /* 2nd Line: (NumAtoms, NumTypes, Float, Float) */
      /* Line-2:  NumAtoms/Type/Float/Float*/
      sscanf( inpBuf, "%u %u %f %f", &totNumAtoms, &nType, &SnapTime, &tmpF);

      /* SRINI Oct 6, 2004: If Ntype > 2, RESET this to 1 */
      if( nType > 2 ) Ntype = 1;
      else            Ntype = nType;

      break;
    case 3:                  /* 3rd line: (Max values of Box A/B/C edges)    */
      sscanf( inpBuf, "%f %f %f", &maxXYZ[0], &maxXYZ[1], &maxXYZ[2] );
      break;
    case 4:                  /* 4th line: (Min values of Box A/B/C edges)    */
      sscanf( inpBuf, "%f %f %f", &minXYZ[0], &minXYZ[1], &minXYZ[2] );
      break;
    default: Error("$$$ReadConBaskes(): Problem Reading File Header\n");
    }
  }  /* for(i) */

  /* This is dynamo file => box origin is at box center */
  MoveOriginToBoxCenterFlag = YES;

  if( BaskesSurfaceConfigFlag == NO ) {
    for( i = 0; i < 3; i++ ) {
      MaxBoxEdges[i] = maxXYZ[i];
      MinBoxEdges[i] = minXYZ[i];
    }
  }

  /* Loop through Num Atom Types and Discard Info */
  if( nType >=  N_TYPE )
    Error("ReadConBaskes(): Too many atomtypes (%d)\n", nType);
  for( i = 1; i <= nType; i++ )  {
    fgets(inpBuf, sizeof(inpBuf), ifp );
      sscanf(inpBuf, "%le %u", &AtomMass[i-1], &AtomicNumber[i-1] );
  }

  /* Initialize Total Atoms Fields */
  FileHeader.TotNumAtoms    = totNumAtoms;
  FileHeader.TotBluAtoms    = totNumAtoms;

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConBaskes(): Time(%lf)/Dim(%u)/TotAtoms(%u)/Blu(%u)/Red(%u)/Grn(%u)\n",
	    FileHeader.TimeOfSnap,  FileHeader.Dimension, 
	    FileHeader.TotNumAtoms, FileHeader.TotBluAtoms,
	    FileHeader.TotRedAtoms, FileHeader.TotGrnAtoms);
    
    fprintf(stderr, 
	    "$$$ ReadConBaskes(): CooFlag(%u)/Vflag(%u)/PEflag(%u)/StrFlag(%u)/NstrComp(%u)\n",
	    FileHeader.CoordFlag, FileHeader.VelFlag, FileHeader.PEFlag,
	    FileHeader.StressFlag, FileHeader.NstressComp);

    fprintf(stderr, 
	    "$$$ ReadConBaskes(): TypeFlag(%u)/Colorflag(%u)/Moveflag(%u)/InitCoordFlag(%u)\n",
	    FileHeader.TypeFlag, FileHeader.ColorFlag, FileHeader.MoveFlag,
	    FileHeader.InitCoordFlag );

    fprintf(stderr, "$$$ ReadConBaskes(): SpinFlag(%u)/padding(%u)\n",
	    FileHeader.SpinFlag, FileHeader.padding );

    fprintf(stderr, "$$$ ReadConBaskes(): MinBoxXYZ(%f %f %f)\n",
	    minXYZ[0], minXYZ[1], minXYZ[2] );

    fprintf(stderr, "$$$ ReadConBaskes(): MaxBoxXYZ(%f %f %f)\n",
	    maxXYZ[0], maxXYZ[1], maxXYZ[2] );

    fflush( stderr );
  }

  if( FileHeader.StressFlag )
    Error("$$$ ReadConBaskes(): CNA code can't handle stress in SN/CK file\n");


  /* Initialize Box Vectors */
  for( i = 0; i < 9; i++ )    {
    float lenB;
    BoxVectors[i] = 0.0;
    BoxVecDeri[i] = 0.0;           /* Read Box Edge Vectors Derivatives */

    switch( i )    {
    case 0: 
      lenB          = (float) (minXYZ[0] - maxXYZ[0]);
      BoxVectors[0] = (float) (ABS(lenB));
      /* BoxVectors[0] = (float) (ABS(minXYZ[0]) + ABS(maxXYZ[0])); */
      break;

    case 4: 
      lenB          = (float) (minXYZ[1] - maxXYZ[1]);
      BoxVectors[4] = (float) (ABS(lenB));
      /* BoxVectors[4] = (float) (ABS(minXYZ[1]) + ABS(maxXYZ[1])); */
      break;

    case 8: 
      lenB          = (float) (minXYZ[2] - maxXYZ[2]);
      BoxVectors[8] = (float) (ABS(lenB));
      /* BoxVectors[8] = (float) (ABS(minXYZ[2]) + ABS(maxXYZ[2])); */
      break;
    }
  }

  /* TRUE box matrix = Transpose of the matrix read in */
  TRANSPOSE_MATRIX( boxP, BoxVectors );

  ASSIGN_MATRICES( boxVec, boxP )

  /* Compute Box matrix Inverse */
  INVERT_MATRIX( boxInv, boxVec, boxVol );

  /* Compute TRUE number of slots */
  NumSlotsXYZ(boxVec, rcut, &nSlotsX, &nSlotsY, &nSlotsZ, verbose );

  /* Initialize true number of slots */
  LCInitializeSlots(nSlotsX, nSlotsY, nSlotsZ, totNumAtoms, verbose );

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConBaskes(): Hvec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
	    (float) boxVec[0], (float) boxVec[1], (float) boxVec[2],
	    (float) boxVec[3], (float) boxVec[4], (float) boxVec[5],
	    (float) boxVec[6], (float) boxVec[7], (float) boxVec[8]);
    fflush( stderr );
  }


  /* Loop, and Read Parameters */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    float    pe     = 0.0;
    Vector_t xyz    = {0., 0., 0.};     /* FLOAT coord of grid points   */
    double   coo[3] = {0., 0., 0.};     /* Tmp     coordinates          */
    Vector_t vel    = {0., 0., 0.};

    uint4    atomicNum = 18;
    uint4    id;                        /* Atom ID                      */
    uint4    arrIn;                     /* Index of atom's position 
					   in AtomP[]                   */
    char     tags[4] = {'M', 'B'};      /* [0] = movetag, [1] = color;  */
    
    numAtomsProcessed++;

    if( numAtomsProcessed%5000 == 0 )
      fprintf( stderr, "$$$ReadConBaskes(): Read %u atoms\n",
	      numAtomsProcessed );

    id = numAtomsProcessed;

    /* A. Read X/Y/Z  */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%lf %lf %lf",  &coo[0], &coo[1], &coo[2] );

    if( BaskesSurfaceConfigFlag == NO ) {
      /* This is dynamo file => box origin is at box center */
      if( ai == 0 ) {
	InitCooAtom01[0] = (float) coo[0];
	InitCooAtom01[1] = (float) coo[1];
	InitCooAtom01[2] = (float) coo[2];
      }
    }

    /* fprintf(stderr, "ReadBaskes(): %u, XYZ Before(%lf %lf %lf)\n",
       id, coo[0], coo[1], coo[2] ); */

    /* Coord values must be between 0 and Box edge => Move origin */
    xyz.x = ((float) -minXYZ[0]) + ((float) coo[0]);
    xyz.y = ((float) -minXYZ[1]) + ((float) coo[1]);
    xyz.z = ((float) -minXYZ[2]) + ((float) coo[2]);

    /* fprintf(stderr, "ReadBaskes(): %u, XYZ After(%lf %lf %lf)\n",
       id, xyz.x, xyz.y, xyz.z );
       fflush( stderr ); */

    /* B. Read VX/Y/Z  */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%f %f %f", 
	   (float*) &vel.x, (float*) &vel.y, (float*) &vel.z );

    /* C. Read UINT Atom Type */
    fgets(inpBuf, sizeof(inpBuf), ifp );
    sscanf( inpBuf, "%u", (uint4*) &atomicNum );

    /* Srini Oct 6, 2004: ap->atomicNum is atomtype ONLY for baskes dynamo format. 
       For ntypes > 2, reset this to 1  */
    /* SRINI Oct 6, 2004: If Ntype > 2, RESET this to 1 */
    if( nType > 2 ) atomicNum = 1;      

    /* Insert into link-cell slots: RETURN value is index of atom's
       position in AtomP[] (ignored for now) */
    arrIn = LCInsertAtomIntoSlot(boxVec, boxInv, &xyz, pe, id, atomicNum,
				 tags, verbose );
  }   /* Loop Over Atoms.               */
  
  fclose(ifp);


  /* if( verbose ) DbgPrintReadInAtoms( totNumAtoms ); */

  printf("$$$ ReadConBaskes(): Read in %u atoms\n", totNumAtoms );

  return( numAtomsProcessed );
} /* int ReadConBaskes(char *ifnameP, float *boxP, float rcut, int4 verbose) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * Write Mike Baskes's Snap file (Text Format). It has following fields: *
 *         ** First Line is Header (String)                                 *
 *         ** 2nd line: (NumAtoms, NumTypes, Float, Float)                  *
 *         ** 3rd line: (Max values of Box A/B/C edges)                     *
 *         ** 4th line: (Min values of Box A/B/C edges)                     *
 *         ** 5th/6th Lines: Read and IGNORE                                *
 *         ** From 7th Line Onwards:                                        *
 *               Line-1: X/Y/Z, Line-2: Vx/Vy/Vz, Line-3: AtomType(int)     *
 ****************************************************************************/
void WriteConBaskes(char *ofnameP, float *boxP, int4 nAtoms, int4 verbose )
{
  char        hdrStr[1024];
  FILE        *ofp              = NULL;
  FILE        *ofp2             = NULL;
  float       boxAngle  = 90.0;
  float       dispXYZ[3]= {0.0, 0.0, 0.0};
  float       maxXYZ[3] = {0.0, 0.0, 0.0};
  float       minXYZ[3] = {0.0, 0.0, 0.0};
  int4        nType     = 5; /* 5 types of atoms: BCC,FCC,HCP,ICO,OTHERS */
  float       atommass[10]       = {1.,1.,1.,1.,1.,1.};  /* Dummy mass   */
  int4        atomnum[10]        = {1,2,3,4,5,6};
  int4        typeTag           = 0;

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0, totNfcc = 0;
  int4        numAtomsProcessed = 0;
  int         atomType          = 1;

  void WriteHcpPairInfo(uint4, int4 verbose );
  void WriteBondInfo(uint4, int4 verbose );

  if( MoveOriginToBoxCenterFlag == YES ) {
    for( i = 0; i < 3; i++ ) {
      maxXYZ[i] = MaxBoxEdges[i];
      minXYZ[i] = MinBoxEdges[i];
    }
  }
  else {
    /*  Max bounds for box edges */
    maxXYZ[0] = boxP[0];
    maxXYZ[1] = boxP[4];
    maxXYZ[2] = boxP[8];
    
    /* SRINI Sept 4, 2002: Add non principal diagonal components */
    maxXYZ[0] += (boxP[3] + boxP[6]);
    maxXYZ[1] += (boxP[1] + boxP[7]);
    maxXYZ[2] += (boxP[2] + boxP[5]);
    /* SRINI Sept 4, 2002: Add non principal diagonal components */
  }

  /* For one component and more than 2 component systems, we'll print
     BCC/FCC/HCP/ICO/UNKNOWN.  For 2 component systems (e.g. Ni-Zr), 
     1/2/3 means FCC/HCP/UNKNOWN for Ni, and 4/5/6 means FCC/HCP/UNKNOWN
     for Zr */
  if( Ntype == 2 ) {
    nType     = 6;  /* 6 types: FCC,HCP,OTHERS for each of the 2 components */
    for(i = 0; i < nType; i++ ) {
      atommass[i]  = 1.;
      atomnum[i]   = i+1;
    }
  } 
  else {
    nType     = 5;  /* 5 types: BCC,FCC,HCP,ICO,OTHERS */
    for(i = 0; i < nType; i++ ) {
      atommass[i]  = 1.;
      atomnum[i]   = i+1;
    }
  }

  if( Ntype > 2 ) 
    Error("WriteConBaskes(): Too many atom types (%d)\n", Ntype );

  /* Count FCC atoms */
  totNfcc = 0;
  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < nAtoms; ai++ )   {
    if( AtomP[ai].cnTag == FCC ) totNfcc++;
  }
  
  /* Write output file in Baskes restart file format */
  if( !(ofp = fopen(ofnameP, "w") ) )
    Error("1. WriteConBaskes(): Unable to Open InFile: %s\n", ofnameP );

  /* May 20, 2003: Write File with No FCC atoms */
  if( !(ofp2 = fopen("NoFccAtoms.CN", "w") ) )
    Error("1b. WriteConBaskes(): Unable to Open InFile: NoFccAtoms.CN\n" );

  /* Write First SIX Lines */
#ifdef DBG
  /* ** First Line is Header (String);                                *
   * ** 2nd line: (NumAtoms, NumTypes, SnapTime, 90)                  *
   * ** 3rd line: (Max values of Box A/B/C edges)                     *
   * ** 4th line: (Min values of Box A/B/C edges)                     *
   * ** 5th/6th Lines: Read and IGNORE                                */
#endif
   /* Line-01: Comment line */
  sprintf(hdrStr, "CN outfile @ t= %.3lf. Types 1,2,3,4,5=BCC,FCC,HCP,ICO,UNKNOWN\0",
	  (double) SnapTime );
  fprintf(ofp, "%s\n", hdrStr );

  /* May 20, 2003: Write File with No FCC atoms */
  fprintf(ofp2, "%s\n", hdrStr );

  /* Line-02: (NumAtoms, NumTypes, time, box-angle) */
  fprintf(ofp, "%10d%10d%15.8lE%10.5lf\n", nAtoms, nType,
	  (double) SnapTime, (double) boxAngle  );
  
  /* May 20, 2003: Write File with No FCC atoms */
  fprintf(ofp2, "%10d%10d%15.8lE%10.5lf\n", nAtoms-totNfcc, nType,
	  (double) SnapTime, (double) boxAngle  );

   /* Line-03: (Max values of Box A/B/C edges) */
   fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", (double) maxXYZ[0],
           (double) maxXYZ[1], (double) maxXYZ[2] );

  /* May 20, 2003: Write File with No FCC atoms */
   fprintf(ofp2, "%25.16lE%25.16lE%25.16lE\n", (double) maxXYZ[0],
           (double) maxXYZ[1], (double) maxXYZ[2] );

   /* Line-04: (Min values of Box A/B/C edges) */
   fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", (double) minXYZ[0],
           (double) minXYZ[1], (double) minXYZ[2] );

  /* May 20, 2003: Write File with No FCC atoms */
   fprintf(ofp2, "%25.16lE%25.16lE%25.16lE\n", (double) minXYZ[0],
           (double) minXYZ[1], (double) minXYZ[2] );

   /* Line-05: (Atomic mass (amu*1.0365e-4) and atomic number */
  /* Loop through Num Atom Types and scan mass/type info */
  for( i = 0; i < nType; i++ )  {
     fprintf( ofp, "%25.16lE%10d\n", (double) atommass[i], atomnum[i]  );

     /* May 20, 2003: Write File with No FCC atoms */
     fprintf( ofp2, "%25.16lE%10d\n", (double) atommass[i], atomnum[i]  );
  }

  /* How much to shift the origin */
  if( MoveOriginToBoxCenterFlag == YES ) {
    dispXYZ[0] =  InitCooAtom01[0] - AtomP[0].coo.x;
    dispXYZ[1] =  InitCooAtom01[1] - AtomP[0].coo.y;
    dispXYZ[2] =  InitCooAtom01[2] - AtomP[0].coo.z;
  }

  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < nAtoms; ai++ )   {
     double x,  y,   z;
     double vx = 0., vy = 0., vz = 0.;

     x = (double)  AtomP[ai].coo.x + dispXYZ[0];
     y = (double)  AtomP[ai].coo.y + dispXYZ[1];
     z = (double)  AtomP[ai].coo.z + dispXYZ[2];

     /* A. Write X/Y/Z  */
     fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", x, y, z );
     
     /* B. Write VX/Y/Z  */
     fprintf(ofp, "%25.16lE%25.16lE%25.16lE\n", vx, vy, vz ); 
     
     /* May 20, 2003: Write File with No FCC atoms */
     if( AtomP[ai].cnTag != FCC ) {
       /* A. Write X/Y/Z;   VX/Y/Z */
       fprintf(ofp2, "%25.16lE%25.16lE%25.16lE\n", x, y, z );
       fprintf(ofp2, "%25.16lE%25.16lE%25.16lE\n", vx, vy, vz ); 

       /* C. Write Atom Type */
       switch( AtomP[ai].cnTag ) {
       case BCC:
	 typeTag = 1;
	 break;
       case FCC:
	 typeTag = 2;
	 break;
       case HCP:
	 typeTag = 3;
	 break;
       case ICO:
	 typeTag = 4;
	 break;
       default:
	 typeTag = 5;
       }

       fprintf(ofp2, "%10d\n", typeTag );
     }


     /* For one component and more than 2 component systems, we'll print
	BCC/FCC/HCP/ICO/UNKNOWN.  For 2 component systems (e.g. Ni-Zr), 
	1/2/3 means FCC/HCP/UNKNOWN for Ni, and 4/5/6 means FCC/HCP/UNKNOWN
	for Zr */
     if( Ntype == 2 ) {
       atomType = (int) AtomP[ai].atomicNum;
       switch( atomType ) {
       case 1:
	 /* C. Write Atom Type */
	 switch( AtomP[ai].cnTag ) {
	 case FCC:
	   typeTag = 1;
	   break;
	 case HCP:
	   typeTag = 2;
	   break;
	 default:
	   typeTag = 3;
	 }
	 break;
       case 2:
	 /* C. Write Atom Type */
	 switch( AtomP[ai].cnTag ) {
	 case FCC:
	   typeTag = 4;
	   break;
	 case HCP:
	   typeTag = 5;
	   break;
	 default:
	   typeTag = 6;
	 }
	 break;
       default: Error("WriteBaskes(): Type must be 1/2. You're reading true atomic number from non-dynamo input & printing here as atomtype?\n");
       }
     } 
     else {
       /* C. Write Atom Type */
       switch( AtomP[ai].cnTag ) {
       case BCC:
	 typeTag = 1;
	 break;
       case FCC:
	 typeTag = 2;
	 break;
       case HCP:
	 typeTag = 3;
	 break;
       case ICO:
	 typeTag = 4;
	 break;
       default:
	 typeTag = 5;
       }
     }
     
     fprintf(ofp, "%10d\n", typeTag );
  }
  
  fclose(ofp);

  /* Print nbor pair infor for HCP atoms */
  WriteHcpPairInfo(nAtoms, verbose );

  /* NOV 21, 2004: Print pair type info for each atom */
  WriteBondInfo(nAtoms, verbose );

} /* int WriteConBaskes(char*ofnameP, float*boxP, int4 verbose) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * Write KAI's Snap file (Text Format). It has following fields:            *
 *         ** each line: (AtomType, X, Y, Z)                                *
 *         ** Last 4 lines have 4.                                          *
 * NOTE: Types 0,1,2,3,4=BCC,FCC,HCP,ICO,UNKNOWN                            *
 ****************************************************************************/
void WriteConVASP(char *ofnameP, int4 nAtoms, int4 verbose )
{
  char        hdrStr[1024];
  FILE        *ofp              = NULL;
  int4        nType     = 5; /* 5 types of atoms: BCC,FCC,HCP,ICO,OTHERS */
  float       atommass[N_TYPE]  = {1.,1.,1.,1.,1.};  /* Dummy mass/number */
  int4        atomnum[N_TYPE]   = {0,1,2,3,4};
  int4        typeTag           = 0;

  int4        ai                = 0;


  /* Write output file in Baskes restart file format */
  if( !(ofp = fopen(ofnameP, "w") ) )
    Error("1. WriteConVASP(): Unable to Open InFile: %s\n", ofnameP );

  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < nAtoms; ai++ )   {
     double x,  y,   z;
     double vx = 0., vy = 0., vz = 0.;

     x = (double)  AtomP[ai].coo.x;
     y = (double)  AtomP[ai].coo.y;
     z = (double)  AtomP[ai].coo.z;

     /* C. Write Atom Type */
     switch( AtomP[ai].cnTag ) {
     case BCC:
       typeTag = 1;
       break;
     case FCC:
       typeTag = 2;
       break;
     case HCP:
       typeTag = 3;
       break;
     case ICO:
       typeTag = 4;
       break;
     default:
       typeTag = 5;
     }
     
     /* A. Write type/X/Y/Z  */
     fprintf(ofp, "%1d %+.4lE%+.4lE%+.4lE\n", typeTag, x, y, z );
  }

  /* Last 4 lines are 4/4/4/4 */
  fprintf(ofp, "4\n4\n4\n4\n");
  
  fclose(ofp);

} /* int WriteConBaskes(char*ofnameP, float*boxP, int4 verbose) */
/*---------------------------------------------------------------------------*/
/*|------------------------------------------------------------------------|
  |   check_for() compares two strings tmp1[]  and  '=' scanned from the   |
  |   file '*fp'  with the two given strings s1 and '='.                   |
  |------------------------------------------------------------------------|*/
static void check_for( FILE *fp, char *s1 ) 
{
  char tmp1[100], tmp2[3];

  fscanf(fp, "%s %s", tmp1, tmp2); 

  if (strcmp(s1, tmp1) != 0)  {
    fprintf(stderr,"Error. %s doesn't match %s in the input data file\n", 
            tmp1, s1);
    exit(1);
  }
  
  if (strcmp(tmp2, "=") != 0)  {
    fprintf(stderr,
            "Error. %s does not match the symbol = in the input data file\n",
            tmp2);
    exit(1);    
  }
}  /*  Close 'check_for()'  */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read Old Adhara Snap File format:                                    *
 *         ** Comments (start with '#')                                     *
 *         ** "Total_Number_of_Atoms           = nAtoms"                    *
 *         ** "Number_of_Blue_Red_Green_Atoms  =       3125 3165 3270"      *
 *         ** "Simulation_Box_Edge_Vectors     ="                           *
 *              Box Vector (A,, B, C) (3 rows)                              *
 *         ** "Box_Edge_Vector_Derivatives     ="                           *
 *              Derivatives (A,, B, C) (3 rows)                             *
 *         ** "Atom_Parameters ="                                           *
 *         ** (R/G/B atomPE X Y X id)                                       *
 ****************************************************************************/
int4 ReadConOldAdhara(char *ifnameP, float *boxP, float rcut, int4 verbose )
{
  float      boxVec[9];
  float      boxInv[9];

  float      h0[9], h0d[9];
  float      boxVol            =  1.0;
  float      tmpDDD[2];
  
  char       linE[1024], buff1[512], buff2[128];
  char       inpBuf[BUF_SIZE];
  FILE       *fpR              = NULL;    /* Turab Input file ptr */ 
  
  int4       ai                = 0;
  int4       i                 = 0;
  int4       j                 = 0;
  int4       numAtomsProcessed = 0;
  
  int4       totNumAtoms       = 0, nTypes = 0, nBlue, nRed, nGreen;
  
  int4       zero         = 0;
  int4       atomicNum    = 18;
  int4       nSlotsX      = 0,
             nSlotsY      = 0,
             nSlotsZ      = 0;
  
  for(i=0;i<9;i++ ) boxVec[i] = boxInv[i] = 0.0;

  /* Read Input file */
  if( !(fpR = fopen(ifnameP, "r") ) )
    Error("0. ReadConOldAdhara(): Unable to Open InFile: %s\n", ifnameP );

  /* Set file header for writing out CK/SN files */
  memset((void*) &FileHeader, zero, (size_t)sizeof(struct FHead) );

  FileHeader.Dimension      = 3;
  FileHeader.CoordFlag      = 1;        /* Single precision coord */
  FileHeader.TypeFlag       = 1;        /* Print CN tag           */

  AtomMass[0]     = 1.0;
  AtomicNumber[0] = 18; 

  /************** (I) Prepare to Read Headers *********************/
  /* a. scan in comments, i.e. first 3 lines beginning with # */
  while(fgets( linE, 1024, fpR ) != NULL)  {
    if(linE[0] == '#' || linE[0] == '\n') continue;
    else if(linE[0] == 'T') break;     /* T (in 'Total_Num....' is 1st char */
    else
      Error("1a. ReadConOldAdhara(): invalid data field\n");
  }

  /* Get Atom Info */
  sscanf( linE,"%s %s %u", buff1, buff2, &totNumAtoms );

  /* Initialize Total Atoms Fields */
  FileHeader.TotNumAtoms    = totNumAtoms;
  FileHeader.TotBluAtoms    = totNumAtoms;

  if(strcmp( buff1,"Total_Number_of_Atoms" ) != 0) 
    Error("1b. ReadConOldAdhara(): Invalid String\n");
  
  if(strcmp( buff2, "=" ) != 0) 
    Error("1c. ReadConOldAdhara(): : invalid field\n");
  
  check_for(fpR, "Number_of_Blue_Red_Green_Atoms");
  fscanf(fpR, "%u %u %u", &nBlue, &nRed, &nGreen);
  
  /* 3 numbers in Row-1 contain components of edge A;
     3 numbers in Row-2 contain components of edge B;
     3 numbers in Row-3 contain components of edge C; */
  check_for(fpR, "Simulation_Box_Edge_Vectors");
  for ( i = 0; i < 9; i++)  {
    fscanf(fpR, "%f", h0 + i);
    BoxVectors[i] = h0[i];
  }

  if( verbose )   {
    fprintf(stderr, 
            "$$$ ReadConOldAdhara(): Hvec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n",
            BoxVectors[0],BoxVectors[1],BoxVectors[2],
            BoxVectors[3],BoxVectors[4],BoxVectors[5],
            BoxVectors[6],BoxVectors[7],BoxVectors[8]);
    fflush( stderr );
  }

  check_for(fpR, "Box_Edge_Vector_Derivatives");
  for ( i = 0; i < 9; i++)  fscanf(fpR, "%f", h0d + i); 
  
  /* Print box edges
  fprintf( stderr, "\n***********\n");
  fprintf(stderr, "$$$ReadConOldAdhara() Ax/Ay/Az    =  ");
  for( i = 0; i < 9; i++) fprintf( stderr, "%.5f ", h0[i] );
  fprintf( stderr, "\n***********\n");
  fflush(stderr);  */

  check_for(fpR, "Atom_Parameters");
  fprintf(stderr, "\n$$$ReadConOldAdhara: READING ATOM PARAMETERS $$$\n");

  /* TRUE box matrix = Transpose of the matrix read in */
  TRANSPOSE_MATRIX( boxP, BoxVectors );
  ASSIGN_MATRICES( boxVec, boxP )

  /* Compute Box matrix Inverse */
  INVERT_MATRIX( boxInv, boxVec, boxVol );

  /* Compute TRUE number of slots */
  NumSlotsXYZ(boxVec, rcut, &nSlotsX, &nSlotsY, &nSlotsZ, verbose );

  /* Initialize true number of slots */
  LCInitializeSlots(nSlotsX, nSlotsY, nSlotsZ, totNumAtoms, verbose );

  /* Loop and Read (R/G/B atomPE X Y X id) for each atom in the system 
     IGNORE all except X/Y/Z */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    Vector_t xyz;                       /* coord of grid points   */
    float    pe = 0.0;
    uint4    atomType;
    uint4    itmp, id;                  /* Atom ID                      */
    uint4    arrIn;                     /* Index of atom's position 
                                           in AtomP[]                   */
    char     tags[4] = {'M', 'R'};      /* [0] = movetag, [1] = color;  */
    numAtomsProcessed++;
    
    id       = numAtomsProcessed;
    atomType = 0;

    if( numAtomsProcessed%5000 == 0 )
      fprintf( stderr, "$$$ReadConOldAdhara(): Read %u atoms\n", id );


    fscanf( fpR, "%s %f %f %f %f %d", 
            buff1, &pe, &xyz.x, &xyz.y, &xyz.z, &itmp );

    /* Insert into link-cell slots: RETURN value is index of atom's
       position in AtomP[] (ignored for now) */
    arrIn = LCInsertAtomIntoSlot(boxVec, boxInv, &xyz, pe, id, atomicNum,
                                 tags, verbose );

  }   /* Loop Over Atoms.               */
  
  fclose(fpR);

  /* if( verbose ) DbgPrintReadInAtoms( totNumAtoms ); */
  fprintf(stderr, "$$$ ReadConOldAdhara():TotAtoms=(%d)\n", totNumAtoms );

  return(totNumAtoms);

} /* void ReadConOldAdhara(char *ifnameP, float *, float rcut, int4 verbose) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * (A) Read VASP POSCAR file. These files have following fields:            *
 *         ** Line-01: Comment Line                                         *
 *         ** Line-02: Alat constant                                        *
 *         ** Line-03 to Line-05: Box A/B/C vectors                         *
 *         ** Line-06: Number of atoms (can handle only one type)           *
 *         ** Line-07: Keyword 'Direct" must be present.                    *
 *         ** Line-08 -- atom coordinates X/Y/Z infractional coords         *
 ****************************************************************************/
int4 ReadConVASP(char *ifnameP, float *boxP, float rcut, int4 verbose )
{
  char       coordSys[12];
  float      alat              = 1.0;
  float      boxVec[9];
  float      boxInv[9];
  float      boxVol            =  1.0;
  float      x[9];
  const float deltaXYZ = 0.000001;


  float      snaptime          = -1.0;

  char        buFF[1024];
  char        inpBuf[BUF_SIZE];
  FILE        *ifp              = NULL;

  int4        ai                = 0;
  int4        i                 = 0;
  int4        j                 = 0;
  int4        numAtomsProcessed = 0;

  int4        totNumAtoms  = 0;
  int4        numBluAtoms  = 0;
  int4        numRedAtoms  = 0;
  int4        numGrnAtoms  = 0;
  int4        zero         = 0;

  int4        nSlotsX      = 0,
              nSlotsY      = 0, 
              nSlotsZ      = 0;

  /* Set file header for writing out CK/SN files */
  memset((void*) &FileHeader, zero, (size_t)sizeof(struct FHead) );

  FileHeader.Dimension      = 3;
  FileHeader.CoordFlag      = 1;        /* Single precision coord */
  FileHeader.TypeFlag       = 1;        /* Print CN tag           */

  /* Read Input file (VASP POSCAR) */
  if( !(ifp = fopen(ifnameP, "r") ) )
    Error("0. ReadConVASP(): Unable to Open InFile: %s\n", ifnameP );


  /* Line-01: is comment line; read and ignore */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("1. ReadConVASP(): Unable to get parameters in %s\n", ifnameP );

  /* Line-02: Scan in alat */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("2a. ReadConVASP(): Unable to get alat in %s\n", ifnameP );
  
    /* Scan into parameters array */
  if( sscanf(buFF, "%f", &alat) != 1 )
    Error("2b. ReadConVASP(): Unable to get alat in %s\n", ifnameP );

  /* Line-03: Scan in Box A-vector */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("3a. ReadConVASP(): Unable to get Box Avect in %s\n", ifnameP );
  
    /* Scan into parameters array */
  if( sscanf(buFF, "%f %f %f", &x[0], &x[1], &x[2]) != 3 )
    Error("3b. ReadConVASP(): Unable to get Box Avect in %s\n", ifnameP );

  /* Line-04: Scan in Box B-vector */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("4a. ReadConVASP(): Unable to get Box Bvect in %s\n", ifnameP );
  
    /* Scan into parameters array */
  if( sscanf(buFF, "%f %f %f", &x[3], &x[4], &x[5]) != 3 )
    Error("4b. ReadConVASP(): Unable to get Box Bvect in %s\n", ifnameP );

  /* Line-05: Scan in Box C-vector */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("5a. ReadConVASP(): Unable to get Box Cvect in %s\n", ifnameP );
  
    /* Scan into parameters array */
  if( sscanf(buFF, "%f %f %f", &x[6], &x[7], &x[8]) != 3 )
    Error("5b. ReadConVASP(): Unable to get Box Cvect in %s\n", ifnameP );

  /* Scale box by alat */
  for( i = 0; i < 9; i++ )  {
    BoxVectors[i] = (float) (alat*x[i]);
    x[i]          = BoxVectors[i];
  }

  /* Line-06: Scan in Natoms */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("6a. ReadConVASP(): Unable to get Natoms in %s\n", ifnameP );
  
    /* Scan into parameters array */
  if( sscanf(buFF, "%d", &totNumAtoms) != 1)
    Error("6b. ReadConVASP(): Unable to get Natoms in %s\n", ifnameP );

  FileHeader.TotNumAtoms    = totNumAtoms;
  FileHeader.TotBluAtoms    = totNumAtoms;

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConVASP(): TotAtoms(%u)\n", totNumAtoms );
    fflush( stderr );
  }

  if( FileHeader.StressFlag )
    Error("$$$ ReadConVASP(): CNA code cannot handle stress in SN/CK file\n");

#ifdef DBGGGG
  ASSIGN_MATRICES( boxVec, BoxVectors )
  ASSIGN_MATRICES( boxP,   boxVec     )
#endif

  /* TRUE box matrix = Transpose of the matrix read in */
  TRANSPOSE_MATRIX( boxP, BoxVectors );
  ASSIGN_MATRICES( boxVec, boxP )



  /* Compute Box matrix Inverse */
  INVERT_MATRIX( boxInv, boxVec, boxVol );

  /* Compute TRUE number of slots */
  NumSlotsXYZ(boxVec, rcut, &nSlotsX, &nSlotsY, &nSlotsZ, verbose );

  /* Initialize true number of slots */
  LCInitializeSlots(nSlotsX, nSlotsY, nSlotsZ, totNumAtoms, verbose );

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ ReadConVASP(): Hvec\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n\t\t(%f %f %f)\n", x[0],x[1],x[2],x[3],x[4],x[5], x[6],x[7],x[8]);
    fflush( stderr );
  }

  /* Line-07: Scan in Coord System Flag */
  if( fgets(buFF, 1024, ifp) == NULL ) 
    Error("7a. ReadConVASP(): Unable to get Natoms in %s\n", ifnameP );

  /* Copy 1st character: can be 'D' or 'C' */
  strncpy(coordSys, buFF, 1);

  switch(buFF[0]) {
  case 'd':
  case 'D':
    ;        /* Direct coordinate. CORRECT => continue */
    break;
  default:   /* ERROR */
    Error("7b. ReadConVASP(): Unable to get Invalid coord system(%s) in %s\n",
	  buFF, ifnameP );
  }
							  
  /* Loop, and Read Parameters */
  for( ai = 0; ai < totNumAtoms; ai++ )   {
    Vector_t scoo;
    float    pe;
    Vector_t xyz;                       /* FLOAT coord of grid points   */
    Vector_t vel;

    uint4    atomicNum;
    uint4    id;                        /* Atom ID                      */
    uint4    arrIn;                     /* Index of atom's position 
					   in AtomP[]                   */
    char     tags[4] = {'M', 'B'};      /* [0] = movetag, [1] = color;  */
    
    numAtomsProcessed++;

    if( numAtomsProcessed%500000 == 0 )
      fprintf( stderr, "$$$ReadConVASP(): Read %u atoms\n",
	      numAtomsProcessed );

    /* A1. Read X/Y/Z in DOUBLE PRECISION */
    switch( FileHeader.CoordFlag )   {
    case 1:
      /* Line-08: Scan in atom's X/Y/Z FRACTIONAL coord */
      if( fgets(buFF, 1024, ifp) == NULL ) 
	Error("8a. ReadConVASP(): Unable to get atom(%d) coord in %s\n",
	      ai+1, ifnameP );
  
      /* Scan into parameters array */
      if( sscanf(buFF, "%f %f %f", &scoo.x, &scoo.y, &scoo.z) != 3 )
	Error("8b. ReadConVASP(): Unable to get atom(%d) coord in %s\n", 
	      ai+1, ifnameP );
      break;
    default: Error("8b. ReadConVASP(): Invalid Precision for Coord\n");
    }
    
    /* ROUND-OFF ERROR FIX:  Add a small deltaXYZ */
    scoo.x += deltaXYZ;
    scoo.y += deltaXYZ;
    scoo.z += deltaXYZ;

    /* Convert to cartesian coord if necessary */
    switch(coordSys[0]) {
    case 'd':
    case 'D':
      /* Direct coordinate. convert to cartesian */
      /* Compute S-coord from cartesian coords */
      MATRIX_3D_VECTOR_MULTIPLY( xyz, boxVec, (scoo) );  
      break;
    default:   /* ERROR */
      Error("8c. ReadConVASP(): Invalid coord system (%c)\n",
	    coordSys[0]);
    }
			
    /* color is BLUE */
    pe        = 0.0;
    atomicNum = 18;
    tags[0]   = 'M';
    tags[1]   = 'B';
    id        = numAtomsProcessed;
    
    /* Insert into link-cell slots: RETURN value is index of atom's
       position in AtomP[] (ignored for now) */
    arrIn = LCInsertAtomIntoSlot(boxVec, boxInv, &xyz, pe, id, atomicNum,
				 tags, verbose );

  }   /* Loop Over Atoms.               */
  
  fclose(ifp);


  /* if( verbose ) DbgPrintReadInAtoms( totNumAtoms ); */

  printf("$$$ ReadConVASP():TotAtoms/Nblu/red/green=(%d %d %d %d)\n",
	 totNumAtoms, numBluAtoms, numRedAtoms, numGrnAtoms );

  return( numAtomsProcessed );
} /* int ReadConVASP(char *ifnameP, float *boxP, float rcut, int4 verbose ) */
/*---------------------------------------------------------------------------*/
/****************************************************************************
 * Write Pair Info for HCP Atoms Only.                                      *
 ****************************************************************************/
void WriteHcpPairInfo(uint4 nAtoms, int4 verbose )
{
  Atom_t      *ap               = NULL;
  FILE        *ofp              = NULL;
  int4        ai                = 0;
  int4        ii                = 0;
  int4        numAtomsProcessed = 0;

  if( !(ofp = fopen("HcpNbors.LOG", "a") ) )
    Error("WriteHcpPairInfo(): Unable to open HcpNbors.LOG file \n");

  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < nAtoms; ai++ )   {
    uint4 nn = 0;

    /* Atom pointer */
    ap  = AtomP + ai;

    /* Check if HCP Atom Type */
    if( (ap->cnTag != HCP) )  continue;

    /* Else, it's a HCP atom */

    /* Number of neighbors */
    nn = (uint4) ap->pairs;

    /* Print I-J and pair info */
    fprintf(ofp, "*** HCP Atom: %6d:\n", ai+1 );

    /* Loop over neighbors */
    for( ii = 0; ii < nn; ii++ ) {
      int4 nborID = 0;

      nborID = AtomP[ap->nbor[ii]].id;

      /* Print I-J and pair info */
      fprintf(ofp, "\tNbor: %6d\tPairType: %1d%1d%1d\n", nborID,
	      (int4) ap->index[ii][0], (int4) ap->index[ii][1],
	      (int4) ap->index[ii][2] );


      /* fprintf(ofp, "\t%6d - %6d\t%1d%1d%1d\n", ai+1, nborID,
	 (int4) ap->index[ii][0], (int4) ap->index[ii][1],
	 (int4) ap->index[ii][2] ); */
    } /* for( ii = 0; ii < nn; ii++ ) */

    fprintf(ofp, "************************\n\n");

    numAtomsProcessed++;

  }
  
  fclose(ofp);

} /* int WriteHcpPairInfo(uint4 nAtoms, int4 verbose ) */
/*---------------------------------------------------------------------------*/
/* SRINI: NOV 21, 2004: Print pairs info for EACH Atom in CNbond.log */
void WriteBondInfo(uint4 nAtoms, int4 verbose )
{
  Atom_t      *ap               = NULL;
  FILE        *ofp              = NULL;
  int4        ai                = 0;
  int4        ii                = 0;
  int4        numAtomsProcessed = 0;

  if( !(ofp = fopen("CNbond.log", "w") ) )
    Error("WriteHcpPairInfo(): Unable to open CNbond.log file \n");

  fprintf(ofp, "\n************ Bond Indices Of EACH Atom ************\n");


  /* Loop and print atom coordinates first: AtomInfo[nAtoms].x */
  for( ai = 0; ai < nAtoms; ai++ )   {
    uint4 nn = 0;

    /* Atom pointer */
    ap  = AtomP + ai;

    /* Number of neighbors */
    nn = (uint4) ap->pairs;

    /* Print I-J and pair info */
    fprintf(ofp, "Atom: %6d (CNtag=%d)\n", ai+1, ap->cnTag );

    /* Loop over neighbors */
    for( ii = 0; ii < nn; ii++ ) {
      int4 nborID = 0;

      nborID = AtomP[ap->nbor[ii]].id;

      /* Print I-J and pair info */
      fprintf(ofp, "\tNbor: %6d\t(%1d%1d%1d) Pair\n", nborID,
	      (int4) ap->index[ii][0], (int4) ap->index[ii][1],
	      (int4) ap->index[ii][2] );
    } /* for( ii = 0; ii < nn; ii++ ) */

    fprintf(ofp, "***\n\n");

    numAtomsProcessed++;

  }
  
  fclose(ofp);

} /* int WriteBondInfo(uint4 nAtoms, int4 verbose ) */
/*---------------------------------------------------------------------------*/
