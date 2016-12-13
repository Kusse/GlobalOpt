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
 * $Date: 2000/01/28 19:58:36 $
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNpairStat.c,v $
 * $Revision: 1.4 $
 * $Log: CNpairStat.c,v $
 * Revision 1.4  2000/01/28 19:58:36  sgsrini
 * Working version with uint2 variables in Atom_t and dynamic mem alloc
 * DON'T DYNAMICALLY allocate MULTI-DIMENSIONAL arrays for uint2/int2/char*
 * Pointersizes are at least 4 bytes (larger than these variable sizes)
 *
 * Revision 1.3  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
 *
 * Revision 1.2  2000/01/27 21:45:31  sgsrini
 * *** empty log message ***
 *
 * Revision 1.1  2000/01/27 21:39:21  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *                  $$$$ CNpairStat.c: CNN pair statistics.                  *
 *****************************************************************************
 * CNpairStat.c: We need to plot the spatial distribution of the different   *
 * pair types (000, 100, 200, 211, 300, 311, 322, 400, 411, 421, 422, 433,   *
 * 444, 532, 533, 544, 555).  To do this bin each atom into a 2D bin based on*
 * their (X, Y) coordinates.  Compute number of pair types in each bin and   *
 * plot this info as a histogram/contour plots for different pair-types. We  *
 * need a bin of size 5500 X 5500 (for a 80000 atom system). This takes up   *
 * lots of memory. so we've only one BIN-ARRAY and compute the distribution  *
 * for each pair-type sequentially. Use a BIN-SIZE of 1.1 (change if needed) *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*---------------------------------------------------------------------------*/
#include "CNmain.h"
#include "CNmacros.h"
#include "CNpairStat.h"
/*---------------------------------------------------------------------------*/
/* Extern Global Functions */
extern void    Error(char*, ... );
/*---------------------------------------------------------------------------*/
/* Extern Global Variables */
extern Atom_t  *AtomP;
extern int     Ntype;                               /* Defined in CNio.c   */
extern int4    InpFileFormatFlag;                   /* Defined in CNmain.c */
/*---------------------------------------------------------------------------*/
/* Initialize CNstat variables */
void InitializeCNvariables( CNstat_t *cnStatTblP )
{
  void  *pOut  = NULL;
  int4  zero   = 0;
  
  /* Zero structure using memset() instead of bzero() */
  pOut = (void*) memset( (void*) cnStatTblP, zero, sizeof(CNstat_t) );
    
  /* Check return value */
  if( pOut != ((void*) cnStatTblP) )  {
    fprintf( stderr, "$$$ InitializeCNvariables(): pOut != cnStatTblP\n" );
    fflush( stderr );
    Error("$$$ InitializeCNvariables(): memset() problem........\n");
  }

} /* void InitializeCNvariables( CNstat_t *cnStatTblP ) */
/*---------------------------------------------------------------------------*/
/*****************************************************************************
 * Compute total number of BCC, FCC, HCP and ICOSAHEDRAL ATOMS. Also, compute*
 * average number of neighbors, numatoms(and %atoms) with each coordination  *
 *****************************************************************************
 *    ASSUME: ONE COMPONENT SYSTEM.                                          *
 *****************************************************************************/
void CountBccFccHcpIcoAtoms( CNstat_t *cnStatTblP, uint4 nAtoms, int4 verbose )
{
  register Atom_t   *ap;                /* Atom_t pointer                    */
  int4              tyI = 0;
  uint4             n421;               /* # 421 pairs for each atom         */
  uint4             n422;               /* # 422 pairs for each atom         */
  uint4             n444;               /* # 444 pairs for each atom         */
  uint4             n555;               /* # 551 pairs for each atom         */
  uint4             n666;               /* # 666 pairs for each atom         */
  uint4             i;
  uint4             ia;
  uint4             isumdist, isum;

  uint2             J, K, L;            /* 3 digit pairs-type indx(421 etc)  */

  /* Initialize structure counting fields */
  for( i = 0; i < N_TYPE; i++ ) {
    cnStatTblP->Nbcc[i]     = 0;
    cnStatTblP->Nfcc[i]     = 0;
    cnStatTblP->Nhcp[i]     = 0;
    cnStatTblP->Nico[i]     = 0;
    cnStatTblP->Nicodist[i] = 0;
  }

  /* Initialize counting array */
  for( J = 0; J < NNBR; J++ )
    for( K = 0; K < NNBR; K++ )
      for( L = 0; L < NNBR; L++ ) cnStatTblP->NumDiffPairs[J][K][L] = 0;

  /* Number of Atoms with Different Coordination number */
  for( ia = 0; ia < nAtoms; ia++ )   {
    uint4  nn = 0;

    /* Get atom pointer */
    ap        = AtomP + ia;
    nn        = (uint4) ap->pairs;
    ap->cnTag = (uint2) ap->pairs;

    /* FCC, HCP, ICO have 12 1st nearest neighbors, 
       BCC has 14 neighbors (1st + 2nd nearest neighbors).
       For atoms that are don't belong to these regular polyhedra, we use
       their coordination number for the tag.  */
    if( nn < 12 )      {
      for( i = 0; i < nn; i++ )   {
	J = (uint2) ap->index[i][0];
	K = (uint2) ap->index[i][1];
	L = (uint2) ap->index[i][2];
 	cnStatTblP->NumDiffPairs[J][K][L]++;   /* count this pair type */
      }
      continue;
    }

    /* if( verbose )  printf("CountBccFccHcpIcoAtoms() = %d neighbors\n", 
       (int4) nn ); */

    /* reinitialize some counters */
    n421 = n422 = n444 = n555 = n666 = 0;

    /* get the 3 digit designation (index) of each pair-type */
    for( i = 0; i < nn; i++ )   { 
      J = (uint2) ap->index[i][0];
      K = (uint2) ap->index[i][1];
      L = (uint2) ap->index[i][2];
      cnStatTblP->NumDiffPairs[J][K][L]++;   /* count this pair type */

      if( (J == 4) && (K == 2) && (L == 1) ) n421++;
      if( (J == 4) && (K == 2) && (L == 2) ) n422++;
      if( (J == 4) && (K == 4) && (L == 4) ) n444++;
      if( (J == 5) && (K == 5) && (L == 5) ) n555++;
      if( (J == 6) && (K == 6) && (L == 6) ) n666++;
    }  /* for(i) */

    /* Get atom type from 0-N-1 */
    if( Ntype == 2 ) {
      if( InpFileFormatFlag != MBASKES_FORMAT )
	Error("CountBccFccHcpIcoAtoms(): ap->atomicNum is atomtype ONLY for baskes dynamo format\n");
      
      /* Baskes's atomtype goes from 1-N => subtract one */
      tyI =  (int4) ap->atomicNum - 1;
    }
    else {
      tyI = 0;
    }

    /* FCC atom */
    if( n421 == 12 && nn == 12)   {
      cnStatTblP->Nfcc[tyI]++;      /* Increment # FCC atoms by 1            */
      ap->cnTag = FCC;         /* Environment of this atom is FCC like shell */
      continue;                /* Skip to next pair                          */
    }

    /* HCP atoms */
    if(n421 == 6 && n422 == 6 && nn == 12) {
      cnStatTblP->Nhcp[tyI]++;
      ap->cnTag = HCP;         /* Environment of this atom is HCP like shell */
      continue;
    }

    /* DIC atoms */
    /* Srini: Nov 5, 2004: Commented and Changed */
    /* if(n555 >= 8 && nn >= 12) { */
    if(n555 >= 8 && (nn != 12)) {
      cnStatTblP->Nicodist[tyI]++;
      ap->cnTag = DIC;         /* Environment is distorted ICOSAHEDRA        */
      
      /* mark the ICOSdistatoms field of this atom and its neighbors */
      ap->tags[2] = 'D';                     /* Part of Distorted Icosahedra */

#ifdef SRINI_COMMENT_NOV_2004
      for( i = 0; i < nn; i++ )
	AtomP[ap->nbor[i]].tags[2] = 'D';    /* Part of Distorted Icosahedra */
#endif

      continue;                  /* skip to next pair */
    }

    /* ICOSAHEDRAL atoms */
    if( n555 == 12 && nn == 12 )   {
      cnStatTblP->Nico[tyI]++;
      ap->cnTag = ICO;         /* Environment is ICOSAHEDRA                  */

      /* mark the ICOSatoms field of this atom and its neighbors */
      ap->tags[2] = 'I';                     /* Part of FULL Icosahedra      */

#ifdef SRINI_COMMENT_NOV_2004
      for( i = 0; i < nn; i++ )
	AtomP[ap->nbor[i]].tags[2] = 'I';    /* Part of FULL Icosahedra      */
#endif

      continue;                  /* skip to next pair */
    }

    /* BCC atoms */
    if( n444 == 6 && n666 == 8 && nn == 14 )    {
      cnStatTblP->Nbcc[tyI]++;
      ap->cnTag = BCC;        /* env is BCC */
    }
  }  /* for( ia = 0; ia < nAtoms; ia++ )  */

  /* Count howmany atoms are in ICOsahedra */
  isumdist = isum = 0;

  /* Number of Atoms that are part of an ICOSAHEDRA/distorted ICOSAHEDRA */
  for( ia = 0; ia < nAtoms; ia++ )   {

    /* Get atom pointer */
    ap        = AtomP + ia;

    /* Srini: Nov 5, 2004: Changed */
    switch( ap->cnTag )   {
    case ICO :  isum++;     break;
    case DIC :  isumdist++; break;
    default  :  break;         /* Do nothing */
    }
#ifdef SRINI_COMMENT_NOV_2004
    switch( ap->tags[2] )   {
    case '\0':              break;         /* Do nothing */
    case 'I' :  isum++;     break;
    case 'D' :  isumdist++; break;
    default  :  Error("$$$ CountBccFccHcpIcoAtoms(): Invalid tags[2]\n");
    }
#endif
  }  /* for( ia = 0; ia < nAtoms; ia++ )  */
  
  /* # atoms in complete and distorted ICOSAHEDRA */
  cnStatTblP->NatmInFullICO = isum;
  cnStatTblP->NatmInDistICO = isumdist;

#ifdef SRINI_COMMENT_NOV_2004
  /* We double-count Nico in Nicodist => make adjustment */
  for( i = 0; i < N_TYPE; i++ )
    cnStatTblP->Nicodist[i] -= cnStatTblP->Nico[i];
#endif
} /* void CountBccFccHcpIcoAtoms( cnStatTblP ) */
/*---------------------------------------------------------------------------*/
/* ASSUME: ONE COMPONENT SYSTEM.  Find total number of neighbors, 
   average num of neighbors, numatoms/%-atoms with each coordination */
void ComputeCoordination(CNstat_t *cnStatTblP, uint4 nAtoms, int4 verbose )
{
  uint4            ia             = 0;
  int4             maxCoordNum    = 0;         /* maximum value of coord num */
  uint4            neighbor_total = 0;

  /* Initialize the counting arrays */
  cnStatTblP->maxCoordNum  = 0;
  cnStatTblP->avneb        = 0.0;

  for( ia = 0; ia < NNBR; ia++ )  cnStatTblP->NcoordnumbM[ia] = 0;

  /* Number of Atoms with Different Coordination number */
  for( ia = 0; ia < nAtoms; ia++ )   {
    register Atom_t *ap = AtomP + ia;
    uint4           nnn = ap->nNbors;

    /* If atom is Fixed, ignore it
    if( ap->tags[0] == 'F' ) continue; */

    if( nnn >= NNBR ) Error("ComputeCoordination(): Too many neighbors (%d)\n",
			    nnn );

    if(ap->cnTag == FCC || ap->cnTag == BCC || ap->cnTag == HCP ||
       ap->cnTag == ICO || ap->cnTag == DIC ) continue;

    cnStatTblP->NcoordnumbM[nnn]++;

    /* compute the maximum value of the coordination number */
    maxCoordNum = MAX( nnn, maxCoordNum );
  }  /* for( ia = 0; ia < nAtoms; ia++ ) */

  /* Find average number of nearest neighbors (irrespective of type) */
  cnStatTblP->maxCoordNum =  maxCoordNum;
  neighbor_total          = 0;
  for( ia = 0; ia <= maxCoordNum; ia++ )   {  /* it is <= */
    neighbor_total += ( (cnStatTblP->NcoordnumbM[ia]) * (uint4) ia );

    /* Compute the percentage of atoms with each coo number */
    cnStatTblP->PercoordnAA[ia] = ( ((float) cnStatTblP->NcoordnumbM[ia])/
				    ((float) nAtoms) );
    cnStatTblP->PercoordnAA[ia] *= ((float)100.0);
  } /* for(ia) */


  /* NOV 12, 2004: This is WRONG => Comment out the print out in CNio.c */
  cnStatTblP->avneb = ((float) neighbor_total)/((float) nAtoms);
} /* static void ComputeCoordination( cnStatTblP ) */
/*---------------------------------------------------------------------------*/










/*---------------------------------------------------------------------------*/
#ifdef NOT_YET_REWRITTEN                                   /* Jan 26, 2000   */
/*---------------------------------------------------------------------------*/
static float  BinSize  = BIN_SIZE;
static int4   NumBinsX = -1;                               /* #bins along X */
static int4   NumBinsY = -1;                               /* #bins along Y */
static int4   TotNumCNPairsOfEachType[N_CN_PTYPES];        /* Total Num pairs 
							      of each type  */
static int4   NumCNPairTypeInBin[NUM_X_BINS][NUM_Y_BINS];  /* Num pairs of each
							       type in bins  */
static int4   NumAtomsInBin[NUM_X_BINS][NUM_Y_BINS];       /* NumAtom in bin */
static cnpType CNpairID[N_CN_PTYPES];
/*---------------------------------------------------------------------------*/
/* Initialize the CN-pair type IDs that we want to count.  The following 17
   CN pair types exist in energy minimized TJ configurations. In other systems 
   we will have other CN-pair types => CHANGE HERE FOR OTHER SYSTEMS */
static void InitCNpairIDs()
{
  int i;

  /* Initialize num pairs of each type to zero */
  for( i = 0; i < N_CN_PTYPES; i++ ) TotNumCNPairsOfEachType[i] = 0;

  CNpairID[0].l = 0; CNpairID[0].m = 0; CNpairID[0].n = 0;  /* 000 Pair */
  CNpairID[1].l = 1; CNpairID[1].m = 0; CNpairID[1].n = 0;  /* 100 Pair */
  CNpairID[2].l = 2; CNpairID[2].m = 0; CNpairID[2].n = 0;  /* 200 Pair */
  CNpairID[3].l = 2; CNpairID[3].m = 1; CNpairID[3].n = 1;  /* 211 Pair */
  CNpairID[4].l = 3; CNpairID[4].m = 0; CNpairID[4].n = 0;  /* 300 Pair */
  CNpairID[5].l = 3; CNpairID[5].m = 1; CNpairID[5].n = 1;  /* 311 Pair */
  CNpairID[6].l = 3; CNpairID[6].m = 2; CNpairID[6].n = 2;  /* 322 Pair */
  CNpairID[7].l = 4; CNpairID[7].m = 0; CNpairID[7].n = 0;  /* 400 Pair */
  CNpairID[8].l = 4; CNpairID[8].m = 1; CNpairID[8].n = 1;  /* 411 Pair */
  CNpairID[9].l = 4; CNpairID[9].m = 2; CNpairID[9].n = 1;  /* 421 Pair */

  CNpairID[10].l = 4; CNpairID[10].m = 2; CNpairID[10].n = 2;  /* 422 Pair */
  CNpairID[11].l = 4; CNpairID[11].m = 3; CNpairID[11].n = 3;  /* 433 Pair */
  CNpairID[12].l = 4; CNpairID[12].m = 4; CNpairID[12].n = 4;  /* 444 Pair */
  CNpairID[13].l = 5; CNpairID[13].m = 3; CNpairID[13].n = 2;  /* 532 Pair */
  CNpairID[14].l = 5; CNpairID[14].m = 3; CNpairID[14].n = 3;  /* 533 Pair */
  CNpairID[15].l = 5; CNpairID[15].m = 4; CNpairID[15].n = 4;  /* 544 Pair */
  CNpairID[16].l = 5; CNpairID[16].m = 5; CNpairID[16].n = 5;  /* 555 Pair */
}  /* static void InitCNpairIDs() */
/*---------------------------------------------------------------------------*/
/* Initialize the 2D array (bins) used to count CN-pair type IDs.
   COLUMNS of hpresent are box vectors */
static void InitBinsUsedToCountCNpairTypes(float *hpresent)
{
  int i, j;

  /* Get the farthest vertex of the box's X-Y face */
  NumBinsX = (hpresent[0]+hpresent[1]+hpresent[2])/BinSize;
  NumBinsY = (hpresent[3]+hpresent[4]+hpresent[5])/BinSize;

  for( i = 0; i < NUM_X_BINS; i++ )
    for( j = 0; j < NUM_Y_BINS; j++ )   {
      NumCNPairTypeInBin[i][j] = 0;
      NumAtomsInBin[i][j]      = 0;        /* June 8, 1998: Added number of 
					      atoms in the bin to normalize 
					      raw CN-pair values */
    }
} /* static void InitBinsUsedToCountCNpairTypes() */
/*---------------------------------------------------------------------------*/
/* July 14, 1998:
   Compute spatial distributioon based purely on atom coordination
   NOTE: COLUMNS OF hpresent are Box-Vectors */
void ComputeOnlyCNSpatialDistr( char *outFile, float *hpresent )
{
  char     bufferArray[384];
  char     fileName[384];

  FILE     *fp;

  nodeType *np, *tnp;          /* temp pointer to head of list     */
  Atom_t   *ap;                /* Atom_t pointer                   */
  int      nbx, nby;           /* Bin coordinates X/Y              */
  int      i, j, k;

  int      len;

  /* Get head of the serial list */
  tnp = (nodeType *) LCGetSerialList();

  /* (A) Initialize counting array */
  InitBinsUsedToCountCNpairTypes( hpresent );

  if( NumBinsX <= 0 || NumBinsY <= 0 ) 
    Error("(00) ComputeOnlyCNSpatialDistr(): Invalid bins");

  /* (B) Loop over atoms  */
  for( np = tnp; np != NULL; np = np->nextInSerialList )   {
    int  nn;

    /* Get atom pointer */
    ap = (Atom_t  *)  &(np->dx);
    
    /* Nov. 28, 1997: count only if it is MOBILE. ADDED movement_tag today */
    if( ap->movement_tag == FIX )  continue;

    /* Get this Atom's bin indices */
    nbx = (int)  (ap->coord.x/BinSize);
    nby = (int)  (ap->coord.y/BinSize);
      
    if( (nbx >= NUM_X_BINS) || (nby >= NUM_Y_BINS) )
      Error("(1a) ComputeOnlyCNSpatialDistr(): Too many bins");
      
    if( (nbx < 0) || (nby < 0) )
      Error("(1b) ComputeOnlyCNSpatialDistr(): invalid num bins");
      
    /* Number of nearest neighbors */
    nn = ap->numNeigbors;
    if( nn > 19 ) 
      Error("(2) ComputeOnlyCNSpatialDistr(): Too many neighbors");;

    /* Increment coord num info in bins */
    NumCNPairTypeInBin[nbx][nby] += nn;
    NumAtomsInBin[nbx][nby]++;
  }  /* for( np = tnp; np != NULL; np = np->nextInSerialList ) */


  /* (C) Open file to write CN spatial distribution data. Remove last 3 
         characters "con"; Copy name (without "con" into tmp[] */
  strcpy( fileName, outFile );

  len = strlen( fileName );
      
  /* Print file name with CNPAIR as suffix into fileName[] */
  sprintf( fileName+len-3, "SpaDisCN\0" );
      
  if( !(fp = fopen(fileName, "w") ) )
    Error("(3) ComputeOnlyCNSpatialDistr(): Unable to open file");

  /* (E) Print File Header */
  fprintf( fp, "#Spatial Distribution of Various Coordination Numbers\n");
  
  fprintf( fp, "#Xcoo\t\tYcoo\t\tAtomsInBin\tRaw/Normalized CoordNum Value\n");
  
  /* (F) Loop through bins and print out the CNpairtype COUNT (RAW DATA) in 
     each bin to the output file */
  for( k = 0; k < NumBinsY; k++ ) {
    for( j = 0; j < NumBinsX; j++ ) {
      
      /* Nov 2, 1997: Print the number even if it's ZERO */
      
      /* All files will have X/Y coordinates also */
      if( NumAtomsInBin[j][k] > 0 ) {
	float nCN;
	nCN = ( ((float)NumCNPairTypeInBin[j][k])/
		     ((float)NumAtomsInBin[j][k]) );

	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t%9.4lf\n",
		 ((float) j) * BinSize, ((float) k) * BinSize, 
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k], nCN );
      }
      else {
	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t0000.0000\n",
		 ((float) j) * BinSize, ((float) k) * BinSize,
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k] );
      }  /* else */
    } /* for(j) */

    fprintf(fp, "\n");     /* To make contour plots (using Gnuplot) we
			      need to leave a blank line between 2 successive
			      column (or row) of data */
  }   /* for(k) */
  /* June 8, 1998 */


  /* Close File */
  fclose( fp );
}  /* void ComputeOnlyCNSpatialDistr() */
/*---------------------------------------------------------------------------*/

/* Compute the spatial distribution of various (17 listed above) CNpair types.
   To limit memory used up by counting array (bins) we compute the distribution
   for each pair separately and write the raw count into 17 different files 
   suffixed by the pair type.  NOTE: COLUMNS OF hpresent are Box-Vectors */
void ComputeCNpairSpatialDistribution( char *outFile, float *hpresent )
{
  char     bufferArray[384];
  char     fileName[384];

  FILE     *fp;

  nodeType *np, *tnp;          /* temp pointer to head of list     */
  Atom_t   *ap;                /* Atom_t pointer                   */
  int      nbx, nby;           /* Bin coordinates X/Y              */
  int      L, M, N;            /* 3 digit pairs-type indx(421 etc) */
  int      i, j, k;

  int      len;

  /* Initialize pair-types */
  InitCNpairIDs();

  /* Get head of the serial list */
  tnp = (nodeType *) LCGetSerialList();

  /* Loop through the various CN pair types */
  for( i = 0; i < N_CN_PTYPES; i++ )   {  /* Loop-1 */
    /* (A) Initialize counting array */
    InitBinsUsedToCountCNpairTypes( hpresent );

    if( NumBinsX <= 0 || NumBinsY <= 0 ) 
      Error("(00) ComputeCNpairSpatialDistribution(): Invalid bins");

    /* (B) Get CN-pairtype into local L, M, N variables */
    L = CNpairID[i].l;
    M = CNpairID[i].m;
    N = CNpairID[i].n;

    /* June 8, 1998: Count only 200/211/300/311/322/411/422/433 pairs */
    if( (i!=LMN_200) && (i!=LMN_211) && (i!=LMN_300) && (i!=LMN_311) && 
	(i!=LMN_322) && (i!=LMN_411) && (i!=LMN_422) && (i!=LMN_433) )
      continue;

    /* (C) Open file (name padded with atmost 2 zeros) for distribution data */
    sprintf(fileName, "%02d.CNPAIR\0", i);
    if( !(fp = fopen(fileName, "w") ) )
      Error("(0) ComputeCNpairSpatialDistribution(): Unable to open file");

    /* (D) Loop over atoms  */
    for( np = tnp; np != NULL; np = np->nextInSerialList )   {
      uint  nn;
      int   icoPairFlag = NO; 

      /* Get atom pointer */
      ap = (Atom_t  *)  &(np->dx);

      /* Nov. 28, 1997: count only if it is MOBILE. ADDED movement_tag today */
      if( ap->movement_tag == FIX )  continue;

      /* Skip if atom is FCC */
      if( ap->cnTag == FCC ) continue;
      
      /* Get this Atom's bin indices */
      nbx = (int)  (ap->coord.x/BinSize);
      nby = (int)  (ap->coord.y/BinSize);
      
      if( (nbx >= NUM_X_BINS) || (nby >= NUM_Y_BINS) )
	Error("(1) ComputeCNpairSpatialDistribution(): Too many bins");
      
      /* Number of nearest neighbors */
      nn = (uint) ap->numNeigbors;
      if( nn > 19 ) 
	Error("(2) ComputeCNpairSpatialDistribution(): Too many neighbors");;

      /* Get the 3 digit designation (index) of each pair-type */
      for( j = 0; j < nn; j++ )   {
	/* Increment this bin value */
	if( (L == ap->pairIndexJ[j]) && (M == ap->pairIndexK[j]) &&
	    (N == ap->pairIndexL[j]) )  {
	  icoPairFlag = YES; 
	  NumCNPairTypeInBin[nbx][nby]++;
	  TotNumCNPairsOfEachType[i]++;
	}  /* if() */
      }  /* for(j) */

      /* June 8, 1998: Atoms in this bin */
      if( icoPairFlag == YES )  {
	NumAtomsInBin[nbx][nby]++;
	icoPairFlag = NO; 
      }
    }  /* for( np = tnp; np != NULL; np = np->nextInSerialList ) */

    /* (E) Print File Header */
    fprintf( fp, "#Total %d%d%dPr= %d;\n", L, M, N,
	     TotNumCNPairsOfEachType[i]  );

    fprintf( fp, "#Xcoo\t\tYcoo\t\tAtomsInBin\tRaw/Normalized %d%d%dPairs\n", 
	     L, M, N );
  
  /* (F) Loop through bins and print out the CNpairtype COUNT (RAW DATA) in 
     each bin to the output file */
  for( k = 0; k < NumBinsY; k++ ) {
    for( j = 0; j < NumBinsX; j++ ) {
      
      /* Nov 2, 1997: Print the number even if it's ZERO */
      
      /* All files will have X/Y coordinates also */
      if( NumAtomsInBin[j][k] > 0 ) {
	float nCNPairs;
	nCNPairs = ( ((float)NumCNPairTypeInBin[j][k])/
		     ((float)NumAtomsInBin[j][k]) );

	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t%9.4lf\n",
		 ((float) j) * BinSize, ((float) k) * BinSize, 
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k], nCNPairs );
      }
      else {
	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t0000.0000\n",
		 ((float) j) * BinSize, ((float) k) * BinSize,
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k] );
      }  /* else */
    } /* for(j) */

    fprintf(fp, "\n");     /* To make contour plots (using Gnuplot) we
			      need to leave a blank line between 2 successive
			      column (or row) of data */
  }   /* for(k) */
  /* June 8, 1998 */


#ifdef PRE_JUNE8_1998
    /* (F) Loop through bins and print out the CNpairtype COUNT (RAW DATA) in 
           each bin to the output file */
    for( k = 0; k < NumBinsY; k++ ) {
      for( j = 0; j < NumBinsX; j++ ) {

	/* Nov 2, 1997: Print the number even if it's ZERO */
	/* If no pair of a particular type; DON'T PRINT
	if( NumCNPairTypeInBin[j][k] == 0 )  continue; */

	/* All files will have X/Y coordinates also */
	fprintf( fp, "%lf\t%lf\t%d\n", ((float) j) * BinSize,
		 ((float) k) * BinSize, NumCNPairTypeInBin[j][k] );
      } /* for(j) */
      fprintf(fp, "\n");     /* To make contour plots (using Gnuplot) we
			   need to leave a blank line between 2 successive
			   column (or row) of data */
    }   /* for(k) */
#endif

    /* Close File */
    fclose( fp );
  } /*   for( i = 0; i < N_CN_PTYPES; i++ ):: Loop-1 */

  /* If a pair of any type is absent remove that file; else rename files
     as OutFile.lmnPAIR */
  for( i = 0; i < N_CN_PTYPES; i++ ) {

    /* June 8, 1998: Count only 200/211/300/311/322/411/422/433 pairs */
    if( (i!=LMN_200) && (i!=LMN_211) && (i!=LMN_300) && (i!=LMN_311) && 
	(i!=LMN_322) && (i!=LMN_411) && (i!=LMN_422) && (i!=LMN_433) )
      continue;

    if( TotNumCNPairsOfEachType[i] == 0 ) {
      sprintf(fileName, "rm -f %02d.CNPAIR\0", i);
      system( fileName );
    }
    else {  /* rename files outFile.lmnPAIR */
      /* Remove last 3 characters "con"; Copy name (without "con" into tmp[] */
      strcpy( fileName, outFile );

      len = strlen( fileName );
      
      /* Print file name with CNPAIR as suffix into fileName[] */
      sprintf( fileName+len-3, "%d%d%dCNPAIR\0",  
	       CNpairID[i].l, CNpairID[i].m, CNpairID[i].n );
      
      /* Concatenate all node sn files and compress the full file
	 name padded with atmost 3 zeros => can use wildcard substitution */
      
      sprintf( bufferArray, "mv %02d.CNPAIR %s; compress %s\0", i, fileName, fileName );
      
      system( bufferArray );
    } /* else */
  } /* for( i = 0; i < N_CN_PTYPES; i++ ) */


  /* The output is a MESS when we paste columns together */
#ifdef OLD_CODE
  
  /******** Join all files (stack them columnwise **********/

  /* If a pair of any type is absent remove that file */
  for( i = 0; i < N_CN_PTYPES; i++ ) {
    if( TotNumCNPairsOfEachType[i] == 0 ) {
      sprintf(fileName, "rm -f %02d.CNPAIR\0", i);
      system( fileName );
    }
  } /* for( i = 0; i < N_CN_PTYPES; i++ ) */


  /* Remove last 3 characters "con"; Copy name (without "con" into tmp[] */
  strcpy( fileName, outFile );

  len = strlen( fileName );

  /* Print file name with CNPAIR as suffix into fileName[] */
  sprintf( fileName+len-3, "CNPAIR\0" );

  /* Concatenate all node sn files and compress the full file
     name padded with atmost 3 zeros => can use wildcard substitution */

  sprintf( bufferArray,
	   "paste [0-9][0-9].CNPAIR >> %s;rm -f [0-9][0-9].CNPAIR;\0",
	   fileName );

  /* sprintf( bufferArray,
     "paste 00.CNPAIR [0-9][1-9].CNPAIR >> %s;rm -f [0-9][0-9].CNPAIR;\0",
     fileName ); */

  system( bufferArray );

#endif
}  /* void ComputeCNpairSpatialDistribution() */
/*---------------------------------------------------------------------------*/
/* Compute the spatial distribution of CN-555 + CN-544 pairs distribution.
   We discretize the X-Y plane into bins and count number of 555+544 CN-pairs
   in each bin.  This value is normalized on a per atom (in each bin) basis.
   We write NORMALIZED count into a file suffixed by the 544P555 pair type.  
   NOTE: COLUMNS OF hpresent are Box-Vectors */
void CN544Plus555PairDistribution( char *outFile, float *hpresent )
{
  char     buf1[384];
  char     buf2[384];
  char     fileName[384];

  FILE     *fp;

  nodeType *np, *tnp;          /* temp pointer to head of list     */
  Atom_t   *ap;                /* Atom_t pointer                   */
  int      nbx, nby;           /* Bin coordinates X/Y              */
  int      L, M, N;            /* 3 digit pairs-type indx 544      */
  int      O, P, Q;            /* 3 digit pairs-type indx 555      */
  int      j, k;

  int      len;

  /* Initialize pair-types */
  InitCNpairIDs();

  /* Get head of the serial list */
  tnp = (nodeType *) LCGetSerialList();

  /* (A) Initialize counting array */
  InitBinsUsedToCountCNpairTypes( hpresent );

  if( NumBinsX <= 0 || NumBinsY <= 0 ) 
    Error("(00) CN544Plus555PairDistribution(): Invalid bins");

  /* (B) Initialize CN-pairtype into local L, M, N variables */
  L = 5;  M = 4;  N = 4;
  O = 5;  P = 5;  Q = 5;
  
  /* (C) Open file (name padded with atmost 2 zeros) for distribution data */
  sprintf(fileName, "544P555.CNPAIR\0");
  if( !(fp = fopen(fileName, "w") ) )
    Error("(0) CN544Plus555PairDistribution(): Unable to open file");

  /* (D) Loop over atoms  */
  for( np = tnp; np != NULL; np = np->nextInSerialList )   {
    uint            nn;
    int             icoPairFlag = NO;   /* atom has 544/555 pair */

    /* Get atom pointer */
    ap = (Atom_t  *)  &(np->dx);
    
    /* Nov. 28, 1997: count only if it is MOBILE. ADDED movement_tag today */
    if( (ap->movement_tag == FIX) || (ap->cnTag == FCC) )  continue;
    
    /* Get this Atom's bin indices */
    nbx = (int)  (ap->coord.x/BinSize);
    nby = (int)  (ap->coord.y/BinSize);
    
    if( (nbx >= NUM_X_BINS) || (nby >= NUM_Y_BINS) )
      Error("(1) CN544Plus555PairDistribution(): Too many bins");
    
    /* Number of nearest neighbors */
    nn = (uint) ap->numNeigbors;
    if( nn > 19 ) 
      Error("(2) CN544Plus555PairDistribution(): Too many neighbors");;
    
    /* Get the 3 digit designation (index) of each pair-type */
    for( j = 0; j < nn; j++ )   {
      /* Increment this bin value */
      if( ((L==ap->pairIndexJ[j]) && (M==ap->pairIndexK[j]) && 
	   (N==ap->pairIndexL[j])) ) {
	/* Turn the flag that says ICO-Pair found in this bin */
	icoPairFlag = YES;

	NumCNPairTypeInBin[nbx][nby]++;
	TotNumCNPairsOfEachType[15]++;  /* increment 544 pairs             */
      } /* if() */
      else if( ((O==ap->pairIndexJ[j]) && (P==ap->pairIndexK[j]) &&
		(Q==ap->pairIndexL[j])) )   {
	/* Turn the flag that says ICO-Pair found in this bin */
	icoPairFlag = YES;

	NumCNPairTypeInBin[nbx][nby]++;
	TotNumCNPairsOfEachType[16]++;  /* increment 555 pairs             */
      }  /* else if() */
    }  /* for(j) */

    /* June 8, 1998: Atoms in this bin */
    if( icoPairFlag == YES )  {
      NumAtomsInBin[nbx][nby]++;
      icoPairFlag = NO; 
    }
  }  /* for( np = tnp; np != NULL; np = np->nextInSerialList ) */

  /* (E) Print File Header */
  fprintf( fp, "#Total %d%d%dPr= %d;\n", L, M, N,
	   TotNumCNPairsOfEachType[15]  );

  fprintf( fp, "#Total %d%d%dPr= %d;\n", O, P, Q,
	   TotNumCNPairsOfEachType[16]  );

  fprintf( fp, "#Xcoo\t\tYcoo\t\tAtomsInBin\tRaw/Normalized (%d%d%d+%d%d%d)Pairs\n", 
	   L, M, N, O, P, Q );
  
  
  /* (F) Loop through bins and print out the CNpairtype COUNT (RAW DATA) in 
     each bin to the output file */
  for( k = 0; k < NumBinsY; k++ ) {
    for( j = 0; j < NumBinsX; j++ ) {
      
      /* Nov 2, 1997: Print the number even if it's ZERO */
      /* If no pair of a particular type; DON'T PRINT
	 if( NumCNPairTypeInBin[j][k] == 0 )  continue; */
      
      /* All files will have X/Y coordinates also */
      if( NumAtomsInBin[j][k] > 0 ) {
	float nCNPairs;
	nCNPairs = ( ((float)NumCNPairTypeInBin[j][k])/
		     ((float)NumAtomsInBin[j][k]) );

	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t%9.4lf\n",
		 ((float) j) * BinSize, ((float) k) * BinSize, 
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k], nCNPairs );
      }
      else {
	fprintf( fp, "%9.4lf\t%9.4lf\t%3d\t%3d\t0000.0000\n",
		 ((float) j) * BinSize, ((float) k) * BinSize,
		 NumAtomsInBin[j][k], NumCNPairTypeInBin[j][k] );
      }  /* else */


    } /* for(j) */
    fprintf(fp, "\n");     /* To make contour plots (using Gnuplot) we
			      need to leave a blank line between 2 successive
			      column (or row) of data */
  }   /* for(k) */
  
  /* Close File */
  fclose( fp );

  /* If a 555+544 pair type is absent remove the file; else rename files
     as OutFile.544P555PAIR */
  if( (TotNumCNPairsOfEachType[15]+TotNumCNPairsOfEachType[16]) == 0 ) {
    sprintf(buf1, "rm -f %s\0", fileName);
    system( buf1 );
  }
  else {  /* rename files outFile.lmnPAIR */
    /* Remove last 3 characters "con"; Copy name (without "con" into tmp[] */
    strcpy( buf1, outFile );
    
    len = strlen( buf1 );
    
    /* Print file name with CNPAIR as suffix into fileName[] */
    sprintf( buf1+len-3, "%d%d%dP%d%d%dCNPAIR\0",  L, M, N, O, P, Q );
    
    sprintf( buf2, "mv %s %s; gzip %s\0", fileName, buf1, buf1 );
    
    system( buf2 );
  } /* else */

}  /* void CN544Plus555PairDistribution() */
/*---------------------------------------------------------------------------*/
#endif                                           /* #ifdef NOT_YET_REWRITTEN */
/*---------------------------------------------------------------------------*/
