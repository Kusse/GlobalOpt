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
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNpairStat.h,v $
 * $Revision: 1.3 $
 * $Log: CNpairStat.h,v $
 * Revision 1.3  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
 *
 * Revision 1.2  2000/01/27 21:39:21  sgsrini
 * Added functions to compute and write CNN pair statistics
 *
 * Revision 1.1  2000/01/27 04:17:09  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *                  $$$$ CNpairStat.h: CNN pair statistics header            *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
/* #define BIN_SIZE    0.05 */                 /* discretizing bin size     */
/*---------------------------------------------------------------------------*/
/* Lattice parameter, alat, is 1.537762 (EQ. FCC has 1.545633 as alat)      */
#define BIN_SIZE       1.1                     /* bin size: 1st peak of RDF
						  used for CN-anal          */
#define N_CN_PTYPES    17                      /* # of diff CN-pair types   */
#define NUM_BINS       500                     /* Max # bins                */
#define NUM_X_BINS     NUM_BINS
#define NUM_Y_BINS     NUM_BINS
/*---------------------------------------------------------------------------*/
#define N_TYPE         8
/*---------------------------------------------------------------------------*/
#define LMN_000        0
#define LMN_100        1
#define LMN_200        2
#define LMN_211        3
#define LMN_300        4
#define LMN_311        5
#define LMN_322        6
#define LMN_400        7
#define LMN_411        8
#define LMN_421        9
#define LMN_422        10
#define LMN_433        11
#define LMN_444        12
#define LMN_532        13
#define LMN_533        14
#define LMN_544        15
#define LMN_555        16
/*---------------------------------------------------------------------------*/
/* 1. CN-pair type */
typedef struct cnpT  {
  int l, m, n;
} cnpType;
/*---------------------------------------------------------------------------*/
/* 2. Coordination number (CN) parameters:  Although we've defined variables
      to handle a two-component system, we've implemented CN anal for ONLY
      ONE-COMPONENT SYSTEM. */
typedef struct cnp    {
  float avneb;                  /* Average number of neighbors              */
  float PercoordnAA[NNBR];      /* % A-atoms with CN 0->maxCoordNum         */

#ifdef NOT_YET
  float AvgEnergyofAnFCCAtom;   /* Energy per atom for a MOBILE FCC atom    */
  float AvgEnergyofANonFCCAtom; /* Energy per atom for MOBILE NON-FCC atom  */
  float SumPairPEforEachCN[FCC+1];/* Sum of pairPEs of atoms with each CN
                                    (0->FCC+1). Sum of the pairPE of all 
				    MOBILE FCC atoms is in 
				    SumPairPEforEachCN[FCC]. At the end of
                                    ComputePairPEforEachCN() this array will
                                    contain pairPE on a PER ATOM BASIS       */
  int4  NumAtomsofEachCN[FCC+1];/* Atoms with each CN used to compute pairPE
				    per atom from SumPairPE...[]. Number of all
                                    NON-FIXED FCC in NumAtomsofEachCN[FCC]   */
#endif

  int4  maxCoordNum;            /* Max value of CN; size of NcoordnumbM[]    */
  uint4 NcoordnumbM[NNBR];      /* # of atoms with CN 0->NNBR - A and B    */

  uint4 Nbcc[N_TYPE];           /* # of bcc        coordinated atoms (CN=14) */
  uint4 Nfcc[N_TYPE];           /* # of fcc        coordinated atoms (CN=12) */
  uint4 Nhcp[N_TYPE];           /* # of hcp        coordinated atoms (CN=12) */
  uint4 Nico[N_TYPE];           /* # of icosahedra coordinated atoms         */
  uint4 Nicodist[N_TYPE];       /* # of distorted  icosahedra coord atoms    */

  uint4 NatmInFullICO;          /* Number of atoms in complete ICOsahedra    */
  uint4 NatmInDistICO;          /* Number of atoms in distorted ICOsahedra   */

  uint4 NumDiffPairs[NNBR][NNBR][NNBR];
                                /* Pairs having index JKL will be counted in
				   this array. Ex: NumDiffPairs[J][K][L] will
				   give us total number of pairs with the 3 
				   digit label JKL */
}  CNstat_t;
/*---------------------------------------------------------------------------*/
