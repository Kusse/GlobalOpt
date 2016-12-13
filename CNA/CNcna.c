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
 * $Date: 2000/01/29 18:58:52 $
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNcna.c,v $
 * $Revision: 1.6 $
 * $Log: CNcna.c,v $
 * Revision 1.6  2000/01/29 18:58:52  sgsrini
 * *** empty log message ***
 *
 * Revision 1.5  2000/01/29 18:54:36  sgsrini
 * Removed a lot of redundant arithmetic
 *
 * Revision 1.4  2000/01/29 18:31:44  sgsrini
 * Removed redundant arithmetic operations
 *
 * Revision 1.3  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
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
 *                   $$$$ CNcna.c: Main CNA functions                        *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*---------------------------------------------------------------------------*/
/* #define _LARGESYSTEM_ */
/*---------------------------------------------------------------------------*/
#undef _LARGESYSTEM_
/*---------------------------------------------------------------------------*/
#include "CNmain.h"
#include "CNmacros.h"
/*---------------------------------------------------------------------------*/
#define	NBND  (NNBR * (NNBR - 1))
/*---------------------------------------------------------------------------*/
/****** extern global function declarations *****************************/
extern void Error( char*, ... );
/*---------------------------------------------------------------------------*/
/****** extern global variable declarations *****************************/
extern Atom_t *AtomP;
/*---------------------------------------------------------------------------*/
/***** static variables ************************************************/
int4 BondTable[NBN2][NBN2];        /* BondTable[i][j]=1=> atom in bondlist[i]
                                      and bond..[j] are nbors of each other 
                                      0 => NOT neighbors; [i][i] = 0        */
/*---------------------------------------------------------------------------*/
/* *** static functions ********************************************** */
/*---------------------------------------------------------------------------*/
/* check if atom 1 is in the neighbor list of atom 2 */
static int4 AreNbrs( int4 ia1, Atom_t *ap2 )
{
  int4	i2;
  int4  nNbors = (int4) ap2->nNbors;

  /* loop over neighbors to a2 */
  for( i2 = 0; i2 < nNbors; i2++)  {
    /* if ap1 (i.e. ia1) is same as nbor to ap2, then return with value 1 */
    if (ia1 == ap2->nbor[i2])   return 1;
  }

  return 0;

}   /* static int4 AreNbrs( int4 ia1, Atom_t *ap2 ) */
/*---------------------------------------------------------------------------*/
/* Form a Nbor table of atoms in cbondlist (list of neighbors given by the
   2nd INDEX in the CNA pair label) */
static void FormBondTblFromSecondIndex(int4 *cbondlist, int4 nbonds, 
                                       int4 bondTable[][NBN2] )
{
  int4 i, j;
  int4 nbonds2 = 2*nbonds;

  for( i = 0; i < nbonds2; i++ )    {

    bondTable[i][i] = 0;                        /* No SELF Bonds     */
    
    for( j = i+1; j < nbonds2; j++ )    {
      bondTable[i][j] = bondTable[j][i] = 0;    /* DEFAULT: NOT NBOR */

      if( AreNbrs( cbondlist[i], (AtomP + cbondlist[j]) ) != 0 )
        bondTable[i][j] = bondTable[j][i] = 1;
    }
  }
} /* static void FormBondTblFromSecondIndex(int4 *cbondlist, int4 cnbBonds, 
                                            int4 **bborTable ) */
/*---------------------------------------------------------------------------*/
/* Form a Nbor table of atoms in cbondlist (list of neighbors given by the
   2nd INDEX in the CNA pair label) */
static int4 GetNborFromSecondIndexNborTbl(int4 bIn1, int4 bIn2 )
{

  return( BondTable[bIn1][bIn2] );                    /* 1 if bonded, else 0 */

}/*static void GetNborFromSecondIndexNborTbl(int4 *cbondlist, int4 cnbBonds) */
/*---------------------------------------------------------------------------*/
/* A RECURSIVE algorithm to compute longest path starting at atom1 
   among bonds in the list */
static int4 LongestContinuousPath(int4 ia, int4 *bondlist, int4 bonds )
{
  int4 bond[2]   = {0, 0};
  int4 ib        = 0,
       ic        = 0,
       id        = 0,
       leng      = 0, 
       max       = 0;
  

  /* fprintf(stderr, "%d {", ia );
     for( ib = 0; ib < bonds; ib++)    {
     if( bondlist[2 * ib] != -1 ) fprintf(stderr, "(%d,%d), ", 
     bondlist[2 * ib], bondlist[2 * ib + 1]);
     }
     fprintf(stderr, "}\n");
   */
  
  max = 0;
  for( ib = 0; ib < bonds; ib++)  {
    ic = 2 * ib;
    id = ic + 1;

    /* connection to first atom bond? */
    if(bondlist[ic] != -1 && 
       AreNbrs( ia, (AtomP + bondlist[ic]) ) != 0 )    {
      bond[0]      = bondlist[ic];
      bond[1]      = bondlist[id];
      bondlist[ic] = -1;
      bondlist[id] = -1;
      if((leng = 1 + 
          LongestContinuousPath(bond[1], bondlist, bonds)) > max )     {
        max = leng;
      }
      bondlist[ic] = bond[0];
      bondlist[id] = bond[1];
    }

    /* connection to second atom in bond? */
    if(bondlist[id] != -1 && 
       AreNbrs( ia, (AtomP + bondlist[id]) ) != 0)  {
      bond[0]      = bondlist[ic];
      bond[1]      = bondlist[id];
      bondlist[ic] = -1;
      bondlist[id] = -1;
      if((leng = 1 + 
          LongestContinuousPath(bond[0], bondlist, bonds)) > max )     {
        max = leng;
      }
      bondlist[ic] = bond[0];
      bondlist[id] = bond[1];
    }
  }
    
  return max;

} /* static int4 LongestContinuousPath(int4 ia, int4 *bondlist, int4 bonds)  */
/*---------------------------------------------------------------------------*/
/* *** global functions ************************************** */
/*---------------------------------------------------------------------------*/
/***************************************************************
   CNAmainLoop(): Performs common neighbor analysis.
 ************************************************************** */
int4 CNAmainLoop(int4 nAtoms, int4 verbose)
{
  float   invNatoms = 0., 
          err       = 0.;
  int4	  cnb       = 0;     /* Common NBor       */
  int4	  cnbBonds  = 0;     /* bonds between cnb */
  int4	  ia1       = 0, 
          ia2       = 0,
          ia12      = 0;     /* nbor of atom ia1 */
  int4	  in1       = 0,
          in2       = 0;
  int4	  leng      = 0, 
          max       = 0,
          min       = 99;

  Atom_t  *ap1      = NULL, 
          *ap2      = NULL;

  const int4 zero   = 0;
  int4    atomsSortedByBondOrder[NNBR];
  int4    numBonds[NNBR];            /* Number of BONDS formed by EACH COMMON 
					NBOR atoms given by first CNA index */

  int4    cbondlist[NBN2];
  int4    cbond[2];
  int4    cnlist[NNBR];          /* CN nbor list stores atom ID via the
				    index to its position in AtomP[] array */
				    
  /* for all atoms, reset number of calculated pairs */
  for( ia1 = 0; ia1 < nAtoms; ia1++ )    AtomP[ia1].pairs = 0; 

  invNatoms = 10.0 / nAtoms;
  err       = 1.1 * invNatoms;

  /* for all atoms, calculate indices */
  for( ia1 = 0; ia1 < nAtoms; ia1++ )        {
    int4 nNbors1;
    int4 tVb = 1;


    if( ia1 % 25000 == 0 )  {
      fprintf( stderr,
	      "$$$CNALoop(): Processing Atom %u\n", ia1 );
    }

    if( !verbose )  {  
      tVb =0;
      verbose=1;
    }

    if(verbose) {
      float tmp = ia1 * invNatoms;
      float diff = fabs( ((double) (tmp - (int)tmp) ));
      if(ia1 && diff < err)
        fprintf(stderr, "%2d%% Bonds Computed\n", (int)tmp * 10);
      /*  else  fprintf(stderr, "%g %g\n", ia1, diff); */
    }
    verbose=tVb;

    ap1     = AtomP + ia1;
    nNbors1 = (int4) ap1->nNbors;

    for( ia2 = 0; ia2 < nNbors1; ia2++ )    {


      /* treat each pair only once, so include only if ia1 < ap1->nbor[ia2] */
      ia12 = ap1->nbor[ia2];
      if( ia1 < ia12  )       { 
	int4 pairs1, pairs2;

	ap2 = AtomP + ia12;

	/* Assign to locals */
	pairs1 = ap1->pairs;
	pairs2 = ap2->pairs;

	/***********************************************
	 *   Compute 1st index: List of common nbr's   *
	 ***********************************************/
	/* loop over neighbors to ap1 */
	cnb = 0;
	for( in1 = 0; in1 < nNbors1; in1++ )	{
	  /* if nbr of ap1 is also a nbr of ap2, then insert in cnlist */
	  if( AreNbrs( ap1->nbor[in1], ap2 ) != 0 )	  {
	    cnlist[cnb++] = ap1->nbor[in1];
	  }
	}

	ap1->index[pairs1][0] = (uint2) cnb;
	ap2->index[pairs2][0] = (uint2) cnb;

	/***********************************************
	 * SPECIAL CASES: 3 or fewer common neighbors  *
	 ***********************************************/
	switch( cnb )  {
	case 0:                             /* NO common neighbors          */
	case 1:                             /* ONE common neighbor          */
	  ap1->index[pairs1][1] = ap2->index[pairs2][1] = 0;
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = 0;

	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;
	case 2:                      /* TWO common neighbors         */
	  /* Default these 2 atoms NOT nbors of each other */
	  ap1->index[pairs1][1] = ap2->index[pairs2][1] = 0;
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = 0;
	  
	  /* These 2 common neighbors also neighbors of each other?? */
	  if( AreNbrs(cnlist[0], (AtomP+cnlist[1])) != 0 )   {
	    ap1->index[pairs1][1] = ap2->index[pairs2][1] = 1;
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	  }

	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;
	case 3:                             /* THREE common neighbors       */
	  cnbBonds = 0;
	  /* Loop over common neighbor pairs; avoid double counting bonds */
	  for( in1 = 0; in1 < cnb; in1++ )
	    for( in2 = in1 + 1; in2 < cnb; in2++ )
	      if( AreNbrs(cnlist[in1], (AtomP+cnlist[in2])) != 0 ) cnbBonds++;

	  /* For three common neighbors case, the SECOND and THIRD 
	     CNA indices are identical */
	  ap1->index[pairs1][1] = ap2->index[pairs2][1] = cnbBonds;
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = cnbBonds;

	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;
	}
	/***********************************************/


	/***********************************************
	 *      FOR MORE THAN 3 COMMON NEIGHBORS       *
	 ***********************************************
	 * Compute 2nd index: Number of bonds between  *
	 * common neighbors; store bonds in bondlist   *
	 ***********************************************/
	/* Aug 28, 2002: Set all array elements (NNBR) to zero?? */
	memset((void*) numBonds, zero, (size_t)(cnb*sizeof(int4)) );
	cnbBonds = 0;
	for( in1 = 0; in1 < cnb; in1++ )   {   /* Loop over common neighbors */

	  /* Loop over remaining common neighbors; 
	     avoid double counting bonds */
	  for( in2 = in1 + 1; in2 < cnb; in2++ )   {

	    /* If nbr1 is also a nbor of nbr2, increment numBonds */
	    if( AreNbrs(cnlist[in1], (AtomP+cnlist[in2])) != 0 )   {
	      int4 cnb2   = 2 * cnbBonds,
	           cnb2P1 = cnb2 + 1;

	      if( cnb2P1 >= NBN2 )
		Error("$$$ CNAmainLoop(): ## Bond list overflow %d\n", 
		      cnb2P1 );
	      cbondlist[cnb2  ] = cnlist[in1];
	      cbondlist[cnb2P1] = cnlist[in2];

	      /* Increment Num Bonds Formed By CN atoms */
	      numBonds[in1]++;
	      numBonds[in2]++;

	      cnbBonds++;
	    }
	  } /* for( in2 = in1 + 1; in2 < cnb; in2++ ) */
	} /* for( in1 = 0; in1 < cnb; in1++ )  */
	
	ap1->index[pairs1][1] = (uint2) cnbBonds;
	ap2->index[pairs2][1] = (uint2) cnbBonds;


	/* Compute Min and Max Bonds per atom */
	memset((void*)atomsSortedByBondOrder,zero,(size_t)(cnb*sizeof(int4)));
	max = 0;
	min = 99;
	for( in1 = 0; in1 < cnb; in1++ )   {
	  min = MIN( min, numBonds[in1] );
	  max = MAX( max, numBonds[in1] );

	  /* Atoms Sorted by their BOND-ORDER */
	  atomsSortedByBondOrder[ numBonds[in1] ]++;
	}  


	/***********************************************
	 *   SPECIAL CASE: 2ND Index is 0, 1, 2, 3, 4  *
	 ***********************************************/
	switch( cnbBonds )  {
	case 0:                             /* NO common Bonds             */
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = 0;
	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;    /* SRINI Added Aug 28, 2002 */
	case 1:                             /* ONE common Bond              */
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;    /* SRINI Added Aug 28, 2002 */
	case 2:                             /* TWO common Bonds             */
	  switch( max )   {
	  case 1:                           /* (cnb, 2, 1) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	    break;    /* SRINI Added Aug 28, 2002 */
	  case 2:                           /* (cnb, 2, 2) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	    break;    /* SRINI Added Aug 28, 2002 */
	  }
	  continue;
	  break;    /* SRINI Added Aug 28, 2002 */
	case 3:                             /* THREE common Bonds           */
	  switch( max )   {
	  case 1:                           /* (cnb, 3, 1) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	    break;    /* SRINI Added Aug 28, 2002 */
	  case 2:                           /* (cnb, 3, 2) PAIR             */
	    switch( atomsSortedByBondOrder[ 2 ] )  {  /* Num Double bonded 
							 atoms              */
	    case 1:                         /* One double bonded atom       */
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;	
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	      break;    /* SRINI Added Aug 28, 2002 */
	    case 2:                         /* TWO double bonded atom       */
	    case 3:                         /* 3   double bonded atom       */
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 3;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	      break;    /* SRINI Added Aug 28, 2002 */
	    default: 
	      Error("$$$CNLoop(): 1a) Longest Bonds: PROBLEM\n");
	    }
	    continue;
	    break;    /* SRINI Added Aug 28, 2002 */
	  default:
	    Error("$$$CNLoop(): 1b) Longest Bonds: PROBLEM\n");
	  } /*  switch( max )  */
	  continue;                         /* NEXT CENTRAL pair (ap1, ap2) */
	  break;    /* SRINI Added Aug 28, 2002 */
	case 4:                             /* FOUR common Bonds            */
	  switch( max )   {
	  case 1:                           /* (cnb, 4, 1) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	  case 2:                           /* (cnb, 4, ??) PAIR            */
	    switch( atomsSortedByBondOrder[ 2 ] )  {  /* Num Double bonded 
							 atoms              */
	    case 1:                         /* 1 double bond atom: (cnb,4,2)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	      
	      /* case 2: With 2 double bond atom configurations have TWO 
		 variants, each with 4 single bond atoms.  But longest path
		 lengths are 3 and 2 respectively => NEED to explicitly
		 compute  the longest path */
	    case 3:                         /* 3 double bond atom: (cnb,4,3)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 3;

	      /* SRINI: 10/18/2002: Reset 3 to 4
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 4;   *//* 3; */

	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    case 4:                         /* 4 double bond atom: (cnb,4,4)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 4;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    } /* switch( atomsSortedByBondOrder[ 2 ] ) */
	    break;                          /* 'break' to handle case 2:    */

	  case 3:                           /* triply bonded atom           */
	    switch( atomsSortedByBondOrder[ 2 ] )  {  /* Num Double bonded 
							 atoms              */
	    case 0:                         /* 0 double bond atom: (cnb,4,2)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    case 2:                         /* 2 double bond atom: (cnb,4,3)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 3;

	      /* SRINI: 10/18/2002: Reset 3 to 4 ???
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 4;  */ /* 3; */

	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    default: 
	      Error("$$$CNLoop(): 2) Longest Bonds: PROBLEM\n");
	    } /* switch( atomsSortedByBondOrder[ 2 ] ) */
	    continue;
	    break;
	  } /* switch(max) */
	  break;                            /* switch( cnbBonds ): case 4:  */

	case 5:                             /* FIVE common Bonds            */
	  switch( max )   {                 /* Max Bond Order               */
	  case 1:                           /* (cnb, 5, 1) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	  case 2:                           /* (cnb, 5, ??) PAIR            */
	    switch( atomsSortedByBondOrder[ 2 ] )  {  /* Num Double bonded 
							 atoms              */
	    case 1:                         /* 1 double bond atom: (cnb,5,2)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	      
	      /* case 2/3/4 from by explicit call to longest path function */
	    case 5:                         /* 5 double bond atom: (cnb,5,5)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 5;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    } /* switch( atomsSortedByBondOrder[ 2 ] ) */
	    break;                          /* 'break' to handle other cases*/
	  } /* switch(max) */

	break;                              /* switch( cnbBonds ): case 5:  */

	case 6:                             /* SIZ common Bonds             */
	  switch( max )   {                 /* Max Bond Order               */
	  case 1:                           /* (cnb, 6, 1) PAIR             */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = 1;
	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	  case 2:                           /* (cnb, 6, ??) PAIR            */
	    switch( atomsSortedByBondOrder[ 2 ] )  {  /* Num Double bonded 
							 atoms              */
	    case 1:                         /* 1 double bond atom: (cnb,6,2)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 2;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	      
	      /* case 2/3/4 from by explicit call to longest path function */
	    case 6:                         /* 6 double bond atom: (cnb,6,6)*/
	      ap1->index[pairs1][2] = ap2->index[pairs2][2] = 6;
	      ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	      continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	    } /* switch( atomsSortedByBondOrder[ 2 ] ) */
	    break;                          /* 'break' to handle other cases*/
	  } /* switch(max) */

	break;                            /* switch( cnbBonds ): case 5:  */

#ifdef WRONG
	default:
	  if( cnbBonds > cnb )    {         /* POSSIBLE Weird Geometry      */
	    ap1->index[pairs1][2] = ap2->index[pairs2][2] = cnb;

	    ap1->pairs++;  ap2->pairs++;    /* ready for next bond          */
	    continue;                       /* NEXT CENTRAL pair (ap1, ap2) */
	  }
	  else
	    break;                          /* Proceed to Longest Path      */
#endif
	}  /* switch( cnbBonds ) */
	/***********************************************/
	  
	/***********************************************
	 * SPECIAL CASE: IF 1st Index<= 5 and 2nd index*
	 *               > 1st index => 3rd INDEX = 1st*
	 ***********************************************/
	if( (cnb <= 5) && (cnbBonds > cnb) )       {
	  ap1->index[pairs1][2] = ap2->index[pairs2][2] = cnb;
	  ap1->pairs++;  ap2->pairs++;  /* ready for next bond          */
	  continue;                     /* NEXT CENTRAL pair (ap1, ap2) */
	}
	/***********************************************/


#ifdef _LARGESYSTEM_
	/***********************************************
	 * SPECIAL CASE: IF 1st and 2nd Index are SAME *
	 *               3rd INDEX is ALSO IDENTICAL   *
	 ***********************************************/
	if( cnb == cnbBonds )    {
	  ap1->index[pairs1][2] = (uint2) cnbBonds;
	  ap2->index[pairs2][2] = (uint2) cnbBonds;
	  
	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* Skip to next (ap1, ap2) pair */
	}
	/***********************************************
	 * SPECIAL CASE: IF 1st INDEX= N and 2nd INDEX *
	 *               = N-1, 3RD INDEX = N-1        *
	 ***********************************************/
	else if( cnbBonds == (cnb - 1) )   {
	  ap1->index[pairs1][2] = (uint2) cnbBonds;
	  ap2->index[pairs2][2] = (uint2) cnbBonds;

	  ap1->pairs++;  ap2->pairs++;	    /* ready for next bond          */
	  continue;                         /* Skip to next (ap1, ap2) pair */
	}
	/***********************************************/
#endif                                              /* LARGESYSTEM           */



	/***********************************************
	 *  FOR OTHER GENERAL CASES COMPUTE 3RD INDEX  *
	 ***********************************************
	 * LongestContinuous chain of bonds among atom *
	 * in the list given by the 2nd Index.         *
	 ***********************************************/
	if( nAtoms < 500 )   {
	  if( verbose )  {
	    int4 izz;
	    fprintf( stderr, "$$$CNALoop(): Pair(%u %u), cnBonds-%d: cNbor: ",
		    ia1, ia12, 2*cnbBonds );
	    for( izz = 0; izz < cnbBonds; izz++ )
	      fprintf( stderr, 
		      "%d %d ", cbondlist[2*izz], cbondlist[2*izz + 1]);
	    fprintf(stderr, "\n");
	  }
	}

	max = 0;
	/* loop over bonds between common neighbors */
	for( in1 = 0; in1 < cnbBonds; in1++ )	{
	  int4 ic = 2 * in1,
	       id = ic + 1;

	  /* save bond temporarily */
	  cbond[0] = cbondlist[ic];
	  cbond[1] = cbondlist[id];

	  /* Remove bond from rest list */
	  cbondlist[ic] = -1;
	  cbondlist[id] = -1;

	  /* Look for path originating in 1st atom */
	  if((leng = 1 + 
	      LongestContinuousPath(cbond[0],cbondlist,cnbBonds)) > max)  {
	    max = leng;
	  }

	  /* Look for path originating in 2nd atom */
	  if((leng = 1 + 
	      LongestContinuousPath(cbond[1],cbondlist,cnbBonds)) > max)  {
	    max = leng;
	  }

	  /* Restore atom to the list */ 
	  cbondlist[ic] = cbond[0];
	  cbondlist[id] = cbond[1];
	} /* for( in1 = 0; in1 < cnbBonds; in1++ ) */

	ap1->index[pairs1][2] = (uint2) max;
	ap2->index[pairs2][2] = (uint2) max;

	/* ready for next bond */
	ap1->pairs++;
	ap2->pairs++;
	
    } /* if( ia1 < ia12 )                                   */
    }   /* Atom-2 Loop: for( ia2 = 0; ia2 < nNbors1; ia2++ )  */
  }     /* Atom-1 Loop: for( ia1 = 0; ia1 < nAtoms; ia1++  )  */

  return 0;

}  /* CNAmainLoop() */ 
/*---------------------------------------------------------------------------*/
/* #undef _LARGESYSTEM_ */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef DBG
/*---------------------------------------------------------------------------*/
/* A RECURSIVE algorithm to compute longest path starting at atom1 
   among bonds in the list */
static int4 LongestContinuousPath(int4 ia, int4 *bondlist, int4 bonds,
				  int4 blIndex )
{
  int4 bond[2]   = {0, 0};
  int4 ib        = 0,
       ic        = 0,
       id        = 0,
       leng      = 0, 
       max       = 0;
  

  /* fprintf(stderr, "%d {", ia );
     for( ib = 0; ib < bonds; ib++)    {
     if( bondlist[2 * ib] != -1 ) fprintf(stderr, "(%d,%d), ", 
     bondlist[2 * ib], bondlist[2 * ib + 1]);
     }
     fprintf(stderr, "}\n");
   */
  
  max = 0;
  for( ib = 0; ib < bonds; ib++)  {
    ic = 2 * ib;
    id = ic + 1;

    /* connection to first atom bond? */
    if(bondlist[ic] != -1 && 
       GetNborFromSecondIndexNborTbl( blIndex, ic ) != 0 )    {
      bond[0]      = bondlist[ic];
      bond[1]      = bondlist[id];
      bondlist[ic] = -1;
      bondlist[id] = -1;
      if((leng = 1 + 
	  LongestContinuousPath(bond[1], bondlist, bonds, id)) > max )     {
	max = leng;
      }
      bondlist[ic] = bond[0];
      bondlist[id] = bond[1];
    }

    /* connection to second atom in bond? */
    if(bondlist[id] != -1 && 
       GetNborFromSecondIndexNborTbl( blIndex, id ) != 0 )    {
      bond[0]      = bondlist[ic];
      bond[1]      = bondlist[id];
      bondlist[ic] = -1;
      bondlist[id] = -1;
      if((leng = 1 + 
	  LongestContinuousPath(bond[0], bondlist, bonds, ic)) > max )     {
	max = leng;
      }
      bondlist[ic] = bond[0];
      bondlist[id] = bond[1];
    }
  }
    
  return max;
} /* static int4 LongestContinuousPath(int4 ia, int4 *bondlist, int4 bonds)  */
/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
