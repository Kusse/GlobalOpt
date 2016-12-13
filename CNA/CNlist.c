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
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNlist.c,v $
 * $Revision: 1.5 $
 * $Log: CNlist.c,v $
 * Revision 1.5  2000/01/28 19:58:36  sgsrini
 * Working version with uint2 variables in Atom_t and dynamic mem alloc
 * DON'T DYNAMICALLY allocate MULTI-DIMENSIONAL arrays for uint2/int2/char*
 * Pointersizes are at least 4 bytes (larger than these variable sizes)
 *
 * Revision 1.4  2000/01/28 16:03:15  sgsrini
 * Replaced all double variables with floats to save memory
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
 *              $$$$ CNlist.c: Link-List based neighbor lists                *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*---------------------------------------------------------------------------*/
#include "CNmain.h"
#include "CNmacros.h"
#include "CNpbc.h"
/*---------------------------------------------------------------------------*/
#define MIN_DIST ((float) 0.5)
/*---------------------------------------------------------------------------*/
/****** extern global function declarations *****************************/
extern void Error( char*, ... );
/*---------------------------------------------------------------------------*/
/****** extern global variable declarations *****************************/
extern Atom_t *AtomP;
/*---------------------------------------------------------------------------*/
/* ***  static variables ************************************************/
/*---------------------------------------------------------------------------*/
static LCType   LinkCell;
static int4     LinkCellCreated = 0;
static int4     Nbrmax          = NNBR;
/*---------------------------------------------------------------------------*/
/* Used to search slots during neighbor forming */
static int4 Direction[13][3] = {   {-1,-1, 1}, { 0,-1, 1}, { 1,-1, 1},
				   {-1, 0, 1}, { 0, 0, 1}, { 1, 0, 1},
				   {-1, 1, 1}, { 0, 1, 1}, { 1, 1, 1},
                                   {-1, 1, 0}, { 0, 1, 0}, { 1, 1, 0},
                                   { 1, 0, 0} };
/*---------------------------------------------------------------------------*/
/*****************************************************************************
 * LCInitialize() Sets Up a Skeletal Frame-Work for Cells. Initially a List  *
 * is Created With Maximum Possible Number of Slots. Next, LCUpdateCellSize()*
 * is Called to Fine-Tune the Slot Sizes by Specifying the Correct Number of *
 * Slots along three box edges.                                              *
 *        CAUTION: LinkCell framework must be created FIRST                  *
 *****************************************************************************/
void LCInitialize( int4 verbose )
{
  void  *pOut  = NULL;
  int4  zero   = 0;
  
  if( LinkCellCreated )
    Error("$$$ InitializeLinkCells(): LinkCells already initialized");
 
  /* Zero structure using memset() instead of bzero() */
  pOut = (void*) memset( (void*) &LinkCell, zero, sizeof(LCType) );
    
  /* Check return value */
  if( pOut != ((void*) &LinkCell) )  {
    fprintf( stderr, "$$$ InitializeLinkCells(): pOut != LinkCell\n" );
    fflush( stderr );
    Error("$$$ InitializeLinkCells(): memset() problem........\n");
  }

  /* We use a UNIT-CUBE simulation box */
  LinkCell.BoxSizeX  = (float) BOX_SIZE_X;     /* size of the simulation box */
  LinkCell.BoxSizeY  = (float) BOX_SIZE_Y;
  LinkCell.BoxSizeZ  = (float) BOX_SIZE_Z;

  /* Slot Size */
  LinkCell.SlotSizeX = LinkCell.BoxSizeX/((double) MAX_CELLS_X);
  LinkCell.SlotSizeY = LinkCell.BoxSizeY/((double) MAX_CELLS_Y);
  LinkCell.SlotSizeZ = LinkCell.BoxSizeZ/((double) MAX_CELLS_Z);

/* #ifdef DOES_NOT_WORK */
  LinkCell.LinkCellArray = NULL;
/* #endif */

  /* Num cells */
  LinkCell.NslotsX   = MAX_CELLS_X;
  LinkCell.NslotsY   = MAX_CELLS_Y;
  LinkCell.NslotsZ   = MAX_CELLS_Z;

  /* Cell structure initialized */
  LinkCellCreated = 1;

  if( verbose )
    fprintf(stderr, "$$$ InitializeLinkCells(): LinkCells Initialized\n");

}  /* void LCInitialize( void ) */
/*---------------------------------------------------------------------------*/
/*****************************************************************************
 * LCInitializeSlots(): initialize ACTUAL number of slots along X/Y/Z Dirs   *
 *****************************************************************************/
void LCInitializeSlots(int4 nSlotsX, int4 nSlotsY, int4 nSlotsZ, uint4 nAtoms,
		       int4 verbose )
{
  unsigned long memAlloc1 = 0;
  unsigned long memAlloc2 = 0;
  unsigned long memAlloc3 = 0;
  int4 i, j;

  if( !LinkCellCreated ) 
    Error("LCInitializeSlots(): LinkCell not yet created");


#ifdef NOW_DYNAMIC_MEM_ALLOC
  if(nSlotsX < 2 || nSlotsX >= MAX_CELLS_X ||
     nSlotsY < 2 || nSlotsY >= MAX_CELLS_Y ||
     nSlotsZ < 2 || nSlotsZ >= MAX_CELLS_Z )
    Error("$$$ LCInitializeSlots): Invalid number of slots (%d %d %d)\n",
	  nSlotsX, nSlotsY, nSlotsZ );
#endif

  if(nSlotsX < 2 || nSlotsY < 2 || nSlotsZ < 2  )
    Error("$$$ LCInitializeSlots): Invalid number of slots (%d %d %d)\n",
	  nSlotsX, nSlotsY, nSlotsZ );

  LinkCell.Natoms          = nAtoms;

/* #ifdef DOES_NOT_WORK */

  /* int **array1 = (int **)malloc(nrows * sizeof(int *));
     for(i = 0; i < nrows; i++)
     array1[i] = (int *)malloc(ncolumns * sizeof(int)); */

  /* Allocate memory for slots** along X-axis */
  LinkCell.LinkCellArray       = ((SlotType***) 
				  calloc( nSlotsX, sizeof(SlotType**) ));

  if( LinkCell.LinkCellArray == NULL )
    Error("$$$ LCInitializeSlots(): Mem Alloc Problem-I\n");

  memAlloc1 = (nSlotsX*sizeof(SlotType**));
  fprintf( stderr, "$$$ LCInitializeSlots(): Allocated %9lu kBytes for ***LinkCellArray\n", (memAlloc1/1000 + 1) );


  /* Allocate memory for slots along Y-axis (columns) */
  for( i = 0; i < nSlotsX; i++ )   {
    LinkCell.LinkCellArray[i]    = ((SlotType**) 
				    calloc(nSlotsY, sizeof(SlotType*)) );

    if( LinkCell.LinkCellArray[i] == NULL )
      Error("$$$ LCInitializeSlots(): Mem Alloc Problem-II\n");

    memAlloc2 += (nSlotsY*sizeof(SlotType*));
  }

  fprintf( stderr, "$$$ LCInitializeSlots(): Allocated %9lu kBytes for **LinkCellArray\n", (memAlloc2/1000 + 1) );

  /* Allocate memory for slots along X/Y-axis (row/columns) */
  for( i = 0; i < nSlotsX; i++ )
    for( j = 0; j < nSlotsY; j++ )   {
      LinkCell.LinkCellArray[i][j] = ((SlotType*) 
				      calloc(nSlotsZ, sizeof(SlotType) ));

      if( LinkCell.LinkCellArray[i][j] == NULL )
	Error("$$$ LCInitializeSlots(): Mem Alloc Problem-III\n");
      
      memAlloc3 += (nSlotsZ*sizeof(SlotType));
    }

  fprintf( stderr, "$$$ LCInitializeSlots(): Allocated %9lu MBytes for *LinkCellArray\n", (memAlloc3/1000000 + 1) );
/* #endif */


  /* Initialize the number of slots to max_nslots in each direction. */
  LinkCell.NslotsX      = nSlotsX; 
  LinkCell.NslotsY      = nSlotsY;
  LinkCell.NslotsZ      = nSlotsZ;
  
  LinkCell.SlotSizeX    = LinkCell.BoxSizeX / (float) nSlotsX;
  LinkCell.SlotSizeY    = LinkCell.BoxSizeY / (float) nSlotsY;
  LinkCell.SlotSizeZ    = LinkCell.BoxSizeZ / (float) nSlotsZ;

  /* Reciprocal of slot size */
  LinkCell.SlotSizeInvX = 1.0/LinkCell.SlotSizeX;
  LinkCell.SlotSizeInvY = 1.0/LinkCell.SlotSizeY;
  LinkCell.SlotSizeInvZ = 1.0/LinkCell.SlotSizeZ;
 
  if( verbose )   {
    fprintf(stderr,
	    "$$$ LCInitializeSlots(): SlotSizeX/Y/Z = (%f, %f, %f)\n", 
	    LinkCell.SlotSizeX, LinkCell.SlotSizeY, LinkCell.SlotSizeZ); 

    fprintf(stderr,
	    "$$$ LCInitializeSlots(): (nsX, nsY, nsZ) = (%d, %d, %d)\n", 
	    LinkCell.NslotsX, LinkCell.NslotsY, LinkCell.NslotsZ);
    fflush( stderr );
  } 

}  /* 'LCInitializeSlots()'  */
/*---------------------------------------------------------------------------*/
/*****************************************************************************
 * ASSUMPTION: Coordinates in the USER Atom struct are Cartesian-coordinates *
 * RETURN    : Index of the atom's position in the AtomP array               *
 *****************************************************************************
 *         boxP     = Pointer to simulation box matrix.                      *
 *         boxInvP  = Pointer to inverse of simulation box matrix.           *
 *         xcooP    = Pointer to atom's cartesian coordinates.               *
 *         pe       = Atom's potential energy.                               *
 *         id       = Atom's numerical id.                                   *
 *         atomicNum= Atom's atomic number (distinguishes atoms in multi-    *
 *                                          component systems)               *
 *         tags     = [0] = F/M; [1] = B/R/G; [2] = not used now.            *
 *         verbose  = Verbose mode flag                                      *
 *****************************************************************************/
int LCInsertAtomIntoSlot(float    *boxP,  float *boxInvP, 
			 Vector_t *xcooP, float pe, uint4 id, int4 atomicNum, 
			 char     *tags,  int4  verbose )
{
  const float  aMinX[3]  = { -1.0, -1.0, -1.0 };
  const float  aX[3]     = {  1.0,  1.0,  1.0 };
  const float  a2X[3]    = {  2.0,  2.0,  2.0 }; 
  Vector_t     scoo;
  static int4  firstTime = 1;
  static uint4 ain       = 0;    /* Index of atom's position in AtomP array */
  Atom_t       *ap       = NULL;
  int4         sx, sy, sz;
  int4         pbcFlag   = INSIDE_THE_BOX;           /* Atom Inside the Box */
  int4         curIn     = 0;
  int4         nSlotAtoms= 0;

  static int4  nOOBatoms = 0;
  static FILE  *oobFP    = NULL;


  /* First call to this function allocate memory etc. */
  if( firstTime )   {
    firstTime = 0;

    if( !LinkCellCreated )  
      Error("LCInsertAtomIntoSlot(): LinkCell not yet created");

    if( LinkCell.Natoms == 0 )
      Error("$$$ LCInsertAtomIntoSlot(): Invalid number of atoms %u\n",
	    LinkCell.Natoms );

    /* Allocate and Initialize memory to zero (IMPORTANT) */
    AtomP = (Atom_t *) calloc( LinkCell.Natoms, sizeof(Atom_t) );

    if( !AtomP )
      Error("$$$ LCInsertAtomIntoSlot(): Cannot Allocate %u MBytes Memory\n",
	    ((LinkCell.Natoms * sizeof(Atom_t))/1000000 + 1) );

    
    fprintf( stderr,
	    "$$$ LCInsertAtomIntoSlot(): Allocated %u MBytes for Atoms\n",
	    ((LinkCell.Natoms * sizeof(Atom_t))/1000000 + 1) );

    /* Open OOB atom listing file */
    if (!( oobFP = fopen("CNoobAtom.log", "a+")))
       Error("LCInsertAtomIntoSlot():\tCan't open 'CNoobAtom.log'\n");

  } /* if( firstTime ) */

  /* Check if valid atom id */
  if( id < 1 || ain >= LinkCell.Natoms )
    Error("$$$ LCInsertAtomIntoSlot(): Invalid atom id/ArrayIndex=(%u, %u)\n",
	  id, ain  ); 

  /* Initialize Atom_t pointer */
  ap = AtomP + ain;

  /* Compute S-coord from cartesian coords */
  MATRIX_3D_VECTOR_MULTIPLY( scoo, boxInvP, (*xcooP) );  

  /* Apply Periodic Boundary Condition */
  APPLY_PBC( scoo, aX, a2X, aMinX, pbcFlag );

  /* Check if atom is way out of the box that it can't be wrapped back in */
  if( pbcFlag == OUT_OF_BOUNDS )  {
     Vector_t xcoo2;

     nOOBatoms++;

     /* Re-Compute Cartesian-coord from s-coords */
     MATRIX_3D_VECTOR_MULTIPLY( xcoo2, boxP, scoo );  

     fprintf(stderr, 
	     "$Insert(%d): OOBatom %u Coo (%f %f %f); Adjusted (%f %f %f)\n", 
	     nOOBatoms, id, (float)xcooP->x, (float)xcooP->y, (float)xcooP->z,
	     (float) xcoo2.x,      (float) xcoo2.y,  (float) xcoo2.z );

     fprintf(oobFP, 
	     "$Insert(%d): OOBatom %u Coo (%f %f %f); Adjusted (%f %f %f)\n", 
	     nOOBatoms, id, (float)xcooP->x, (float)xcooP->y, (float)xcooP->z,
	     (float) xcoo2.x,      (float) xcoo2.y,  (float) xcoo2.z );

     /* Error("$$$ LCInsertAtom(): Atom %u Coord (%f %f %f) out of bounds\n", 
	id, (float) xcooP->x, (float) xcooP->y, (float) xcooP->z ); */
 }
  /* Find the slot for this element  */
  sx = (int4) ( scoo.x * LinkCell.SlotSizeInvX );
  sy = (int4) ( scoo.y * LinkCell.SlotSizeInvY );
  sz = (int4) ( scoo.z * LinkCell.SlotSizeInvZ ); 


  /***** Print slot index *
  fprintf(stderr, "(%2d) scoo(%+.6f  %+.6f %+.6f)\n", id, scoo.x,scoo.y,scoo.z);
  fprintf(stderr, "(%2d) sxyz(%2d  %2d %2d)\n", id, sx,sy,sz);
  ***** Print slot index */



  /* Initialize slot info: Number of atoms in this slot */
  nSlotAtoms = LinkCell.LinkCellArray[sx][sy][sz].numAtomsInThisCell;
  LinkCell.LinkCellArray[sx][sy][sz].numAtomsInThisCell++;

  /* Check for array overflow */
  if( nSlotAtoms >= MAX_ATOMS_IN_SLOT )
    Error("$$$ LCInsertAtom(): Too many (%d) atoms in slot (%d %d %d)\n",
	  nSlotAtoms, sx, sy, sz );

#ifdef DBG
  if( verbose )  {
    fprintf( stderr,
	    "$$$ LCInsertAtomIntoSlot(): Insert: Id-%u, arIn-%u, Sxyz(%d %d %d), NatomsInCell-%d\n", id, ain, sx, sy, sz,
	    LinkCell.LinkCellArray[sx][sy][sz].numAtomsInThisCell );
  }
#endif

  /* Initialize slot table with atom array index "ain" */
  LinkCell.LinkCellArray[sx][sy][sz].indexInAtomArray[nSlotAtoms] = ain;

  /* Re-Compute Cartesian-coord from s-coords */
  MATRIX_3D_VECTOR_MULTIPLY( (ap->coo), boxP, scoo );  

  ap->atomicNum = (Uchar) atomicNum;
  ap->tags[0]   = tags[0];                          /* [0] = Fix(F)/Move('M')*/
  ap->tags[1]   = tags[1];                          /* [1] = 'B'/'R'/'G'     */
  ap->id        = id;                               /* Atom's numerical ID   */


#ifdef SMALL_SYSTEM
  /* Initialize other fields of ap */
  ap->potE      = pe;                               /* Atom Potential Energy */
#endif

  /* Increment array index */ 
  curIn = ain;   
  ain++;

  /* Clode OOB atom file at the end */
  if( id == LinkCell.Natoms ) fclose(oobFP);

  /* Return index of the position in array for this atom */
  return( curIn );

}    /* 'LCInsertAtom()'  */
/*---------------------------------------------------------------------------*/
/* CARTESIAN COORDINATE ARGUMENTS */
#define FORM_PAIR( xcoo1, xcoo2, in1, in2 )                                   \
{                                                                             \
  float dX, dY, dZ;                                                           \
									      \
  dX  = xcoo1.x - xcoo2.x; dY = xcoo1.y - xcoo2.y; dZ = xcoo1.z - xcoo2.z;    \
  rIJ = dX * dX + dY * dY + dZ * dZ;                                          \
									      \
  pairUp = 1;                                  /* Default is pairup */        \
									      \
  if( rIJ > rcutSQ ) pairUp = 0;                                              \
									      \
  if( pairUp )      {                                                         \
   /* Fill neighbor Table of Atom_t; increment num neighbors field            \
      num neighbors must start from zero as array index starts from 0.        \
      => post increment num neighbors field; Save pair-separation */          \
   if((ap1->nNbors >= ap1->pairs) || (ap2->nNbors >= ap2->pairs) )  {         \
      printf("FORM_PAIR(): Too many nbors; PairID(%u %u), #nNbors(%d %d)\n",  \
	     ap1->id, ap2->id, (int4)ap1->nNbors, (int4)ap2->nNbors );        \
      Error("FORM_PAIR(): Too many neighbors; nbor array overflow\n");        \
   }                                                                          \
                							      \
   /* Instead of saving,pointers to atoms we store the INDEX                  \
      of their position in the AtomP array */                                 \
   ap1->nbor[ (ap1->nNbors) ] = in2;    ap2->nbor[ (ap2->nNbors) ] = in1;     \
									      \
   ap1->nNbors++;   ap2->nNbors++;    /*nPairs++;   */                        \
									      \
   if(rIJ < MIN_DIST) fprintf(stderr,"Atoms(%5d,%5d): rIJ=%.4lf\n",           \
                              ap1->id, ap2->id, rIJ);                         \
  }                                                                           \
}
/*---------------------------------------------------------------------------*/
/* ASSUMPTION: Coordinates in the Atom struct are CARTESIAN coordinates.
           boxP      = array of 9 elements (simulation box vectors)
	   rcut      = upper bound for the pair cut-off distance 
	   verbose   = Print info in verbose mode flag */
void LCFormNbrTable( float *boxP, float rcut, int verbose )
{
  Vector_t   scoo2       = { 0., 0., 0. };
  Vector_t   xcoo2       = { 0., 0., 0. };
  float      boxXYZ[3]   = { 0., 0., 0. };
  float      cellVol     = 0.0;
  float      hVect[9]    = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float      hinvVect[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float      minRijSQ    = 9999.0;
  float      maxRijSQ    =    0.0;
  float      rIJ         = 0.0;                          /* pair separation  */
  float      rcutSQ      = 0.0;                          /* Square of cutOff */
  Atom_t     *ap1        = NULL,
             *ap2        = NULL;
  int4       in1, in2    = 0;
  int4       i           = 0;
  int4       ia          = 0;
  int4       ib          = 0;
  int4       nAtomsInThisCell= 0;
  int4       nAtomsInNborCell= 0;
  int4       dir         = 0,
             ni          = 0, 
             nj          = 0,
             nk          = 0;
  int4       nx          = 0, 
             ny          = 0,
             nz          = 0;
  int4       nSlotsX     = 0, 
             nSlotsY     = 0,
             nSlotsZ     = 0;                           /* numSlots in X/Y/Z */
  int4       pairUp      = 0;
  uint4      nPairs      = 0;
  uint4      nAtoms      = 0;

  uint2      minNbors    = 9999;
  uint2      maxNbors    = 0;

  /* extern void DbgPrintReadInAtoms( uint4 );
     if( verbose ) DbgPrintReadInAtoms( LinkCell.Natoms ); */


  if( !LinkCellCreated ) Error("LCFormNbrTable():LinkCell not yet created");
  

  /* Assign to local and Invert hVector matrix */
  for( i = 0; i < 9; i++ ) hVect[i] = boxP[i];
  INVERT_MATRIX( hinvVect, hVect, cellVol)

  /* MOTTO: SAVE CPU, AVOID sqrt() => DON'T compute distance => USE rIJ^2 */
  rcutSQ    = rcut * rcut;
  nSlotsX   = LinkCell.NslotsX; 
  nSlotsY   = LinkCell.NslotsY; 
  nSlotsZ   = LinkCell.NslotsZ;

  boxXYZ[0] = LinkCell.BoxSizeX;
  boxXYZ[1] = LinkCell.BoxSizeY;
  boxXYZ[2] = LinkCell.BoxSizeZ;

  if( verbose )   {
    fprintf( stderr,
	    "$$$ LCFormNbrTable(): rcutSQ=%f, Nslots(%d %d %d); BoxSz(%f %f %f)\n",
	    (float) rcutSQ, nSlotsX, nSlotsY, nSlotsZ, 
	    (float) boxXYZ[0], (float) boxXYZ[1], (float) boxXYZ[2] );
    fflush(stderr);
  }

  /* Initialize num nbors to zero */
  nAtoms = LinkCell.Natoms;
  for( ia = 0; ia < nAtoms; ia++ ) AtomP[ia].nNbors = 0;

  /* Loop over all the slots */
  for( nx = 0; nx < nSlotsX; nx++ )
    for( ny = 0; ny < nSlotsY; ny++ )
      for( nz = 0; nz < nSlotsZ; nz++ )       {

	/* Number of atoms in this slot */
	nAtomsInThisCell = 
	  LinkCell.LinkCellArray[nx][ny][nz].numAtomsInThisCell;

	/* if( verbose )   {
	   fprintf( stderr,
	   "$$$ LCFormNbrTbl(): nx/y/z(%d %d %d); nAtomsInCell(%d)\n",
	   nx, ny, nz, nAtomsInThisCell );
	   fflush(stderr);
	   } */

	/* If cell is empty; skip to next cell */
	if( !nAtomsInThisCell ) continue;
	
        /* 1. Pair each element in this slot with other atoms 
	      also within this SAME CELL */
	for( ia = 0; ia < nAtomsInThisCell; ia++ )     {

	  /* Array Index of Atom1 */
	  in1 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ia];


#ifdef DBG
	  if( verbose )   {
	    fprintf( stderr,
		    "$$$ LCFormNbrTbl(): nx/y/z(%d %d %d); nAtomsInCell(%d): ia-%d, IndexInAtomArray-%u\n",
		    nx, ny, nz, nAtomsInThisCell, ia, in1 );
	    fflush(stderr);
	  }
#endif


	  /* Pointer to atom-1 in pair */
	  ap1 = AtomP + in1;

	  for( ib = ia+1; ib < nAtomsInThisCell; ib++ )    {

	    /* Array Index of Atom1 */
	    in2 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ib];
	    
	    /* Pointer to atom-2 in pair */
	    ap2 = AtomP + in2;

            /* Not a neighbor pair if rIJ > rcut */
            FORM_PAIR((ap1->coo), (ap2->coo), in1, in2);
	  }  /* for( ib = ia+1; ; )                          */
	}    /* for( ia = 0; ; ): Pair up atoms in same slot */


        /* 2. For slot (nx, ny, nz), consider 13 of its neighboring slots.
              Find pairs (p,q) such that 'p' is in (nx, ny, nz) and 'q'
              is in its neighboring slot (ni, nj, nk).  */
        for( dir = 0; dir < 13; dir++ ) {
          float sx = 0., sy = 0., sz = 0.;            /* Box edge for PBC*/

          ni = nx + Direction[dir][0];
          nj = ny + Direction[dir][1];
          nk = nz + Direction[dir][2];

          /* handle periodic boundary conditions. */
          if( ni < 0 )  {                              /* PBC along X     */
            sx = - boxXYZ[0];                          /* boxXYZ[0] = 1.0 */
            ni = (ni+nSlotsX) % nSlotsX;
          } 
          else if( ni >= nSlotsX )   {
            sx = boxXYZ[0];
            ni = (ni+nSlotsX) % nSlotsX;
          }
          
          if( nj < 0 )  {                              /* PBC along Y     */
            sy = - boxXYZ[1];                          /* boxXYZ[1] = 1.0 */
            nj = (nj+nSlotsY) % nSlotsY;
          }
          else if( nj >= nSlotsY )  {
            sy = boxXYZ[1];
            nj = (nj+nSlotsY) % nSlotsY;
          }
          
          if( nk < 0 )  {                              /* PBC along Z     */
            sz = - boxXYZ[2];                          /* boxXYZ[2] = 1.0 */
            nk = (nk+nSlotsZ) % nSlotsZ;
          } 
          else if( nk >= nSlotsZ ) {
            sz = boxXYZ[2];
            nk = (nk+nSlotsZ) % nSlotsZ;
          }


	  if(ni < 0 || ni >= nSlotsX ||
	     nj < 0 || nj >= nSlotsY ||
	     nk < 0 || nk >= nSlotsZ )
	    Error("$$$ LCFormNbrTable(): This Slot(%d %d %d); Invalid Nbor Slot(%d %d %d)\n", nx, ny, nz, ni, nj, nk);
	  

	  /* Pair each element in this slot (nx,ny,nz) with  atoms 
	     in neighboring slot (ni, nj, nk)  */
	  /* Loop over this slot (nx,ny,nz) */
	  for( ia = 0; ia < nAtomsInThisCell; ia++ )     {

	    /* Array Index of Atom1 */
	    in1 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ia];
	    
	    /* Pointer to atom-1 in pair */
	    ap1 = AtomP + in1;
	
	    /* Number of atoms in nbor cell */
	    nAtomsInNborCell =
	      LinkCell.LinkCellArray[ni][nj][nk].numAtomsInThisCell;
	    
	    /* Loop atoms in neighboring slot (ni, nj, nk) */
	    for( ib = 0; ib < nAtomsInNborCell; ib++ )    {
	      
	      /* Array Index of Atom1 */
	      in2 = LinkCell.LinkCellArray[ni][nj][nk].indexInAtomArray[ib];
	      
	      /* Pointer to atom-2 in pair */
	      ap2 = AtomP + in2;

              /* assign ap2 cartesian coord to xcoo */
              ASSIGN_3D_VECTORS( xcoo2, (ap2->coo) );

              /* Compute S-Coord for ap2 */
              MATRIX_3D_VECTOR_MULTIPLY( scoo2, hinvVect, xcoo2 );  

              /* Apply PBC */
              scoo2.x += sx;
              scoo2.y += sy;
              scoo2.z += sz;

              /* Convert to CARTESIAN COORD */
              MATRIX_3D_VECTOR_MULTIPLY( xcoo2, hVect, scoo2 ); 

	      /* Not a neighbor pair if rIJ > rcut */
	      FORM_PAIR((ap1->coo), (xcoo2), in1, in2);

	    }  /* for( ib = 0; ; )                             */
	  }    /* for( ia = 0; ; ): Pair up atoms in NBOR slot */
        } /* for ( dir = 0; dir < 13; dir++ ) */
      } /* for each slot (i,j,k) -- three nested loops */

  /******* DEBUG *******/
  printf("$$$ LCFormNbrTable(): Natoms=%u\n", LinkCell.Natoms );
  /******* DEBUG *******/


#ifdef NOT_NOW
  /* Min/Max num neighbors */
  nAtoms = LinkCell.Natoms;
  for( ia = 0; ia < nAtoms; ia++ )   {
    uint2 nnnn;
    ap1 = AtomP + ia;

    nnnn = ap1->nNbors;
    minNbors = MIN( nnnn, minNbors );
    maxNbors = MAX( nnnn, maxNbors );
  }

  /* Minimum and Maximum nbors/separation */
  LinkCell.MinRij   = sqrt( (double) minRijSQ );
  LinkCell.MaxRij   = sqrt( (double) maxRijSQ );
  LinkCell.MinNbors = (int4) minNbors;
  LinkCell.MaxNbors = (int4) maxNbors;

  /* if( verbose )
     fprintf( stderr, "$$$ LCFormNbrTable(): Number of Pairs = %u\n",
     nPairs ); */

  if( verbose )   {
    fprintf( stderr,
	    "$$$ LCFormNbrTable(): min/maxNbors=(%d %d), Nslots(%d %d %d)\n",
	    (int) minNbors, (int) maxNbors, nSlotsX, nSlotsY, nSlotsZ );

    fprintf( stderr,
	    "$$$ LCFormNbrTable(): rcutSQ/min/maxRijSQ=(%f %f %f), BoxSz(%f %f %f)\n",
	    (float) rcutSQ, (float) minRijSQ, (float) maxRijSQ, 
	    (float) boxXYZ[0], (float)  boxXYZ[1], (float)  boxXYZ[2] );
    fflush(stderr);
  }
#endif

}   /* LCFormNbrTable()  */
/*---------------------------------------------------------------------------*/
LCType *LCgetLinkCellP( void )
{
  return( (LCType *) &LinkCell );
}
/*---------------------------------------------------------------------------*/
/* CARTESIAN COORDINATE ARGUMENTS */
#define CHECK_IF_NBOR( xcoo1, xcoo2)                                          \
{                                                                             \
  float dX, dY, dZ;                                                           \
									      \
  npairs++;                                                                   \
  dX  = xcoo1.x - xcoo2.x; dY = xcoo1.y - xcoo2.y; dZ = xcoo1.z - xcoo2.z;    \
  rIJ = dX * dX + dY * dY + dZ * dZ;                                          \
									      \
  pairUp = 1;                                  /* Default is pairup */        \
									      \
  if( rIJ > rcutSQ ) pairUp = 0;                                              \
									      \
  if( pairUp )      {                                                         \
   ap1->nNbors++;   ap2->nNbors++;                                            \
									      \
   minRijSQ = MIN( rIJ, minRijSQ );                                           \
   maxRijSQ = MAX( rIJ, maxRijSQ );                                           \
									      \
  }                                                                           \
}
/*---------------------------------------------------------------------------*/
/* ASSUMPTION: Coordinates in the Atom struct are CARTESIAN coordinates.
           boxP      = array of 9 elements (simulation box vectors)
	   rcut      = upper bound for the pair cut-off distance 
	   verbose   = Print info in verbose mode flag */
void LCCountNumNbors( float *boxP, float rcut, int verbose )
{
  Vector_t   scoo2       = { 0., 0., 0. };
  Vector_t   xcoo2       = { 0., 0., 0. };
  float      boxXYZ[3]   = { 0., 0., 0. };
  float      cellVol     = 0.0;
  float      hVect[9]    = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float      hinvVect[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float      minRijSQ    = 9999.0;
  float      maxRijSQ    =    0.0;
  float      rIJ         = 0.0;                          /* pair separation  */
  float      rcutSQ      = 0.0;                          /* Square of cutOff */
  Atom_t     *ap1        = NULL,
             *ap2        = NULL;
  unsigned long long npairs = 0;
  int4       in1, in2    = 0;
  int4       i           = 0;
  int4       ia          = 0;
  int4       ib          = 0;
  int4       nAtomsInThisCell= 0;
  int4       nAtomsInNborCell= 0;
  int4       dir         = 0,
             ni          = 0, 
             nj          = 0,
             nk          = 0;
  int4       nx          = 0, 
             ny          = 0,
             nz          = 0;
  int4       nSlotsX     = 0, 
             nSlotsY     = 0,
             nSlotsZ     = 0;                           /* numSlots in X/Y/Z */
  int4       pairUp      = 0;
  uint4      nPairs      = 0;
  uint4      nAtoms      = 0;

  uint2      minNbors    = 9999;
  uint2      maxNbors    = 0;

  uint4      tblMemAlloc = 0;

  if( !LinkCellCreated ) Error("LCCountNumNbors():LinkCell not yet created");
  
  /* Initialize num nbors to zero */
  nAtoms = LinkCell.Natoms;
  for( ia = 0; ia < nAtoms; ia++ ) AtomP[ia].nNbors = 0;

  /* Assign to local and Invert hVector matrix */
  for( i = 0; i < 9; i++ ) hVect[i] = boxP[i];
  INVERT_MATRIX( hinvVect, hVect, cellVol)

  /* MOTTO: SAVE CPU, AVOID sqrt() => DON'T compute distance => USE rIJ^2 */
  rcutSQ    = rcut * rcut;
  nSlotsX   = LinkCell.NslotsX; 
  nSlotsY   = LinkCell.NslotsY; 
  nSlotsZ   = LinkCell.NslotsZ;

  boxXYZ[0] = LinkCell.BoxSizeX;
  boxXYZ[1] = LinkCell.BoxSizeY;
  boxXYZ[2] = LinkCell.BoxSizeZ;

  if( verbose )   {
    fprintf( stderr,
	    "$$$ LCCountNumNbors(): rcutSQ=%f, Nslots(%d %d %d); BoxSz(%f %f %f)\n",
	    (float) rcutSQ, nSlotsX, nSlotsY, nSlotsZ, 
	    (float) boxXYZ[0], (float) boxXYZ[1], (float) boxXYZ[2] );
    fflush(stderr);
  }

  /* Loop over all the slots */
  for( nx = 0; nx < nSlotsX; nx++ )
    for( ny = 0; ny < nSlotsY; ny++ )
      for( nz = 0; nz < nSlotsZ; nz++ )       {

	/* Number of atoms in this slot */
	nAtomsInThisCell = 
	  LinkCell.LinkCellArray[nx][ny][nz].numAtomsInThisCell;

	/* if( verbose )   {
	   fprintf( stderr,
	   "$$$ LCFormNbrTbl(): nx/y/z(%d %d %d); nAtomsInCell(%d)\n",
	   nx, ny, nz, nAtomsInThisCell );
	   fflush(stderr);
	   } */

	/* If cell is empty; skip to next cell */
	if( !nAtomsInThisCell ) continue;
	
        /* 1. Pair each element in this slot with other atoms 
	      also within this SAME CELL */
	for( ia = 0; ia < nAtomsInThisCell; ia++ )     {

	  /* Array Index of Atom1 */
	  in1 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ia];

	  /* Pointer to atom-1 in pair */
	  ap1 = AtomP + in1;

	  for( ib = ia+1; ib < nAtomsInThisCell; ib++ )    {

	    /* Array Index of Atom1 */
	    in2 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ib];
	    
	    /* Pointer to atom-2 in pair */
	    ap2 = AtomP + in2;

            /* Not a neighbor pair if rIJ > rcut */
            CHECK_IF_NBOR((ap1->coo), (ap2->coo));
	  }  /* for( ib = ia+1; ; )                          */
	}    /* for( ia = 0; ; ): Pair up atoms in same slot */


        /* 2. For slot (nx, ny, nz), consider 13 of its neighboring slots.
              Find pairs (p,q) such that 'p' is in (nx, ny, nz) and 'q'
              is in its neighboring slot (ni, nj, nk).  */
        for( dir = 0; dir < 13; dir++ ) {
          float sx = 0., sy = 0., sz = 0.;            /* Box edge for PBC*/

          ni = nx + Direction[dir][0];
          nj = ny + Direction[dir][1];
          nk = nz + Direction[dir][2];

          /* handle periodic boundary conditions. */
          if( ni < 0 )  {                              /* PBC along X     */
            sx = - boxXYZ[0];                          /* boxXYZ[0] = 1.0 */
            ni = (ni+nSlotsX) % nSlotsX;
          } 
          else if( ni >= nSlotsX )   {
            sx = boxXYZ[0];
            ni = (ni+nSlotsX) % nSlotsX;
          }
          
          if( nj < 0 )  {                              /* PBC along Y     */
            sy = - boxXYZ[1];                          /* boxXYZ[1] = 1.0 */
            nj = (nj+nSlotsY) % nSlotsY;
          }
          else if( nj >= nSlotsY )  {
            sy = boxXYZ[1];
            nj = (nj+nSlotsY) % nSlotsY;
          }
          
          if( nk < 0 )  {                              /* PBC along Z     */
            sz = - boxXYZ[2];                          /* boxXYZ[2] = 1.0 */
            nk = (nk+nSlotsZ) % nSlotsZ;
          } 
          else if( nk >= nSlotsZ ) {
            sz = boxXYZ[2];
            nk = (nk+nSlotsZ) % nSlotsZ;
          }


	  if(ni < 0 || ni >= nSlotsX ||
	     nj < 0 || nj >= nSlotsY ||
	     nk < 0 || nk >= nSlotsZ )
	    Error("$$$ LCCountNumNbors(): This Slot(%d %d %d); Invalid Nbor Slot(%d %d %d)\n", nx, ny, nz, ni, nj, nk);
	  

	  /* Pair each element in this slot (nx,ny,nz) with  atoms 
	     in neighboring slot (ni, nj, nk)  */
	  /* Loop over this slot (nx,ny,nz) */
	  for( ia = 0; ia < nAtomsInThisCell; ia++ )     {

	    /* Array Index of Atom1 */
	    in1 = LinkCell.LinkCellArray[nx][ny][nz].indexInAtomArray[ia];
	    
	    /* Pointer to atom-1 in pair */
	    ap1 = AtomP + in1;
	
	    /* Number of atoms in nbor cell */
	    nAtomsInNborCell =
	      LinkCell.LinkCellArray[ni][nj][nk].numAtomsInThisCell;
	    
	    /* Loop atoms in neighboring slot (ni, nj, nk) */
	    for( ib = 0; ib < nAtomsInNborCell; ib++ )    {
	      
	      /* Array Index of Atom1 */
	      in2 = LinkCell.LinkCellArray[ni][nj][nk].indexInAtomArray[ib];
	      
	      /* Pointer to atom-2 in pair */
	      ap2 = AtomP + in2;

              /* assign ap2 cartesian coord to xcoo */
              ASSIGN_3D_VECTORS( xcoo2, (ap2->coo) );

              /* Compute S-Coord for ap2 */
              MATRIX_3D_VECTOR_MULTIPLY( scoo2, hinvVect, xcoo2 );  

              /* Apply PBC */
              scoo2.x += sx;
              scoo2.y += sy;
              scoo2.z += sz;

              /* Convert to CARTESIAN COORD */
              MATRIX_3D_VECTOR_MULTIPLY( xcoo2, hVect, scoo2 ); 

	      /* Not a neighbor pair if rIJ > rcut */
	      CHECK_IF_NBOR((ap1->coo), (xcoo2));

	    }  /* for( ib = 0; ; )                             */
	  }    /* for( ia = 0; ; ): Pair up atoms in NBOR slot */
        } /* for ( dir = 0; dir < 13; dir++ ) */
      } /* for each slot (i,j,k) -- three nested loops */

  /******* DEBUG *******/
  printf("$$$ LCCountNumNbors(): Natoms=%u\n", LinkCell.Natoms );
  /******* DEBUG *******/

  /* Min/Max num neighbors */
  nAtoms = LinkCell.Natoms;
  for( ia = 0; ia < nAtoms; ia++ )   {
    uint2 nnnn;
    ap1 = AtomP + ia;

    nnnn       = ap1->nNbors;
    ap1->pairs = nnnn;
    ap1->nNbors= 0;

    ap1->nbor     = (uint4*) calloc( (size_t) nnnn, (size_t) sizeof(uint4) );

#ifdef DO_NOT_USE
    ap1->index = (uint2**) calloc( (size_t) nnnn, (size_t) sizeof(uint2*) );

    for( i = 0; i < nnnn; i++ )
      ap1->index[i] = (uint2*) calloc( (size_t) 3, (size_t) sizeof(uint2) );
    tblMemAlloc += 3*( ((size_t) nnnn) * ((size_t) sizeof(uint2)) );
#endif

    tblMemAlloc += ( ((size_t) nnnn) * ((size_t) sizeof(uint4)) );

    minNbors   = MIN( nnnn, minNbors );
    maxNbors   = MAX( nnnn, maxNbors );
  }


  fprintf( stderr, "$$$ LCCountNumNbors(): Allocate %u MB for NborTable\n",
	  (tblMemAlloc/1000000 + 1) );

  /* Minimum and Maximum nbors/separation */
  LinkCell.MinRij   = sqrt( (double) minRijSQ );
  LinkCell.MaxRij   = sqrt( (double) maxRijSQ );
  LinkCell.MinNbors = (int4) minNbors;
  LinkCell.MaxNbors = (int4) maxNbors;

  if( verbose )   {
    fprintf( stderr,
	    "$$$ LCCountNumNbors(): min/maxNbors=(%d %d), Nslots(%d %d %d)\n",
	    (int) minNbors, (int) maxNbors, nSlotsX, nSlotsY, nSlotsZ );

    fprintf( stderr,
	    "$$$ LCCountNumNbors(): rcutSQ/min/maxRijSQ=(%f %f %f), BoxSz(%f %f %f)\n",
	    (float) rcutSQ, (float) minRijSQ, (float) maxRijSQ, 
	    (float) boxXYZ[0], (float)  boxXYZ[1], (float)  boxXYZ[2] );
  }

}   /* LCCountNumNbors()  */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef DBG
/*---------------------------------------------------------------------------*/
/* CARTESIAN COORDINATE ARGUMENTS */
#define FORM_PAIR( xcoo1, xcoo2, in1, in2 )                                   \
{                                                                             \
  float dX, dY, dZ;                                                           \
									      \
  dX  = xcoo1.x - xcoo2.x; dY = xcoo1.y - xcoo2.y; dZ = xcoo1.z - xcoo2.z;    \
  rIJ = dX * dX + dY * dY + dZ * dZ;                                          \
									      \
  pairUp = 1;                                  /* Default is pairup */        \
									      \
  if( rIJ > rcutSQ ) pairUp = 0;                                              \
									      \
  if( pairUp )      {                                                         \
   /* Fill neighbor Table of Atom_t; increment num neighbors field            \
      num neighbors must start from zero as array index starts from 0.        \
      => post increment num neighbors field; Save pair-separation */          \
   if((ap1->nNbors+1 >= NNBR) || (ap2->nNbors+1 >= NNBR) )  {                 \
      printf("FORM_PAIR(): Too many nbors; PairID(%u %u), #nNbors(%d %d)\n",  \
	     ap1->id, ap2->id, (int4)ap1->nNbors, (int4)ap2->nNbors );        \
      Error("FORM_PAIR(): Too many neighbors; nbor array overflow\n");        \
   }                                                                          \
                							      \
   /* Instead of saving,pointers to atoms we store the INDEX                  \
      of their position in the AtomP array */                                 \
   ap1->nbor[ (ap1->nNbors) ] = in2;    ap2->nbor[ (ap2->nNbors) ] = in1;     \
									      \
   ap1->nNbors++;   ap2->nNbors++;    /*nPairs++;   */                        \
									      \
   minRijSQ = MIN( rIJ, minRijSQ );                                           \
   maxRijSQ = MAX( rIJ, maxRijSQ );                                           \
									      \
  }                                                                           \
}
/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
