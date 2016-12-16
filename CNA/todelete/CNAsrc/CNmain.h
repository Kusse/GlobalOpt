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
 * $Source: /home/srini/LJ_EAM/CN2000/RCS/CNmain.h,v $
 * $Revision: 1.1 $
 * $Log: CNmain.h,v $
 * Revision 1.1  2002/08/13 20:36:38  sgsrini
 * Initial revision
 *
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
 * Revision 1.2  2000/01/27 04:17:09  sgsrini
 * Changed CN index table to ushort
 *
 * Revision 1.1  2000/01/27 00:30:58  sgsrini
 * Initial revision
 *
 * Revision 1.1  2000/01/27 00:25:31  sgsrini
 * Initial revision
 *****************************************************************************
 *                  $$$$ CNmain.h: Defs for CNA code                         *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
typedef char           *string;
typedef unsigned char  Uchar;
/*---------------------------------------------------------------------------*/
typedef short          int2;
typedef unsigned short uint2;
/*---------------------------------------------------------------------------*/
typedef int            int4;
typedef unsigned int   uint4;
/*---------------------------------------------------------------------------*/
typedef long           int8;
typedef unsigned long long uint8;
/*---------------------------------------------------------------------------*/
typedef float          Real4;
typedef double         Real8;
/*---------------------------------------------------------------------------*/
#define MAGICNUMBER     0x1743
/*---------------------------------------------------------------------------*/
#define NO              0
#define YES             1
/*---------------------------------------------------------------------------*/
#define LASSFIL_FORMAT         0
#define SPASM_1_FORMAT         1
#define SPASM_2_FORMAT         2
#define MBASKES_FORMAT         3
#define VASP_FORMAT            4
#define ADHARA_OLD_SNAP_FORMAT 5
/*---------------------------------------------------------------------------*/
/* Atom environment: BCC, FCC, HCP, ICOSAHEDRAL, DISTORTED ICOSAHEDRA (DIC).
   We will use an atoms coordination number as a tag if an atom does NOT
   belong to one of the above REGULAR POLYHEDRA  */
#define BCC                     30
#define HCP                     35
#define ICO                     40
#define DIC                     45
#define FCC                     50
/*---------------------------------------------------------------------------*/
/********* constants *********************************************************/
#define	PI	           3.14159265
#define NCHR	           1024	   /* max number of characters in a record   */
#define NNBR	             26	   /* max number of neighbors                */
#define NBN2	             94	   /* 2 times max number bonds between nbors */
#define NDAT	              0	   /* number of extra data columns in output */
/*---------------------------------------------------------------------------*/
#define BUF_SIZE            512
/* #define MAX_ATOMS_IN_SLOT   120 */
#define MAX_ATOMS_IN_SLOT   60
/*---------------------------------------------------------------------------*/
#define MAX_CELLS_X        110              /* Max linkcells along X/Y/Z     */
#define MAX_CELLS_Y        110
#define MAX_CELLS_Z        250
/*---------------------------------------------------------------------------*/
/* Max size of box along X, Y, Z */
#define BOX_SIZE_X         1.0
#define BOX_SIZE_Y         1.0
#define BOX_SIZE_Z         1.0
/*---------------------------------------------------------------------------*/
/******** data types *********************************************************/
/*---------------------------------------------------------------------------*/
/* Structure with info on atoms in each link-cell. indexInAtomArray[] has
   index of AtomP[] array and hence has atom ID/Coord info of the atom in this
   link-cell. indexInAtomArray[] is accessed by looping through 
   numAtomsInThisCell */
typedef struct slots  {
  int4 numAtomsInThisCell;
  int4 indexInAtomArray[MAX_ATOMS_IN_SLOT];
} SlotType;
/*---------------------------------------------------------------------------*/
/* LinkCell Framework */
typedef struct lcn {
  float         MinRij;             /* Minimum neighbor separation          */
  float         MaxRij;             /* Minimum neighbor separation          */

  float         BoxSizeX;           /* size of the simulation box.          */
  float         BoxSizeY;
  float         BoxSizeZ;

  float         SlotSizeX;          /* size of slot.                        */
  float         SlotSizeY;
  float         SlotSizeZ;

  float         SlotSizeInvX;       /* 1/(size of slot)                     */
  float         SlotSizeInvY;
  float         SlotSizeInvZ;
                                    /* LinkCell Slots along X/Y/Z           */
  /* SlotType      LinkCellArray[MAX_CELLS_X][MAX_CELLS_Y][MAX_CELLS_Z];
     SlotType      *LinkCellArray; */

  SlotType      ***LinkCellArray;

  uint4         Natoms;             /* Number of Atoms                      */
  int4          NslotsX;            /* current no of slots in each dim.     */
  int4          NslotsY;
  int4          NslotsZ;

  int4          MinNbors;           /* Atoms with Minimum number of nbors   */
  int4          MaxNbors;           /* Atoms with Maximum number of nbors   */

} LCType;
/*---------------------------------------------------------------------------*/
typedef struct vec3 {
  float x, y, z;
} Vector_t;
/*---------------------------------------------------------------------------*/
/* Atom information */
typedef struct atom   {
  Vector_t      coo;             /* coordinate: S-system or Cartesian        */
  uint4         id;

#ifdef STATIC_ALLOC
  uint4         nbor[NNBR];      /* Neighbor ID (index in the AtomP[] array) */

  uint2         atomicNum;
  uint2         cnTag;           /* BCC, FCC, HCP, ICO, DIC, or coord-num    */
  uint2         nNbors;          /* Number of neighbors in nbor[] table      */
  uint2		pairs;
  uint2		index[NNBR][3];  /* CNA (l,m,n) index                        */

#endif

  uint4         *nbor;           /* Neighbor ID (index in the AtomP[] array) */

  uint2         nNbors;          /* Number of neighbors in nbor[] table      */
  uint2		pairs;
  Uchar         atomicNum;       /* For Baskes Dynamo format, atomtype is
				    stored here                              */
  Uchar         cnTag;           /* BCC, FCC, HCP, ICO, DIC, or coord-num    */
  Uchar		index[NNBR][3];  /* CNA (l,m,n) index                        */
  char          tags[4];         /* [0] = Fix('F')/Move('M')
				    [1] = Blue('B')/Red ('R')/Green('G')
				    [2] = 'I' if atom part of an ICOSAHEDRA
				          'D' if atom in DISTORTED ICOSAHEDRA*/
#ifdef SMALL_SYSTEM
  float         potE;            /* Atom Potential Energy                    */
#endif

#ifdef NOT_NEEDED_NOW
  uint4         indexInAtomPArray;
#endif

} Atom_t;
/*---------------------------------------------------------------------------*/
