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
 * change(s). Do NOT USE or DISTRIBUTE is without our written consent.       *
 *****************************************************************************
 * Details of CNA is in:                                                     *
 *          1. Clarke and Jonsson, PRE 47, 3975 (1991)                       *
 *          2. Faken and Jonsson,  Comp. Mater. Sci., (1994)                 *
 *****************************************************************************
 * $Author: sgsrini $ : Srinivasan G. Srivilliputhur (T-11, LANL)
 * $Date: 2000/01/28 19:58:36 $
 * $Source: /g12/sgsrini/MPI_MD98/CNcodes/CN2000/RCS/CNmain.c,v $
 * $Revision: 1.5 $
 * $Log: CNmain.c,v $
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
 * Usage  : "cnA.exc -R rcut -I [infile] -O [outfile] -v"                    *
 * Options:                                                                  *
 *          -C  X/Y/Z:MinCoo:MaxCoo        (Carve a region based on X/Y/Z)   *
 *          -R	cut-off for neighbor atoms (units: Angstrom).                *
 *          -v 	verbose mode                                                 *
 *          -I  infile                                                       *
 *          -F  infileformat   = (lass, spasm1, spasm2, baskes,adsnold)              *
 *          -S  NATOMS:A:B:C   = For SPaSM only
 *          -O  output atom file name                                        *
 *****************************************************************************
 * Description: The program reads a configuration from a atom-file and       *
 *              calculates for each 'bond' (i.e., pair of atoms at distance  *
 *              < rcut) three indices as described in above REFERENCES.      *
 *              An outfile is written with essentially same info as infile.  *
 *              For each atom, the Outfile in addition contains, numberof    *
 *              neighbors and the 3 index label for each bond.               *
 *****************************************************************************/
/*---------------------------------------------------------------------------*/
/* $$$ main.c */
/*---------------------------------------------------------------------------*/
#include <stdarg.h>                        /* variable arg handling */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>                        /* size_t                */
#include <string.h>
/*---------------------------------------------------------------------------*/
#include "CNmain.h"
#include "CNpairStat.h"
/*---------------------------------------------------------------------------*/
extern int strcasecmp(const char *, const char *);
/*---------------------------------------------------------------------------*/
/****** extern global function declarations *****************************/
extern void LCInitialize( int4 verbose );
extern void LCFormNbrTable( float *boxP, float rcut, int4 verbose );
extern void LCCountNumNbors( float *boxP, float rcut, int verbose );
extern int4 ReadConLASS(char*, float *boxP, float rcut, int4 verbose );
extern int4 ReadConBaskes(char*, float *boxP, float rcut, int4 verbose );
extern int4 ReadConOldAdhara(char*, float *boxP, float rcut, int4 verbose );
extern void ReadConSPaSM(char*, float*, float, uint4, int4,  
			 char*, float, float, int4 );
extern int4 CNAmainLoop(uint4 nAtoms, int4 verbose);
extern void WriteConLASS( char *outFil, int4 verbose );
extern void WriteConBaskes(char *outFil,float *boxP,int4 nAtoms,int4 verbose);
extern void WritePairsInfo( CNstat_t*, char*, float*, float, uint4, int4 );
extern void InitializeCNvariables( CNstat_t* );
extern void ComputeCoordination(CNstat_t*, uint4, int4 );
extern void CountBccFccHcpIcoAtoms( CNstat_t*, uint4, int4 );

extern void WriteConVASP( char*, int4, int4 );
extern int4 ReadConVASP(char*, float *boxP, float rcut, int4 verbose );

/*---------------------------------------------------------------------------*/
/* Global Variables here */
Atom_t           *AtomP = NULL;                  /* Global Atom Array        */
static CNstat_t  CNstatTbl;                      /* CN statistics saved here */
/*---------------------------------------------------------------------------*/
char             InputFileName[BUF_SIZE];
int4             InpFileFormatFlag = MBASKES_FORMAT;
/*---------------------------------------------------------------------------*/
void Error(char *msg, ...)
{
  va_list args;
  char    buf1[1024];

  /* Check if error message fits in buf[]; If the array bounds are not 
     checked there may be overflow and coredump.   */
  if( (strlen(msg) ) >=  1024 )
    strcat( buf1," Error Buffer OverFlow. Can't EVEN Print Correct Msg\0" );

  /* Flush all buffers */
  fflush( stderr );
  
  /* handle variable arguments */
  va_start(args, msg);
  vsprintf(buf1, msg, args);
  va_end(args);

  /* Print the message */
  fprintf(stderr, "\n***********\n%s\n***********\n", buf1 );
  fflush( stderr );

  exit( -1 );
} /* void Error(char *msg, ...) */
/*--------------------------------------------------------------------------*/
/******************************************************
 * main(): Main function for Common Neighbor Analysis *
 ******************************************************/
int main(int argc, char *argv[])
{
  float   tmpFlot   = 0.0;
  float   minCoo    = -1.0;
  float   maxCoo    = -1.0;
  float   simBox[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  float   rcut      = 5.0;            /* Default: 5.0                        */
  char    ifformat[BUF_SIZE];         /* lass/spasm1/spasm2  formats         */
  char    ifname[BUF_SIZE];
  char    ofname[BUF_SIZE];
  char    *coordToCarveSubRegion = {"N"};/* Don't CARVE              */
  FILE    *ifp      = NULL;           /* input  file ptr                     */
  FILE    *ofp      = NULL;           /* output file ptr                     */
  int4    conFil    = LASSFIL_FORMAT; /* Default is LASS file                */
  uint4   nAtoms    = 0;
  int4    ifFlag    = 0;              /* File read flag                      */
  int4    ofFlag    = 0;
  int4    rcutFlag  = 0;
  int4    spasmFlag = 0;
  int4    verbose   = 0;              /* Default: NOT verbose mode           */
  int4    c;

  extern char *optarg;                /* For getopt()                        */
  extern int4 optind;                 /* For getopt()                        */

  /* Handle commanline args */
  while( (c = getopt( argc, argv, "C:F:I:O:R:S:v" )) != EOF )
    switch (c) {
    case 'v':                         /* Verbose mode                        */
      verbose  = 1;
      break;

    case 'F':                         /* Input  File Format                  */

      /* File format */
      sprintf( ifformat, "%s", optarg );


      fprintf(stderr, "conFil(%d) fileNm(%s)\n",
	      conFil, ifformat );

      /* Atom file format */
      conFil = LASSFIL_FORMAT;
      if( strcasecmp(ifformat, "lass") == 0 )         
	InpFileFormatFlag = conFil = LASSFIL_FORMAT;
      else if( strcasecmp(ifformat, "spasm1") == 0 )  
	InpFileFormatFlag = conFil = SPASM_1_FORMAT;
      else if( strcasecmp(ifformat, "spasm2") == 0 )  
	InpFileFormatFlag = conFil = SPASM_2_FORMAT;
      else if( strcasecmp(ifformat, "baskes") == 0 )  
	InpFileFormatFlag = conFil = MBASKES_FORMAT;
      else if( strcasecmp(ifformat, "vasp"  ) == 0 )  
	InpFileFormatFlag = conFil = VASP_FORMAT;
      else if( strcasecmp(ifformat, "adoldsn" )== 0 ) {
	InpFileFormatFlag = conFil = ADHARA_OLD_SNAP_FORMAT;
	fprintf(stderr, "conFil(%d) fileNm(%s)\n",
		conFil, ifformat );
      }
      else
	Error("main(1): %s is unknown file format\n", ifformat );

      break;

    case 'I':                         /* Input  File Name                    */
      /* In file name */
      sprintf( ifname,"%s", optarg );
      sprintf( InputFileName, "%s", ifname );

      ifFlag++;                       /* Input file read                     */
      break;

    case 'O':                         /* Output File Name                    */
      sprintf( ofname,"%s", optarg );
      ofFlag++;                       /* output file read                    */
      break;

    case 'R':                         /* rcut                                */
      rcut     = (float) atof(optarg);
      rcutFlag++;
      break;

    case 'S':                         /* Box Vector and Atoms in SPaSM file  */
      /* Read number of atoms  in SPaSM atom file */
      nAtoms   = (uint4) atoi(optarg);
      spasmFlag++;

      /* SPaSM uses an orthogonal box => get Ax/By/Cz box vectors */
      while (*optarg++ != ':');
      simBox[0]     = (float) atof(optarg);

      while (*optarg++ != ':');
      simBox[4]     = (float) atof(optarg);

      while (*optarg++ != ':');
      simBox[8]     = (float) atof(optarg);
      break;

    case 'C':                         /* Carve out a sub region based on 
					 min/max values of either x/y/z      */
      /* X/Y/Z coord to carve out region */
      sprintf( coordToCarveSubRegion, "%s", optarg );

      /* Read Min/Max value of region's X/Y/Z coord */
      while (*optarg++ != ':');
      minCoo   = (float) atof(optarg);

      /* SPaSM uses an orthogonal box => get Ax/By/Cz box vectors */
      while (*optarg++ != ':');
      maxCoo   = (float) atof(optarg);

      if( maxCoo < minCoo )    {
	tmpFlot = maxCoo;
	maxCoo  = minCoo;
	minCoo  = tmpFlot;
      }

      switch( coordToCarveSubRegion[0] )   {
      case 'N':
      case 'X':
      case 'x':
      case 'Y':
      case 'y':
      case 'Z':
      case 'z':
	break;
      default:
	Error("main(): %s is Unknown option for -C\n", coordToCarveSubRegion);
      }
      break;

    default:
      (void) fprintf(stderr, "## ERROR in USAGE\n");
      Error("main(): Usage: '%s -C X/Y/Z:minCoo:maxCoo -R rcut -I infile -F infileformat -O outFil -S nAtoms:Ax:By:Cz -v'\n", argv[0]);
    } /* while{switch()} */


  if( ifFlag   != 1  )    Error("main(): Input  file NOT read correctly\n");
  if( ofFlag   != 1  )    Error("main(): Output file NOT read correctly\n");
  if( rcutFlag != 1  )    Error("main(): Rcut NOT read correctly\n");

  if( conFil == SPASM_1_FORMAT && conFil != SPASM_2_FORMAT )
    /* if( spasmFlag != 1 ) */
      Error("main(): Provide num atoms and box vectors for SPaSM file\n");

  if( verbose )
    (void) fprintf(stderr, "Cut-off for neighbours is %.4lf Angstrom\n", rcut);

  if(verbose) fprintf(stderr, "1. Enter LCInitialize()\n");

  /* Initialize CNstat variables */
  InitializeCNvariables( &CNstatTbl );

  /* Initialize Link-Cell Variables */
  LCInitialize( verbose );
  
  if(verbose) fprintf(stderr, "2. Enter ReadCon()\n");

  /* Read Input File */
  switch( conFil )   {
  case LASSFIL_FORMAT:
    /* Read Atom Configuration File */
    nAtoms = ReadConLASS( ifname, simBox, rcut, verbose );
    break;

  case ADHARA_OLD_SNAP_FORMAT:
    /* Read Atom Configuration File */
    nAtoms = ReadConOldAdhara( ifname, simBox, rcut, verbose );
    break;

  case VASP_FORMAT:
    /* Read Atom Configuration File */
    nAtoms = ReadConVASP( ifname, simBox, rcut, verbose );
    break;

  case SPASM_1_FORMAT:
  case SPASM_2_FORMAT:
    ReadConSPaSM(ifname, simBox, rcut, nAtoms, conFil, 
		 coordToCarveSubRegion, minCoo, maxCoo, verbose );
    break;
  case MBASKES_FORMAT:
    /* Read Atom Configuration File */
    nAtoms = ReadConBaskes( ifname, simBox, rcut, verbose );
    break;

  default:
    Error("main(): Unknown file format\n");
  }

  if( verbose )   {
    fprintf(stderr, 
	    "$$$ main(): Hvec\n\t\t(%lf %lf %lf)\n\t\t(%lf %lf %lf)\n\t\t(%lf %lf %lf)\n",
	    simBox[0],simBox[1],simBox[2],simBox[3],simBox[4],simBox[5],
	    simBox[6],simBox[7],simBox[8]);
    fflush( stderr );
  }

  /* Count neighbors to allocate table memory */
  LCCountNumNbors( simBox, rcut, verbose );

  /* Form neighbor table */
  if(verbose) fprintf(stderr, "3. Enter LCFormNbrTable()\n");
  LCFormNbrTable( simBox, rcut, verbose );

  /* Perform common neighbor analysis */
  if(verbose) fprintf(stderr, "4. Enter CNAmainLoop()\n");
  CNAmainLoop( nAtoms, verbose );

  /* Compute Number of FCC, BCC, HCP atoms */
  if(verbose) fprintf(stderr, "5. Enter CountBccFccHcpIcoAtoms()\n");
  CountBccFccHcpIcoAtoms( &CNstatTbl, nAtoms, verbose );

  /* Compute Average Coordination Number */
  if(verbose) fprintf(stderr, "6. Enter ComputeCoordination()\n");
  ComputeCoordination( &CNstatTbl, nAtoms, verbose );

  /* Write out CNN pairs related info */
  if(verbose) fprintf(stderr, "7. Enter WritePairsInfo()\n");
  WritePairsInfo( &CNstatTbl, InputFileName, simBox, rcut, nAtoms, verbose );

  /* Write out output configuration */
  if(verbose) fprintf(stderr, "8. Enter WriteConfig()\n");

  switch( conFil )   {
  case LASSFIL_FORMAT:
  case SPASM_1_FORMAT:
  case SPASM_2_FORMAT:
    WriteConLASS( ofname, verbose );
    break;
  case VASP_FORMAT:
    WriteConVASP( ofname, nAtoms, verbose );
    break;
  case ADHARA_OLD_SNAP_FORMAT:
  case MBASKES_FORMAT:
    WriteConBaskes(ofname, simBox, nAtoms, verbose );
    break;
  default:
    Error("main(): Unknown file format\n");
  }

  return(EXIT_SUCCESS);
}  /* main() */
/*---------------------------------------------------------------------------*/
