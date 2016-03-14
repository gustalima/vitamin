/*
 * Neighbouring Contact Locator
 *
 * This program is designed to scan Brookhaven Protein Databank files for
 * neighbouring atom contacts.
 */

/* Version 0.9     by David T. Jones, August 1990 */

/* Version 0.9a    by Dorica Naylor, May 1991
                   modified to take account of closest neighbours on a
                   per residue basis. Currently only intra-chain contacts
                   are reported, and for any given residue-residue pair only
                   the closest atom-atom contact is reported. */

/* Version 0.9b    by Dorica Naylor, 9th August 1991  */

/* Version 0.9c    by Dorica Naylor, 28th August 1991
                   modified to change the way CA coordinates are handled
                   and to accommodate insertion codes.  */

/* Version 0.9d    by Dorica Naylor, 25th October 1991
                   modified to take account of insertion residues and
                   missing residues in calculating nbgap. */

/* Version 0.9e    by Dorica Naylor, 19th Noveber 1991
                   fixed bug that caused meaningless values for GAP
                   in proteins with NO chain breaks. */

/* Version 0.9f    by Dorica Naylor, 13th March 1992
                   modified to take account of non-standard amino acids
                   and ligand molecules. */

/* Version 0.9g    by Dorica Naylor, 7th April 1992
                   modified to take account of non-standard amino acids,
                   missing atoms (non-standard chain brakes), and ligand
                   molecules. */

/* Version 0.9h    by Dorica Naylor, 5th May 1992
                   MAXCNCOLS changed from 11 to 5-reads covalent bonds 
                   only */

/* Compile with   "cc  -o nb  nb.c  /usr/lib/libm.a"    */

/* Amendments made for the PROCHECK suite are identified by CHECK v.m.n,
   where m.n is the version number of the PROCHECK */

/* v.2.0.3  -     Amendments to correct two minor bugs:
                     1. Last residue always excluded from nearest-neighbour 
                        checks
                     2. For multimeric proteins, with chains identified by
                        chain identifiers, HETATM with a blank chain
                        identifier were not being compared against the
                        protein atoms - simply because the chain identifiers
                        were not equal in this case
                                            Dorica Naylor 21 October 1992 */

/* v.2.1    -     Removal of Ctrl-G characters in Error messages
                                          Roman Laskowski 13 January 1993 */

/* v.2.1.3  -     Addition of asterisks to error prints so that can be plucked
                  out of the log files and displayed by the script file.
                                            Roman Laskowski 26 March 1993 */

/* v.2.1.4  -     Amendment of '\/' to '/' which is unnecessary and causes
                  problems with certain compilers
                                              Roman Laskowski 12 May 1993 */

/* v.2.3   -      Tiny change so that VAX .COM routines don't show error if
                  no asterisks found in log file
                                              Roman Laskowski 14 Nov 1993 */

/* v.3.2   -      Amendments suggested by Dave Love of Daresbury.
                                  Dave Love / Roman Laskowski 18 Oct 1994 */

/* v.3.4   -      Amendments suggested by Pete Dunten to fix the problem
                  of the last peptide bond being picked up as a bad contact
                  whenever no HET groups or other residues follow.
                                  Dave Love / Roman Laskowski 18 Oct 1994 */

/* v.3.6.2  -     Bug-fixes to deal with very large files, or files with  
                  large number of waters/ligands.
                                              Roman Laskowski 16 Feb 2004 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
/* CHECK v.3.2 */
#include <stdlib.h>
/* CHECK v.3.2 */

#define max(x,y) ((x)>(y) ? (x) : (y))
#define min(x,y) ((x)<(y) ? (x) : (y))

#define FALSE              0
#define TRUE               1

#ifdef __MSDOS__
#define MAXNATM 1000
#define MAXNRES  150
#else
/* CHECK v.3.6.2 */
/* #define MAXNATM 50000
#define MAXNRES 10000
#define MAXBRKS  1000 */
#define MAXNATM 100000
#define MAXNRES  10000
#define MAXBRKS   5000
/* CHECK v.3.6.2 */
#endif

#define STDAA       20   /* no of standard amino acids recognized by neighbour*/
#define NONSTDAA    35   /* no of non-standard amino acids recognized by nb   */
#define TOTNAA      55   /* total no of amino acids recognized by neighbour   */
#define CHNBRKFLAG   4   /* mode for determination of chain breaks (1-4)      */
#define MAXCON       8   /* max no of bonds allowed for any 1 atom            */
#define MAXCNRECS 5000   /* max no of CONECT records allowed in a pdbfile     */
#define MAXCNCOLS    5   /* max no of columns of data in a CONECT record      */
#define IGNLB        1   /* ignore bonds between ligands and proteins?        */

#define SQR(X)  ((X)*(X))

#define NTOKENS            25

/* A list of common PDB record types... */

#define HEADER             1
#define COMPND             2
#define SOURCE             3
#define AUTHOR             4
#define REVDAT             5
#define REMARK             6
#define SEQRES             7
#define FTNOTE             8
#define HET                9
#define FORMUL             10
#define HELIX              11
#define CRYST1             12
#define ORIGX1             13
#define ORIGX2             14
#define ORIGX3             15
#define SCALE1             16
#define SCALE2             17
#define SCALE3             18
#define ATOM               19
#define TER                20
#define HETATM             21
#define CONECT             22
#define ENDENT             23
#define JRNL               24
#define TURN               25

struct pdbatm
{
    float       x, y, z;
    int         aanum, aacode, atmtyp, caindex, atmnum, ncon;
    char        altcode, inscode, chnid, hetflg, ssflg, strucsum,
                resnam[4], atmnam[5];
}   atom[MAXNATM];

float   caxyz[MAXNRES][3], cxyz[MAXNRES][3], nxyz[MAXNRES][3];

short   cafound[MAXNRES], brake[MAXBRKS][2], nfound[MAXNRES], firstatom[MAXNRES],
        cfound[MAXNRES];

char    csdfn[160], logfn[160], keyword[40], buf[160];
char    brcode[5], tablefn[160];

/* Record names for decoding record types */
char   *tokstr[25] =
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET",    "FORMUL",
 "HELIX",  "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM",   "TER",
 "HETATM", "CONECT", "END",    "JRNL",   "TURN"
};

/* Residue name to allow conversion of a.a. name into numeric code */
char   *rnames[TOTNAA] =
{
 "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
 "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
 "AIB", "PHL", "SEC", "ALM", "MPR", "FRD", "LYM", "GLM", "PPH", "PGL",
 "OLE", "ABA", "NLE", "B2V", "B2I", "B1F", "BNO", "B2A", "B2F", "IVA",
 "LOV", "STA", "PVL", "CAL", "PHA", "DCI", "AHS", "CHS", "MSE", "ETA",
 "PCA", "ASX", "GLX", "INI", "LOL"
};

enum RESCODE
{
    Ala, Cys, Asp, Glu, Phe, Gly, His, Ile, Lys, Leu,
    Met, Asn, Pro, Gln, Arg, Ser, Thr, Val, Trp, Tyr,
    Aib, Phl, Sec, Alm, Mpr, Frd, Lym, Glm, Pph, Pgl,
    Ole, Aba, Nle, B2v, B2i, B1f, Bno, B2a, B2f, Iva,
    Lov, Sta, Pvl, Cal, Pha, Dci, Ahs, Chs, Mse, Eta,
    Pca, Asx, Glx, Ini, Lol
};

/* Bonding Atom Dictionary - this is the list of atom names, arranged on a
   residue by residue basis in RESCODE order, of those atoms capable of
   forming neighbouring contacts (i.e. all atoms!). Only atoms in the input
   file (.NEW or .PDB) whose names match those in NECATM will be extracted
   from the file.
*/
char *necatm[TOTNAA] =
{
      " N   CA  C   O   CB  OXT",
      " N   CA  C   O   CB  SG  OXT",
      " N   CA  C   O   CB  CG  OD1 OD2 OXT",
      " N   CA  C   O   CB  CG  CD  OE1 OE2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   C   OXT",
      " N   CA  C   O   CB  CG  ND1 CD2 CE1 NE2 OXT",
      " N   CA  C   O   CB  CG1 CG2 CD1 OXT",
      " N   CA  C   O   CB  CG  CD  CE  NZ  OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 OXT",
      " N   CA  C   O   CB  CG  SD  CE  OXT",
      " N   CA  C   O   CB  CG  OD1 ND2 OXT",
      " N   CA  C   O   CB  CG  CD  OXT",
      " N   CA  C   O   CB  CG  CD  OE1 NE2 OXT",
      " N   CA  C   O   CB  CG  CD  NE  CZ  NH1 NH2 OXT",
      " N   CA  C   O   CB  OG  OXT",
      " N   CA  C   O   CB  OG1 CG2 OXT",
      " N   CA  C   O   CB  CG1 CG2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OH  OXT",
      " N   CA  C   O   CB1 CB2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   CB SEG  OD1 OD2 OXT",
      " N   CA  C   O   CB  CM  OXT",
      " ??? CA  C   O   CB  SG ",
      " N   CA  C   ??? CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   CB  CG  CD  CE  NZ  CM  OXT",
      " N   CA  C   O   CM  OXT",
      " N   CA  P   OP1 OP2 CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  P   O1P O2P OXT",
      " ON  CA  C   O   CB  CG  CD1 CD2 OXT",
      " N   CA  C   O   CB  CG  OXT",
      " N   CA  C   O   CB  CG  CD  CE  OXT",
      " N   CA  B   O1  O2  CB  CG1 CG2 OXT",
      " N   CA  B   O1  O2  CB  CG1 CG2 CD1 OXT",
      " N   CA  B   O1  O2  CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  B   O1  O2  CB  CG  CD  CE  OXT",
      " N   CA  B   O1  O2  CB  OXT",
      " N   CA  B   O1  O2  CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " ??? CA  C   O   CB  CG1 CG2 OXT",
      " N   CA  C   O   CB  CG1 CG2 CD1 CD2 C1G C1B C1A CS  OS  CT  OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CH  OH  CM  OXT",
      " ON  CA  C   O   CB ",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  CA2 CB2 CG2 CD3 CD4 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  ??? ??? CB  CG1 CG2 CD1",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  N1  CB2 CG2 CD3 CD4 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  OXT",
      " N   CA  C   O   CB  CG SED  CE  OXT",
      " N   CA  ??? O   CB ",
      " N   CA  C   O   CB  CG  CD  OE  OXT",
      " N   CA  C   O   CB  CG  AD1 AD2 OXT",
      " N   CA  C   O   CB  CG  CD  AE1 AE2 OXT",
      " N   CA  C   O   CB  CG  CD  CE  CZ  CI1 CI2 CI3 CI4 NI2 CI5 CI6 OXT",
      " N   CA  CH  OH  CB  CG  CD1 CD2 OXT"
};

/* NBONDS contains the number of bonds contained by each of the STDAA amino
   acids, arranged in RNAMES/RESCODE order.
   Each subsequent array (e.g.alabonds) contains the names of the two atoms
   that make up each of the bonds for that particular residue.
*/

int   nbonds[STDAA]={5,6,8,9,12,4,11,8,9,8,8,8,8,9,11,6,7,7,16,13};

char *alabonds[5]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB " };

char *cysbonds[6]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  SG " };

char *aspbonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  OD1", " CG  OD2" };

char *glubonds[9]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD ", " CD  OE1",
                      " CD  OE2" };

char *phebonds[12] = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD1", " CG  CD2",
                      " CD1 CE1", " CD2 CE2", " CE1 CZ ", " CE2 CZ " };

char *glybonds[4]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT" };

char *hisbonds[11] = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  ND1", " CG  CD2",
                      " ND1 CE1", " CD2 NE2", " CE1 NE2" };

char *ilebonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG1", " CB  CG2", " CG1 CD1" };

char *lysbonds[9]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD ", " CD  CE ",
                      " CE  NZ " };

char *leubonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD1", " CG  CD2" };

char *metbonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  SD ", " SD  CE " };

char *asnbonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  OD1", " CG  ND2" };

char *probonds[8]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD ", " CD  N  " };

char *glnbonds[9]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD ", " CD  OE1",
                      " CD  NE2" };

char *argbonds[11] = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD ", " CD  NE ",
                      " NE  CZ ", " CZ  NH1", " CZ  NH2" };

char *serbonds[6]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  OG " };

char *thrbonds[7]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  OG1", " CB  CG2" };

char *valbonds[7]  = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG1", " CB  CG2" };

char *trpbonds[16] = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD1", " CG  CD2",
                      " CD1 NE1", " NE1 CE3", " CD2 CE3", " CD2 CE2",
                      " CE3 CZ3", " CE2 CZ2", " CZ2 CH2", " CZ3 CH2" };

char *tyrbonds[13] = {" N   CA ", " CA  C  ", " C   O  ", " C   OXT",
                      " CA  CB ", " CB  CG ", " CG  CD1", " CG  CD2",
                      " CD1 CE1", " CD2 CE2", " CE1 CZ ", " CE2 CZ ",
                      " CZ  OH " };

#define  NBDIST    4.0
#define  PEPBND    2.5
#define  SSDIST    3.0
#define  CAWARN   12.0
#define  CADBND    5.0

FILE    *ofp;
FILE    *ifp;
FILE    *dbgfp;

int      natoms, debug;
int      strtat[MAXNRES], stopat[MAXNRES], nres, nbrakes;
int      icon[MAXNATM][MAXCON], bonds[MAXCNRECS][MAXCNCOLS], nconrecs;


/*============================================================================*/

/* CHECK v.2.0 */
/* getnam  -  Peel off the directory path and extension from the full name
              of the .pdb file */

int getnam(pdbfil,istart,iend)
char  *pdbfil;
int   *istart, *iend;
{
     char pchar;
     int  finish, gotdot, ipos, istate;


/*-- Initialise variables --*/
     finish = FALSE;
     *iend = -1;
     *istart = 0;
     istate = 1;
     ipos = strlen(pdbfil) - 1;
     gotdot = FALSE;

/*-- Check through the filename from right to left --*/
     while (finish == FALSE && ipos > -1)
     {

/*--     Pick off next character --*/
         pchar = pdbfil[ipos];

/*--     State 1: Searching for first non-blank character --*/
         if (istate == 1)
         {
/* CHECK v.2.1.4
             if (pchar == '\/' || pchar == '\\' || pchar == ']')
CHECK v.2.1.4 */
             if (pchar == '/' || pchar == '\\' || pchar == ']')
/* CHECK v.2.1.4 */
             {
                 printf("*** ERROR in supplied name of file: [%s\n",
                    pdbfil);
                 return(-1);
             }

             if (pchar != ' ' && pchar != '.')
             {
                 *iend = ipos;
                 istate = 2;
             }
         }

/*--     State 2: Searching for end of extension, or end of directory path --*/
         else if (istate == 2)
         {

/*--        If character is a dot, and is the first dot, then note position--*/
            if (pchar == '.' && gotdot == FALSE)
            {
                *iend = ipos - 1;
                gotdot = TRUE;
            }

/*--        If character signifies the end of a directory path, note pstn --*/
/* CHECK v.2.1.4
            else if (pchar == '\/' || pchar == '\\' || pchar == ']')
CHECK v.2.1.4 */
            else if (pchar == '/' || pchar == '\\' || pchar == ']')
/* CHECK v.2.1.4 */
            {
                *istart = ipos + 1;
                finish = TRUE;
            }
        }

/*--    Step back a character --*/
        ipos = ipos - 1;
     }

/*-- Check whether file name is sensible --*/
     if (*istart > *iend)
     {
         printf("*** ERROR in supplied name of file: %s\n", pdbfil);
         printf("    No name found\n");
         return(-1);
     }
     return(0);
}
/* CHECK v.2.0 */

void fail(errstr)
char     *errstr;
{
/* CHECK v.2.1 */
    printf(" ");
/* CHECK v.2.1 */
    printf("\n*** %s\n\n", errstr);
    exit(-1);
}



/*============================================================================*/


static char *instr(ct, cs)
char *ct, *cs;
{
    int  l = strlen(cs);

    /* Input to INSTR consists of two strings CT and CS. It returns the address
       of the position in CT of the occurrence of the first character of the
       substring CS (if it is present) or zero. */

    /* printf("'%s' '%s' %d\n", ct, cs, l); */

    for (; *ct; ct++)
        if (!strncmp(ct, cs, l))
            return (ct);
    return (NULL);
}



/*============================================================================*/


void getcoord(x, y, z, chain, n, aacode, resnam, atmnam, atmnum)
float        *x, *y, *z;
char         *chain;
int          *n, *aacode, *atmnum;
char         *resnam, *atmnam;
{
    int       i, itemp;

    /* Extracts the following info from the line of input held in BUF[] :
       x,y,z     - the XYZ coordinates
       chain     - the one-letter chain code (if any)
       n         - the amino acid residue number
       resnam    - the three-letter residue name
       atmnam    - the four-letter atom name
       aacode    - an index into RNAMES for this RESNAM
                   (0 to STDAA-1 if a standard AA)
                   (STDAA to TOTNAA-1 if a non-standard AA)
                   (= TOTNAA if not-recognized)
       atmnum    - the atom id number
    */
    strncpy(atmnam, buf+12, 4);
    atmnam[4] = '\0';
    if (atmnam[2] == ' ')
        atmnam[3] = ' ';

    sscanf(buf+6, "%d", atmnum);

    sscanf(buf+17, "%c%c%c", resnam, resnam+1, resnam+2);
    resnam[3] = '\0';

    *chain = buf[21];

    sscanf(buf+22, "%4d", n);

    sscanf(buf+27, "%f%f%f", x, y, z);

    itemp = TOTNAA;
    for (i=0; i<TOTNAA; i++)
    {
        if (!strncmp(rnames[i], resnam,3))
        {
            itemp = i;
            break;
        }
    }
    *aacode = itemp;
}

/*============================================================================*/

int  alreadybonded(atom1, atom2)
int  atom1, atom2;

/* returns 1(true) if atoms 1+2 are already bonded (according to ICON), else 0*/

{
    int  i;

    for (i=0; i<MAXCON; i++)
    {
        if (icon[atom1][i]==atom2 || icon[atom2][i]==atom1)
            return(1);
    }
    return(0);
}

/*============================================================================*/

int nearlybonded(atom1, atom2)
int atom1, atom2;
/* returns 1 (true) if atoms 1+2 are 1-3 or 1-4 atoms (according to ICON), or
   else 0
*/

{
     int  i, j, k, l, itemp;

     itemp = debug;
     debug = 1;
/* CHECK v.1.0 */
     debug = 0;
/* CHECK v.1.0 */

     for (i=0; i<atom[atom1].ncon; i++)
     {
         j = icon[atom1][i];
         if (j == -1)
             continue;
         if (alreadybonded(j,atom2))
         {
             if (debug == 1)
             {
/* CHECK v.1.0
                 fprintf(dbgfp, "Atoms %5d and %5d are 1-3 neighbours\n",
                         atom1+1, atom2+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
                 printf("Atoms %5d and %5d are 1-3 neighbours\n",
                         atom1+1, atom2+1);
/* CHECK v.1.0 */
             }
             debug = itemp;
             return(1);
         } /* end if */
     } /* end i loop */
     for (i=0; i<atom[atom1].ncon; i++)
     {
         k = icon[atom1][i];
         if (k == -1)
             continue;
         for (j=0; j<atom[atom2].ncon; j++)
         {
             l = icon[atom2][j];
             if (l == -1)
                 continue;
             if (alreadybonded(k,l))
             {
                 if (debug == 1)
                 {
/* CHECK v.1.0
                     fprintf(dbgfp, "Atoms %5d and %5d are 1-4 neighbours\n",
                             atom1+1, atom2+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
                     printf("Atoms %5d and %5d are 1-4 neighbours\n",
                             atom1+1, atom2+1);
/* CHECK v.1.0 */
                 }
                 debug = itemp;
                 return(1);
             } /* end if */
         }  /* end of j loop */
     } /* end of i loop */
     debug = itemp;
     return(0);
} /* end of procedure */


/*============================================================================*/

int  atomindex(atomid)
int  atomid;
{
    int  i;

    for (i=0; i<natoms; i++)
        if (atom[i].atmnum == atomid)
            return(i);
    return(-1);
}

/*============================================================================*/

void printicon(string)
char *string;
{
    int i, j, k;

/* CHECK v.1.0
    fprintf(dbgfp, "\n%s\n", string);
CHECK v.1.0 */
/* CHECK v.1.0 */
    printf("\n%s\n", string);
/* CHECK v.1.0 */
    for (i=0; i<natoms; i++)
    {
/* CHECK v.1.0
        fprintf(dbgfp, "  atom %5d : ", i+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
        printf("  atom %5d : ", i+1);
/* CHECK v.1.0 */
        for (j=0; j<MAXCON; j++)
        {
            k = icon[i][j];
            if (k < 0)
/* CHECK v.1.0
                fprintf(dbgfp, "     ");
CHECK v.1.0 */
/* CHECK v.1.0 */
                printf("     ");
/* CHECK v.1.0 */
            else
/* CHECK v.1.0
                fprintf(dbgfp, "%5d", k+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
                printf("%5d", k+1);
/* CHECK v.1.0 */
        }
/* CHECK v.1.0
        fprintf(dbgfp, "\n");
CHECK v.1.0 */
/* CHECK v.1.0 */
        printf("\n");
/* CHECK v.1.0 */
    }
}

/*============================================================================*/

void printbonds()
{
    int i, j;

/* CHECK v.1.0
    fprintf(dbgfp, "\nContents of BONDS array:\n");
CHECK v.1.0 */
/* CHECK v.1.0 */
    printf("\nContents of BONDS array:\n");
/* CHECK v.1.0 */

    for (i=0; i<nconrecs; i++)
    {
/* CHECK v.1.0
        fprintf(dbgfp, "  CONECT %5d : ", bonds[i][0]);
CHECK v.1.0 */
/* CHECK v.1.0 */
        printf("  CONECT %5d : ", bonds[i][0]);
/* CHECK v.1.0 */
        for (j=1; j<MAXCNCOLS; j++)
        {
/* CHECK v.1.0
            fprintf(dbgfp, "%5d", bonds[i][j]);
CHECK v.1.0 */
/* CHECK v.1.0 */
            printf("%5d", bonds[i][j]);
/* CHECK v.1.0 */
        }
/* CHECK v.1.0
        fprintf(dbgfp, "\n");
CHECK v.1.0 */
/* CHECK v.1.0 */
        printf("\n");
/* CHECK v.1.0 */
    }
}

/*============================================================================*/

void readpdb(fname)
char      *fname;

{
    FILE  *pdbfp;
    char  line[160], keywd[7], stemp[6], pdbfn[128];
    int   i, iatom, jatom[10], offset, numfnd, itemp, token, dbgtemp;

    nconrecs = 0;
    dbgtemp  = debug;
    debug = 1;
/* CHECK v.1.0 */
    debug = 0;
/* CHECK v.1.0 */

    pdbfp = fopen(fname,"r");
    if (!pdbfp)
    {
/* CHECK v.2.1.3
        printf("Failed to open pdb file %s - no CONECT records used.\n", pdbfn);
*/
        printf("**** Failed to open pdb file %s - no CONECT records used.\n",
            pdbfn);
/* CHECK v.2.1.3 */
        return;
    }
    while (!feof(pdbfp))
    {
        if (!fgets(line, 160, pdbfp))
            break;
        sscanf (line, "%s", keywd);
        if (!keywd[0])
            break;
        token = 0;
        for (i=1; i<=NTOKENS; i++)
            if (!strcmp(keywd, tokstr[i-1]))
                token = i;
        switch (token)
        {
        case CONECT:
            if (nconrecs>=MAXCNRECS)
                fail("Too many CONECT records!");
            for (i=0; i<MAXCNCOLS; i++)
                bonds[nconrecs][i] = 0;
            offset = 6;
            strncpy (stemp, line+offset, 5);
            stemp[5] ='\0';
            numfnd = sscanf (stemp, "%d", &iatom);
            if (numfnd != 1)
                break;
            bonds[nconrecs][0] = iatom;
            for (i=1; i<MAXCNCOLS; i++)
            {
                offset = offset + 5;
                strncpy (stemp, line+offset, 5);
                stemp[5] = '\0';
                numfnd = sscanf(stemp, "%d", &itemp);
                if (numfnd == 1)
                    bonds[nconrecs][i] = itemp;
            }
            nconrecs++;
            break;
        default:
            break;
        }          /* end of switch statement */
    }              /* end of while statement */
    fclose (pdbfp);
    printf ("PDB file contained %5d CONECT records \n", nconrecs);
    if (debug == 1)
        printbonds();
    debug = dbgtemp;
}                  /* end of procedure */

/*============================================================================*/

void load_resbonds(restyp, ires)
int  restyp, ires;

/* Procedure to load into the icon[][] array all the bonds belonging to
   residue ires (numbered from 0 to (nres-1)). The residue is of the type
  "restyp", numbered from 0 to (STDAA -1). Only applies to standard amino
   acids! */
{
    char *p1, *p2, tempstr[9];
    int  iatom, jatom, ibond, kount, nconi, nconj;

    kount =nbonds[restyp];
    for (ibond=0; ibond<kount; ibond++)
    {
        switch (restyp)
        {
             case Ala:
                      strcpy(tempstr, alabonds[ibond]);
                      break;
             case Cys:
                      strcpy(tempstr, cysbonds[ibond]);
                      break;
             case Asp:
                      strcpy(tempstr, aspbonds[ibond]);
                      break;
             case Glu:
                      strcpy(tempstr, glubonds[ibond]);
                      break;
             case Phe:
                      strcpy(tempstr, phebonds[ibond]);
                      break;
             case Gly:
                      strcpy(tempstr, glybonds[ibond]);
                      break;
             case His:
                      strcpy(tempstr, hisbonds[ibond]);
                      break;
             case Ile:
                      strcpy(tempstr, ilebonds[ibond]);
                      break;
             case Lys:
                      strcpy(tempstr, lysbonds[ibond]);
                      break;
             case Leu:
                      strcpy(tempstr, leubonds[ibond]);
                      break;
             case Met:
                      strcpy(tempstr, metbonds[ibond]);
                      break;
             case Asn:
                      strcpy(tempstr, asnbonds[ibond]);
                      break;
             case Pro:
                      strcpy(tempstr, probonds[ibond]);
                      break;
             case Gln:
                      strcpy(tempstr, glnbonds[ibond]);
                      break;
             case Arg:
                      strcpy(tempstr, argbonds[ibond]);
                      break;
             case Ser:
                      strcpy(tempstr, serbonds[ibond]);
                      break;
             case Thr:
                      strcpy(tempstr, thrbonds[ibond]);
                      break;
             case Val:
                      strcpy(tempstr, valbonds[ibond]);
                      break;
             case Trp:
                      strcpy(tempstr, trpbonds[ibond]);
                      break;
             case Tyr:
                      strcpy(tempstr, tyrbonds[ibond]);
                      break;
             default:
/* CHECK v.2.1.3
                      printf("This message should never appear! \n");
*/
                      printf("**** This message should never appear! \n");
/* CHECK v.2.1.3 */
                      break;

        }   /* end of switch block */
        tempstr[8] = '\0'; /* this may be superfluos, but does no harm*/
        for (iatom = strtat[ires]; iatom<stopat[ires]; iatom++)
        {
            p1 = instr(tempstr, atom[iatom].atmnam);
            if (p1== NULL)
                continue;  /* atom "iatom" is NOT involved in bond "ibond: */
            for (jatom = iatom + 1; jatom<=stopat[ires]; jatom++)
            {
                p2 = instr(tempstr, atom[jatom].atmnam);
                if (p2 == NULL)
                    continue; /* atom "jatom" is NOT involved in bond "ibond" */
                nconi = atom[iatom].ncon;
                nconj = atom[jatom].ncon;
                if (nconi >= MAXCON)
                {
                    printf("Too many bonds to atom %5d", iatom);
                    printf("- excess bonds ignored. \n");
                    break; /* exit jatom loop, try next iatom */
                }
                if (nconj >= MAXCON)
                {
                    printf("Too many bonds to atom %5d", jatom);
                    printf("- excess bonds ignored. \n");
                    continue; /* try next jatom */
                }
                if (alreadybonded(iatom,jatom))
                    break;
                icon[iatom][nconi] = jatom;
                icon[jatom][nconj] = iatom;
                atom[iatom].ncon++;
                atom[jatom].ncon++;
            }  /* end of jatom loop */
        }  /* end of iatom loop */
    }  /* end of ibond loop */
}  /* end of procedure */

/*============================================================================*/

void load_icon()
{
    int    i, j, bondatom1, bondatom2, iatom1, iatom2, ncon1, ncon2,
           ires, restyp, jres, ityp, jtyp, itemp, jtemp, dbgtemp, katom;
    float  dist;

    dbgtemp = debug;
    debug   = 1;
/* CHECK v.1.0 */
    debug = 0;
/* CHECK v.1.0 */

    for (i=0; i<natoms; i++)
    {
        atom[i].ncon = 0;
        for (j=0; j<MAXCON; j++)
            icon[i][j] = -1;

    }

    /* Begin by loading the contents of the BONDS array into ICON */

    for (i=0; i<nconrecs; i++)
    {
        bondatom1 = bonds[i][0];
        if (bondatom1 == 0)
            continue;   /* increment i - try next row of BONDS */
        iatom1 = atomindex(bondatom1);
        if (iatom1 < 0)
        {
            printf("Could not find atom with id number %5d in ", bondatom1);
            printf("atom[].atmnum array.\n");
            printf("  atom %5d was defined by CONECT record %4d \n",
                    bondatom1, i+1);
            continue;   /* increment i - try next row of BONDS */
        }
        for (j=1; j<MAXCNCOLS; j++)
        {
            bondatom2 = bonds[i][j];
            if (bondatom2 == 0)
                continue;              /* next j */
            iatom2 = atomindex(bondatom2);
            if (iatom2 < 0)
            {
                printf("Could not find atom with id number %5d in ", bondatom2);
                printf("atom[].atmnum array.\n");
                printf("  atom %5d was defined by CONECT record %4d\n",
                        bondatom2, i+1);
                continue;   /* this must be continue not break */
            }
            if (iatom2 == iatom1)
            {
                printf("CONECT record %4d: ATOM2=ATOM1! \n", i+1);
                continue;     /* next j */
            }
            if (IGNLB == 1)
            {
                if (atom[iatom1].hetflg && (atom[iatom1].aacode == TOTNAA)&&
                                           (atom[iatom2].aacode < TOTNAA))
                {
                /*  printf("Ignoring CONECT bond from record %d", i+1);
                    printf(" - atom ids are %d and %d\n", bondatom1, bondatom2); */
                    continue; /* try next column j */
                }
                if (atom[iatom2].hetflg && (atom[iatom2].aacode == TOTNAA) &&
                                           (atom[iatom1].aacode < TOTNAA))
                {
                /*  printf("Ignoring CONECT bond from record %d", i+1);
                    printf(" - atom ids are %d and %d\n", bondatom2, bondatom1); */
                    continue;  /* try next column j */
                }
            }
            if (atom[iatom1].ncon >= MAXCON)
            {
                printf("Warning: too many connectivities for atom id %5d\n",
                        bondatom1);
                break;    /* goto end j loop, next conect record */
            }
            if (atom[iatom2].ncon >= MAXCON)
            {
                printf("Warning: too many connectivities for atom id %5d\n",
                        bondatom2);
                continue;   /* try next j column */
            }
            if (alreadybonded(iatom1,iatom2))
                continue;
            ncon1 = atom[iatom1].ncon;
            ncon2 = atom[iatom2].ncon;
            icon[iatom1][ncon1] = iatom2;
            icon[iatom2][ncon2] = iatom1;
            atom[iatom1].ncon++;
            atom[iatom2].ncon++;
        }    /* end of j loop */
    }        /* end of i loop */
    if (debug == 1)
        printicon("  ICON array after phase 1:");

    /* Explicit connectivities from CONECT records now loaded into ICON[][]. We
       are now ready to load implicit connectivities based on atom and residue
       names! (for the standard residues only at present) */

/* CHECK v.3.4
    for (ires=0; ires<nres; ires++)
CHECK v.3.4 */
    for (ires=0; ires<=nres; ires++)
/* CHECK v.3.4 */
    {
        itemp = strtat[ires];
        restyp = atom[itemp].aacode;
        if (restyp < STDAA)
            load_resbonds(restyp, ires);
    }
    if (debug == 1)
        printicon("  ICON array after phase 2:");

    /* Now load into ICON[][] all peptide bonds, from C of 1 residue to N of the
       next provided they are (a) from the same chain, and (b) close enough
       together! */

/* CHECK v.3.4
    for (ires=0; ires < nres-1; ires++)
CHECK v.3.4 */
    for (ires=0; ires <= nres-1; ires++)
/* CHECK v.3.4 */
    {
        jres = ires + 1;
        itemp = strtat[ires];
        jtemp = strtat[jres];
        ityp = atom[itemp].aacode;
        jtyp = atom[jtemp].aacode;
        if (ityp < TOTNAA && jtyp < TOTNAA)
        {
            iatom1 = -1;
            for (i=strtat[ires]; i<=stopat[ires]; i++)
            {
                if (atom[i].atmtyp == 2)    /* C atom */
                {
                    iatom1 = i;
                    break;
                }                    /* end if */
            }                        /* end of i loop */
            if (iatom1 == -1)
            {
                printf ("Cannot find atom C for residue %5d ", ires);
                katom = firstatom[ires];
                printf("(%3s_%c_%4d%c)\n", atom[katom].resnam,
                                           atom[katom].chnid,
                                           atom[katom].aanum,
                                           atom[katom].inscode);
                continue;   /* try next ires */
            }
            iatom2 = -1;
            for (j=strtat[jres]; j<=stopat[jres]; j++)
            {
                if (atom[j].atmtyp == 0)  /* N atom */
                {
                    iatom2 = j;
                    break;
                }                   /* end if */
            }                       /* end of j loop */
            if (iatom2 == -1)
            {
                printf ("Cannot find atom N for residue %5d ", jres);
                katom = firstatom[jres];
                printf("(%3s %c %4d%c)\n", atom[katom].resnam,
                                           atom[katom].chnid,
                                           atom[katom].aanum,
                                           atom[katom].inscode);
                continue;   /* try next ires */
            }
            ncon1 = atom[iatom1].ncon;
            ncon2 = atom[iatom2].ncon;
            if (atom[iatom1].chnid != atom[iatom2].chnid)
                continue;  /* try next ires */
            dist = SQR(atom[iatom1].x - atom[iatom2].x) +
                   SQR(atom[iatom1].y - atom[iatom2].y) +
                   SQR(atom[iatom1].z - atom[iatom2].z);
            if (dist > SQR(PEPBND))
                continue;  /* try next ires */
            if (alreadybonded(iatom1, iatom2))
                continue;  /* try next ires */

            if (ncon1 >= MAXCON)
            {
                printf("C atom from residue %5d has too many connectivities!\n",
                       ires);
                continue;   /* try next ires */
            }
            if (ncon2 >= MAXCON)
            {
                printf("N atom from residue %5d has too many connectivities!\n",
                        jres);
                continue;   /* try next ires */
            }
            icon[iatom1][ncon1] = iatom2;
            icon[iatom2][ncon2] = iatom1;
            atom[iatom1].ncon++;
            atom[iatom2].ncon++;
        }   /* end if */
    }   /* end of ires loop */
    if (debug == 1)
        printicon("  ICON array after phase 3:");
    debug = dbgtemp;
}   /* end of procedure */

/*============================================================================*/

void find_brakes()
{
    float  dist;
    int    ires, jres, iatom, jatom;

    debug   =  1;
/* CHECK v.1.0 */
    debug = 0;
/* CHECK v.1.0 */

    nbrakes = -1;
    for (ires = 0; ires<(nres -1); ires++)
    {
        jres = ires + 1;
        if (CHNBRKFLAG == 1 || CHNBRKFLAG == 4)
        {
            dist = 999.9;
            if (cafound[ires] && cafound[jres])
            {
                dist = SQR(caxyz[ires][0] - caxyz[jres][0]) +
                       SQR(caxyz[ires][1] - caxyz[jres][1]) +
                       SQR(caxyz[ires][2] - caxyz[jres][2]);
            }
            if (dist>SQR(CADBND))
            {
/* CHECK v.3.6.2 */
	      if (nbrakes < MAXBRKS - 1)
		{
/* CHECK v.3.6.2 */
                nbrakes = nbrakes + 1;
                brake[nbrakes][0] = ires;
                brake[nbrakes][1] = jres;
/* CHECK v.3.6.2 */
		}
/* CHECK v.3.6.2 */
                if (debug == 1)
                {
                    if (nbrakes == 0)
/* CHECK v.1.0
                        fprintf(dbgfp, "\n");
CHECK v.1.0 */
/* CHECK v.1.0 */
                        printf("\n");
/* CHECK v.1.0 */
/* CHECK v.1.0
                    fprintf(dbgfp,"Chain break detected between residues");
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("Chain break detected between residues");
/* CHECK v.1.0 */
/* CHECK v.1.0
                    fprintf(dbgfp, "%5d and %5d", ires, jres);
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("%5d and %5d", ires, jres);
/* CHECK v.1.0 */
                    iatom = firstatom[ires];
                    jatom = firstatom[jres];
/* CHECK v.1.0
                    fprintf(dbgfp,"%3s%c%4d%c and %3%c%4d%c\n", 
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
CHECK v.1.0 */
/* CHECK v.1.0 */
/* CHECK v.3.2
                    printf("%3s%c%4d%c and %3%c%4d%c\n", 
CHECK v.3.2 */
/* CHECK v.3.2 */
                    printf("%3s%c%4d%c and %3s%c%4d%c\n", 
/* CHECK v.3.2 */
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
/* CHECK v.1.0 */
                }
                continue;  /* look at next residue pair */
            }   /* end if */
        }   /* end if */
        if (CHNBRKFLAG == 1)
            continue;

        if (CHNBRKFLAG == 2 || CHNBRKFLAG == 4)
        {
            dist = 999.9;
            if (cfound[ires] && nfound[jres])
            {
                dist = SQR(cxyz[ires][0] - nxyz[jres][0]) +
                       SQR(cxyz[ires][1] - nxyz[jres][1]) +
                       SQR(cxyz[ires][2] - nxyz[jres][2]);
            }
            if (dist > SQR(PEPBND))
            {
/* CHECK v.3.6.2 */
	      if (nbrakes < MAXBRKS - 1)
		{
/* CHECK v.3.6.2 */
                nbrakes = nbrakes + 1;
                brake[nbrakes][0] = ires;
                brake[nbrakes][1] = jres;
/* CHECK v.3.6.2 */
		}
/* CHECK v.3.6.2 */
                if (debug == 1)
                {
                    if (nbrakes == 0)
/* CHECK v.1.0
                        fprintf(dbgfp,"\n");
CHECK v.1.0 */
/* CHECK v.1.0 */
                        printf("\n");
/* CHECK v.1.0 */
/* CHECK v.1.0
                    fprintf(dbgfp,"Chain break detected between residues");
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("Chain break detected between residues");
/* CHECK v.1.0 */
/* CHECK v.1.0
                    fprintf(dbgfp, "%5d and %5d", ires, jres);
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("%5d and %5d", ires, jres);
/* CHECK v.1.0 */
                    iatom = firstatom[ires];
                    jatom = firstatom[jres];
/* CHECK v.1.0
                    fprintf(dbgfp,"%3s%c%4d%c and %3%c%4d%c\n", 
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
CHECK v.1.0 */
/* CHECK v.1.0 */
/* CHECK v.3.2
                    printf("%3s%c%4d%c and %3%c%4d%c\n", 
CHECK v.3.2 */
/* CHECK v.3.2 */
                    printf("%3s%c%4d%c and %3s%c%4d%c\n", 
/* CHECK v.3.2 */
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
/* CHECK v.1.0 */
                }
                continue;  /* look at next residue pair */
            }   /* end if */
        }   /* end if */
        if (CHNBRKFLAG == 2)
            continue;

        if (CHNBRKFLAG == 3 || CHNBRKFLAG == 4)
        {
            if (!cafound[ires] || !cfound[ires] ||
                 !nfound[jres] || !cafound[jres])
            {
/* CHECK v.3.6.2 */
	      if (nbrakes < MAXBRKS - 1)
		{
/* CHECK v.3.6.2 */
                nbrakes = nbrakes + 1;
                brake[nbrakes][0] = ires;
                brake[nbrakes][1] = jres;
/* CHECK v.3.6.2 */
		}
/* CHECK v.3.6.2 */
                if (debug == 1)
                {
                    if (nbrakes == 0)
                        printf("\n");
/* CHECK v.1.0
                    fprintf(dbgfp,"Chain break detected between residues");
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("Chain break detected between residues");
/* CHECK v.1.0 */
/* CHECK v.1.0
                    fprintf(dbgfp,"%5d and %5d", ires, jres);
CHECK v.1.0 */
/* CHECK v.1.0 */
                    printf("%5d and %5d", ires, jres);
/* CHECK v.1.0 */
                    iatom = firstatom[ires];
                    jatom = firstatom[jres];
/* CHECK v.1.0
                    fprintf(dbgfp,"%3s%c%4d%c and %3%c%4d%c\n", 
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
CHECK v.1.0 */
/* CHECK v.1.0 */
/* CHECK v.3.2
                    printf("%3s%c%4d%c and %3%c%4d%c\n", 
CHECK v.3.2 */
/* CHECK v.3.2 */
                    printf("%3s%c%4d%c and %3s%c%4d%c\n", 
/* CHECK v.3.2 */
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
/* CHECK v.1.0 */
                }
                continue; /* next ires-jres pair */
            } /* end if */
        } /* end if */
    } /* end loop */
} /* end procedure */


/*============================================================================*/

void load_ststan()
{
    int     i, stres, curres;
    char    stchn, curchn, stinscode, curinscode;

    /* Load start and stop atom numbers for each residue */

    debug = 1;
/* CHECK v.1.0 */
    debug = 0;
/* CHECK v.1.0 */

    stchn     = atom[0].chnid;
    stres     = atom[0].aanum;
    stinscode = atom[0].inscode;
    nres      = 0;
    strtat[nres] = 0;

    for (i=0; i < natoms; i++)
    {
        curchn     = atom[i].chnid;
        curres     = atom[i].aanum;
        curinscode = atom[i].inscode;
        if (curchn     != stchn ||
            curres     != stres ||
            curinscode != stinscode)
        {
            stopat[nres] = i-1;
            if (debug == 1)
            {
/* CHECK v.1.0
            fprintf(dbgfp,
                "Residue %4d starts at atom %5d and stops at atom %5d \n",
                    nres+1, strtat[nres]+1, stopat[nres]+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
            printf("Residue %4d starts at atom %5d and stops at atom %5d \n",
                    nres+1, strtat[nres]+1, stopat[nres]+1);
/* CHECK v.1.0 */
            }
            nres++;
            strtat[nres] = i;
            stchn     = curchn;
            stres     = curres;
            stinscode = curinscode;
        }                              /* end if */
    }                                  /* end of i loop */
    stopat[nres] = natoms-1;
    if (debug == 1)
    {
/* CHECK v.1.0
        fprintf(dbgfp,
                  "Residue %4d starts at atom %5d and stops at atom %5d \n",
                       nres+1, strtat[nres]+1, stopat[nres]+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
        printf("Residue %4d starts at atom %5d and stops at atom %5d \n",
                       nres+1, strtat[nres]+1, stopat[nres]+1);
/* CHECK v.1.0 */
    } /* end if */
} /* end of procedure load_ststan */


/*============================================================================*/


void find_nb()
{
    int     i, ii, j, jj, k, nnb, ri, rj, ires, jres, gap, ibrk;
    int     iatom, jatom;
    char    bndtyp[3], acctyp[4], dontyp[4], space;
    float   ca_d, d;
    float   mindist;

    debug = 0;
    nnb   = 0;
    space = ' ';
/*CHECK*/
    strncpy(brcode,"    ",4);
	brcode[4] = '\0';
/*CHECK*/

    /* First Pass : Check for disulphide bridges */

    printf("Checking for disulphide bridges . . . ");
    for (i = 0; i < natoms; i++)
    {
     /* printf("%4d", i); */
        for (j = i+1; j < natoms; j++)
        {
            if (atom[i].aacode == Cys && atom[i].atmtyp == 5 &&
                atom[j].aacode == Cys && atom[j].atmtyp == 5 &&
                fabs(atom[i].x - atom[j].x) < SSDIST &&
                fabs(atom[i].y - atom[j].y) < SSDIST &&
                fabs(atom[i].z - atom[j].z) < SSDIST)
            {
                d = SQR(atom[i].x - atom[j].x) +
                    SQR(atom[i].y - atom[j].y) +
                    SQR(atom[i].z - atom[j].z);

                if (d <= SQR(SSDIST))
                {
                    atom[i].ssflg = j;
                    atom[j].ssflg = i;
                }
            }
        } /* end of j loop */
    } /* end of i loop */
    printf("done.\n");

    /* Second Pass : check for neighbours */

    puts("Checking for neighbours . . . ");

/* CHECK v.2.0.3
    for (ires=0; ires<nres-1; ires++) */
/* CHECK v.2.0.3 */
    for (ires=0; ires<=nres-1; ires++)
/* CHECK v.2.0.3 */
    {
/* CHECK v.2.0.3
        for (jres=ires+1; jres<nres; jres++) */
/* CHECK v.2.0.3 */
        for (jres=ires+1; jres<=nres; jres++)
/* CHECK v.2.0.3 */
        {
            /* Check that residues IRES and JRES are from the same chain. */
            iatom = strtat[ires];
            jatom = strtat[jres];
/* CHECK v.2.0.3 
            if (atom[iatom].chnid != atom[jatom].chnid)
                continue; */
/* CHECK v.2.0.3 */
            if ((atom[iatom].chnid != atom[jatom].chnid)&&
               (atom[iatom].chnid != '-')&&
               (atom[jatom].chnid != '-'))
                continue;
/* CHECK v.2.0.3 */

            if (debug == 1)
            {
                printf("checking residue %4d against residue %4d\n",
                        ires+1, jres+1);
            }
            mindist = 999.9;

            for (i=strtat[ires]; i<=stopat[ires]; i++)
            {
                for (j=strtat[jres]; j<=stopat[jres]; j++)
                {
                    /* Check that atoms I and J are within a "reasonable"
                       distance of each other. */
                    if (fabs(atom[i].x - atom[j].x) > NBDIST ||
                        fabs(atom[i].y - atom[j].y) > NBDIST ||
                        fabs(atom[i].z - atom[j].z) > NBDIST)
                        continue;

                    /* Check that atoms I and J are not S-S bridged. */
                    if (atom[i].ssflg == j)
                        continue;

                    /* Check that atoms I and J are not bonded to each other! */
                    if (alreadybonded(i,j))
                        continue;

                    /* Check that they are not 1-3 or 1-4 contacts from
                       adjacent residues */
                    if (atom[i].aacode<TOTNAA && atom[j].aacode<TOTNAA
                                              && nearlybonded(i,j))
                        continue;

                    /* Go ahead and determine the distance IJ - if the current
                       distance is the closest so far (for this residue pair)
                       retain its value and those of I and J for later. */

                    d = SQR(atom[i].x - atom[j].x) +
                        SQR(atom[i].y - atom[j].y) +
                        SQR(atom[i].z - atom[j].z);

                    if (d <= mindist)
                    {
                        ii      = i;
                        jj      = j;
                        mindist = d;
                    }                                               /* end if */
                }                                            /* end of j loop */
            }                                                /* end of i loop */
            if (mindist <= SQR(NBDIST))
            {
                nnb++;
                bndtyp[0] = bndtyp[1] = '?';
                acctyp[3] = dontyp[3] = bndtyp[2] = '\0';

                if (atom[ii].aacode == TOTNAA)
                    bndtyp[0] = 'H';
                else
                {
                    if (atom[ii].aacode <STDAA)
                    {
                        bndtyp[0] = (atom[ii].atmtyp <=3) ? 'M' : 'S';
                    }
                    else
                    {
                        bndtyp[0] = (atom[ii].atmtyp <=3) ? 'm' : 's';
                    }
                }
                if (atom[jj].aacode == TOTNAA)
                    bndtyp[1] = 'H';
                else
                {
                    if (atom[jj].aacode <STDAA)
                    {
                        bndtyp[1] = (atom[jj].atmtyp <=3) ? 'M' : 'S';
                    }
                    else
                    {
                        bndtyp[1] = (atom[jj].atmtyp <=3) ? 'm' : 's';
                    }
                }
             /* if (atom[ii].hetflg)
                    strcpy(dontyp, atom[ii].resnam);
                else
                    strcpy(dontyp, "A-A");

                if (atom[jj].hetflg)
                    strcpy(acctyp, atom[jj].resnam);
                else
                    strcpy(acctyp, "A-A"); */

                ri = atom[ii].caindex;
                rj = atom[jj].caindex;

                if (atom[ii].aacode==TOTNAA || (ri == -99) || !cafound[ri] ||
                    atom[jj].aacode==TOTNAA || (rj == -99) || !cafound[rj]  )
                    gap = -1;
                else
                {
                    if (nbrakes < 0)
                    {
                       gap = abs (ri -rj);
                    }
                    else
                    {
                       for (ibrk = 0; ibrk <= nbrakes; ibrk++)
                       {
                         /*printf("ibrk=%4d  start=%4d  stop=%4d  ",      */
                         /*        ibrk, brake[ibrk][0], brake[ibrk][1]); */
                         /*printf("ri=%4d  rj=%4d  min=%4d  max=%4d\n",   */
                         /*        ri, rj, min(ri,rj), max(ri,rj) );      */
                           if (brake[ibrk][0] >= min (ri,rj) &&
                               brake[ibrk][1] <= max (ri,rj))
                           {
                             /*printf("Break detected\n"); */
                               gap = -1;
                               break;
                           }
                           else
                               gap = abs (ri - rj);
                       }     /* end for */
                    }
                }

                if (atom[ii].aacode==TOTNAA || (ri == -99) || !cafound[ri] ||
                    atom[jj].aacode==TOTNAA || (rj == -99) || !cafound[rj]  )
                    ca_d = -1.0;
                else
                {
                    ca_d = SQR(caxyz[ri][0] - caxyz[rj][0]) +
                           SQR(caxyz[ri][1] - caxyz[rj][1]) +
                           SQR(caxyz[ri][2] - caxyz[rj][2]) ;
                    if (ca_d > SQR(CAWARN))
                    {
                        printf("Warning: unusual CA-CA distance for ");
                        printf("contact number %5d", nnb);
                        printf(" : distance = %5.1f\n", sqrt(ca_d));
/* CHECK v.3.6.2 */
			iatom = firstatom[ri];
			printf("   Residues   %s %4d%c %c  <-> ",
			       atom[iatom].resnam,
			       atom[iatom].aanum,
			       atom[iatom].inscode,
			       atom[iatom].chnid);
			iatom = firstatom[rj];
			printf(" %s %4d%c %c",
			       atom[iatom].resnam,
			       atom[iatom].aanum,
			       atom[iatom].inscode,
			       atom[iatom].chnid);
                        printf("   (ri = %5d,  rj = %5d)\n", ri, rj);
/* CHECK v.3.6.2 */
                    }
                }
#if 0
                printf("%4s",   brcode);
                printf("%c",    (atom[ii].chnid=='-') ? space : atom[ii].chnid);
                printf("%4s/",  brcode);
                printf("%c",    atom[ii].chnid);
                printf("%04d",  atom[ii].aanum);
                printf("%c",    atom[ii].inscode);
             /* printf("%3s",   dontyp);               */
                printf("%3s",   atom[ii].resnam);
                printf("%4s",   atom[ii].atmnam);
                printf("%c",    atom[ii].strucsum);
                printf("%c",    (atom[jj].chnid=='-') ? space : atom[jj].chnid);
                printf("%4s/",  brcode);
                printf("%c",    atom[jj].chnid);
                printf("%04d",  atom[jj].aanum);
                printf("%c",    atom[jj].inscode);
             /* printf("%3s",   acctyp);               */
                printf("%3s",   atom[jj].resnam);
                printf("%4s",   atom[jj].atmnam);
                printf("%c",    atom[jj].strucsum);
                printf("%4.1f", sqrt(mindist));
                printf(" %2s",  bndtyp);
                printf("%4d",   gap);
                printf("%4.1f", (ca_d <= 0.0) ? -1.0 : sqrt(ca_d));
                printf("%6d\n", nnb);
#endif
                fprintf(ofp, "%4s",   brcode);
                fprintf(ofp, "%c",    (atom[ii].chnid == '-') ? space :
                                       atom[ii].chnid);
                fprintf(ofp, "%4s/",  brcode);
                fprintf(ofp, "%c",    atom[ii].chnid);
                fprintf(ofp, "%04d",  atom[ii].aanum);
                fprintf(ofp, "%c",    atom[ii].inscode);
             /* fprintf(ofp, "%3s",   dontyp);               */
                fprintf(ofp, "%3s",   atom[ii].resnam);
                fprintf(ofp, "%4s",   atom[ii].atmnam);
                fprintf(ofp, "%c",    atom[ii].strucsum);
                fprintf(ofp, "%c",    (atom[jj].chnid == '-') ? space :
                                       atom[jj].chnid);
                fprintf(ofp, "%4s/",  brcode);
                fprintf(ofp, "%c",    atom[jj].chnid);
                fprintf(ofp, "%04d",  atom[jj].aanum);
                fprintf(ofp, "%c",    atom[jj].inscode);
             /* fprintf(ofp, "%3s",   acctyp);               */
                fprintf(ofp, "%3s",   atom[jj].resnam);
                fprintf(ofp, "%4s",   atom[jj].atmnam);
                fprintf(ofp, "%c",    atom[jj].strucsum);
                fprintf(ofp, "%4.1f", sqrt(mindist));
                fprintf(ofp, " %2s",  bndtyp);
                fprintf(ofp, "%4d",   gap);
                fprintf(ofp, "%4.1f", (ca_d <= 0.0) ? -1.0 : sqrt(ca_d));
                fprintf(ofp, "%6d\n", nnb);
            }                                                       /* end if */
        }                                                 /* end of jres loop */
    }                                                     /* end of ires loop */
    printf("%5d neighbouring contacts found.\n", nnb);
}



/*============================================================================*/


void parsefn(fname, rootnam)
char        *fname, *rootnam;
{
    char    *p;

    p = fname + strlen(fname);
    while (*--p != '.');
    while (isalnum(*--p))
        if (p == fname)
        {
            p--;
            break;
        }
    strcpy(rootnam, p + 1);
    p = rootnam + strlen(rootnam);
    while (*--p != '.');
    *(p + 1) = '\0';
}



/*============================================================================*/


void procfile(fname)
char         *fname;
{
    int       i, token, namino, aac, oldresnum, reskount, atmnum;
    short     caonly;
    char     *p, atmnam[5], resnam[4], chain, strucsum, inscode, contchar,
              oldchain, oldinscode, space, altcode;
    float     x, y, z, dist;
    char      outfn[128], sstfn[128], sstbuf[160], pdbfn[128];
    FILE     *ifp, *sstfp;

    printf("Processing : %s\n", fname);
    ifp = fopen(fname, "r");
    if (!ifp)
    {
/* CHECK v.2.1 */
/* CHECK v.2.1.3
        printf("Failed to open specified input file!\n");
*/
        printf("**** Failed to open specified input file!\n");
/* CHECK v.2.1.3 */
/* CHECK v.2.1 */
        return;
    }

    parsefn(fname, outfn);
    strcat(outfn, "nb");

    if (!tablefn[0])
    {
        ofp = fopen(outfn, "w");
        if (!ofp)
        {
/* CHECK v.2.1 */
/* CHECK v.2.1.3
            printf("Failed to open output file %s\n", outfn);
*/
            printf("**** Failed to open output file %s\n", outfn);
/* CHECK v.2.1.3 */
/* CHECK v.2.1 */
            return;
        }
    }

    debug      = 1;
/* CHECK v.1.0 */
    debug = 0;
/* CHECK v.1.0 */
    natoms     = 0;
    space      = ' ';
    caonly     = TRUE;
    chain      = '?';
    oldchain   = '?';
    oldinscode = '?';
    oldresnum  = -999;
    reskount   = -1;

    while (!feof(ifp))
    {
        if (!fgets(buf, 160, ifp))
            break;
        sscanf(buf, "%s", keyword);                   /* Read the record name */

        if (!keyword[0])            /* Odd - there isn't a record name! Exit. */
            break;

        token = 0;
        for (i = 1; i <= NTOKENS; i++)                  /* Decode record type */
            if (!strcmp(keyword, tokstr[i - 1]))
                token = i;

        if (token == ENDENT)
            break;

        switch (token)
        {
/* CHECK v.1.0
        case HEADER:
            contchar = buf[9];
            if (contchar != space)
            {
                printf("HEADER continuation line identified and ignored.\n");
                break;
            }
            strncpy(brcode, buf+62, 4);
            brcode[4] = '\0';
            strcpy(pdbfn, "/data/pdb/p");
            strcat(pdbfn, brcode);
            strcat(pdbfn, ".pdb");
            for (p = pdbfn; *p; p++)
                if (isupper(*p))
                    *p = tolower(*p);
            readpdb(pdbfn);
            strcpy(sstfn, "/oracle/sst/p");
            strcat(sstfn, brcode);
            strcat(sstfn, ".sst");
            for (p=sstfn; *p; p++)
                if (isupper(*p))
                    *p = tolower(*p);
            sstfp = fopen(sstfn, "r");
            if (!sstfp)
            {
                printf("Failed to open corresponding SST file %s\n", sstfn);
                fclose(ifp);
                if (!tablefn[0])
                    fclose(ofp);
                return;
            }
            for (i=0; i<7; i++)
                fgets(sstbuf,160,sstfp);
            break;
CHECK v.1.0 */

        case HETATM:
/* CHECK v.1.0 */
            strucsum = space;
/* CHECK v.1.0 */

        case ATOM:
            if (natoms >= MAXNATM)
/* CHECK v.2.1.3
                fail("Too many atoms! Increase MAXNATM.");
*/
                fail("**** Too many atoms! Increase MAXNATM.");
/* CHECK v.2.1.3 */

            altcode = buf[16];
            if (altcode != space && altcode != 'A')
                continue;

            inscode = buf[26];
            if (inscode == space)
                inscode = '-';

            getcoord(&x,&y,&z,&chain,&namino,&aac,resnam,atmnam,&atmnum);
            p = NULL;
            if (token == ATOM && aac == TOTNAA)
            {
                printf("Residue %3s is not recognized by NEIGHBOUR \n", resnam);
                break;
                /* ignore unrecognized residues */
            }
            /* Arrival here implies the atom is either an ATOM or a HETATM from
               a recognized residue, OR a HETATM from a ligand. Check first
               whether or not this atom is the first atom of a new residue. If
               it is, update the residue counters, and initialize the XYZ
               coordinate arrays. Additionally, if it is NOT from a ligand,
               read in its SSSUM from the SST file. */

            if (chain != oldchain   ||
                namino != oldresnum  ||
                inscode != oldinscode)
            {
                oldchain   = chain;
                oldresnum  = namino;
                oldinscode = inscode;
                reskount   = reskount + 1;
                firstatom[reskount] = natoms;
                if (aac < TOTNAA)
                {
/* CHECK v.1.0
                    strucsum = '?';
                    fgets(sstbuf, 160, sstfp);
                    strucsum = sstbuf[25];
CHECK v.1.0 */
/* CHECK v.1.0 */
                    strucsum = space;
/* CHECK v.1.0 */
                }
                else
                    strucsum = space;
/* CHECK v.1.0
                fprintf(dbgfp, "Residue %5d is %3s %c %4d%c %c\n",
                                reskount+1, resnam, chain, namino, inscode,
                                strucsum);
CHECK v.1.0 */
/* CHECK v.1.0 */
                printf("Residue %5d is %3s %c %4d%c %c\n",
                                reskount+1, resnam, chain, namino, inscode,
                                strucsum);
/* CHECK v.1.0 */
                if (reskount >= MAXNRES)
/* CHECK v.2.1.3
                    fail("Too many residues! Increase MAXNRES.");
*/
                    fail("**** Too many residues! Increase MAXNRES.");
/* CHECK v.2.1.3 */

                caxyz[reskount][0] = -999.9;
                caxyz[reskount][1] = -999.9;
                caxyz[reskount][2] = -999.9;
                cxyz[reskount][0]  = -999.9;
                cxyz[reskount][1]  = -999.9;
                cxyz[reskount][2]  = -999.9;
                nxyz[reskount][0]  = -999.9;
                nxyz[reskount][1]  = -999.9;
                nxyz[reskount][2]  = -999.9;
                cafound[reskount]  =  FALSE;
                cfound[reskount]   =  FALSE;
                nfound[reskount]   =  FALSE;
            }
            if ((token == ATOM) || (token == HETATM && aac < TOTNAA))
            {
                if ((p = instr(necatm[aac], atmnam)) != NULL)
                {
                    atom[natoms].caindex = reskount;
                    atom[natoms].atmtyp  = (p-necatm[aac])/4;
                    if (atom[natoms].atmtyp == 1)
                    {
                        caxyz[reskount][0] = x;
                        caxyz[reskount][1] = y;
                        caxyz[reskount][2] = z;
                        cafound[reskount]  = TRUE;
                    }
                    else
                    {
                        caonly = FALSE;
                    }
                    if (atom[natoms].atmtyp == 0)
                    {
                        nxyz[reskount][0] = x;
                        nxyz[reskount][1] = y;
                        nxyz[reskount][2] = z;
                        nfound[reskount]  = TRUE;
                    }
                    if (atom[natoms].atmtyp == 2)
                    {
                        cxyz[reskount][0] = x;
                        cxyz[reskount][1] = y;
                        cxyz[reskount][2] = z;
                        cfound[reskount]  = TRUE;
                    }
                }
                else
                {
                    printf("Unrecognized atom name %4s from residue %3s \n",
                            atmnam, resnam);
                    break; /* ignore unrecognized atoms */
                }
            }
            else  /* atom is from a ligand */
            {
                atom[natoms].atmtyp  = -99;
                atom[natoms].caindex = -99;
            }
            atom[natoms].atmnum   = atmnum;
            atom[natoms].x        = x;
            atom[natoms].y        = y;
            atom[natoms].z        = z;
            atom[natoms].aanum    = namino;
            atom[natoms].aacode   = aac;
            atom[natoms].strucsum = strucsum;
            atom[natoms].chnid    = (chain   == space) ? '-' : chain;
            atom[natoms].altcode  = (altcode == space) ? '-' : altcode;
            atom[natoms].inscode  = inscode;
            atom[natoms].hetflg   = (token == HETATM);
            atom[natoms].ssflg    = -1;
            strcpy(atom[natoms].resnam, resnam);
            strcpy(atom[natoms].atmnam, atmnam);

            if (debug == 1)
            {
/* CHECK v.1.0
                fprintf(dbgfp,
                     "%5d %c %5d %c %3s %4s (%3d) %c %7.2f%7.2f%7.2f %5d\n",
                       natoms+1, chain, namino, inscode, resnam, atmnam,
                       atom[natoms].atmtyp, strucsum, x, y, z, reskount+1);
CHECK v.1.0 */
/* CHECK v.1.0 */
                printf("%5d %c %5d %c %3s %4s (%3d) %c %7.2f%7.2f%7.2f %5d\n",
                       natoms+1, chain, namino, inscode, resnam, atmnam,
                       atom[natoms].atmtyp, strucsum, x, y, z, reskount+1);
/* CHECK v.1.0 */
            }
            natoms++;
            break;

        default:                    /* Ignore all other types in this version */
            break;

        } /* end switch */
    } /* end while */

    printf("%5d potential contact atoms selected from %5d residues.\n\n",
            natoms, reskount+1);


    if (natoms && !caonly)
    {
        load_ststan();
        load_icon();
        find_brakes();
        find_nb();
    }

    if (caonly)
        puts("\n**** C-Alpha only file! ****\n");

/* CHECK v.1.0
    fclose(sstfp);
CHECK v.1.0 */
    fclose(ifp);
/* CHECK v.1.0
    if (!tablefn[0])
        fclose(ofp);
CHECK v.1.0 */
} /* end of procedure */


/*============================================================================*/


main()
{
    char     tmpstr[128];
    FILE    *bfp;
/* CHECK v.2.0 */
    char     tmpfn[160];
    int      iend, ierror, inew, ipos, istart;
/* CHECK v.2.0 */

    printf("\nNeighbouring Contact Locator - By David Jones, August 1990\n");

    printf("\nProgram is configured for %5d atoms and %5d residues.\n",
            MAXNATM, MAXNRES);

/* CHECK v.1.0
    dbgfp = fopen ("nbdebug.dat", "w");
CHECK v.1.0 */

/* CHECK v.2.0
    printf("\nEnter table output filename (or '@' for multiple files):");

    fgets(tmpstr, sizeof(tmpstr), stdin);
    sscanf(tmpstr, "%s", tablefn);
    if (tablefn[0] != '@')
    {
        ofp = fopen(tablefn, "w");
        if (!ofp)
            fail("Failed to open table output file!");
    }
    else
        tablefn[0] = '\0';

    for (;;)
    {
        printf("\nEnter the name of the next file to be processed:\n");
        printf("  just a file name for a Brookhaven Format (.NEW) file or\n");
        printf("  @file name for a file holding names of Brookhaven ");
        printf("format files.\n");
        printf("If you've had enough, a blank line lets you escape!\n\n");
        fgets(tmpstr, sizeof(tmpstr), stdin);
        if (tmpstr[0] == '\0' || tmpstr[0] == '\n')
            break;
        if (tmpstr[0] == '@')
        {
            bfp = fopen(tmpstr + 1, "r");
            if (!bfp)
            {
                printf("Failed to open file of file names as input.\n");
                continue;
            }
            while (!feof(bfp))
            {
                if (!fgets(tmpstr, 128, bfp))
                {
                    fclose(bfp);
                    break;
                }
                tmpstr[strlen(tmpstr) - 1] = '\0';
                procfile(tmpstr);
            }
            fclose(bfp);
        }
        else
            procfile(tmpstr);
    }
CHECK v.2.0 */

/* CHECK v.2.0 */
    printf("\nEnter the name of original coordinates file\n");
    fgets(tmpstr, sizeof(tmpstr), stdin);
    if (tmpstr[0] != '\0' && tmpstr[0] != '\n')
    {
        ierror = getnam(tmpstr,&istart,&iend);
        if (ierror != -1)
        {

/*--        Derive name of input file --*/
            for (ipos = istart, inew = 0; ipos < iend + 1; ipos++, inew++)
                tmpfn[inew] = tmpstr[ipos];
            tmpfn[inew] = '\0';
            strcat(tmpfn,".new");
            printf("Input file name %s\n",tmpfn);            

/*--        Derive name of output file --*/
            for (ipos = istart, inew = 0; ipos < iend + 1; ipos++, inew++)
                tablefn[inew] = tmpstr[ipos];
            tablefn[inew] = '\0';
            strcat(tablefn,".nb");
            printf("Output file name %s\n",tablefn);            

/*--        Open output file --*/
            ofp = fopen(tablefn, "w");
            if (!ofp)
/* CHECK v.2.1.3
                fail("Failed to open table output file!"); 
*/
                fail("**** Failed to open table output file!"); 
/* CHECK v.2.1.3 */

/*--        Process input file --*/
            procfile(tmpfn);
        }
    }
/* CHECK v.2.0 */

/* CHECK v.2.3 */
    printf("* Program completed\n");
/* CHECK v.2.3 */

/* CHECK v.3.2 */
    exit(0);
/* CHECK v.3.2 */
}

/*============================================================================*/
