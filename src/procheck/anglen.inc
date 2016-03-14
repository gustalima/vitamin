C**************************************************************************
C
C  ANGLEN.INC   -   Include file for ANGLEN.FOR
C
C----------------------------------------------------------------------+--- 

CVAX      IMPLICIT NONE

      INTEGER       MAXEXC, MCHAIN, MCHANG, MCHLEN, MXAMIN, MXITEM,
     -              MXONE, MXTYPE, NMAIN, NMATOM, NTYPE, SCHAIN
      PARAMETER    (
     -              MAXEXC =      6,
     -              MCHAIN =      1,
     -              MCHANG =     20,
     -              MCHLEN =     11,
     -              MXAMIN =     20,
     -              MXITEM =   6000,
CHECK v.3.4.4-->
C     -              MXONE  =   6000,
     -              MXONE  =  10000,
CHECK v.3.4.4<--
     -              MXTYPE =    100,
     -              NMAIN  =  MCHANG + MCHLEN,
     -              NMATOM =      5,
     -              NTYPE  =     10,
     -              SCHAIN =      2
     -             )

CHECK v.2.2-->
      INTEGER       NPLATM, NPLANE
      PARAMETER    (
     -              NPLATM = 10,
     -              NPLANE = 9
     -             )
CHECK v.2.2<--

CHECK v.2.2-->
C      CHARACTER*3   BONDAT(3,NMAIN), EXTYPE(MAXEXC,NTYPE),
C     -              MATOM(NMATOM), RTYPE(MXTYPE)
      CHARACTER*3   ANAME(NPLATM+1,NPLANE), BONDAT(3,NMAIN),
     -              EXTYPE(MAXEXC,NTYPE), MATOM(NMATOM), RTYPE(MXTYPE)
CHECK v.2.2<--
      CHARACTER*4   BRCODE
      CHARACTER*13  ENGNAM(NMAIN)
      CHARACTER*16  EXDESC(NTYPE)
CHECK v.2.2-->
C      CHARACTER*78  PDBFIL
C      CHARACTER*80  FILLAN, FILNEW
C      COMMON/CHARS/ BONDAT, EXTYPE, MATOM, RTYPE,
CHECK v.3.4.3-->
C      CHARACTER*78  FILLAN, FILNEW, FILPLN, PDBFIL
      CHARACTER*78  FILLAN, FILNEW, FILNRV, FILPLN, PDBFIL
CHECK v.3.4.3<--
      COMMON/CHARS/ ANAME, BONDAT,
     -              EXTYPE, MATOM, RTYPE,
CHECK v.2.2<--

     -              BRCODE,

     -              ENGNAM,

     -              EXDESC,

CHECK v.2.2-->
C     -              PDBFIL,
C     -              FILLAN, FILNEW
CHECK v.3.4.3-->
C     -              FILLAN, FILNEW, FILPLN, PDBFIL
     -              FILLAN, FILNEW, FILNRV, FILPLN, PDBFIL
CHECK v.3.4.3<--
CHECK v.2.2<--

      INTEGER       BATNUM(3,NMAIN), EXCEPT(NMAIN), FSTLET(26),
     -              INAME(NMAIN), LSTLET(26), MCHNUM(NMAIN),
CHECK v.2.2-->
C     -              NBAT(NMAIN), NRESID, NXTLET(MXTYPE), NXTPOS,
C     -              RELPOS(3,NMAIN)
     -              NBAT(NMAIN), NPATOM(NPLANE), NRESID, NXTLET(MXTYPE),
     -              NXTPOS, RELPOS(3,NMAIN)
CHECK v.2.2<--
      COMMON/INTS/  BATNUM, EXCEPT, FSTLET,
     -              INAME, LSTLET, MCHNUM,
CHECK v.2.2-->
C     -              NBAT, NRESID, NXTLET, NXTPOS,
C     -              RELPOS
     -              NBAT, NPATOM, NRESID, NXTLET,
     -              NXTPOS, RELPOS
CHECK v.2.2<--

      LOGICAL       BANGLE(NMAIN), IFAIL
      COMMON/LOGS/  BANGLE, IFAIL

      REAL          ENGMEA(NMAIN), ENGSTD(NMAIN), MCHMAX(NMAIN),
     -              MCHMIN(NMAIN), MCHMEA(NMAIN), MCHSTD(NMAIN), PI,
     -              RADDEG
      COMMON/REALS/ ENGMEA, ENGSTD, MCHMAX,
     -              MCHMIN, MCHMEA, MCHSTD, PI,
     -              RADDEG

CHECK v.3.2-->
      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/
CHECK v.3.2<--

C--------------------------------------------------------------------------
C
C VARIABLE DESCRIPTIONS:-
C ---------------------
C
C Parameters:-
C ----------
C ANAME(NPLATM+1,NPLANE) - Residue-by-residue atom names defining planar
C          groups
C BANGLE(NMAIN) - Identifier which is true if the corresponding elements
C          of BONDAT define an angle rather than a length
C BATNUM(3,NMAIN) - Numerical array corresponding to BONDAT array
C BONDAT(3,NMAIN) - Array defining the atoms involved in each of the
C          main-chain bond-lengths and bond angles to be calculated
C BRCODE - 4-letter Brookhaven code for the protein being processed
C ENGNAM(NMAIN) - X-plor/Engh & Huber name for bond length/angle
C ENGMEA(NMAIN) - Engh & Huber (small molec) mean for bond length/angle
C ENGSTD(NMAIN) - Engh & Huber standard deviation for bond length/angle
C EXCEPT(NMAIN) - Pointer to exception type, EXTYPE, for this bond length
C          or angle
C EXDESC(NTYPE) - Description corresponding to EXTYPE
C EXTYPE(MAXEXC,NTYPE) - Exception categories for the bond lengths and angles.
C          The first element of each condition is one of:
C              ALL  - all residues
C             space - only residues listed
C              NOT  - excluding the residues listed
C FILNEW - Input file, <filename>.new, being the cleaned up version of the
C          .pdb file
C FILLAN - Output file, <filename>.lan, containing means and standard
C          deviations
C FILNRV - Output file, <filename>.nrv, containing the "restraint
C          violations" for use with viol2pdb
C FILPLN - Output file, <filename>.pln, containing transformed coords of all
C          planar groups
C FSTLET(26) - Pointer to RTYPE array for each letter of the alphabet. The
C          residue name stored in RTYPE at the given position is the first
C          residue starting with that particular letter of the alphabet
C INAME(NMAIN) - Marker defining which of the spanned residues is the
C          one whose name is used here
C LSTLET(26) - As for FSTLET, but pointing to the last occurrence of the
C          given letter in the RTYPE array
C MATOM(NMATOM) - Names of the main-chain atoms
C MAXEXC - Number of exception conditions held in the EXCEPT table for each
C          exception type
C MXAMIN - Number of residue types (ie 20 amino acids. Can be reduced for
C          testing purposes)
C MCHAIN - Parameter identifying main chain
C MCHANG - Number of main-chain bond angle types
C MCHLEN - Number of main-chain bond length types
C MCHMAX(NMAIN) - Maximum value for each type of main-chain bond length and
C          angle
C MCHMIN(NMAIN) - Minimum value for each type of main-chain bond length and
C          angle
C MCHMEA(NMAIN) - Mean value for each type of main-chain bond length and
C          angle
C MCHNUM(NMAIN) - Number of values for each  type of main-chain bond
C          length and angle
C MCHSTD(NMAIN) - Standard deviation for each type of main-chain bond
C          length and angle
C MXITEM - Maximum number of length/angles that can be stored for
C          Calculating standard deviations
C MXONE  - Maximum number of residues allowed in single protein
C MXTYPE - Maximum number of residue types
C NBAT(NMAIN) - Number of atoms in current bond length (2) or angle (3)
C NMAIN  - Number of main-chain bond length and bond angle types
C NMATOM - Number of main-chain atoms
C NPATOM(NPLANE) - Number of atoms defining each planar group
C NPLANE - Number of defined planar groups
C NPLATM - Maximum number of atoms defining a planar group
C NRESID - Number of residues in protein
C NTYPE  - Number of exception types defining which residues which lengths
C          and angles refer to
C NXTLET(MXTYPE) - Pointer to next entry in residue types list, RTYPE, for
C          next amino acid starting with the same first letter
C NXTPOS - Next free slot in the RTYPE array
C PDBFIL - Name of original .pdb file holding coords of structure
C PI     - Value of pi
C RADDEG - Conversion factor between radians and degrees
C RELPOS(3,NMAIN) - Relative residue positions of the 2nd and 3rd atoms
C          in current bond length/angle
C RTYPE(MXTYPE) - Residue names, starting with the common 20, and adding
C          any new ones encountered in the given structure
C SCHAIN - Parameter identifying side chain
C
C--------------------------------------------------------------------------
