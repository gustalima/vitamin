C----------------------------------------------------------------------+--- 
C GFAC2PDB.INC - Include variables for GFAC2PDB.FOR
C----------------------------------------------------------------------+--- 

CVAX      IMPLICIT NONE

      INTEGER       MXATOM, MXRES
      PARAMETER    (
     -              MXATOM  = 100000,
     -              MXRES   =  10000
     -             )

      INTEGER       NAMINO, NCOVAL, NDIHED, NDISTR
      PARAMETER    (
     -              NAMINO  =   20,
     -              NCOVAL =     2,
     -              NDIHED =     5,
     -              NDISTR =    NCOVAL + NDIHED + 3
     -             )

CHECK v.3.6.3-->
      INTEGER       NFSTEP, NGSTEP
      PARAMETER    (
     -              NFSTEP  =    4,
     -              NGSTEP  =    6
     -             )

      REAL          FMIN, FSIZE, FSTEP, GMIN, GSTEP
      PARAMETER    (
     -              FMIN    =  0.2,
     -              FSTEP   =  0.5,
     -              FSIZE   =  0.1,
     -              GMIN    = -2.0,
     -              GSTEP   =  0.2
     -             )
CHECK v.3.6.3<--
      

      CHARACTER*1   CHAIN
      CHARACTER*3   AMINO(NAMINO+1), VALRES(MXRES)
      CHARACTER*6   RESDET(MXRES), RESNO(MXATOM)
      CHARACTER*9   VALSNO(MXRES)
      CHARACTER*13  ATOMID(MXATOM)
      CHARACTER*40  TITLE
CHECK v.3.6.3-->
C      CHARACTER*80  FILOUT, FILNEW, FILSDH, PDBFIL
      CHARACTER*80  FILOUT, FILNEW, FILRAS, FILSDH, PDBFIL
CHECK v.3.6.3<--
      COMMON/CHARS/
     -              CHAIN,

     -              AMINO, VALRES,

     -              RESDET, RESNO,

     -              VALSNO,

     -              ATOMID,

     -              TITLE,

CHECK v.3.6.3-->
C     -              FILOUT, FILNEW, FILSDH, PDBFIL
     -              FILOUT, FILNEW, FILRAS, FILSDH, PDBFIL
CHECK v.3.6.3<--

CHECK v.3.6.3-->
C      INTEGER       ACONST, ENASTR, FSTATM(MXRES), FSTRES, INTRLX,
C     -              MODEL, NATOMS, NCOUNT(NAMINO+1), NRESID,
C     -              RESNUM(MXATOM)
C
C      COMMON/INTS/
C     -              ACONST, ENASTR, FSTATM, FSTRES, INTRLX,
C     -              MODEL, NATOMS, NCOUNT, NRESID,
C     -              RESNUM
C----------------------------------------------------------------------+---
      INTEGER       ACONST, ATOMNO(MXATOM), ENASTR, FSTATM(MXRES),
     -              FSTRES, INTRLX, MODEL, NATOMS, NCOUNT(NAMINO+1),
     -              NRESID, RESNUM(MXATOM)

      COMMON/INTS/
     -              ACONST, ATOMNO, ENASTR, FSTATM,
     -              FSTRES, INTRLX, MODEL, NATOMS, NCOUNT,
     -              NRESID, RESNUM
CHECK v.3.6.3<--

      LOGICAL       IFAIL
      COMMON/LOGS/
     -              IFAIL

      REAL          ATMXYZ(3,MXATOM), BVALUE(MXATOM), GFACTR(NDISTR),
     -              SFACT
      COMMON/REALS/
     -              ATMXYZ, BVALUE, GFACTR,
     -              SFACT

      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/

C----------------------------------------------------------------------+--- 
C Parameters
C ----------
C
C ACONST - Total number of dihedral angle constraints
C AMINO(NAMINO) - Three-letter residue code
C ATOMID(MXATOM) - Details uniquely identifying each atom in the structure
C          (ie atom name, chain and residue number)
C ATOMNO(MXATOM) - Stored atom numbers
C ATMXYZ(3,MXATOM) - Atomic coords of each atom
C BOTHND - Flag indicating that restraints are only to be included if
C          both ends are on residues that are included in the
C          residue-selection
C BVALUE(MXATOM) - Atomic temperature factors
C CHAIN  - Selected chain (not used here)
C FILOUT - Output PDB format file containing the G-factors translated into
C          atomic G-factors
C FILNEW - File containing cleaned up version of PDB file
C FILRAS - Output file containing RasMol script colouring structure by G-factor
C FILSDH - File containing all the residue-by-residue G-factor values
C FMIN   - Minimum thickness of G-factor bonds
C FSTATM(MXRES) - Starting atom in each residue
C FSIZE  - Thickness increment for G-factor bonds
C FSTEP  - G-factor range size for bond thickness increments
C FSTRES - Pointer to the first residue of the sequence
C GFACTR(NDISTR) - G-factors giving measures of the stereochemical quality
C          of different aspects of the protein's structure
C GMIN   - Minimum G-factor below which torsion angle G-factors are
C          highlighted
C GSTEP  - Range size for G-factors below minimum corresponding to
C          different highlight colours
C IFAIL  - Failure flag causing program to abort due to some error
C          condition
C MODALL - Flag indicating whether all models in the ensemble are required
C MODEL  - Selected model number, where only a single model to be
C          processed
C MXATOM - Maximum number of atoms allowed for the structure
C MXRCON - Maximum number of distance restraints per residue
C MXRES  - Maximum number of residues allowed for the protein structure
C          being plotted
C MXTORS - Maximum number of torsion-angle restraints per residue
C NATOMS - Number of atoms in the protein
C NDISTR - Total number of log-odds scores in the .sdh file
C NFSTEP - Number of bond thickness increments
C NGSTEP - Number of G-factors increments
C NRESID - Number of residues in the protein
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C RESDET(MXRES) - Residue details for each residue
C RESNO(MXATOM) - Residue-number of each atom
C RESNUM(MXATOM) - Sequential residue-number of each atom
C RESOPT - Restraint selection option
C SFACT  - Stereochemical R-factor, as calculated by program TPLOT and
C          held at the end of the .sdh file
C VALRES(MXRES) - Residue name for each residue in the sequence
C VALSNO(MXRES) - Sequence number for each residue in the sequence
C
C---------------------------------------------------------------------------
