C----------------------------------------------------------------------+--- 
C VIOL2PDB.INC - Include variables for VIOL2PDB.FOR
C----------------------------------------------------------------------+--- 

CVAX      IMPLICIT NONE

      INTEGER       MAXARR, MXATOM, MAXRNG, MXCONS, MXMODL, MXRES,
     -              MXTORS
      PARAMETER    (
     -              MAXARR  = 100000,
     -              MXATOM  = 100000,
     -              MAXRNG  =    200,
     -              MXCONS  =  10000,
     -              MXMODL  =     60,
     -              MXRES   =   1000,
     -              MXTORS  =      4
     -             )

      INTEGER       NAMINO, NCTYPE
      PARAMETER    (
     -              NAMINO  =   20,
     -              NCTYPE  =    3
     -             )

      CHARACTER*1   RESTYP(MXCONS)
      CHARACTER*3   AMINO(NAMINO+1), VALRES(MXRES)
CHECK v.3.4.3-->
      CHARACTER*6   RESFRM(MAXRNG), RESNO(MXATOM), RESTO(MAXRNG)
CHECK v.3.4.3<--
      CHARACTER*9   VALSNO(MXRES)
      CHARACTER*13  ATOMID(MXATOM), ATOMRS(2,MXCONS)
      CHARACTER*40  TITLE
      CHARACTER*80  BRCALL, BRCODE(MXMODL), FILOUT(2,NCTYPE),
CHECK v.3.4.3-->
C     -              FILPSU, FILTOR, FILVIO(NCTYPE), PDBFIL
     -              FILPSU, FILRNG, FILTOR, FILVIO(NCTYPE), PDBFIL
CHECK v.3.4.3<--
      COMMON/CHARS/
     -              RESTYP,

     -              AMINO, VALRES,

CHECK v.3.4.3-->
     -              RESFRM, RESNO, RESTO,
CHECK v.3.4.3<--

     -              VALSNO,

     -              ATOMID, ATOMRS,

     -              TITLE,

     -              BRCALL, BRCODE, FILOUT,
CHECK v.3.4.3-->
C     -              FILPSU, FILTOR, FILVIO, PDBFIL
     -              FILPSU, FILRNG, FILTOR, FILVIO, PDBFIL
CHECK v.3.4.3<--

CHECK v.3.4.3-->
C      INTEGER       ACONST, CONNEC(3,MXCONS),
      INTEGER       ACONST, ACTNUM(MXMODL), CONNEC(3,MXCONS),
CHECK v.3.4.3<--
     -              DVIOL(MXTORS,MXRES), ENASTR,
     -              FSTCON(MXRES), FSTRES,
     -              INTRLX, LSTATM(MXMODL), MODEL,
     -              MODFRM(MAXRNG), MODNUM(MXATOM), MODTO(MAXRNG),
     -              MRANGE, NATOMS, NCONEC, NCONST(NCTYPE),
     -              NCOUNT(NAMINO+1),
CHECK v.3.4.3-->
C     -              NMODEL, NMRMOD(MXMODL), NRES,
     -              NMODEL, NMRMOD(MXMODL), NRANGE, NRES,
CHECK v.3.4.3<--
     -              NXTRES(MXRES), RESNUM(MXATOM),
CHECK v.3.4.3-->
C     -              SSTCNT(3,MXRES), TCONST, TLEN, VIOATM(2,MXCONS)
     -              SSTCNT(3,MXRES), TCONST, TLEN, TOPMOD,
     -              VIOATM(2,MXCONS)
CHECK v.3.4.3<--

      COMMON/INTS/
CHECK v.3.4.3-->
C     -              ACONST, CONNEC,
     -              ACONST, ACTNUM, CONNEC,
CHECK v.3.4.3<--
     -              DVIOL, ENASTR,
     -              FSTCON, FSTRES,
     -              INTRLX, LSTATM, MODEL,
     -              MODFRM, MODNUM, MODTO, MRANGE,
     -              NATOMS, NCONEC, NCONST,
     -              NCOUNT,
CHECK v.3.4.3-->
C     -              NMODEL, NMRMOD, NRES,
     -              NMODEL, NMRMOD, NRANGE, NRES,
CHECK v.3.4.3<--
     -              NXTRES, RESNUM,
CHECK v.3.4.3-->
C     -              SSTCNT, TCONST, TLEN, VIOATM
     -              SSTCNT, TCONST, TLEN, TOPMOD,
     -              VIOATM
CHECK v.3.4.3<--

CHECK v.3.4.3-->
C      LOGICAL       ENSEMB, IFAIL, INSCRP, MODALL, MODFST, NMR, WLONG,
C     -              WMEDIU, WSHORT
C      COMMON/LOGS/
C     -              ENSEMB, IFAIL, INSCRP, MODALL, MODFST, NMR, WLONG,
C     -              WMEDIU, WSHORT
      LOGICAL       BOTHND, ENSEMB, HAVRAN, IFAIL, INSCRP, MODALL,
     -              MODFST, MWANT(MXMODL), NMR, RSELEC, WLONG, WMEDIU,
     -              WSHORT
      COMMON/LOGS/
     -              BOTHND, ENSEMB, HAVRAN, IFAIL, INSCRP, MODALL,
     -              MODFST, MWANT, NMR, RSELEC, WLONG, WMEDIU,
     -              WSHORT
CHECK v.3.4.3<--

      REAL          ATMXYZ(3,MXATOM), DVIOLM(MXTORS,MXRES,MXMODL)
      COMMON/REALS/
     -              ATMXYZ, DVIOLM

      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/

C----------------------------------------------------------------------+--- 
C Parameters
C ----------
C
C ACONST - Total number of dihedral angle constraints
C ACTNUM(MXMODL) - Actual model number corresponding to every entry in
C          MODNUM(MXMODL)
C AMINO(NAMINO) - Three-letter residue code
C ATOMID(MXATOM) - Details uniquely identifying each atom in the structure
C          (ie atom name, chain and residue number)
C ATMXYZ(3,MXATOM) - Atomic coords of each atom
C ATOMRS(2,MXCONS) - Atoms involved in each distance restraint
C BOTHND - Flag indicating that restraints are only to be included if
C          both ends are on residues that are included in the
C          residue-selection
C BRCALL - Brookhaven code for ensemble
C CONNEC(3,MXCONS) - Connectivity records for the bonds representing
C          the restraints and restraint violations
C DESC(NQUALP) - Descriptions of the stereochemical parameters plotted on
C          the Quality Plots
C DVIOL(MXTORS,MXRES) - Number of constraint violations for the dihedral
C          angle constraints
C DVIOLM(MXTORS,MXRES,MXMODL) - Dihedral-angle violation for each residue
C          in each model
C FILOUT(2,NCTYPE) - Output filenames of various output PDB files
C FILPSU - Filename for pseudo-atoms file, created by AQUA.
C FILTOR - Filename for torsion angle constraints file, created by AQUA.
C FILVIO(NCTYPE) - Input violation files generated by AQUA
C FSTCON(MXRES) - Pointer to first constraint for the given residue
C FSTRES - Pointer to the first residue of the sequence
C HAVRAN - Flag indicating whether have a ranges file
C IFAIL  - Failure flag causing program to abort due to some error
C          condition
C INSCRP - Flag indicating that program running within script file, so
C          requests for input data not required
C LSTATM(MXMODL) - Last atom-number for each NMR model
C MAXARR - Maximum size of constraints array
C MAXRNG - Maximum number of model-number ranges for selection
C MODALL - Flag indicating whether all models in the ensemble are required
C MODEL  - Selected model number, where only a single model to be
C          processed
C MODFRM(MAXRNG), MODTO(MAXRNG) - Starting and ending model-number for each
C          selected model-number range
C MODFST - Flag indicating that only the first model of the ensemble is
C          required
C MODNUM(MXATOM) - Model to which atom belongs (true model numbers, as in
C          PDB file, not the sequential model-count)
C MRANGE - Number of selected model-number ranges
C MWANT(MXMODL) - Flags indicating which models are required from the
C          original list in mplot.in
C MXATOM - Maximum number of atoms allowed for the structure
C MXCONS - Maximum number of restraints that can be stored per model
C MXRCON - Maximum number of distance restraints per residue
C MXMODL - Maximum number of models allowed in the ensemble
C MXRES  - Maximum number of residues allowed for the protein structure
C          being plotted
C MXTORS - Maximum number of torsion-angle restraints per residue
C NATOMS - Number of atoms in the protein
C NCONEC - Number of connectivity records stored
C NCONRS(MXRES) - Number of constraints per residue
C NCONST(NCTYPE) - Total number of distance constraints of each type picked
C          up from the .rvio file
C NMODEL - Number of model structures in the ensemble
C NMR    - Flag indicating that the ensemble is a set of NMR models coming
C          from a single file
C NMRMOD(MXMODL) - Model numbers of models as defined by the MODEL record
C NRANGE - Number of residue-number ranges in the ranges file
C NRES   - Number of residues in the protein
C NVIOL(MXCONS) - Number of violations of each constraint
C NXTRES(MXRES) - Pointer to next residue in sequence (necessary to cope
C          with situation where have slight differences in the sequences
C          of the .rin files
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C RESFRM(MAXRNG), RESTO(MAXRNG) - Residue-number ranges read in from the
C          ranges file
C RESMOD(MXRES,MXMODL) - Flags indicating which residues are present in
C          each model (necessary for insertions and deletions)
C RESNO(MXATOM) - Residue-number of each atom
C RESNUM(MXATOM) - Sequential residue-number of each atom
C RESOPT - Restraint selection option
C RESTYP(MXCONS) - Type of each restraint (L=Lower,U=Upper,B=Both)
C TCONST - Total number of restraint violations read in
C TOPMOD - Highest model-number in the ensemble
C VALRES(MXRES) - Residue name for each residue in the sequence
C VALSNO(MXRES) - Sequence number for each residue in the sequence
C VIOATM(2,MXCONS) - Atom-pair making up each restraint
C WLONG  - Flag indicating that long-range restraints required
C WMEDIU - Flag indicating that medium-range restraints required
C WSHORT - Flag indicating that short-range restraints required
C
C---------------------------------------------------------------------------
