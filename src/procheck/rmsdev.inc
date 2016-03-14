C**************************************************************************
C
C  RMSDEV.INC   -   Include file for RMSDEV.FOR
C
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INTEGER       MXAMIN, MXATOM, MXCONS, MXRES, MXTORS, NCTYPE,
     -              NMATOM
      PARAMETER    (
     -              MXAMIN =     20,
CHECK v.3.5.1-->
C     -              MXATOM =   8000,
     -              MXATOM =  50000,
CHECK v.3.5.1<--
     -              MXCONS =     10,
CHECK v.3.5.1-->
C     -              MXRES  =    800,
     -              MXRES  =   2000,
CHECK v.3.5.1<--
     -              MXTORS =      3,
     -              NCTYPE =      6,
     -              NMATOM =      4
     -             )

      INTEGER       EQUAL, LBOUND, UBOUND
      PARAMETER    (
     -              EQUAL  =      0,
     -              LBOUND =     -1,
     -              UBOUND =      1
     -             )

      CHARACTER*1   HAMARK(MXATOM)
      CHARACTER*4   BRCODE
      CHARACTER*9   RESDET(MXRES)
      CHARACTER*14  ATOMID(MXATOM)
      CHARACTER*78  FILAVE, FILCNS, FILCON, FILNEW, FILRMS, PDBFIL
      COMMON/CHARS/ HAMARK,

     -              BRCODE,

     -              RESDET,

     -              ATOMID,

     -              FILAVE, FILCNS, FILCON, FILNEW, FILRMS, PDBFIL

      INTEGER       ATRESN(MXATOM), CONATM(2,MXCONS), CONBOU(MXCONS),
     -              CONTYP(MXCONS), FSTATM, LSTCON(2,MXRES), MAXLEN,
     -              MINLEN, MVIOL, NATALL, NATOMS, NCONST,
     -              NCOORD(MXATOM), NDIHED, NINTRA, NMALL, NMODEL,
     -              NRESID, NRMS(3,MXRES), NRMALL(3,MXRES), NSALL,
     -              NUMCON, NVIOL(MXCONS), NXTATM(MXATOM),
     -              NXTCON(2,MXCONS), TVIOL

      COMMON/INTS/  ATRESN, CONATM, CONBOU,
     -              CONTYP, FSTATM, LSTCON, MAXLEN,
     -              MINLEN, MVIOL, NATALL, NATOMS, NCONST,
     -              NCOORD, NDIHED, NINTRA, NMALL, NMODEL,
     -              NRESID, NRMS, NRMALL, NSALL,
     -              NUMCON, NVIOL, NXTATM,
     -              NXTCON, TVIOL

      LOGICAL       DOUBL(MXATOM), ENDFIL, ENSEMB, GOTFST(MXATOM),
     -              HAVATM(MXATOM), HAVCON, IFAIL
      COMMON/LOGS/  DOUBL, ENDFIL, ENSEMB, GOTFST,
     -              HAVATM, HAVCON, IFAIL

      REAL          ACDIST(MXCONS), CNSTRN(MXCONS),
     -              DIHCON(2,MXTORS,MXRES), FCOORD(3,MXATOM), MALLDV,
     -              MCOORD(3,MXATOM), RMS(3,MXRES), RMSALL(3,MXRES),
     -              SALLDV, SUMALL, XYZ(3,MXATOM)

      COMMON/REALS/ ACDIST, CNSTRN,
     -              DIHCON, FCOORD, MALLDV,
     -              MCOORD, RMS, RMSALL,
     -              SALLDV, SUMALL, XYZ

C----------------------------------------------------------------------+---
C
C VARIABLE DESCRIPTIONS:-
C ---------------------
C ACDIST(MXCONS) - Actual distance corresponding to each constraint
C ATOMID(MXATOM) - Details uniquely identifying each atom in the structure
C          (ie atom name, chain and residue number)
C ATRESN(MXATOM) - Sequential residue number of residue to which atom
C          belongs
C BRCODE - 4-letter Brookhaven code for the protein being processed
C CNSTRN(MXCONS) - Value of the distance constraint
C CONATM(2,MXCONS) - Two atoms defining the distance constraint
C CONBOU(MXCONS) - Marker indicating the constraint type - whether it
C          corresponds to an upper- or lower-bound value
C CONTYP(MXCONS) - Identifier defining the type of restraint this is
C DIHCON(2,MXTORS,MXRES) - Dihedral angle constraints (lower and upper
C          bounds)
C DOUBL(MXATOM) - Flags indicating apparently common atoms in different
C          residue types in the different models
C ENDFIL - Flag indicating whether end of file-list file has been reached
C ENSEMB - Flag indicating whether list of filenames contained in
C          supplied filename
C FCOORD(3,MXATOM) - Stored coordinates of first model
C FILAVE - Output PDB file, <filename>.ave, giving the mean coords of each
C          atom in the structure
C FILCNS - Output file of the distance constraints, together with the
C          actual values for each model
C FILCON - Input NMR constraints file, <filename>.mr
C FILNEW - Input file, <filename>.new, being the cleaned up version of the
C          .pdb file
C FILRMS - Output file, <filename>.rms, containing the residue-by-residue
C          RMS deviations from the mean coords
C FSTRES - Pointer to the first atom of the ensemble sequence
C GOTFST(MXATOM) - Flag indicating which of all the atoms are in the
C          first model (and hence whose coordinates are stored in FCOORD,
C          used when there are only 2 models and pairwise RMS deviations
C          are computed, rather than deviations from the mean coords)
C HAMARK(MXATOM) - Marker showing whether atom comes from an ATOM or a
C          HETATM record
C HAVATM(MXATOM) - Flag indicating whether a given atom is present in the
C          current model
C HAVCON - Flag indicating whether the constraints file exists
C IFAIL  - Flag indicating an error encountered in the program
C LSTCON(2,MXRES) - Position of last constraint involving this residue
C          (for 1st and 2nd atoms of the constraint)
C MALLDV - RMS dev for all main-chain atoms in all structures
C MAXLEN - Longest sequence in the ensemble in terms of the number of
C          atoms it has
C MCOORD(3,MXATOM) - Mean coordinates for each atom in the structure
C MINLEN - Shortest sequence in the ensemble in terms of the number of
C          atoms it has
C MVIOL  - Number of constraint violations in the current model
C MXATOM - Maximum number of atoms allowed for the structure
C MXCONS - Maximum number of NMR constraints that can be stored
C MXRES  - Maximum number of residues in structure
C MXTORS - Maximum number of dihedral angle constraints
C NATALL - Number of all atoms (main-chain plus sidechain) in all structures
C NATOMS - Number of atoms in the protein
C NCONST - Number of distance constraints read in from the constraints
C          file
C NCOORD(MXATOM) - Number of instances of each atom encountered in the
C          ensemble
C NCTYPE - Number of constraint types
C NDIHED - Number of dihedral angle constraints read in
C NINTRA - Number of intra-residue constraints
C NMALL  - Number of all main-chain atoms in all structures
C NMATOM - Number of main-chain atoms
C NMODEL - Number of models making up the ensemble, encountered in the
C          .new file
C NRESID - Number of residues in protein
C NRMS(3,MXRES) - Counts of data items contributing to RMS(3,MXRES)
C NRMALL(3,MXRES) - Counts of data items contributing to RMSALL(3,MXRES)
C NSALL  - Number of all sidechain atoms in all structures
C NUMCON - Sequential numbering of constraints, in residue order
C NVIOL(MXCONS)  - Number of models in which each constraint is violations
C NXTATM(MXATOM) - Pointer to next atom in sequence (necessary to cope
C          with situation where have slight differences in the different
C          models in the ensemble
C NXTCON(2,MXRES) - Next constraint involving this residue
C          (for 1st and 2nd atoms of the constraint)
C PDBFIL - Name of original .pdb file holding coords of structure
C RESDET(MXRES) - Identifying details for each residue
C RMS(3,MXRES) - RMS deviation from mean coordinates of each residue
C          in a given model for main-chain, side-chain and overall atoms
C RMSALL(3,MXRES) - As RMS above, but for the whole ensemble of structures
C SALLDV - RMS dev for all sidechain atoms in all structures
C SUMALL - RMS dev for all atoms in all structures
C TVIOL  - Total number of constraint violations
C XYZ(3,MXATOM) - Atom coords for current model - used for calculating
C          constraint violations
C
C--------------------------------------------------------------------------
