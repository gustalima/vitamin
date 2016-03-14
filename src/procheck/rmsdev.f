C**************************************************************************
C
C  RMSDEV.FOR   -   Program to calculate RMS deviations between NMR models
C                   of a protein in a PDB file
C
C                   The program reads the protein structure in from the
C                   user-defined .new file and outputs the results to a
C                   .rms file.
C
C  Written by R A Laskowski, University College, London, May 1994.
C
C  Original version was part of v.3.2 of the PROCHECK suite of programs.
C  Subsequent amendments have been labelled by CHECK v.m.n--> and
C  CHECK v.m.n<-- where m.n is the version number corresponding to the
C  change
C
C  v.3.2.1 - Change wherein if an ensemble comprises only 2 model
C            structures, the RMS differences between the two structures
C            are computed, rather than the RMS differences from the mean
C            coordinates.
C            Amendment so that error message is not printed if there's a
C            blank line in the filelist file.
C                                             Roman Laskowski (30 Nov 1994)
C  v.3.4   - Bug-fix on atom-search when any model contains the last two
C            atoms the other way round relative to the first model.
C                                             Roman Laskowski (11 Apr 1996)
C            Error message altered to a warning message when no MODEL
C            records found in the file
C  v.3.5.1 - Increase of MXRES from 800 residues to 2000, and MXATOM
C            from 8000 to 50000.
C                                              Roman Laskowski ( 2 Mar 1998)
C
C----------------------------------------------------------------------+---
C
C Files
C -----
C
C 1  <filename>.new - Cleaned-up version of the .pdb file holding the protein
C                     structure in Brookhaven format by-residue information
C                     on the given structure
C 2  <filename>.mr  - Constraints file holding the NMR distance constraints
c 3  <filelist>     - For ensembles, name of file holding list of PDB files
C                     to be processed
C 7  <filename>.rms - Output file holding the rms deviation of each residue
C                     in each model from the mean position for that residue.
C                     The deviations are shown for the main-chain, side-chain
C                     and total atoms in the residue.
C 8  <filename>.ave - Output file holding the mean coordinates for each atom
C                     in the structure
C 9  <filename>.cns - Output file holding the distance constraints and actual
C                     distances for each model
C
C--------------------------------------------------------------------------
C
C Subroutine calling tree
C -----------------------
C
C MAIN    --> GETCOD  --> GETNAM
C         --> OPNFIL
C         --> INITS
C         --> OPENXT
C         --> READPR
C         --> GETRES
C         --> WRIAVE
C         --> GETCON
C         --> OPENXT
C         --> CALRMS  --> VIOLA
C                     --> PUTVIO
C         --> ALLMOD
C         --> PUTCON
C
C----------------------------------------------------------------------+---


      PROGRAM RMSDEV

      INCLUDE 'rmsdev.inc'

      LOGICAL       FIRST

C---- Read in the filename of the Brookhaven file
      CALL GETCOD
      IF (IFAIL) GO TO 990

C---- Open data files
      CALL OPNFIL
      IF (IFAIL) GO TO 990

C---- Initialise variables
      CALL INITS
      FIRST = .TRUE.

C---- If have an ensemble, loop through the filename to be processed
 100  CONTINUE

C----     For ensemble, read in the next filename and open
          IF (ENSEMB) THEN
              CALL OPENXT(FIRST)
          ENDIF

C----     Read through the .new file to accumulate mean coords
          IF (.NOT.ENDFIL .AND. .NOT.IFAIL) THEN
              CALL READPR
CHECK v.3.4-->
C              IF (IFAIL) GO TO 990
              IF (IFAIL .OR. NMODEL.EQ.0) GO TO 990
CHECK v.3.4<--
          ENDIF
      IF (ENSEMB .AND. .NOT.ENDFIL) GO TO 100
      IF (ENSEMB) THEN
          PRINT*, 'Total no. of unique atoms   ', NATOMS
      ENDIF

C---- Assign a sequential residue number to each atom record
      CALL GETRES
      IF (IFAIL) GO TO 990
      IF (ENSEMB) THEN
          PRINT*, 'Total no. of unique residues', NRESID
      ENDIF

C---- Calculate and write out the averaged atomic coordinates
C     of the structure
      CALL WRIAVE

C---- Read in the distance restraints from the .mr file
      IF (HAVCON) THEN
          CALL GETCON
          IF (IFAIL) THEN
              HAVCON = .FALSE.
              IFAIL = .FALSE.
          ENDIF
      ENDIF

C---- If have an ensemble, rewind the file-list file
      IF (ENSEMB) THEN
          REWIND(3)
          ENDFIL = .FALSE.
      ENDIF

C---- Read through the .new file(s) a second time to calculate the RMS
C     deviations from the mean coordinates

C---- If have an ensemble, loop through the filename to be processed
 200  CONTINUE

C----     For ensemble, read in the next filename and open
          IF (ENSEMB) THEN
              CALL OPENXT(FIRST)
          ENDIF

C----     Read through the .new file to calculate RMS deviations
          IF (.NOT.ENDFIL) THEN
              CALL CALRMS
          ENDIF
      IF (ENSEMB .AND. .NOT.ENDFIL) GO TO 200

C---- Write out the overall RMS differences
      CALL ALLMOD

C---- Write out the constraints, and counts of constraint violations,
C     to the output .cns file
      IF (HAVCON) THEN
          CALL PUTCON
      ENDIF

C---- Close any open files
      CLOSE(1)
      IF (ENSEMB) CLOSE(3)
      CLOSE(7)
      CLOSE(8)

 990  CONTINUE
      IF (IFAIL) THEN
          PRINT*, '* Program RMSDEV terminated with error'
      ELSE
          PRINT*, '* Program completed'
      ENDIF
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETCOD  -  Read in the filename of the PDB file holding the
C                        ensemble of NMR structures
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCOD

      INCLUDE 'rmsdev.inc'

      INTEGER       IEND, ILEN, ISTART
      LOGICAL       IERROR

C---- Initialise variables
      ENSEMB = .FALSE.
      IFAIL = .FALSE.

C---- Accept name of original .pdb file holding the structure
      PRINT*, 'Enter filename containing coordinates of structure'
      PRINT*, '  (or, for a list of files, enter %filelist, where',
     -    ' filelist'
      PRINT*, '   contains the PDB filenames)'
      READ(*,10) PDBFIL
10    FORMAT(A)
      IF (PDBFIL(1:1).EQ.'%') THEN
          ENSEMB = .TRUE.
          PDBFIL = PDBFIL(2:)
      ENDIF

C---- Peel off directory path and extension
      CALL GETNAM(PDBFIL,ISTART,IEND,IERROR)
      IF (IERROR) GO TO 990

C---- Form names of other files that will be required in default directory
      ILEN = IEND - ISTART + 1
      FILAVE = PDBFIL(ISTART:IEND) // '.ave'
      FILCNS = PDBFIL(ISTART:IEND) // '.cns'
      FILCON = PDBFIL(ISTART:IEND) // '.mr'
      FILNEW = PDBFIL(ISTART:IEND) // '.new'
      FILRMS = PDBFIL(ISTART:IEND) // '.rms'

      GO TO 999

C---- Fatal errors
990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE GETNAM  -  Peel off the directory path and extension from the
C                        full name of the .pdb file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETNAM(PDBFIL,ISTART,IEND,IERROR)

      CHARACTER*1   PCHAR
      CHARACTER*78  PDBFIL
      INTEGER       IEND, IPOS, ISTART, ISTATE
      LOGICAL       FINISH, GOTDOT, IERROR

C---- Initialise variables
      FINISH = .FALSE.
      IEND = 0
      IERROR = .FALSE.
      ISTART = 1
      ISTATE = 1
      IPOS = 78
      GOTDOT = .FALSE.

C---- Check through the filename from right to left
100   CONTINUE

C----     Pick off next character
          PCHAR = PDBFIL(IPOS:IPOS)

C----     State 1: Searching for first non-blank character
          IF (ISTATE.EQ.1) THEN
              IF (PCHAR.EQ.'/' .OR. PCHAR.EQ.'\\' .OR.
     -            PCHAR.EQ.']') THEN
                  GO TO 900
              ENDIF
              IF (PCHAR.NE.' ' .AND. PCHAR.NE.'.') THEN
                  IEND = IPOS
                  ISTATE = 2
              ENDIF

C----     State 2: Searching for end of extension, or end of directory path
          ELSE IF (ISTATE.EQ.2) THEN

C----         If character is a dot, and is the first dot, then note position
              IF (PCHAR.EQ.'.' .AND. .NOT.GOTDOT) THEN
                  IEND = IPOS - 1
                  GOTDOT = .TRUE.

C----         If character signifies the end of a directory path, note pstn
              ELSE IF (PCHAR.EQ.'/' .OR. PCHAR.EQ.'\\'
     -            .OR. PCHAR.EQ.']') THEN
                  ISTART = IPOS + 1
                  FINISH = .TRUE.
              ENDIF
          ENDIF

C----     Step back a character
          IPOS = IPOS - 1

C---- Loop back for next character
      IF (.NOT.FINISH .AND. IPOS.GT.0) GO TO 100

C---- Check whether file name is sensible
      IF (ISTART.GT.IEND) GO TO 900

      GO TO 999

C---- Error in file name
900   CONTINUE
      IEND = 40
      IF (PDBFIL(41:78).NE.' ') IEND = 78
      PRINT*,' *** ERROR in supplied name of file: [', PDBFIL(1:IEND),
     -    ']'
      IERROR = .TRUE.

999   CONTINUE
      RETURN
      END
C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE OPNFIL  -   Open data files
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPNFIL

      INCLUDE 'rmsdev.inc'

C---- Initialise variables
      HAVCON = .TRUE.
      HAVCON = .FALSE.
      IFAIL = .FALSE.

C---- If ensemble, then open file-list
      IF (ENSEMB) THEN
          OPEN(UNIT=3,FILE=PDBFIL,STATUS='OLD',ERR=912,
CVAX     -     READONLY,
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')

C---- Otherwise, open input .new file
      ELSE
          OPEN(UNIT=1,FILE=FILNEW,STATUS='OLD',ERR=900,FORM='FORMATTED',
CVAX     -     READONLY,
     -     ACCESS='SEQUENTIAL')
      ENDIF

C---- Open output file, <filename>.rms
      OPEN(UNIT=7,FILE=FILRMS,STATUS='UNKNOWN',ERR=904,
CVAX     -     CARRIAGECONTROL='LIST',
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')

C---- Open output file, <filename>.ave
      OPEN(UNIT=8,FILE=FILAVE,STATUS='UNKNOWN',ERR=906,
CVAX     -     CARRIAGECONTROL='LIST',
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')

C---- Open output constraints file
CXXX      OPEN(UNIT=9,FILE=FILCNS,STATUS='UNKNOWN',ERR=910,
CXXXCVAX     -     CARRIAGECONTROL='LIST',
CXXX     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')

C---- Open input constraints file
CXXX      OPEN(UNIT=2,FILE=FILCON,STATUS='OLD',ERR=908,FORM='FORMATTED',
CXXXCVAX     -     READONLY,
CXXX     -     ACCESS='SEQUENTIAL')

      GO TO 999

C---- Fatal errors
900   CONTINUE
      PRINT*, '*** ERROR. Unable to open data file: '
      PRINT*, FILNEW, '*'
      PRINT*, '***        Run program CLEAN to create it'
      GO TO 990

904   CONTINUE
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILRMS, '*'
      GO TO 990

906   CONTINUE
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILAVE, '*'
      GO TO 990

 912  CONTINUE
      PRINT*, '*** ERROR. Unable to open file-list:'
      PRINT*, PDBFIL, '*'
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE INITS  -  Initialise variables
C
C----------------------------------------------------------------------+---

      SUBROUTINE INITS

      INCLUDE 'rmsdev.inc'

      INTEGER       I, IATOM, ICONST, IRESID, ITORS

C---- Initialise variables
      ENDFIL = .FALSE.
      IFAIL = .FALSE.
      NATOMS = 0
      NDIHED = 0
      NRESID = 0
      NMODEL = 0

C---- Initialise the pointers to the atoms making up the averaged
C     ensemble structure, and zero all the mean coordinates
      FSTATM = 1
      DO 300, IATOM = 1, MXATOM
          ATOMID(IATOM) = ' '
          ATRESN(IATOM) = 0
          NXTATM(IATOM) = 0
          NCOORD(IATOM) = 0
          FCOORD(1,IATOM) = 0.0
          FCOORD(2,IATOM) = 0.0
          FCOORD(3,IATOM) = 0.0
          GOTFST(IATOM) = .FALSE.
          MCOORD(1,IATOM) = 0.0
          MCOORD(2,IATOM) = 0.0
          MCOORD(3,IATOM) = 0.0
300   CONTINUE
      TVIOL = 0

C---- Zero all the RMS arrays
      DO 500, IRESID = 1, MXRES
          DO 400, I = 1, 3
              NRMS(I,IRESID) = 0
              NRMALL(I,IRESID) = 0
              RMS(I,IRESID) = 0.0
              RMSALL(I,IRESID) = 0.0
400       CONTINUE
          DO 450, ITORS = 1, MXTORS
              DIHCON(1,ITORS,IRESID) = 999.9
              DIHCON(2,ITORS,IRESID) = 999.9
 450      CONTINUE
          RESDET(IRESID) = ' '
500   CONTINUE
      NATALL = 0
      NMALL = 0
      NSALL = 0
      MALLDV = 0.0
      SALLDV = 0.0
      SUMALL = 0.0

C---- Zero the counts of constraint violations
      DO 600, ICONST = 1, MXCONS
          NVIOL(ICONST) = 0
 600  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE OPENXT  -   Get the next .new file in the list and open in
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPENXT(FIRST)

      INCLUDE 'rmsdev.inc'

      CHARACTER*80  IREC
      INTEGER       NAMLEN
      LOGICAL       FIRST

C---- Initialise variables
      IFAIL = .FALSE.

C---- Close the previous file, if open
      IF (.NOT.FIRST) THEN
          CLOSE(1)
      ENDIF

C---- Read in the next record from the file-list
      READ(3,110,ERR=902,END=904) IREC
 110  FORMAT(A)
      FILNEW = IREC

C---- Open input .new file
      NAMLEN = INDEX(FILNEW,' ')
CHECK v.3.2.1-->
      IF (NAMLEN.LT.1 .OR. FILNEW.EQ.' ') GO TO 904
CHECK v.3.2.1<--
      PRINT*
      PRINT*, 'Processing ', FILNEW(1:NAMLEN), '  ...'
      OPEN(UNIT=1,FILE=FILNEW,STATUS='OLD',ERR=900,FORM='FORMATTED',
CVAX     -     READONLY,
     -     ACCESS='SEQUENTIAL')

      GO TO 999

C---- Fatal errors
900   CONTINUE
      PRINT*, '*** ERROR. Unable to open data file: '
      PRINT*, FILNEW, '*'
      IF (.NOT.ENSEMB) THEN
          PRINT*, '***        Run program CLEAN to create it'
      ENDIF
      GO TO 990

 902  CONTINUE
      PRINT*, '*** DATA ERROR reading file-list'
      PRINT*, FILNEW, '*'
      GO TO 990

 904  CONTINUE
      ENDFIL = .TRUE.
      GO TO 999

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE READPR  -  Read through appropriate .new file and
C                        accumulate coords to calculate mean positions
C
C----------------------------------------------------------------------+---

      SUBROUTINE READPR

      SAVE

      INCLUDE 'rmsdev.inc'

      CHARACTER*1   CHAIN, HETMRK
      CHARACTER*3   RESDUE
      CHARACTER*4   ATNAME
      CHARACTER*5   SEQNO
      CHARACTER*6   IDENT, INSEQ
      CHARACTER*14  ATMKEY
      CHARACTER*80  IREC
      INTEGER       I, IATOM, IPOS, JATOM, JPOS, LATOM, LINE, LPOS,
     -              LSTPOS, NMRMOD
      LOGICAL       DUPLIC, FIRST, FIRSTM, GOTPOS, HATOM, NMR, READIN,
     -              WANTED
      REAL          COORDS(3)

      DATA FIRST  / .TRUE. /

C---- Initialise variables
      FIRSTM = .TRUE.
      IFAIL = .FALSE.
      LINE = 0
      IPOS = 0
      LPOS = 0
      IATOM = 0
      NMR = .FALSE.
      NMRMOD = 0
      IF (ENSEMB) THEN
          READIN = .TRUE.
          NMODEL = NMODEL + 1
      ELSE
          NMR = .TRUE.
          READIN = .FALSE.
      ENDIF

C---- Loop through the .new file, looking for all the ATOM and HETATM
C     records
100   CONTINUE
          READ(1,150,END=800,ERR=904) IREC
150       FORMAT(A)
          LINE = LINE + 1
          IDENT = IREC(1:6)

C----     Determine whether this record is required for processing
          WANTED = .FALSE.
          IF (READIN) THEN
              IF (IDENT.EQ.'ATOM  ' .OR. IDENT.EQ.'HETATM') THEN
                  WANTED = .TRUE.
                  IF (IREC(18:20).EQ.'HOH' .OR.
     -                IREC(18:20).EQ.'WAT') WANTED = .FALSE.
              ENDIF
          ENDIF

C----     Process record if required
          IF (WANTED) THEN

C----         Pick off atom name, residue name, chain and sequence number
              ATNAME = IREC(13:16)
              RESDUE = IREC(18:20)
              CHAIN = IREC(22:22)
              INSEQ = IREC(22:27)
              SEQNO = IREC(23:27)

C----         Check if this is a hydrogen atom or pseudo-atom
              IF (ATNAME(2:2).EQ.'H' .OR. ATNAME(2:2).EQ.'Q' .OR.
     -            ATNAME(1:1).EQ.'H' .OR. ATNAME(1:1).EQ.'Q') THEN
                  HATOM = .TRUE.
              ELSE
                  HATOM = .FALSE.
              ENDIF

C----         Retrieve the atomic coordinates
              READ(IREC,180,ERR=902) (COORDS(I), I = 1, 3)
180           FORMAT(30X,3F8.0)

C----         Find the atom's location in the sequence and add its
C             coords to the mean

C----         Set HETATM marker if this is a HETATM record
              IF (IDENT.EQ.'ATOM  ') THEN
                  HETMRK = ' '
              ELSE
                  HETMRK = 'H'
              ENDIF

C----         Form the unique, identifying key for this atom
              ATMKEY = HETMRK // INSEQ // RESDUE // ATNAME
              DUPLIC = .FALSE.
              IPOS = IPOS + 1
              IF (IPOS.GT.MXATOM) GO TO 906

C----         If this is the first instance of this atom, then
C             store it and update the appropriate pointer
              IF (FIRST) THEN
                  IATOM = IPOS
                  IF (LPOS.GT.0) THEN
                      NXTATM(LPOS) = IPOS
                  ENDIF
                  LPOS = IPOS

C----         Otherwise, check whether have met this atom in a
C             previous model
              ELSE
                  GOTPOS = .FALSE.
                  IATOM = 0
                  JATOM = 0
                  LPOS = -1
                  DO 400, JPOS = 1, NATOMS

C----                 Save previous atom
                      LATOM = JATOM

C----                 Get the pointer to the next atom
                      IF (JPOS.EQ.1) THEN
                          JATOM = FSTATM
                      ELSE
                          JATOM = NXTATM(JATOM)
                      ENDIF
                      LSTPOS = JATOM

C----                 Check whether current atom fits in here
                      IF (ATMKEY.EQ.ATOMID(JATOM)) THEN
                          IATOM = JATOM
                          LPOS = JATOM
                          GO TO 410

C----                 If atom doesn't match exactly, see whether the
C                     residue number is right, but just the residue type
C                     that is different
                      ELSE
                          IF (ATMKEY(1:7).EQ.ATOMID(JATOM)(1:7)) THEN
                              LPOS = JATOM
                              DUPLIC = .TRUE.
                              DOUBL(JATOM) = .TRUE.
                              GOTPOS = .TRUE.
                          ELSE IF (ATMKEY(1:7).LT.ATOMID(JATOM)(1:7)
     -                        .AND. .NOT.GOTPOS) THEN
                              LPOS = LATOM
                              GOTPOS = .TRUE.
                          ELSE IF (GOTPOS) THEN
                              GO TO 410
                          ENDIF
                      ENDIF
 400              CONTINUE
 410              CONTINUE

C----             If existing slot has not been been found, then put
C                 atom at end of list and adjust pointers accordingly
C                 so that it "slots in" between the last matched atom
C                 and the one that follows it
                  IF (IATOM.EQ.0) THEN
                      NATOMS = NATOMS + 1
                      IF (NATOMS.GT.MXATOM) GO TO 906
                      IATOM = NATOMS
                      IF (LPOS.EQ.-1) LPOS = LSTPOS
                      IF (LPOS.EQ.0) THEN
                          NXTATM(IATOM) = FSTATM
                          FSTATM = IATOM
                      ELSE
                          NXTATM(IATOM) = NXTATM(LPOS)
                          NXTATM(LPOS) = IATOM
                      ENDIF
                      LPOS = IATOM
                      IF (DUPLIC) DOUBL(IATOM) = .TRUE.
                  ENDIF
              ENDIF

C----         Add the coordinates of this atom to the mean coords
              IF (IATOM.GT.0) THEN
                  MCOORD(1,IATOM) = MCOORD(1,IATOM) + COORDS(1)
                  MCOORD(2,IATOM) = MCOORD(2,IATOM) + COORDS(2)
                  MCOORD(3,IATOM) = MCOORD(3,IATOM) + COORDS(3)
                  NCOORD(IATOM) = NCOORD(IATOM) + 1

CHECK v.3.2.1-->
C----             Store the atom's coords if this is the first file
                  IF (NMODEL.EQ.1) THEN
                      FCOORD(1,IATOM) = COORDS(1)
                      FCOORD(2,IATOM) = COORDS(2)
                      FCOORD(3,IATOM) = COORDS(3)
                      GOTFST(IATOM) = .TRUE.
                  ENDIF
CHECK v.3.2.1<--

C----             Store marker identifying whether this atom is a CA (C),
C                 a hydrogen (h), other heavy atom (a), or HETATM hydrogen
C                 (H), or HETATM heavy atom (A)
                  IF (IDENT.EQ.'ATOM  ') THEN
                      IF (HATOM) THEN
                          HAMARK(IATOM) = 'h'
                      ELSE
                          HAMARK(IATOM) = 'a'
                      ENDIF
                  ELSE
                      IF (HATOM) THEN
                          HAMARK(IATOM) = 'H'
                      ELSE
                          HAMARK(IATOM) = 'A'
                      ENDIF
                  ENDIF

C----             Save the ATOM ID
                  ATOMID(IATOM) = ATMKEY
              ENDIF

C----     If MODEL record, then increment count of models processed and
C         reinitialise variables
          ELSE IF (IDENT.EQ.'MODEL ') THEN
              IF (ENSEMB .AND. FIRSTM) THEN
                  NMR = .TRUE.
                  NMRMOD = NMRMOD + 1
                  FIRSTM = .FALSE.
              ELSE
                  NMODEL = NMODEL + 1
                  NMRMOD = NMRMOD + 1
              ENDIF
              READIN = .TRUE.
              IPOS = 0
              IATOM = 0
              LPOS = 0

C----     If ENDMDL record, then adjust sequence length
          ELSE IF (IDENT.EQ.'ENDMDL') THEN

C----         Adjust the maximum and minimum sequence length
              IF (NATOMS.EQ.0) THEN
                  NATOMS = IATOM
              ENDIF
              IF (IPOS.GT.0) THEN
                  IF (MAXLEN.EQ.0) THEN
                      MINLEN = IPOS
                      MAXLEN = IPOS
                  ELSE
                      MINLEN = MIN(IPOS,MINLEN)
                      MAXLEN = MAX(IPOS,MAXLEN)
                  ENDIF
              ENDIF
              FIRST = .FALSE.
              READIN = .FALSE.
          ENDIF
      GO TO 100

C---- End of Brookhaven file reached
 800  CONTINUE

C---- Show statistics on the models read in
      IF (FIRST .AND. .NOT.NMR) NATOMS = IATOM
      FIRST = .FALSE.
      IF (NMODEL.EQ.0) GO TO 900
      IF (NATOMS.EQ.0) GO TO 908
      IF (NMR) THEN
          PRINT*
          PRINT*, 'NMR structure'
          PRINT*, '-------------'
          PRINT 810, 'Number of models read in      ', NMRMOD
 810      FORMAT(1X,A,I6)
          PRINT*
          PRINT 810, 'Minimum atoms in sequence     ', MINLEN
          PRINT 810, 'Maximum atoms in sequence     ', MAXLEN
          IF (.NOT.ENSEMB) THEN
              PRINT*
              PRINT 810, 'Total no. of unique atoms     ', NATOMS
          ENDIF
          PRINT*
      ENDIF

      GO TO 999

C---- Fatal errors
900   CONTINUE
CHECK v.3.4-->
C      PRINT*, '**** ERROR. No MODEL records found in .new file!'
C      GO TO 990
      PRINT*, '* Warning. No MODEL records found in PDB file'
      PRINT*, '*          Single structure assumed; no RMS deviations ',
     -    'to be calculated'
      GO TO 999
CHECK v.3.4<--

902   CONTINUE
      PRINT*, '**** Error in coords for atom: ', ATNAME, '-', RESDUE,
     -    '-', INSEQ
      GO TO 100

904   CONTINUE
      PRINT*, '**** Data error reading file:  [', FILNEW, ']'
      PRINT*, '    Line: ', LINE
      GO TO 990

906   CONTINUE
      PRINT*, '**** Maximum number of atoms per protein exceeded:',
     -    MXATOM
      GO TO 990

908   CONTINUE
      PRINT*, '**** ERROR. No coordinates found in .new file!'
      PRINT*, '**** May be due to missing MODEL or ENDMDL records'
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETRES  -  Assign a sequential number to each residue for
C                        calculating residue RMS deviations later
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETRES

      INCLUDE 'rmsdev.inc'

      CHARACTER*10  LSTRES, RESID
      INTEGER       IATOM, IPOS, IRESID

C---- Initialise variables
      IATOM = FSTATM
      IRESID = 0
      LSTRES = ' '

C---- Loop through all the atoms in the structure
      DO 500, IPOS = 1, NATOMS

C----     Get the residue identifier for this atom
          RESID = ATOMID(IATOM)(1:10)

C----     If residue has changed, then increment residue number
          IF (RESID.NE.LSTRES) THEN
              IRESID = IRESID + 1
              IF (IRESID.GT.MXRES) GO TO 900
              LSTRES = RESID
          ENDIF

C----     Save the sequential residue number for this atom
          ATRESN(IATOM) = IRESID
          RESDET(IRESID) = RESID(2:)

C----     Get the next atom in the linked list
          IATOM = NXTATM(IATOM)
500   CONTINUE

C---- Save the number of residues encountered
      NRESID = IRESID
      IF (NRESID.EQ.0) GO TO 902

      GO TO 999

C---- Error messages
900   CONTINUE
      PRINT*, '**** Maximum number of residues per protein exceeded:',
     -    MXRES
      GO TO 990

902   CONTINUE
      IF (ENSEMB) THEN
          PRINT*, '**** ERROR. No residues found in PDB files!'
      ELSE
          PRINT*, '**** ERROR. No residues found in .new file!'
      ENDIF
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE WRIAVE  -  Write out the averaged coordinates of the
C                        structure
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRIAVE

      INCLUDE 'rmsdev.inc'

      CHARACTER*3   RESDUE
      CHARACTER*4   ATNAME
      CHARACTER*6   IDENT, INSEQ
      INTEGER       I, IATOM, IPOS

C---- Initialise variables
      IATOM = FSTATM

C---- Loop through all the atoms in the structure
      DO 800, IPOS = 1, NATOMS

C----     Calculate the mean coordinates for this atom
          IF (NCOORD(IATOM).GT.0) THEN
              MCOORD(1,IATOM) = MCOORD(1,IATOM) / NCOORD(IATOM)
              MCOORD(2,IATOM) = MCOORD(2,IATOM) / NCOORD(IATOM)
              MCOORD(3,IATOM) = MCOORD(3,IATOM) / NCOORD(IATOM)
          ENDIF

C----     Write out the mean coordinates
          IF (HAMARK(IATOM).EQ.'h' .OR. HAMARK(IATOM).EQ.'a') THEN
              IDENT = 'ATOM'
          ELSE
              IDENT = 'HETATM'
          ENDIF
          ATNAME = ATOMID(IATOM)(11:14)
          RESDUE = ATOMID(IATOM)(8:10)
          INSEQ = ATOMID(IATOM)(2:7)
          WRITE(8,100) IDENT, IPOS, ATNAME, RESDUE, INSEQ,
     -        (MCOORD(I,IATOM), I = 1, 3), REAL(NCOORD(IATOM)),
     -        ATRESN(IATOM)
100       FORMAT(A6,I5,1X,A4,1X,A3,1X,A6,3X,3F8.3,'  1.00',F6.2,I6)

CHECK v.3.2.1-->
C----     If there are only two models, then store first model's
C         mean values (if there are any)
          IF (NMODEL.EQ.2 .AND. GOTFST(IATOM)) THEN
              MCOORD(1,IATOM) = FCOORD(1,IATOM)
              MCOORD(2,IATOM) = FCOORD(2,IATOM)
              MCOORD(3,IATOM) = FCOORD(3,IATOM)
          ENDIF
CHECK v.3.2.1<--

C----     Get the pointer to the next atom
          IATOM = NXTATM(IATOM)
800   CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETCON  -  Read in the distance constraints from the .mr file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCON

      INCLUDE 'rmsdev.inc'

      CHARACTER*1   CHAIN
      CHARACTER*3   RESDUE(2,2)
      CHARACTER*4   ATNAME(2,2)
      CHARACTER*5   DIHEDR(2), RNUMB(2,2)
      CHARACTER*6   CONSID(NCTYPE), IDENT
      CHARACTER*14  ATMKEY
      CHARACTER*80  IREC
      INTEGER       ATPOS(2), BOUND(NCTYPE), CCOUNT(NCTYPE), CTYPE,
     -              IATOM, ICONST, IEND, INCON, IPOINT, IRESID, ITORS,
     -              ITYPE, LINE, NEND, NODD, NOTFND, ONLINE
      LOGICAL       HAVEIT
      REAL          LOWER(2), UPPER(2), VALUE(2)

      DATA CONSID / 'NOELOW', 'NOEUPP', 'HBLOW', 'HBUPP', 'SSLOW',
     -              'SSUPP' /

      DATA BOUND  /       -1,        1,      -1,       1,      -1,
     -                     1 /

C---- Initialise variables
      CHAIN = ' '
      IATOM = 0
      LINE = 0
      NINTRA = 0
      NODD = 0
      NOTFND = 0
      DO 50, CTYPE = 1, NCTYPE
          CCOUNT(CTYPE) = 0
 50   CONTINUE

C---- Initialise the pointers giving the restraints involving each
C     residue
      DO 60, IRESID = 1, MXRES
          LSTCON(1,IRESID) = 0
          LSTCON(2,IRESID) = 0
 60   CONTINUE
      DO 80, ICONST = 1, MXCONS
          NXTCON(1,ICONST) = 0
          NXTCON(2,ICONST) = 0
 80   CONTINUE
      ICONST = 0

C---- Loop through the constraint file, storing all the distance constraints
100   CONTINUE
          READ(2,150,END=1200,ERR=900) IREC
150       FORMAT(A)
          LINE = LINE + 1
          IDENT = IREC(1:6)
          CTYPE = 0

C----     Determine whether this is one of the constraint records
          DO 200, ITYPE = 1, NCTYPE
              IF (IDENT.EQ.CONSID(ITYPE)) CTYPE = ITYPE
 200      CONTINUE

C----     If record is a constraint record, determine which atoms are
C         involved prior to storing the constraint
          IF (CTYPE.NE.0) THEN

C----         Pick up the details of the one (or two) distance constraints
C             on this record
              READ(IREC,250,ERR=902) RESDUE(1,1), RNUMB(1,1),
     -            ATNAME(1,1), RESDUE(2,1), RNUMB(2,1), ATNAME(2,1),
     -            VALUE(1), 
     -            ATNAME(1,2), RESDUE(2,2), RNUMB(2,2), ATNAME(2,2),
     -            VALUE(2)
 250          FORMAT(8X,A3,1X,A5,2(A4,1X,A3,1X,A5,A4,F6.2,5X))
              RESDUE(1,2) = RESDUE(1,1)
              RNUMB(1,2) = RNUMB(1,1)

C----         Determine whether there is only one constraint or two on
C             this line
              ONLINE = 2
              IF (ATNAME(1,2).EQ.' ') ONLINE = 1

C----         Loop for each constraint picked up from the line
              DO 600, INCON = 1, ONLINE

C----             Initialise search for this constraint
                  ATPOS(1) = 0
                  ATPOS(2) = 0

C----             Loop for the two atoms involved in the constraint
                  DO 400, IEND = 1, 2

C----                 If atom name is ' HN ', replace it by ' H  '
                      IF (ATNAME(IEND,INCON).EQ.' HN ') THEN
                          ATNAME(IEND,INCON) = ' H  '
                      ENDIF

C----                 Initialise the search key for this atom
                      ATMKEY = ' ' // CHAIN // RNUMB(IEND,INCON) //
     -                    RESDUE(IEND,INCON) // ATNAME(IEND,INCON)
                      HAVEIT = .FALSE.
                      IATOM = 0

C----                 Loop until match found, or all atoms searched
300                   CONTINUE
                          IATOM = IATOM + 1

C----                     Check whether current atom matches here
                          IF (ATMKEY.EQ.ATOMID(IATOM)) HAVEIT = .TRUE.

C----                 Loop back, if necessary
                      IF (.NOT.HAVEIT .AND. IATOM.LT.NATOMS) GO TO 300

C----                 If the atom has been found, store its location
                      IF (HAVEIT) THEN
                          ATPOS(IEND) = IATOM
                      ENDIF
 400              CONTINUE

C----             If have located both atoms defining the ends of this
C                 constraint, then store the constraint
                  IF (ATPOS(1).GT.0 .AND. ATPOS(2).GT.0) THEN

C----                 Increment constraints count
                      ICONST = ICONST + 1
                      IF (ICONST.GT.MXCONS) GO TO 904
                      CCOUNT(CTYPE) = CCOUNT(CTYPE) + 1

C----                 Store the constraint
                      CONATM(1,ICONST) = ATPOS(1)
                      CONATM(2,ICONST) = ATPOS(2)
                      CONBOU(ICONST) = BOUND(CTYPE)
                      CONTYP(ICONST) = CTYPE
                      CNSTRN(ICONST) = VALUE(INCON)

C----                 Check whether the constraint is an intra-residue
C                     one, in which case only need to count it once for
C                     this residue
                      IF (ATRESN(ATPOS(1)).EQ.ATRESN(ATPOS(2))) THEN
                          NEND = 1
                          NINTRA = NINTRA + 1
                      ELSE
                          NEND = 2
                      ENDIF

C----                 Update the linked-list of residue-pointers giving
C                     list of constraints involving any given residue
                      DO 550, IEND = 1, NEND

C----                     Get the residue for this atom
                          IRESID = ATRESN(ATPOS(IEND))

C----                     Find the last constraint involving this residue
                          IPOINT = LSTCON(IEND,IRESID)

C----                     Update the last-constraint pointers
                          NXTCON(IEND,ICONST) = IPOINT
                          LSTCON(IEND,IRESID) = ICONST
 550                  CONTINUE

C----             If one or other atom not found, then increment count of
C                 missed constraints
                  ELSE
                      NOTFND = NOTFND + 1
                  ENDIF
 600          CONTINUE

C----     Otherwise, if this is a dihedral angle constraint, read in and
C         store the constraint details
          ELSE IF (IDENT.EQ.'ANGLE') THEN

C----         Pick up the details of the one (or two) residue dihedrals
C             on this record
              READ(IREC,620,ERR=902) RESDUE(1,1), RNUMB(1,1),
     -            DIHEDR(1), LOWER(1), UPPER(1),
     -            DIHEDR(2), LOWER(2), UPPER(2)
 620          FORMAT(8X,A3,1X,A5,2(A5,2F8.2,8X))

C----         Determine whether there is only one constraint or two on
C             this line
              ONLINE = 2
              IF (DIHEDR(2).EQ.'    ') ONLINE = 1

C----         Loop for each dihedral constraint picked up from the line
              DO 800, INCON = 1, ONLINE

C----             Search for this residue

C----             Initialise the search key for this atom
                  ATMKEY = ' ' // CHAIN // RNUMB(1,1) // RESDUE(1,1)
                  HAVEIT = .FALSE.
                  IATOM = 0

C----             Loop until match found, or all atoms searched
700               CONTINUE
                      IATOM = IATOM + 1

C----                 Check whether current atom matches here
                      IF (ATMKEY(1:10).EQ.ATOMID(IATOM)(1:10)) THEN
                          HAVEIT = .TRUE.
                      ENDIF

C----             Loop back, if necessary
                  IF (.NOT.HAVEIT .AND. IATOM.LT.NATOMS) GO TO 700

C----             If the atom has been found, get its corresponding
C                 sequential residue number
                  IF (HAVEIT) THEN
                      IRESID = ATRESN(IATOM)

C----                 Store the constraint according to its type
                      ITORS = 0
                      IF (DIHEDR(INCON).EQ.'PHI ') ITORS = 1
                      IF (DIHEDR(INCON).EQ.'PSI ') ITORS = 2
                      IF (DIHEDR(INCON).EQ.'CHI1') ITORS = 3
                      IF (ITORS.NE.0) THEN
                          DIHCON(1,ITORS,IRESID) = LOWER(INCON)
                          DIHCON(2,ITORS,IRESID) = UPPER(INCON)
                          NDIHED = NDIHED + 1
                      ELSE
                          NODD = NODD + 1
                      ENDIF
                  ELSE
                      NODD = NODD + 1
                  ENDIF

 800          CONTINUE
          ENDIF

      GO TO 100

C---- End of constraints file reached
 1200 CONTINUE

C---- Store the number of distance constraints read in
      NCONST = ICONST
      IF (NCONST.EQ.0) GO TO 906
      PRINT*, '* Distance constraints read in:-'
      PRINT*, '* ------------------------------'
      DO 850, CTYPE = 1, NCTYPE
          PRINT 840, '*    ', CONSID(CTYPE), CCOUNT(CTYPE)
 840      FORMAT(A,A,5X,I7,A)
 850  CONTINUE
      PRINT 840, '*    ', '            ------'
      PRINT 840, '*    ', 'Total ', NCONST
      IF (NOTFND.GT.0) THEN
          PRINT*
          PRINT 840, '**** ', 'Missed', NOTFND, '   (One or both ' //
     -        'atoms defining the constraint not found)'
      ENDIF
      PRINT*
      PRINT*, '* Dihedral angle constraints:', NDIHED
      PRINT*
      IF (NODD.GT.0) THEN
          PRINT*, '*** Unknown angle restraints found:', NODD
      ENDIF

      GO TO 999

C---- Fatal errors
 900  CONTINUE
      PRINT*, '**** Error in constraint record at line', LINE
      GO TO 990

 902  CONTINUE
      PRINT*, '**** Data error reading file:  [', FILCON, ']'
      PRINT*, '*    Line: ', LINE
      GO TO 990

 904  CONTINUE
      PRINT*, '**** Maximum number of distance restraints exceeded:',
     -    MXCONS
      GO TO 990

 906  CONTINUE
      PRINT*, '**** WARNING. No constraints found in constraints file!'
      HAVCON = .FALSE.
      GO TO 999

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE CALRMS  -  Read through the .new file a second time to
C                        calculate the atomic RMS deviations from the mean
C                        coordinates
C
C----------------------------------------------------------------------+---

      SUBROUTINE CALRMS

      SAVE

      INCLUDE 'rmsdev.inc'

      CHARACTER*1   CHAIN, HETMRK
      CHARACTER*3   RESDUE
      CHARACTER*4   ATNAME, MATOM(NMATOM)
      CHARACTER*5   SEQNO
      CHARACTER*6   IDENT, INSEQ
      CHARACTER*14  ATMKEY
      CHARACTER*78  IREC
      INTEGER       I, IATOM, IMAT, IMODEL, IPOS, IRESID, LINE,
     -              NATS, NMATS, NSATS
      LOGICAL       ENDNEW, FIRST, FIRSTM, HATOM, HAVEIT, MCHAIN, NMR,
     -              READIN, WANTED
      REAL          ATMDEV, COORDS(3), DEV1, DEV2, DEV3, MDEV, SDEV,
     -              SUMDEV

      DATA FIRST  / .TRUE. /
      DATA  MATOM / ' C  ', ' N  ', ' O  ', ' CA ' /

C---- Rewind the .new file
      REWIND(1)
      ENDNEW = .FALSE.

C---- Blank out the atom coordinates for the new model
      DO 50, IATOM = 1, MXATOM
          HAVATM(IATOM) = .FALSE.
          XYZ(1,IATOM) = 0.0
          XYZ(2,IATOM) = 0.0
          XYZ(3,IATOM) = 0.0
 50   CONTINUE

C---- Initialise variables
      FIRSTM = .TRUE.
      IFAIL = .FALSE.
      IATOM = FSTATM
      MVIOL = 0
      NATS = 0
      NSATS = 0
      NMATS = 0
      MDEV = 0.0
      SDEV = 0.0
      SUMDEV = 0.0
      IF (FIRST) THEN
          IMODEL = 0
          LINE = 0
      ENDIF
      FIRST = .FALSE.
      NMR = .FALSE.
      IF (ENSEMB) THEN
          READIN = .TRUE.
          IMODEL = IMODEL + 1
          WRITE(IREC,80) IMODEL
 80       FORMAT('MODEL',I9)
          WRITE(7,150) IREC
          IF (HAVCON) WRITE(9,150) IREC
      ELSE
          NMR = .TRUE.
          READIN = .FALSE.
      ENDIF

C---- Loop through the .new file, reading in all the ATOM and HETATM
C     records
100   CONTINUE
          READ(1,150,END=160,ERR=904) IREC
150       FORMAT(A)
          LINE = LINE + 1
          IDENT = IREC(1:6)
          GO TO 170

C----     End of file reached (equivalent to reaching end of current model)
 160      CONTINUE
          IF (READIN) THEN
              IDENT = 'ENDMDL'
          ENDIF
          ENDNEW = .TRUE.

C----     Determine whether this record is required for processing
 170      CONTINUE
          WANTED = .FALSE.
          IF (READIN) THEN
              IF (IDENT.EQ.'ATOM  ' .OR. IDENT.EQ.'HETATM') THEN
                  WANTED = .TRUE.
                  IF (IREC(18:20).EQ.'HOH' .OR.
     -                IREC(18:20).EQ.'WAT') WANTED = .FALSE.
              ENDIF
          ENDIF

C----     Process record if required
          IF (WANTED) THEN

C----         Pick off atom name, residue name, chain and sequence number
              ATNAME = IREC(13:16)
              RESDUE = IREC(18:20)
              CHAIN = IREC(22:22)
              INSEQ = IREC(22:27)
              SEQNO = IREC(23:27)

C----         Check if this is a hydrogen atom or a pseudo-atom
              IF (ATNAME(2:2).EQ.'H' .OR. ATNAME(2:2).EQ.'Q' .OR.
     -            ATNAME(1:1).EQ.'H' .OR. ATNAME(1:1).EQ.'Q') THEN
                  HATOM = .TRUE.
              ELSE
                  HATOM = .FALSE.
              ENDIF

C----         Retrieve the atomic coordinates
              READ(IREC,180,ERR=902) (COORDS(I), I = 1, 3)
180           FORMAT(30X,3F8.0)

C----         Find the atom's location in the sequence and compute the
C             deviation of its coordinates from the mean

C----         Set HETATM marker if this is a HETATM record
              IF (IDENT.EQ.'ATOM  ') THEN
                  HETMRK = ' '
              ELSE
                  HETMRK = 'H'
              ENDIF

C----         Form the unique, identifying key for this atom
              ATMKEY = HETMRK // INSEQ // RESDUE // ATNAME

C----         Check whether key matches the current atom we're looking at
              HAVEIT = .FALSE.
              IF (ATMKEY.EQ.ATOMID(IATOM)) THEN
                  HAVEIT = .TRUE.
              ENDIF

C----         If a match not found, start the search at the beginning
              IF (.NOT.HAVEIT) THEN
                  IPOS = 0

C----             Loop until match found, or all atoms searched
200               CONTINUE
                      IPOS = IPOS + 1

C----                 Get the pointer to the next atom
                      IF (IPOS.EQ.1) THEN
                          IATOM = FSTATM
                      ELSE
                          IATOM = NXTATM(IATOM)
                      ENDIF

C----                 Check whether current atom matches here
                      IF (ATMKEY.EQ.ATOMID(IATOM)) HAVEIT = .TRUE.

C----             Loop back, if necessary
                  IF (.NOT.HAVEIT .AND. IPOS.LT.NATOMS) GO TO 200
              ENDIF

CHECK v.3.2.1-->
C----         If there are only 2 models, calculate RMS's only for
C             the second model
              IF (NMODEL.EQ.2 .AND. IMODEL.EQ.1) HAVEIT = .FALSE.
              IF (NMODEL.EQ.2 .AND. IMODEL.EQ.2 .AND.
     -            .NOT.GOTFST(IATOM)) HAVEIT = .FALSE.
CHECK v.3.2.1<--

C----         Calculate RMS deviation of coords for this atom from the
C             mean
              IF (HAVEIT) THEN
                  DEV1 = MCOORD(1,IATOM) - COORDS(1)
                  DEV2 = MCOORD(2,IATOM) - COORDS(2)
                  DEV3 = MCOORD(3,IATOM) - COORDS(3)
                  ATMDEV = DEV1 * DEV1 + DEV2 * DEV2 + DEV3 * DEV3

C----             Determine whether this is a main-chain or side-chain
C                 atom
                  MCHAIN = .FALSE.
                  DO 250, IMAT = 1, NMATOM
                      IF (ATMKEY(11:14).EQ.MATOM(IMAT)) MCHAIN = .TRUE.
250               CONTINUE

C----             Get the corresponding residue number
                  IRESID = ATRESN(IATOM)

C----             Store the atom's coordinates
                  HAVATM(IATOM) = .TRUE.
                  XYZ(1,IATOM) = COORDS(1)
                  XYZ(2,IATOM) = COORDS(2)
                  XYZ(3,IATOM) = COORDS(3)

C----             If this is a heavy atom increment accumulators for this
C                 model
                  IF (.NOT.HATOM) THEN

C----                 Main-chain atoms
                      IF (MCHAIN) THEN
                          RMS(1,IRESID) = RMS(1,IRESID) + ATMDEV
                          NRMS(1,IRESID) = NRMS(1,IRESID) + 1
                          MDEV = MDEV + ATMDEV
                          NMATS = NMATS + 1

C----                 Side-chain atoms
                      ELSE
                          RMS(2,IRESID) = RMS(2,IRESID) + ATMDEV
                          NRMS(2,IRESID) = NRMS(2,IRESID) + 1
                          SDEV = SDEV + ATMDEV
                          NSATS = NSATS + 1
                      ENDIF

C----                 All heavy atoms
                      RMS(3,IRESID) = RMS(3,IRESID) + ATMDEV
                      NRMS(3,IRESID) = NRMS(3,IRESID) + 1
                      SUMDEV = SUMDEV + ATMDEV
                      NATS = NATS + 1
                  ENDIF

C----             Prepare pointer to the next atom
                  IATOM = NXTATM(IATOM)
CHECK v.3.4-->
                  IF (IATOM.EQ.0) IATOM = FSTATM
CHECK v.3.4<--
              ELSE
                  IATOM = FSTATM
              ENDIF

C----     If ENDMDL record, then calculate RMS deviation for the
C         model and write out
          ELSE IF (IDENT.EQ.'ENDMDL') THEN

C----         Calculate the constraint violations for this model
              IF (HAVCON) THEN
                  CALL VIOLA
              ENDIF
              NUMCON = 0

C----         Write out the residue-by-residue RMS deviations and
C             constraints for this model
              DO 400, IRESID = 1, NRESID
                  DO 300, I = 1, 3
                      NRMALL(I,IRESID)
     -                    = NRMALL(I,IRESID) + NRMS(I,IRESID)
                      RMSALL(I,IRESID)
     -                    = RMSALL(I,IRESID) + RMS(I,IRESID)
                      IF (NRMS(I,IRESID).GT.0) THEN
                          RMS(I,IRESID)
     -                         = SQRT(RMS(I,IRESID) / NRMS(I,IRESID))
                      ENDIF
300               CONTINUE

C----             Write the details out to the .rms file
CHECK v.3.2.1-->
                  IF (NMODEL.NE.2) THEN
CHECK v.3.2.1<--
                      WRITE(7,320) IRESID, RESDET(IRESID),
     -                    (RMS(I,IRESID), I = 1, 3),
     -                    (NRMS(I,IRESID), I = 1, 3)
320                   FORMAT(I6,A9,3F10.3,3I6)
CHECK v.3.2.1-->
                  ENDIF
CHECK v.3.2.1<--

C----             Write out all the constraints involving this residue
                  IF (HAVCON) CALL PUTVIO(IRESID)
400           CONTINUE

C----         Write out the number of constraint violations for this
C             model
              IF (HAVCON) THEN
                  PRINT*, 'Distance constraint violations for model',
     -                IMODEL, ':   ', MVIOL
                  TVIOL = TVIOL + MVIOL
              ENDIF

C----         Write out the ENDMDL record to the .rms file
CHECK v.3.2.1-->
                  IF (NMODEL.NE.2) THEN
CHECK v.3.2.1<--
                  WRITE(7,150) 'ENDMDL'
CHECK v.3.2.1-->
                  ENDIF
CHECK v.3.2.1<--
              IF (HAVCON) WRITE(9,150) 'ENDMDL'
              READIN = .FALSE.

C----         Add RMS deviation to overall RMS deviation between models
              MALLDV = MALLDV + MDEV
              NMALL = NMALL + NMATS
              SALLDV = SALLDV + SDEV
              NSALL = NSALL + NSATS
              SUMALL = SUMALL + SUMDEV
              NATALL = NATALL + NATS

C----         Calculate RMS deviation for model as a whole and write out
              IF (NMATS.GT.0) THEN
                  MDEV = SQRT(MDEV / NMATS)
              ENDIF
              IF (NSATS.GT.0) THEN
                  SDEV = SQRT(SDEV / NSATS)
              ENDIF
              IF (NATS.GT.0) THEN
                  SUMDEV = SQRT(SUMDEV / NATS)
              ENDIF
CHECK v.3.2.1-->
                  IF (NMODEL.NE.2) THEN
CHECK v.3.2.1<--
                  WRITE(7,450) MDEV, SDEV, SUMDEV, NMATS, NSATS, NATS
450               FORMAT('Model RMS deviations:',3F10.4,3I8)
CHECK v.3.2.1-->
                  ENDIF
CHECK v.3.2.1<--

C----         Reinitialise for the next model
              DO 600, IRESID = 1, NRESID
                  DO 500, I = 1, 3
                      NRMS(I,IRESID) = 0
                      RMS(I,IRESID) = 0.0
500               CONTINUE
600           CONTINUE

C----     If MODEL record, then reinitialise variables
          ELSE IF (IDENT.EQ.'MODEL ') THEN

C----         Blank out the atom coordinates for the new model
              DO 650, IATOM = 1, MXATOM
                  HAVATM(IATOM) = .FALSE.
                  XYZ(1,IATOM) = 0.0
                  XYZ(2,IATOM) = 0.0
                  XYZ(3,IATOM) = 0.0
 650          CONTINUE

C----         Initialise the other variables
              IATOM = FSTATM
              NATS = 0
              NMATS = 0
              NSATS = 0
              MDEV = 0.0
              SDEV = 0.0
              SUMDEV = 0.0
              MVIOL = 0
              IF (ENSEMB .AND. FIRSTM) THEN
                  NMR = .TRUE.
                  FIRSTM = .FALSE.
              ELSE

C----             Write out the MODEL record to the .rms and .cns files
                  IMODEL = IMODEL + 1
                  WRITE(IREC,80) IMODEL
CHECK v.3.2.1-->
                  IF (NMODEL.NE.2) THEN
CHECK v.3.2.1<--
                  WRITE(7,150) IREC
CHECK v.3.2.1-->
                  ENDIF
CHECK v.3.2.1<--
                  IF (HAVCON) WRITE(9,150) IREC
              ENDIF
              READIN = .TRUE.
          ENDIF
      IF (.NOT.ENDNEW) GO TO 100

      GO TO 999

C---- Fatal errors
902   CONTINUE
      PRINT*, '**** Error in coords for atom: ', ATNAME, '-', RESDUE,
     -    '-', INSEQ
      GO TO 100

904   CONTINUE
      PRINT*, '**** Data error reading file:  [', FILNEW, ']'
      PRINT*, '    Line: ', LINE
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE VIOLA  -  Calculate the constraint violations for the
C                       current model  
C
C----------------------------------------------------------------------+---

      SUBROUTINE VIOLA

      INCLUDE 'rmsdev.inc'

      INTEGER       IATOM1, IATOM2, ICONST
      REAL          DIST, DX, DY, DZ

C---- Loop through all the stored constraints
      DO 800, ICONST = 1, NCONST

C----     Get the two atoms defining this distance restraint
          IATOM1 = CONATM(1,ICONST)
          IATOM2 = CONATM(2,ICONST)
          ACDIST(ICONST) = -999.9

C----     Proceed if both atoms are present in the current model
          IF (HAVATM(IATOM1) .AND. HAVATM(IATOM2)) THEN

C----         Calculate the distance between the two atoms
              DX = XYZ(1,IATOM1) - XYZ(1,IATOM2)
              DY = XYZ(2,IATOM1) - XYZ(2,IATOM2)
              DZ = XYZ(3,IATOM1) - XYZ(3,IATOM2)
              DIST = SQRT(DX * DX + DY * DY + DZ * DZ)

C----         Store the actual distance with this constraint
              ACDIST(ICONST) = DIST
          ENDIF
 800  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PUTVIO  -  Write out the constraints involving the current
C                        residue
C
C----------------------------------------------------------------------+---

      SUBROUTINE PUTVIO(IRESID)

      INCLUDE 'rmsdev.inc'

      CHARACTER*3   OUTSID
      INTEGER       IATOM, ICONST, IEND, IPOINT, IRESID, JATOM, JEND,
     -              JRESID

C---- Loop through the two sets of constraint linked-lists
      DO 500, IEND = 1, 2

C----     Get the index of the other end
          JEND = 3 - IEND

C----     Get the pointer to the last-read constraint for this residue
          IPOINT = LSTCON(IEND,IRESID)

C----     Loop until get to the end of the chain
 100      CONTINUE

C----         If this is a valid constraint position, extract the constraint
              IF (IPOINT.NE.0) THEN
                  ICONST = IPOINT

C----             Increment sequential constraint number
                  NUMCON = NUMCON + 1

C----             Proceed only if actual distance was calculated
                  IF (ACDIST(ICONST).GE.0.0) THEN

C----                 Get the two atoms involved in this constraint
                      IATOM = CONATM(IEND,ICONST)
                      JATOM = CONATM(JEND,ICONST)

C----                 Get the other residue involved in this constraint
                      JRESID = ATRESN(JATOM)

C----                 Determine whether the actual distance violates the
C                     given constraint
                      OUTSID = ' '
                      IF (CONBOU(ICONST).EQ.LBOUND) THEN
                          IF (ACDIST(ICONST).LT.CNSTRN(ICONST)) THEN
                              OUTSID = '*<*'
                          ENDIF
                      ENDIF
                      IF (CONBOU(ICONST).EQ.UBOUND) THEN
                          IF (ACDIST(ICONST).GT.CNSTRN(ICONST)) THEN
                              OUTSID = '*>*'
                          ENDIF
                      ENDIF

C----                 If this is a constraint violation (and the first
C                     time it's been encountered) increment the count
C                     of violations
                      IF (OUTSID(1:1).EQ.'*' .AND.
     -                    IRESID.LE.JRESID) THEN
                          MVIOL = MVIOL + 1
                          NVIOL(ICONST) = NVIOL(ICONST) + 1
                      ENDIF

C----                 Write the constraint details out to file
                      WRITE(9,120) IRESID, RESDET(IRESID), NUMCON,
     -                    CNSTRN(ICONST), ACDIST(ICONST), OUTSID
 120                  FORMAT(I6,A9,I6,2F7.2,2X,A3)
                  ENDIF

C----             Get pointer to next constraint involving this residue
C                 and loop back
                  IPOINT = NXTCON(IEND,IPOINT)
                  GO TO 100
              ENDIF
 500  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE ALLMOD  -  Write out the overall RMS deviations for all
C                        models read in
C
C----------------------------------------------------------------------+---

      SUBROUTINE ALLMOD

      SAVE

      INCLUDE 'rmsdev.inc'

      INTEGER       I, IRESID

C---- Write out overall residue-by-residue RMS deviations
      WRITE(7,150) 'ALL MODELS'
 150  FORMAT(A)

C---- Loop through all the residues
      DO 200, IRESID = 1, NRESID
          DO 100, I = 1, 3
              IF (NRMALL(I,IRESID).GT.0) THEN
                  RMSALL(I,IRESID)
     -                 = SQRT(RMSALL(I,IRESID) / NRMALL(I,IRESID))
              ENDIF
 100      CONTINUE

C----     Write the details out to the .rms file
          WRITE(7,320) IRESID, RESDET(IRESID),
     -        (RMSALL(I,IRESID), I = 1, 3),
     -        (NRMALL(I,IRESID), I = 1, 3)
 320      FORMAT(I6,A9,3F10.3,3I6)
 200  CONTINUE

C---- Write out the ENDMDL record to the .rms file
      WRITE(7,150) 'ENDMDL'

C---- Calculate overall RMS deviation and write out
      IF (NMALL.GT.0) THEN
          MALLDV = SQRT(MALLDV / NMALL)
      ENDIF
      IF (NSALL.GT.0) THEN
          SALLDV = SQRT(SALLDV / NSALL)
      ENDIF
      IF (NATALL.GT.0) THEN
          SUMALL = SQRT(SUMALL / NATALL)
      ENDIF
      WRITE(7,250) MALLDV, SALLDV, SUMALL, NMALL, NSALL, NATALL
 250  FORMAT('Overall RMS deviation:',3F10.4,3I8)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PUTCON  -  Write out the constraints to the output .cns
C                        file in residue order
C
C----------------------------------------------------------------------+---

      SUBROUTINE PUTCON

      INCLUDE 'rmsdev.inc'

      CHARACTER*1   AST
      CHARACTER*4   CNUMB, CTYPE(MXTORS)
      INTEGER       IATOM, ICONST, IEND, IPOINT, IRESID, ITORS, JATOM,
     -              JEND, JRESID

      DATA CTYPE  / 'PHI ',  'PSI ', 'CHI1' /

C---- Write out the header to the .cns constraints output file
      WRITE(9,*) 'Distance constraints'
      NUMCON = 0

C---- Loop through all the residues
      DO 600, IRESID = 1, NRESID

C----     Loop through the two sets of constraint linked-lists
          DO 500, IEND = 1, 2

C----         Get the index of the other end
              JEND = 3 - IEND

C----         Get the pointer to the last-read constraint for this residue
              IPOINT = LSTCON(IEND,IRESID)

C----         Loop until get to the end of the chain
 100          CONTINUE

C----             If this is a valid constraint position, extract the
C                 constraint
                  IF (IPOINT.NE.0) THEN
                      NUMCON = NUMCON + 1
                      ICONST = IPOINT

C----                 Get the two atoms involved in this constraint
                      IATOM = CONATM(IEND,ICONST)
                      JATOM = CONATM(JEND,ICONST)

C----                 Get the other residue involved in this constraint
                      JRESID = ATRESN(JATOM)

C----                 Write the constraint out to the output .cns file
                      WRITE(CNUMB,320) NVIOL(ICONST)
 320                  FORMAT(I4)
                      IF (NVIOL(ICONST).EQ.0) THEN
                          AST = ' '
                          CNUMB = ' '
                      ELSE
                          AST = '*'
                      ENDIF
                      WRITE(9,420) IRESID, RESDET(IRESID),
     -                    ATOMID(IATOM)(11:14),
     -                    JRESID, RESDET(JRESID),
     -                    ATOMID(JATOM)(11:14),
     -                    CNSTRN(ICONST), CONTYP(ICONST), NUMCON,
     -                    CNUMB, AST, ICONST
 420                  FORMAT(2(I6,A9,A4),F7.2,I2,I6,3X,A4,A1,I6)

C----                 Get pointer to next constraint involving this residue
C                     and loop back
                      IPOINT = NXTCON(IEND,IPOINT)
                      GO TO 100
                  ENDIF
 500      CONTINUE
 600  CONTINUE

C---- Write out the dihedral angle constraints to the output file
      WRITE(9,*) 'Dihedral angle constraints'

C---- Loop through all the residues to adjust the Chi-1 constraints
C     if necessaryu
      DO 650, IRESID = 1, NRESID
          IF (DIHCON(1,3,IRESID).LT.0.0) THEN
              DIHCON(1,3,IRESID) = DIHCON(1,3,IRESID) + 360.0
              DIHCON(2,3,IRESID) = DIHCON(2,3,IRESID) + 360.0
          ENDIF
 650  CONTINUE

C---- Loop through all the residues
      DO 800, IRESID = 1, NRESID
          DO 700, ITORS = 1, MXTORS
              IF (DIHCON(1,ITORS,IRESID).LT.900.0) THEN
                  WRITE(9,620) IRESID, RESDET(IRESID), CTYPE(ITORS),
     -                ITORS, DIHCON(1,ITORS,IRESID),
     -                DIHCON(2,ITORS,IRESID)
 620              FORMAT(I6,A9,1X,A4,I3,2F8.2)
              ENDIF
 700      CONTINUE
 800  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------

