C**************************************************************************
C
C  ANGLEN.FOR   -   Program to calculate means and standard deviations
C                   of main-chain bond lengths and bond angles
C
C                   The program reads the protein structure in from the
C                   user-defined .new file and outputs the results to a
C                   .lan file.
C
C  Written by R A Laskowski, University College, London, May 1992.
C
C  Original version was part of v.1.0 of the PROCHECK suite of programs.
C  Subsequent amendments have been labelled by CHECK v.m.n--> and
C  CHECK v.m.n<-- where m.n is the version number corresponding to the
C  change
C
C  v.2.0.1 - Bug-fix. Problem with certain compilers not liking uninitialised
C            and unSAVEd variables              Roman Laskowski 12 Aug 1992
C  v.2.1   - Initialisation of uninitialised variables
C                                               Roman Laskowski  8 Jan 1993
C  v.2.1.3 - Addition of asterisks to error prints so that can be plucked
C            out of the log files and displayed by the script file.
C                                               Roman Laskowski 26 Mar 1993
C  v.2.2   - Calculation of best-fit planes for planar groups, and RMS
C            deviations of the atoms from the plane.
C                                               Roman Laskowski 19 Oct 1993
C            Write out transformed coords of each planar group to .pln file
C                                               Roman Laskowski  5 Nov 1993
C
C v.2.3   - Tiny change so that VAX .COM routines don't show error if no
C           asterisks found in log file
C                                               Roman Laskowski 14 Nov 1993
C
C
C v.3.0.1 - Bug-fix. Problem picked up on SG compiler, due to 2D array
C           incorrectly defined as a 1D array.
C           Commented out unreferenced label (producing a compiler warning
C           when compiling on Convex).
C                                            Roman Laskowski 29/30 Mar 1994
C
C v.3.2   - Write out both residue names to .lan file for all bond lengths
C           and angles that span two residues so that both are shown in
C           the Distorted Geometry plots of program bplot.f.
C                                                Roman Laskowski 3 May 1994
C           Minor amendments to various statements to make them
C           acceptable to f2c, and to deal with various uninitialised
C           variables.(Amendments supplied by Dave Love at Daresbury).
C                                    David Love/Roman Laskowski 13 Oct 1994
C
C v.3.4.3 - Write out of "restraint violations" files for use with viol2pdb.
C                                               Roman Laskowski 11 Jul 1996
C
C v.3.4.4 - Increase in maximum number of residues (MXONE) to 10000.
C                                            Roman Laskowski (15 Oct 1996)
C
C v.3.5.5 - Minor fix for some compilers.
C                                            Roman Laskowski (25 Apr 2001)
C
C----------------------------------------------------------------------+---
C
C Files
C -----
C
C 1  <filename>.new - Cleaned-up version of the .pdb file holding the protein
C                     structure in Brookhaven format by-residue information 
C                     on the given structure
C 2  <filename>.lan - Output file holding the means and standard deviations
C                     of the main-chain bond lengths and bond angles
C 3  <filename>.pln - Output file holding the transformed coordinates of
C                     each planar group in the structure
C 7  <filename>.nrv - Output file holding the "restraint violations" for
C                     mainchain bond lengths
C
C--------------------------------------------------------------------------
C
C Subroutine calling tree
C -----------------------
C
C MAIN    --> GETCOD  --> GETNAM
C         --> OPNFIL
C         --> INITS
C         --> READPR  --> STORAT  --> CALC
C                     --> STORPL  --> PLANE   --> GETCOM
C                                             --> SCPMAT
C                                             --> JACOBI
C                                             --> EIGSRT
C                                 --> GETRMS
C                                 --> WRIPLN
C
C----------------------------------------------------------------------+---


      PROGRAM ANGLEN

      INCLUDE 'anglen.inc'


C---- Read in the code of the Brookhaven file and the structure's resolution
      CALL GETCOD
      IF (IFAIL) GO TO 990

C---- Open data files
      CALL OPNFIL
      IF (IFAIL) GO TO 990

C---- Initialise variables
      CALL INITS

C---- Read through the Brookhaven file
      CALL READPR

C---- Close any open files
      CLOSE(1)
      CLOSE(2)
CHECK v.2.2-->
      CLOSE(3)
CHECK v.2.2<--

 990  CONTINUE
      IF (IFAIL) THEN
CHECK v.2.3-->
C          PRINT*, 'Program terminated with error'
C      ELSE
C          PRINT*, 'Program completed'
          PRINT*, '* Program anglen terminated with error'
      ELSE
          PRINT*, '* Program completed'
CHECK v.2.3<--
      ENDIF
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETCOD  -  Read in the 4-letter Brokhaven code of the protein
C                        in question
C
C----------------------------------------------------------------------+---
 
      SUBROUTINE GETCOD
 
      INCLUDE 'anglen.inc'
 
CHECK v.2.0-->
C      INTEGER       I 
CHECK v.2.0<--
CHECK v.2.0-->
      INTEGER       IEND, ILEN, ISTART
      LOGICAL       IERROR
CHECK v.2.0<--

C---- Initialise variables

CHECK v.2.0-->
C---- Accept 4-letter Brookhaven code
C      PRINT*, 'Enter 4-letter Brookhaven code of protein'
C      READ(*,10) BRCODE
C 10   FORMAT(A4)

C---- Check that the code doesn't contain any blanks
C      DO 100, I = 1, 4
C          IF (BRCODE(I:I).EQ.' ') GO TO 900
C 100  CONTINUE
C      FILNEW = 'p' // BRCODE // '.new'
C      FILLAN = 'p' // BRCODE // '.lan'
CHECK v.2.0<--

CHECK v.2.0-->
C---- Accept name of original .pdb file holding the structure
      PRINT*, 'Enter filename containing coordinates of structure'
      READ(*,10) PDBFIL
10    FORMAT(A)

C---- Peel off directory path and extension
      CALL GETNAM(PDBFIL,ISTART,IEND,IERROR)
      IF (IERROR) GO TO 990

C---- Form names of other files that will be required in default directory
      ILEN = IEND - ISTART + 1
      FILNEW = PDBFIL(ISTART:IEND) // '.new'
      FILLAN = PDBFIL(ISTART:IEND) // '.lan'
CHECK v.3.4.3<--
      FILNRV = PDBFIL(ISTART:IEND) // '.nrv'
CHECK v.3.4.3-->
CHECK v.2.2-->
      FILPLN = PDBFIL(ISTART:IEND) // '.pln'
CHECK v.2.2<--
CHECK v.2.0<--

      GO TO 999

C---- Fatal errors
990   CONTINUE
      IFAIL = .TRUE.
 
999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
 
CHECK v.2.0-->
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
CHECK v.2.0<--
C**************************************************************************
C
C  SUBROUTINE OPNFIL  -   Open data files
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPNFIL

      INCLUDE 'anglen.inc'


C---- Initialise variables
      IFAIL = .FALSE.

C---- Open input .new file
      OPEN(UNIT=1,FILE=FILNEW,STATUS='OLD',ERR=900,FORM='FORMATTED',
CVAX     -     READONLY,
     -     ACCESS='SEQUENTIAL')

C---- Open output file, <filename>.lan
      OPEN(UNIT=2,FILE=FILLAN,STATUS='UNKNOWN',ERR=904,
CVAX     -     CARRIAGECONTROL='LIST',
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')

CHECK v.2.2-->
C---- Open output file, <filename>.pln
      OPEN(UNIT=3,FILE=FILPLN,STATUS='UNKNOWN',ERR=906,
CVAX     -     CARRIAGECONTROL='LIST',
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')
CHECK v.2.2<--

CHECK v.3.4.3-->
CC---- Open output file, <filename>.nrv and write out header records
C      OPEN(UNIT=7,FILE=FILNRV,STATUS='UNKNOWN',ERR=908,
CCVAX     -     CARRIAGECONTROL='LIST',
C     -     FORM='FORMATTED',ACCESS='SEQUENTIAL')
C      WRITE(7,100)
C 100  FORMAT(
C     -    '$ DATA  DISTANCES',//,
C     -    'Chn Resnm Res Atom  Chn Resnm Res Atom  Type  LoBnd UpBnd',
C     -    ' Str_01')
CHECK v.3.4.3<--

      GO TO 999

C---- Fatal errors
900   CONTINUE
CHECK v.2.2-->
C      PRINT*, '*** ERROR. Unable to open data file: ', FILNEW
      PRINT*, '*** ERROR. Unable to open data file: '
      PRINT*, FILNEW, '*'
CHECK v.2.2<--
      PRINT*, '***        Run program CLEAN to create it'
      GO TO 990

904   CONTINUE
CHECK v.2.2-->
C      PRINT*, '*** ERROR. Unable to open output file: ', FILLAN
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILLAN, '*'
CHECK v.2.2<--
      GO TO 990

CHECK v.2.2-->
906   CONTINUE
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILPLN, '*'
      GO TO 990
CHECK v.2.2<--

CHECK v.3.4.3-->
 908  CONTINUE
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILNRV, '*'
      GO TO 990
CHECK v.3.4.3<--

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

      INCLUDE 'anglen.inc'

      CHARACTER*1   A, AL
      CHARACTER*13  BONLEN(NMAIN)
CHECK v.2.2-->
C      INTEGER       ADJUST, I, IBAT, ILET, IMAIN, LEN, LEN1, LEN2, LEN3
      INTEGER       ADJUST, I, IATOM, IBAT, ILET, IMAIN, IPLANE, LEN,
     -              LEN1, LEN2, LEN3
CHECK v.2.2<--

C---- Initialise variables
      A = 'A'
      ADJUST = ICHAR(A) - 1
      IFAIL = .FALSE.
      NRESID = 0
      NXTPOS = 0
      PI = 4.0 *  ATAN(1.0)
      RADDEG = 180.0 / PI

C---- Set up pointers system for residue names
      DO 100, I = 1, MXTYPE
          NXTLET(I) = 0
100   CONTINUE
      DO 200, I = 1, 26
          FSTLET(I) = 0
          LSTLET(I) = 0
200   CONTINUE
      DO 300, I = 1, MXTYPE
          IF (RTYPE(I).EQ.'   ') GO TO 500
          ILET = ICHAR(RTYPE(I)(1:1)) - ADJUST
          IF (ILET.GT.26) ILET = 26
          NXTPOS = NXTPOS + 1
          IF (LSTLET(ILET).EQ.0) THEN
              FSTLET(ILET) = NXTPOS
              LSTLET(ILET) = NXTPOS
          ELSE
              NXTLET(LSTLET(ILET)) = NXTPOS
              LSTLET(ILET) = NXTPOS
          ENDIF
300   CONTINUE

C---- Initialise parameters defining each of the bond lengths and bond
C     angles
 500  CONTINUE

C---- Loop through all the bond lengths and bond angles
      DO 1000, IMAIN = 1, NMAIN

C----     Initialise valriables giving the relative residues from
C         which the atoms defining each bond length/angle come (all
C         zero means they come from the same residue).
          INAME(IMAIN) = 0
          RELPOS(1,IMAIN) = 0
          RELPOS(2,IMAIN) = 0
          RELPOS(3,IMAIN) = 0

C----     Set up array of integer values corresponding to the atom types
C         in BONDAT

C----     Loop through the triplet of atoms defining this bond length/angle
          DO 700, I = 1, 3
              BATNUM(I,IMAIN) = 0

C----         Loop through all the main-chain atoms (as defined by the
C             MATOM array).
              DO 600, IBAT = 1, NMATOM

C----             If this main-chain atom matches that in the BONDAT
C                 array, store its position-number in the BATNUM aray
                  IF (MATOM(IBAT).EQ.BONDAT(I,IMAIN)) THEN
                      BATNUM(I,IMAIN) = IBAT
                  ENDIF
 600          CONTINUE
 700      CONTINUE

C----     Check whether dealing with a bond length or a bond angle (ie if
C         the third atom-name in the current triplet is not blank, it is
C         a bond angle)
          IF (BONDAT(3,IMAIN).NE.' ') THEN
              BANGLE(IMAIN) = .TRUE.
              NBAT(IMAIN) = 3

C----         Determine whether all atoms are in the same residue, or whether
C             they span adjacent residues

C----         If the 3rd atom is N, then first 2 atoms must come from the
C             previous residue, and that is the residue associated with
C             this angle
              IF (BONDAT(3,IMAIN).EQ.'N  ') THEN
                  RELPOS(1,IMAIN) = -1
                  RELPOS(2,IMAIN) = -1
                  INAME(IMAIN) = -1

C----         If 3rd atom is CA, then 1st atom must come from previous
C             residue
              ELSE IF (BONDAT(3,IMAIN).EQ.'CA ') THEN
                  RELPOS(1,IMAIN) = -1

C----         Otherwise, all 3 atoms come from the same residue
              ELSE
                  RELPOS(1,IMAIN) = 0
                  RELPOS(2,IMAIN) = 0
              ENDIF

C----     If dealing with a bond length
          ELSE
              BANGLE(IMAIN) = .FALSE.
              NBAT(IMAIN) = 2

C----         Determine whether the two atoms are in the same residue, or
C             whether they span adjacent residues

C----         If the 2nd atom is N, then first atom must come from the
C             previous residue
              IF (BONDAT(2,IMAIN).EQ.'N  ') THEN
                  RELPOS(1,IMAIN) = -1
              ENDIF
          ENDIF

C----     Form description of this bond length/angle from its component
C         atom names
          LEN1 = INDEX(BONDAT(1,IMAIN),' ') - 1
          LEN2 = INDEX(BONDAT(2,IMAIN),' ') - 1
          LEN3 = INDEX(BONDAT(3,IMAIN),' ') - 1
          IF (LEN1.LT.1) LEN1 = 3
          IF (LEN2.LT.1) LEN2 = 3
          IF (LEN3.LT.1) LEN3 = 3
          BONLEN(IMAIN) = BONDAT(1,IMAIN)
          BONLEN(IMAIN)(LEN1 + 1:LEN1 + 1) = '-'
          BONLEN(IMAIN)(LEN1 + 2:) = BONDAT(2,IMAIN)
          LEN = LEN1 + LEN2 + 1
          IF (BANGLE(IMAIN)) THEN
              BONLEN(IMAIN)(LEN + 1: LEN + 1) = '-'
              BONLEN(IMAIN)(LEN +2:) = BONDAT(3,IMAIN)
              LEN = LEN + LEN3 + 1
          ENDIF

 1000 CONTINUE

C---- Write out the Engh & Huber means and standard deviations for each of
C     the main-chain bond lengths and angles
      DO 1200, IMAIN = 1, NMAIN
          IF (BANGLE(IMAIN)) THEN
              AL = 'A'
          ELSE
              AL = 'L'
          ENDIF
          WRITE(2,1100) AL, BONLEN(IMAIN), ENGNAM(IMAIN),
     -        EXDESC(EXCEPT(IMAIN)), ENGMEA(IMAIN), ENGSTD(IMAIN)
 1100     FORMAT(A1,1X,A10,A17,4X,A16,2F8.3)
 1200 CONTINUE

CHECK v.2.2-->
C---- Count the number of planar atoms in each sidechain
      DO 1280, IPLANE = 1, NPLANE
          DO 1260, IATOM = 1, NPLATM
              IF (ANAME(IATOM + 1,IPLANE).NE.'   ')
     -            NPATOM(IPLANE) = IATOM
 1260     CONTINUE
 1280 CONTINUE

C---- Write out the descriptions of the planar groups
      DO 1300, IPLANE = 1, NPLANE
          WRITE(2,1290) 'P', NPATOM(IPLANE), 'Planar group     ',
     -        ANAME(1,IPLANE)
 1290     FORMAT(A1,1X,I4,6X,A17,4X,A3)
 1300 CONTINUE
CHECK v.2.2<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE READPR  -  Read through appropriate Brookhaven file and
C                        perform required distance and angle calculations
C
C----------------------------------------------------------------------+---

      SUBROUTINE READPR

      INCLUDE 'anglen.inc'

      CHARACTER*1   CHAIN
      CHARACTER*3   ATNAME, RSNAME
      CHARACTER*5   SEQNO
      CHARACTER*6   IDENT, INSEQ, LSTSEQ
      CHARACTER*80  IREC
      INTEGER       FRESID, I, IBAT, IRESID, LINE, RESIDC, RESIDN
      LOGICAL       NEWCHN
      REAL          COORDC(3), COORDN(3), COORDS(3), DIST2, PEPBN2

      PARAMETER    (PEPBN2 = 2.5 * 2.5)

C---- Initialise variables
      IFAIL = .FALSE.
      IRESID = 0
CHECK v.2.1-->
      LINE = 0
CHECK v.2.1<--
      LSTSEQ = '^^^^^^'
      NEWCHN = .TRUE.
      RESIDC = 0
      RESIDN = 0

C---- Loop through Brookhaven file, searching for main-chain atoms
100   CONTINUE
          READ(1,150,END=800,ERR=904) IREC
150       FORMAT(A)
          LINE = LINE + 1
          IDENT = IREC(1:6)

C----     If record is an atom record, check its atom type
          IF (IDENT.EQ.'ATOM  ') THEN

C----         Pick off atom name, residue names, sequence number, and B-value
              ATNAME = IREC(14:16)
              RSNAME = IREC(18:20)
              CHAIN = IREC(22:22)
              INSEQ = IREC(22:27)
              SEQNO = IREC(23:27)

C----         If new chain encountered, then flag a chain-break
CHECK v.3.0.1-->
C170           CONTINUE
CHECK v.3.0.1<--
              IF (INSEQ(1:1).NE.LSTSEQ(1:1)) THEN
                  NEWCHN = .TRUE.
                  LSTSEQ = '^^^^^^'
              ENDIF

C----         If sequence number has changed, increment residue-count
              IF (INSEQ.NE.LSTSEQ) THEN
                  IRESID = IRESID + 1
                  IF (IRESID.GT.MXONE) GO TO 906
              ENDIF
              LSTSEQ = INSEQ

C             Check for chain break by looking for excessively long peptide
C             bonds C(i-1)->N(i)

C----         If atom is a main-chain C or N, store its coords if not
C             already done so
              IF (IRESID.NE.RESIDC .AND. (ATNAME.EQ.'C  ' .OR.
     -            ATNAME.EQ.'H  ')) THEN
                  READ(IREC,180,ERR=902) (COORDC(I), I = 1, 3)
180               FORMAT(30X,3F8.0)
                  RESIDC = IRESID
              ELSE IF (IRESID.NE.RESIDN .AND. (ATNAME.EQ.'N  ' .OR.
     -            ATNAME.EQ.'D  ')) THEN
                  READ(IREC,180,ERR=902) (COORDN(I), I = 1, 3)
                  RESIDN = IRESID

C----             If have coords for C atom from previous residue, check
C                 whether distance between it and the current N indicates
C                 a chain break here
                  IF (.NOT.NEWCHN .AND. RESIDN.EQ.RESIDC + 1 .AND.
     -                RESIDC.NE.0) THEN
                      DIST2 = (COORDN(1) - COORDC(1)) ** 2
     -                    + (COORDN(2) - COORDC(2)) ** 2
     -                    + (COORDN(3) - COORDC(3)) ** 2
                      IF (DIST2.GT.PEPBN2) THEN
                          NEWCHN = .TRUE.
                      ENDIF
                  ELSE
                      NEWCHN = .TRUE.
                  ENDIF
              ENDIF

C----         If have a chain-break here, store the number of the first
C             residue in the new chain
              IF (NEWCHN) FRESID = IRESID

CHECK v.2.2-->
C----         Retrieve the atomic coordinates
              READ(IREC,180,ERR=902) (COORDS(I), I = 1, 3)
CHECK v.2.2<--

C----         Determine whether the atom is a main-chain atom
              IBAT = 0
              DO 400, I = 1, NMATOM
                  IF (ATNAME.EQ.MATOM(I)) IBAT = I
 400          CONTINUE

C----         If atom is a main-chain atom call the calculation routine
              IF (IBAT.NE.0) THEN
CHECK v.2.2-->
C                  READ(IREC,180,ERR=902) (COORDS(I), I = 1, 3)
CHECK v.2.2<--
                  CALL STORAT(IRESID,RSNAME,CHAIN,SEQNO,IBAT,COORDS,
     -                NEWCHN)
              ENDIF

CHECK v.2.2-->
C             See if atom is part of a planar group
              CALL STORPL(IRESID,RSNAME,ATNAME,CHAIN,SEQNO,COORDS,
     -                NEWCHN)
CHECK v.2.2<--

C----     If TER record, then have encountered end of chain
          ELSE IF (IDENT.EQ.'TER   ') THEN
              LSTSEQ = '^^^^^^'
          ENDIF
      GO TO 100

C---- End of Brookhaven file reached, calculate last set of lengths and angles
800   CONTINUE
      IBAT = 0
      NEWCHN = .TRUE.
      CALL STORAT(IRESID,RSNAME,CHAIN,SEQNO,IBAT,COORDS,NEWCHN)

CHECK v.2.2-->
C---- Calculate planarity for last residue, if required
CHECK v.3.2-->
C      CALL STORPL(-999,' ',' ',' ',SEQNO,COORDS,NEWCHN)
      CALL STORPL(-999,'   ','   ',' ',SEQNO,COORDS,NEWCHN)
CHECK v.3.2<--
CHECK v.2.2<--

      GO TO 999

C---- Fatal errors
902   CONTINUE
CHECK v.2.1.3-->
C      PRINT*, 'Error in coords for atom: ', ATNAME, '-', RSNAME, '-',
C     -    INSEQ
      PRINT*, '**** Error in coords for atom: ', ATNAME, '-', RSNAME,
     -    '-', INSEQ
CHECK v.2.1.3<--
      GO TO 100

904   CONTINUE
CHECK v.2.1.3-->
C      PRINT*, 'Data error reading file:  [', FILNEW, ']'
C      PRINT*, '    Line: ', LINE
      PRINT*, '**** Data error reading file:  [', FILNEW, ']'
      PRINT*, '    Line: ', LINE
CHECK v.2.1.3<--
      GO TO 990

906   CONTINUE
CHECK v.2.1.3-->
C      PRINT*, 'Maximum number of residues per protein exceeded:',
C     -    MXONE
      PRINT*, '**** Maximum number of residues per protein exceeded:',
     -    MXONE
CHECK v.2.1.3<--
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE STORAT  -  Store atom coordinates for current atom and, when
C                        have required bond length or angle, calculate its
C                        value
C
C----------------------------------------------------------------------+---

      SUBROUTINE STORAT(IRESID,RSNAME,INCHN,INSEQ,IBAT,COORDS,NEWCHN)

CHECK v.2.0.1-->
      SAVE
CHECK v.2.0.1<--

      INCLUDE 'anglen.inc'

      CHARACTER*1   CHAIN(2), INCHN
      CHARACTER*3   RNAME(2), RSNAME, RTEST
CHECK v.3.4.3-->
      CHARACTER*4   ATNAME(3)
CHECK v.3.4.3<--
      CHARACTER*5   INSEQ, SEQNO(2)
      INTEGER       CURRES, CURRNT, FOUND, I, IBAT, IMAIN, IPOINT,
     -              IPOS, IRESID, ISLOT, LAST, NONSP, SRESID(2)
CHECK v.3.2-->
C      LOGICAL       FILLED(2,NMATOM), FIRST, MATCH, NEWCHN, THISUN
CHECK v.3.5.1-->
C      LOGICAL       FILLED(2,NMATOM), FIRST, MATCH, NEWCHN, SPANER,
C     -              THISUN
      LOGICAL       FILLED(2,NMATOM), FIRST, INVALD, MATCH, NEWCHN,
     -              SPANER, THISUN
CHECK v.3.5.1<--
CHECK v.3.2<--
CHECK v.3.4.3-->
C      REAL          COORDS(3), SAVXYZ(3,3), VALUE, XYZ(3,2,NMATOM)
      REAL          COORDS(3), LOWER, SAVXYZ(3,3), STDEVS,
     -              UPPER, VALUE, XYZ(3,2,NMATOM)
CHECK v.3.4.3<--

CHECK v.2.0.1-->
      DATA CURRES, CURRNT, LAST / 0, 1, 2 /
CHECK v.2.0.1<--
      DATA FIRST  / .TRUE. /

CHECK v.2.1-->
C---- Initialise variables on first entry into routine
      IF (FIRST) THEN
          CHAIN(1) = ' '
          CHAIN(2) = ' '
          DO 100, I = 1, NMATOM      
              FILLED(1,I) = .FALSE.
              FILLED(2,I) = .FALSE.
              DO 50, IPOS = 1, 3
                  XYZ(IPOS,1,I) = 0.0
                  XYZ(IPOS,2,I) = 0.0
 50           CONTINUE
 100      CONTINUE
          RNAME(1) = ' '
          RNAME(2) = ' '
          SEQNO(1) = ' '
          SEQNO(2) = ' '
          SRESID(1) = 0
          SRESID(2) = 0
      ENDIF
CHECK v.2.1<--

C---- If residue has changed, calculate all possible angles for previous
C     residue and then re-initialise pointers
      IF (.NOT.FIRST .AND. (NEWCHN .OR. IRESID.NE.CURRES)) THEN

C----     Loop through all bond lengths and angles to be tried
          DO 200, IMAIN = 1, NMAIN

C----         Test whether the current residue is appropriate for the
C             current bond length/angle

C----         Initialise variables
              THISUN = .TRUE.
              NONSP = 0

C----         Store the first exception definition (which will be one
C             of: ALL, NOT, or blank)
              RTEST = EXTYPE(1,EXCEPT(IMAIN))

C----         If the appropriate residue name belongs to the previous
C             residue, set the name-pointer to that
              IF (INAME(IMAIN).EQ.-1) THEN
                  IPOINT = LAST

C----         Otherwise, set the name-pointer to the current residue
              ELSE
                  IPOINT = CURRNT
              ENDIF

C----         If the name of the residue we're looking at is blank
C             then skip the rest of the routines altogether
              IF (RNAME(IPOINT).EQ.'   ') THEN
                  THISUN = .FALSE.
                  GO TO 140
              ENDIF

C----         If the exception test matches ALL residues, then keep
              IF (RTEST.EQ.'ALL') GO TO 140

C----         If the excpetion test is blank, then the residue names
C             that follow in the EXTYPE table are the required ones
              IF (RTEST.EQ.'   ') MATCH = .TRUE.

C----         If the excpetion test is NOT, then the residue names
C             that follow in the EXTYPE table are the exceptions
              IF (RTEST.EQ.'NOT') MATCH = .FALSE.

C----         Check through the residue names in the remainder of the
C             current exception definition to see if the current
C             residue appears
              DO 120, I = 2, MAXEXC
                  RTEST = EXTYPE(I,EXCEPT(IMAIN))

C----             If the residue name is in the table, then can
C                 proceed if we're looking for an exact match, or can
C                 reject the current bond length/angle if we're looking
C                 for a NOT match
                  IF (RTEST.EQ.RNAME(IPOINT)) THEN
                      IF (MATCH) THEN
                          GO TO 140
                      ELSE
                          THISUN = .FALSE.
                          GO TO 140
                      ENDIF
                  ENDIF
 120          CONTINUE

C----         If looking for the given residue in the table and have got
C             to here, then a match was not found
              IF (MATCH) THISUN = .FALSE.

C----         If residue is appropriate, see if have all the required
C             atoms
 140          CONTINUE
              IF (THISUN) THEN

C----             Loop through the 2 or 3 atoms defining current bond
C                 length/angle to determine if all are present
CHECK v.3.2-->
                  SPANER = .FALSE.
CHECK v.3.2<--
                  FOUND = 0
                  DO 150, I = 1, NBAT(IMAIN)

C----                 Check whether the atom at this position in the
C                     pair/triplet comes from the previous residue or
C                     from the current one
                      IF (RELPOS(I,IMAIN).EQ.-1) THEN
                          IPOS = LAST
CHECK v.3.2-->
                          SPANER = .TRUE.
CHECK v.3.2<--
                      ELSE
                          IPOS = CURRNT
                      ENDIF

C----                 If have the appropriate atom, save its coordinates
                      IF (FILLED(IPOS,BATNUM(I,IMAIN))) THEN
                          FOUND = FOUND + 1
                          SAVXYZ(1,I) = XYZ(1,IPOS,BATNUM(I,IMAIN))
                          SAVXYZ(2,I) = XYZ(2,IPOS,BATNUM(I,IMAIN))
                          SAVXYZ(3,I) = XYZ(3,IPOS,BATNUM(I,IMAIN))
CHECK v.3.4.3-->
                          ATNAME(I) = MATOM(BATNUM(I,IMAIN))
CHECK v.3.4.3<--
                      ENDIF
 150              CONTINUE

C----             If have all the atoms, calculate distance/angle
                  IF (FOUND.EQ.NBAT(IMAIN)) THEN
CHECK v.3.5.1-->
C                      CALL CALC(BANGLE(IMAIN),SAVXYZ,RADDEG,VALUE)
                      CALL CALC(BANGLE(IMAIN),SAVXYZ,RADDEG,VALUE,
     -                    INVALD)

C----                 If VERY bizarre value calculated, show bond/angle
C                     involved
                      IF (INVALD) THEN

C----                     Show residue(s) involved
                          IF (SPANER) THEN
                              PRINT*, '***          ',
     -                            RNAME(LAST), SEQNO(LAST), CHAIN(LAST),
     -                            RNAME(CURRNT), SEQNO(CURRNT),
     -                            CHAIN(CURRNT)
                          ELSE
                              PRINT*, '***          ',
     -                            RNAME(CURRNT), SEQNO(CURRNT),
     -                            CHAIN(CURRNT)
                          ENDIF

C----                     Show atom names involved
                          IF (BANGLE(IMAIN)) THEN
                              PRINT*, '***    Atoms ',
     -                            ATNAME(1), ATNAME(2), ATNAME(3)
                          ELSE
                              PRINT*, '***    Atoms ',
     -                            ATNAME(1), ATNAME(2)
                          ENDIF
                      ENDIF
CHECK v.3.5.1<--

CHECK v.3.4.3-->
C----                 Calculate difference from mean value
                      STDEVS = 2.0
                      LOWER = ENGMEA(IMAIN) - STDEVS * ENGSTD(IMAIN)
                      UPPER = ENGMEA(IMAIN) + STDEVS * ENGSTD(IMAIN)
CHECK v.3.4.3<--

CHECK v.3.2-->
CC----                 Write value to output file
C                      WRITE(2,180) SRESID(IPOINT), CHAIN(IPOINT),
C     -                    SEQNO(IPOINT), RNAME(IPOINT), IMAIN, VALUE
C 180                  FORMAT(I6,A1,A5,A3,I2,F9.4)
C
C----                 If the bond length/angle spans two residues, write
C                     out the value, showing both residues
                      IF (SPANER) THEN
                          WRITE(2,180) SRESID(IPOINT), CHAIN(IPOINT),
     -                        SEQNO(IPOINT), RNAME(IPOINT), IMAIN,
     -                        VALUE, CHAIN(LAST), SEQNO(LAST),
     -                        RNAME(LAST), CHAIN(CURRNT),
     -                        SEQNO(CURRNT), RNAME(CURRNT)
 180                      FORMAT(I6,A1,A5,A3,I2,F9.4,2(1X,A1,A5,A3))
CHECK v.3.4.3-->
C----                     If bond length, write out to "restraint
C                         violations" file
C                          IF (.NOT.BANGLE(IMAIN)) THEN
C                              WRITE(7,184) CHAIN(IPOINT),
C     -                            RNAME(IPOINT), SEQNO(IPOINT),
C     -                            ATNAME(2), CHAIN(LAST), RNAME(LAST),
C     -                            SEQNO(LAST), ATNAME(1), LOWER, UPPER,
C     -                            VALUE
C 184                          FORMAT(2('_',A1,1X,A3,2X,A5,1X,A4,2X),
C     -                            '  0 ',3F7.3)
C                          ENDIF
CHECK v.3.4.3<--

C----                 Otherwise, just write out the current residue
                      ELSE
                          WRITE(2,180) SRESID(IPOINT), CHAIN(IPOINT),
     -                        SEQNO(IPOINT), RNAME(IPOINT), IMAIN,
     -                        VALUE
CHECK v.3.4.3-->
C----                     If bond length, write out to "restraint
C                         violations" file
C                          IF (.NOT.BANGLE(IMAIN)) THEN
C                              WRITE(7,184) CHAIN(IPOINT),
C     -                            RNAME(IPOINT), SEQNO(IPOINT),
C     -                            ATNAME(1), CHAIN(IPOINT),
C     -                            RNAME(IPOINT), SEQNO(IPOINT), 
C     -                            ATNAME(2), LOWER, UPPER, VALUE
C                          ENDIF
CHECK v.3.4.3<--
                      ENDIF
CHECK v.3.2<--
                  ENDIF
              ENDIF
 200      CONTINUE

C----     Set current residue as previous residue and reinitialise pointers
          CURRNT = 3 - CURRNT
          LAST = 3 - LAST
          DO 300, ISLOT = 1, NMATOM
              FILLED(CURRNT,ISLOT) = .FALSE.
 300      CONTINUE
      ENDIF

C---- If at start of new chain or at break in chain, then initialise variables
      IF (FIRST .OR. NEWCHN) THEN

C----     Initialise flags showing whether have the 5 main-chain coords
C         for current and prior residues
          DO 400, ISLOT = 1, NMATOM
              FILLED(1,ISLOT) = .FALSE.
              FILLED(2,ISLOT) = .FALSE.
 400      CONTINUE
          CURRNT = 1
          LAST = 2
          CHAIN(LAST) = ' '
          CHAIN(CURRNT) = ' '
          FIRST = .FALSE.
          NEWCHN = .FALSE.
          CURRES = IRESID
          RNAME(LAST) = '   '
          RNAME(CURRNT) = '   '
          SEQNO(LAST) = '   '
          SEQNO(CURRNT) = '   '
      ENDIF

C---- Store coordinates of current atom
      IF (IBAT.GT.0) THEN
          CURRES = IRESID
          CHAIN(CURRNT) = INCHN
          FILLED(CURRNT,IBAT) = .TRUE.
          RNAME(CURRNT) = RSNAME
          SEQNO(CURRNT) = INSEQ
          SRESID(CURRNT) = IRESID
          XYZ(1,CURRNT,IBAT) = COORDS(1)
          XYZ(2,CURRNT,IBAT) = COORDS(2)
          XYZ(3,CURRNT,IBAT) = COORDS(3)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE CALC  -  Calculate the required bond length/angle and
C                      return the value
C
C----------------------------------------------------------------------+---

CHECK v.3.5.1-->
C      SUBROUTINE CALC(BANGLE,XYZ,RADDEG,VALUE)
      SUBROUTINE CALC(BANGLE,XYZ,RADDEG,VALUE,INVALD)
CHECK v.3.5.1<--

CVAX      IMPLICIT NONE

      INTEGER       ICOORD
CHECK v.3.5.1-->
C      LOGICAL       BANGLE
      LOGICAL       BANGLE, INVALD
CHECK v.3.5.1<--
      REAL          COSANG, DIFF1, DIFF2, DIFF3, DIST1, DIST2, DIST3,
     -              RADDEG, XYZ(3,3), VALUE

CHECK v.3.5.1-->
C---- Initialise variables
      INVALD = .FALSE.
      VALUE = 0.0
CHECK v.3.5.1<--

C---- Calculate the value of the angle
      IF (BANGLE) THEN
          
C----     Calculate distances between atoms 1->2, 2->3, 3->1
          DIST1 = 0
          DIST2 = 0
          DIST3 = 0
          DO 300, ICOORD = 1, 3
              DIFF1 = XYZ(ICOORD,2) - XYZ(ICOORD,1)
              DIFF2 = XYZ(ICOORD,3) - XYZ(ICOORD,2)
              DIFF3 = XYZ(ICOORD,1) - XYZ(ICOORD,3)
              DIST1 = DIST1 + DIFF1 * DIFF1
              DIST2 = DIST2 + DIFF2 * DIFF2
              DIST3 = DIST3 + DIFF3 * DIFF3
300       CONTINUE

C----     Calculate angle
          IF (DIST1.GT.0.0 .AND. DIST2.GT.0.0 .AND.
     -        DIST3.GE.0.0) THEN
              DIST1 = SQRT(DIST1)
              DIST2 = SQRT(DIST2)
              DIST3 = SQRT(DIST3)
              COSANG = (DIST3 * DIST3 - DIST1 * DIST1
     -            - DIST2 * DIST2) / (-2.0 * DIST1 * DIST2)
              IF (ABS(COSANG).LE.1.0) THEN
                  VALUE = RADDEG * ACOS(COSANG)
              ENDIF
          ENDIF

C---- Calculate value of bond length
      ELSE
          DIST1 = 0
          DO 400, ICOORD = 1, 3
              DIFF1 = XYZ(ICOORD,2) - XYZ(ICOORD,1)
              DIST1 = DIST1 + DIFF1 * DIFF1
400       CONTINUE
          VALUE = SQRT(DIST1)
      ENDIF

CHECK v.3.5.1-->
C---- If value outside printable range, then set to zero
      IF (VALUE.GT.9999.99 .OR. VALUE.LT.-999.99) THEN
          IF (BANGLE) THEN
              PRINT*, '*** Warning. Very strange angle calculated:',
     -            VALUE
          ELSE
              PRINT*, '*** Warning. Very strange length calculated:',
     -            VALUE
          ENDIF

C----     Set flag that value is invalid
          INVALD = .TRUE.
          VALUE = 0.0
      ENDIF
CHECK v.3.5.1<--

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.2.2-->
C**************************************************************************
C
C  SUBROUTINE STORPL  -  Store atom coordinates for current atom and, when
C                        have required set for a planar group, calculate
C                        the best-fit plane and atomic deviation from it
C
C----------------------------------------------------------------------+---

      SUBROUTINE STORPL(IRESID,RSNAME,ATNAM,INCHN,INSEQ,COORDS,NEWCHN)

      SAVE

      INCLUDE 'anglen.inc'

      CHARACTER*1   CHAIN, INCHN
      CHARACTER*3   ATNAM, RSNAME, SNAME
      CHARACTER*5   INSEQ, SEQNO
CHECK v.3.2-->
C      INTEGER       CURRES, CURRNT, IATOM, IP, IRESID, LAST, NGOT,
C     -              ORDVEC(3), PLATOM, PLRES, SRESID
      INTEGER       IATOM, IP, IRESID, NGOT, ORDVEC(3), PLATOM, PLRES,
     -              SRESID
CHECK v.3.2<--
      LOGICAL       FILLED(NPLATM), FIRST, NEWCHN, PERROR
      REAL          CENTRE(3), COORDS(3), EIGVEC(3,3), NORMAL(3),
CHECK v.3.0.1-->
C     -              RMSDIF, TCOORD(NPLATM), XYZ(3,NPLATM)
     -              RMSDIF, TCOORD(3,NPLATM), XYZ(3,NPLATM)
CHECK v.3.0.1<--

CHECK v.3.2-->
C      DATA CURRES, CURRNT, LAST / 0, 1, 2 /
CHECK v.3.2<--
      DATA FIRST  / .TRUE. /

C---- Initialise variables on first entry into routine
      PERROR = .FALSE.
      IF (FIRST) THEN

C----     Initialise the residue details
          CHAIN = ' '
          FIRST = .FALSE.
CHECK v.3.0.1-->
C          NEWCHN = .TRUE.
CHECK v.3.0.1<--
          NGOT = 0
          DO 300, IATOM = 1, NPLATM
              FILLED(IATOM) = .FALSE.
              XYZ(1,IATOM) = 0.0
              XYZ(2,IATOM) = 0.0
              XYZ(3,IATOM) = 0.0
 300      CONTINUE
CHECK v.3.5.5-->
          PLRES = 0
CHECK v.3.5.5<--
          SEQNO = ' '
          SNAME = ' '
          SRESID = -9999
      ENDIF

C---- If residue has changed, may need to perform planarity calculations
C     on previous residue
      IF (NEWCHN .OR. IRESID.NE.SRESID) THEN

C----     Check that have at least 3 atoms
          IF (NGOT.GT.3) THEN

C----         Check that have all the atoms required
              NGOT = 0
              DO 400, IATOM = 1, NPATOM(PLRES)
                  IF (FILLED(IATOM)) NGOT = NGOT + 1
 400          CONTINUE

C----         If have all the planar atoms for this residue, then
C             calculate a best-fit plane
              IF (NGOT.EQ.NPATOM(PLRES)) THEN

C----             Calculate best-fit plane
                  CALL PLANE(XYZ,NGOT,CENTRE,NORMAL,EIGVEC,ORDVEC,
     -                 PERROR)
                  IF (.NOT.PERROR) THEN

C----                 Calculate RMS distance of all points from best-fit
C                     plane
CHECK v.3.2-->
C                      CALL GETRMS(XYZ,NGOT,CENTRE,NORMAL,RMSDIF)
                      CALL GETRMS(XYZ,NGOT,NORMAL,RMSDIF)
CHECK v.3.2<--

C----                 Transform coordinates of atoms and write out
                      CALL WRIPLN(XYZ,NGOT,EIGVEC,ORDVEC,SRESID,SEQNO,
     -                    SNAME,CHAIN,RMSDIF,TCOORD,ANAME(2,PLRES),
     -                    PDBFIL,PLRES)
                  ENDIF
              ENDIF
          ENDIF

C----     Reinitalise all the variables
          NGOT = 0
          DO 500, IATOM = 1, NPLATM
              FILLED(IATOM) = .FALSE.
              XYZ(1,IATOM) = 0.0
              XYZ(2,IATOM) = 0.0
              XYZ(3,IATOM) = 0.0
 500      CONTINUE
      ENDIF

C---- Check whether the residue is one that has a planar group
      IP = 0
      PLRES = 0
 600  CONTINUE
          IP = IP + 1
          IF (RSNAME.EQ.ANAME(1,IP)) PLRES = IP
      IF (PLRES.EQ.0 .AND. IP.LT.NPLANE) GO TO 600

C---- If not a planar group, then return
      IF (PLRES.EQ.0) GO TO 999

C---- Check whether this is one of the atoms in the planar group
      IP = 0
      PLATOM = 0
 700  CONTINUE
          IP = IP + 1
          IF (ATNAM.EQ.ANAME(IP + 1,PLRES)) PLATOM = IP
CHECK v.3.0.1-->
C      IF (PLATOM.EQ.0 .AND. IP.LT.NPLATM .AND.
C     -    ANAME(IP + 1,PLRES).NE.'   ') GO TO 700
      IF (PLATOM.EQ.0 .AND. IP.LT.NPLATM) GO TO 700
CHECK v.3.0.1<--

C---- If not one of the atoms in the planar group, then return
      IF (PLATOM.EQ.0) GO TO 999

C---- Store coordinates of current atom
      CHAIN = INCHN
      FILLED(PLATOM) = .TRUE.
CHECK v.3.0.1-->
C      NEWCHN = .FALSE.
CHECK v.3.0.1<--
      NGOT = NGOT + 1
      SEQNO = INSEQ
      SNAME = RSNAME
      SRESID = IRESID
      XYZ(1,PLATOM) = COORDS(1)
      XYZ(2,PLATOM) = COORDS(2)
      XYZ(3,PLATOM) = COORDS(3)

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE PLANE  -  Calculate a best-fit plane through the supplied
C                       data points
C
C----------------------------------------------------------------------+---

      SUBROUTINE PLANE(COORD,NPOINT,CENTRE,NORMAL,EIGVEC,ORDVEC,IFAIL)

      INTEGER       NPOINT, NROT, ORDVEC(3)
      LOGICAL       IFAIL
      REAL          CENTRE(3), COORD(3,NPOINT), D, EIGVEC(3,3),
     -              EIGVAL(3), MATRIX(3,3), NORMAL(3)

C---- Calculate centre of mass of points and adjust coords such that
C     CofM is at the origin
      CALL GETCOM(COORD,NPOINT,CENTRE)

C---- Calculate scalar products matrix
      CALL SCPMAT(COORD,NPOINT,MATRIX)

C---- Calculate eigen values and eigenvectors of matrix
      CALL JACOBI(MATRIX,3,3,EIGVAL,EIGVEC,IFAIL,NROT)
      IF (IFAIL) GO TO 900

C---- Sort the eigen vectors into descending order of eigen value
      CALL EIGSRT(EIGVAL,ORDVEC)

C---- Compute parameters of best-fit plane
      NORMAL(1) = EIGVEC(1,ORDVEC(3))
      NORMAL(2) = EIGVEC(2,ORDVEC(3))
      NORMAL(3) = EIGVEC(3,ORDVEC(3))
      D = NORMAL(1) * CENTRE(1) + NORMAL(2) * CENTRE(2)
     -    + NORMAL(3) * CENTRE(3)

      GO TO 999

C---- Error messages
 900  CONTINUE
      PRINT*, '**** Error. Failed to find eigenvalues of matrix'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE GETCOM  -  Calculate centre of mass of supplied data points
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCOM(COORD,NPOINT,CENTRE)

      INTEGER       IPT, NPOINT
      REAL          CENTRE(3), COORD(3,NPOINT), XSUM, YSUM, ZSUM

C---- Initialise values
      XSUM = 0.0
      YSUM = 0.0
      ZSUM = 0.0

C---- Loop through the data points
      DO 300, IPT = 1, NPOINT
          XSUM = XSUM + COORD(1,IPT)
          YSUM = YSUM + COORD(2,IPT)
          ZSUM = ZSUM + COORD(3,IPT)
300   CONTINUE

C---- Calculate centre of mass and adjust all coordinates accordingly
      CENTRE(1) = XSUM / REAL (NPOINT)
      CENTRE(2) = YSUM / REAL (NPOINT)
      CENTRE(3) = ZSUM / REAL (NPOINT)
      DO 400, IPT = 1, NPOINT
          COORD(1,IPT) = COORD(1,IPT) - CENTRE(1)
          COORD(2,IPT) = COORD(2,IPT) - CENTRE(2)
          COORD(3,IPT) = COORD(3,IPT) - CENTRE(3)
400   CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE SCPMAT  -  Calculate scalar products matrix from supplied
C                        coordinates
C
C----------------------------------------------------------------------+---

      SUBROUTINE SCPMAT(COORD,NPOINT,MATRIX)

      INTEGER       IPT, NPOINT
      REAL          COORD(3,NPOINT), MATRIX(3,3), XXSUM, XYSUM, XZSUM,
     -              YYSUM, YZSUM, ZZSUM


C---- Initialise values
      XXSUM = 0.0
      XYSUM = 0.0
      XZSUM = 0.0
      YYSUM = 0.0
      YZSUM = 0.0
      ZZSUM = 0.0

C---- Loop through all the data points
      DO 500, IPT = 1, NPOINT
          XXSUM = XXSUM + COORD(1,IPT) * COORD(1,IPT)
          XYSUM = XYSUM + COORD(1,IPT) * COORD(2,IPT)
          XZSUM = XZSUM + COORD(1,IPT) * COORD(3,IPT)
          YYSUM = YYSUM + COORD(2,IPT) * COORD(2,IPT)
          YZSUM = YZSUM + COORD(2,IPT) * COORD(3,IPT)
          ZZSUM = ZZSUM + COORD(3,IPT) * COORD(3,IPT)
500   CONTINUE

C---- Form the scalar products matrix
      MATRIX(1,1) = XXSUM
      MATRIX(1,2) = XYSUM
      MATRIX(1,3) = XZSUM
      MATRIX(2,1) = XYSUM
      MATRIX(2,2) = YYSUM
      MATRIX(2,3) = YZSUM
      MATRIX(3,1) = XZSUM
      MATRIX(3,2) = YZSUM
      MATRIX(3,3) = ZZSUM

      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE EIGSRT  -  Sort the eigen vectors into descending order
C                        of eigen value
C
C----------------------------------------------------------------------+---

      SUBROUTINE EIGSRT(EIGVAL,ORDVEC)

      INTEGER       ORDVEC(3), I, ISWAP, J, K
      REAL          EIGVAL(3)

C---- Bubble sort the eigenvalues
      ORDVEC(1) = 1
      ORDVEC(2) = 2
      ORDVEC(3) = 3

C---- Repeate the sort twice
      DO 580, K = 1, 2
          DO 560, J = 1, 2
              DO 540, I = J + 1, 3
                  IF (EIGVAL(ORDVEC(I)).GT.EIGVAL(ORDVEC(J))) THEN
                      ISWAP = ORDVEC(I)
                      ORDVEC(I) = ORDVEC(J)
                      ORDVEC(J) = ISWAP
                  ENDIF
 540          CONTINUE
 560      CONTINUE
 580  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C     ******************************************************************
      SUBROUTINE JACOBI(A,N,NP,D,V,IFAIL,NROT)
C     ******************************************************************
C
C     THIS SUBROUTINE COMPUTES ALL EIGENVALUES AND VECTORS OF A REAL
C     SYMMETRIC MATRIX A USING THE JACOBI METHOD (Taken from Numerical
C     Recipes - Press et al.)
C
C     PARAMETERS:
C
C     A     = REAL SYMMETRIC MATRIX, WHICH IS OF SIZE N BY N, STORED IN
C             A PHYSICAL NP BY NP ARRAY. ON EXIT, ELEMENTS OF A ABOVE
C             THE DIAGONAL ARE DESTROYED.
C     N     = ON ENTRY, ORDER OF MATRIX A FOR WHICH COMPUTATIONS ARE
C             DONE. UNCHANGED ON EXIT.
C     NP    = PHYSICAL ARRAY SUBSCRIPTS OF MATRIX A. UNCHANGED ON EXIT.
C     D     = REAL ARRAY OF PHYSICAL DIMENSION NP WHICH, ON EXIT,
C             CONTAINS EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
C     V     = REAL MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS
C             AS A WHOSE COLUMNS CONTAIN, ON EXIT, THE NORMALIZED
C             EIGENVECTORS OF A.
C     IFAIL = FAILURE INDICATOR. THIS IS SET TO 0 FOR NORMAL RETURN, AND
C             TO 1 WHEN MORE THAN 50 * N(N-1)/2 ROTATIONS ARE NEEDED.
C     NROT  = COUNTER WHICH, ON EXIT, CONTAINS THE NUMBER OF ROTATIONS.
C
C     DURING COMPUTATION:
C
C     Z    = VECTOR WHICH ACCUMULATES TERMS OF FORM TA(PQ)
C     B    = VECTOR WHICH CONTAINS CURRENT EIGENVALUES OF A
C
C     OPERATION COUNT:
C
C     TYPICAL MATRICES REQUIRE 6 TO 10 SWEEPS OF N*(N-1)/2 JACOBI
C     ROTATIONS TO ACHIEVE CONVERGENCE, OR 3*N**2 TO 5*N**2
C     ROTATIONS. EACH ROTATION REQUIRES OF ORDER 6*N OPERATIONS, SO
C     THE TOTAL LABOR IS OF THE ORDER 18*N**3 TO 30*N**3.
C
C
      REAL              EMIN
      PARAMETER        (EMIN=1.0E-8)

      INTEGER          NMAX
      PARAMETER        (NMAX=100)

      INTEGER          I, IP, IQ, J, N, NP, NROT
      LOGICAL          IFAIL
      REAL             A(NP,NP), B(NMAX), C, D(NP), G, H, V(NP,NP), S,
     -                 SM, T, TAU, THETA, TRESH, Z(NMAX)

      IFAIL=.FALSE.

C
C   JACOBI METHOD NOT MOST EFFICIENT METHOD FOR MATRICES OF ORDER
C   GREATER THAN 10 ; MAXIMUM ALLOWED PRESET BY NMAX
C
      IF (N.GT.NMAX) THEN
          WRITE(6,10)
10        FORMAT(/' *** MAXIMUM NUMBER OF ELEMENTS EXCEEDED IN ',
     -    'SUBROUTINE JACOBI; CHANGE NMAX - RUN TERMINATED')
          STOP
      ENDIF
C
C   INITIALISE V TO THE IDENTITY MATRIX
C
      DO 12 IP=1,N
          DO 11 IQ=1,N
              V(IP,IQ)=0.0
11        CONTINUE
          V(IP,IP)=1.0
12    CONTINUE
C
C   INITIALISE B AND D TO THE DIAGONAL OF A AND Z TO 0.0
C
      DO 13 IP=1,N
          B(IP)=A(IP,IP)
          D(IP)=B(IP)
          Z(IP)=0.0
13    CONTINUE
C
C   START OF SWEEPS
C
      NROT=0
      DO 24 I=1,50
C
C---SUM OFF-DIAGONAL ELEMENTS OF MATRIX A
C
          SM=0.0
          DO 15 IP=1,N-1
              DO 14 IQ=IP+1,N
                  SM=SM+ABS(A(IP,IQ))
 14           CONTINUE
 15       CONTINUE
C
C---NORMAL RETURN ON QUADRATIC CONVERGENCE TO MACHINE UNDERFLOW
C
          IF (SM.EQ.0.0) RETURN
C
C---SET TRESHOLD
C
C---ON THE FIRST THREE SWEEPS
          IF (I.LT.4) THEN
              TRESH=0.2*SM/N**2
          ELSE

C---THEREAFTER
              TRESH=0.0
          ENDIF

C---PERFORM PQ ROTATIONS
C
          DO 22 IP=1,N-1
              DO 21 IQ=IP+1,N
                  G=100.0*ABS(A(IP,IQ))

C---AFTER 4 SWEEPS, SKIP ROTATION IF THE OFF-DIAGONAL ELEMENT IS SMALL
                  IF ( (I.GT.4) .AND. (ABS(D(IP))+G.EQ.ABS(D(IP)))
     -                .AND. (ABS(D(IQ))+G.EQ.ABS(D(IQ))) ) THEN
                      A(IP,IQ)=0.0

C---OTHERWISE ROTATE DEPENDENT ON THE TRESHOLD
                  ELSE IF (ABS(A(IP,IQ)).GT.TRESH) THEN
                      H=D(IQ)-D(IP)
                      IF (ABS(H).LT.EMIN) H = 0.0
                      IF (ABS(H)+G.EQ.ABS(H)) THEN

C   T=1/(2THETA)
                          T=A(IP,IQ)/H
                      ELSE

C   T=SIGN(THETA) / (ABS(THETA)+ SQRT(THETA**2 + 1))
                          THETA=0.5*H/A(IP,IQ)
                          T=1.0/(ABS(THETA)+SQRT(1.0+THETA**2))
                          IF (THETA.LT.0.0) T=-T
                      ENDIF
                      C=1.0/SQRT(1+T**2)
                      IF (ABS(C).LT.EMIN) C = 0.0
                      S=T*C
                      IF (ABS(S).LT.EMIN) S = 0.0
                      TAU=S/(1.0+C)
                      IF (ABS(TAU).LT.EMIN) TAU = 0.0
                      H=T*A(IP,IQ)
                      IF (ABS(H).LT.EMIN) H = 0.0
                      Z(IP)=Z(IP)-H
                      IF (ABS(Z(IP)).LT.EMIN) Z(IP) = 0.0
                      Z(IQ)=Z(IQ)+H
                      IF (ABS(Z(IQ)).LT.EMIN) Z(IQ) = 0.0
                      D(IP)=D(IP)-H
                      IF (ABS(D(IP)).LT.EMIN) D(IP) = 0.0
                      D(IQ)=D(IQ)+H
                      IF (ABS(D(IQ)).LT.EMIN) D(IQ) = 0.0
                      A(IP,IQ)=0.0

C   CASE OF ROTATIONS 1.LE.J.LT.P
                      DO 16 J=1,IP-1
                          G=A(J,IP)
                          IF (ABS(G).LT.EMIN) G = 0.0
                          H=A(J,IQ)
                          IF (ABS(H).LT.EMIN) H = 0.0
                          A(J,IP)=G-S*(H+G*TAU)
                          IF (ABS(A(J,IP)).LT.EMIN) A(J,IP) = 0.0
                          A(J,IQ)=H+S*(G-H*TAU)
                          IF (ABS(A(J,IQ)).LT.EMIN) A(J,IQ) = 0.0
16                    CONTINUE

C   CASE OF ROTATIONS P.LT.J.LT.Q
                      DO 17 J=IP+1,IQ-1
                          G=A(IP,J)
                          H=A(J,IQ)
                          A(IP,J)=G-S*(H+G*TAU)
                          IF (ABS(A(IP,J)).LT.EMIN) A(IP,J) = 0.0
                          A(J,IQ)=H+S*(G-H*TAU)
                          IF (ABS(A(J,IQ)).LT.EMIN) A(J,IQ) = 0.0
17                    CONTINUE

C   CASE OF ROTATIONS Q.LT.J.LE.N
                      DO 18 J=IQ+1,N
                          G=A(IP,J)
                          H=A(IQ,J)
                          A(IP,J)=G-S*(H+G*TAU)
                          IF (ABS(A(IP,J)).LT.EMIN) A(IP,J) = 0.0
                          A(IQ,J)=H+S*(G-H*TAU)
                          IF (ABS(A(IQ,J)).LT.EMIN) A(IQ,J) = 0.0
18                    CONTINUE
                      DO 19 J=1,N
                          G=V(J,IP)
                          H=V(J,IQ)
                          V(J,IP)=G-S*(H+G*TAU)
                          IF (ABS(V(J,IP)).LT.EMIN) V(J,IP) = 0.0
                          V(J,IQ)=H+S*(G-H*TAU)
                          IF (ABS(V(J,IQ)).LT.EMIN) V(J,IQ) = 0.0
19                    CONTINUE
                      NROT=NROT+1
                  ENDIF
21            CONTINUE
22        CONTINUE

C---END OF PQ ROTATIONS
C

C---UPDATE D WITH THE SUM OF TA(PQ) AND REINITIALIZE Z TO 0.0
C
          DO 23 IP=1,N
              B(IP)=B(IP)+Z(IP)
              D(IP)=B(IP)
              Z(IP)=0.0
23        CONTINUE
24    CONTINUE

C   END OF SWEEPS
C
C---ABNORMAL RETURN; 50 SWEEPS IN JACOBI SHOULD NEVER HAPPEN
C
      IFAIL=.TRUE.
      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE GETRMS  -  Calculate the RMS distance of the atoms from the
C                        best-fit plane
C
C----------------------------------------------------------------------+---

CHECK v.3.2-->
C      SUBROUTINE GETRMS(COORD,NPOINT,CENTRE,NORMAL,RMSDIF)
      SUBROUTINE GETRMS(COORD,NPOINT,NORMAL,RMSDIF)
CHECK v.3.2<--

      INTEGER       IPT, NPOINT
CHECK v.3.2-->
C      REAL          CENTRE(3), COORD(3,NPOINT), DOTPRD, NORMAL(3),
C     -              RMSDIF, VECTPQ(3)
      REAL          COORD(3,NPOINT), DOTPRD, NORMAL(3), RMSDIF,
     -              VECTPQ(3)
CHECK v.3.2<--

C---- Initialise values
      RMSDIF = 0.0

C---- Loop through all the data points
      DO 1500, IPT = 1, NPOINT

C----     Calculate distance of point from plane

C----     Form vector between point and an arbitrary point in the plane
C         (here the CofG is used as the arbitrary point)
          VECTPQ(1) = COORD(1,IPT)
          VECTPQ(2) = COORD(2,IPT)
          VECTPQ(3) = COORD(3,IPT)

C----     Calculate dot product PQ.n
          DOTPRD = VECTPQ(1) * NORMAL(1) + VECTPQ(2) * NORMAL(2)
     -        + VECTPQ(3) * NORMAL(3)

C----     Add distance squared to RMS-diff summation
          RMSDIF = RMSDIF + DOTPRD * DOTPRD
1500  CONTINUE

C---- Calculate and print RMS difference
      RMSDIF = SQRT(RMSDIF / REAL(NPOINT))

      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE WRIPLN  -  Write out the RMS distance of atoms from the
C                        best-fit plane, and transform so that can picture
C                        looking edge-on to plane
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRIPLN(COORD,NPOINT,EIGVEC,ORDVEC,SRESID,SEQNO,SNAME,
     -    CHAIN,RMSDIF,TCOORD,ATNAME,PDBFIL,PLRES)

CHECK v.3.2-->
      INTEGER       NPOINT
CHECK v.3.2<--
      CHARACTER*1   CHAIN
      CHARACTER*3   ATNAME(NPOINT), SNAME
      CHARACTER*5   SEQNO
      CHARACTER*78  PDBFIL
      INTEGER       ICOORD, IPT, IVEC
CHECK v.3.2-->
C      INTEGER       ILEN, NPOINT, ORDVEC(3), PLRES, SRESID
      INTEGER       ILEN, ORDVEC(3), PLRES, SRESID
CHECK v.3.2<--
      REAL          COORD(3,NPOINT), EIGVEC(3,3), RMSDIF,
     -              TCOORD(3,NPOINT)

C---- Write out the RMS distance
      ILEN = INDEX(PDBFIL,' ')
      IF (ILEN.GT.40) ILEN = 40
      WRITE(2,100) SRESID, CHAIN, SEQNO, SNAME, PLRES, RMSDIF,
     -    PDBFIL(1:ILEN)
 100  FORMAT(I6,A1,A5,A3,I2,F9.4,'  PLANE ',A)

C---- Tranform the coordinates to the new axes
      DO 400, IPT = 1, NPOINT
          DO 300, ICOORD = 1, 3

C----         Calculate eigen vector required (chosen s.t. longest axis
C             is along y, middle is in z-direction, and shortest in x)
              IF (ICOORD.EQ.1) IVEC = ORDVEC(3)
              IF (ICOORD.EQ.2) IVEC = ORDVEC(1)
              IF (ICOORD.EQ.3) IVEC = ORDVEC(2)

C----         Calculate transformed coordinates
              TCOORD(ICOORD,IPT)
     -            = COORD(1,IPT) * EIGVEC(1,IVEC)
     -            + COORD(2,IPT) * EIGVEC(2,IVEC)
     -            + COORD(3,IPT) * EIGVEC(3,IVEC)
 300      CONTINUE
 400  CONTINUE

C---- Write the transformed coordinates out to file
      DO 600, IPT = 1, NPOINT
          WRITE(3,420) PLRES, ATNAME(IPT), SNAME, CHAIN, SEQNO,
     -        (TCOORD(ICOORD,IPT), ICOORD = 1, 3), RMSDIF
 420      FORMAT('ATOM     ',I2,2X,A3,1X,A3,1X,A1,A5,3X,3F8.3,6X,F8.4)
 600  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.2.2<--
C**************************************************************************
C
C  BLOCK DATA
C
C----------------------------------------------------------------------+---

      BLOCK DATA

      INCLUDE 'anglen.inc'

CHECK v.2.2-->
C---- Definitions of planar groups, given by residue name and the
C     atoms defining the group
      DATA ANAME  / 'ARG', 'NE ', 'CZ ', 'NH1', 'NH2', '   ',
     -                     '   ', '   ', '   ', '   ', '   ',
     -              'ASN', 'CB ', 'CG ', 'OD1', 'ND2', '   ',
     -                     '   ', '   ', '   ', '   ', '   ',
     -              'ASP', 'CB ', 'CG ', 'OD1', 'OD2', '   ',
     -                     '   ', '   ', '   ', '   ', '   ',
     -              'GLN', 'CG ', 'CD ', 'OE1', 'NE2', '   ',
     -                     '   ', '   ', '   ', '   ', '   ',
     -              'GLU', 'CG ', 'CD ', 'OE1', 'OE2', '   ',
     -                     '   ', '   ', '   ', '   ', '   ',
     -              'HIS', 'CB ', 'CG ', 'ND1', 'CD2', 'CE1',
     -                     'NE2', '   ', '   ', '   ', '   ',
     -              'PHE', 'CB ', 'CG ', 'CD1', 'CD2', 'CE1',
     -                     'CE2', 'CZ ', '   ', '   ', '   ',
C     -              'PRO', 'N  ', 'CA ', 'CB ', 'CG ', 'CD ',
C     -                     '   ', '   ', '   ', '   ', '   ',
     -              'TRP', 'CB ', 'CG ', 'CD1', 'NE1', 'CE2',
     -                     'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2',
     -              'TYR', 'CB ', 'CG ', 'CD1', 'CD2', 'CE1',
     -                     'CE2', 'CZ ', 'OH ', '   ', '   ' /
CHECK v.2.2<--

C---- Definitions of the 20 different types of main-chain bond-lengths
C     and the 11 different types of main-chain bond angles to
C     be calculated. Each triplet represents a bond angle unless the
C     3rd position is blank, in which case it is a bond length
      DATA BONDAT / 
     -              'CA ', 'C  ', 'N  ', 'CA ', 'C  ', 'N  ',
     -              'CA ', 'C  ', 'N  ',
     -              'O  ', 'C  ', 'N  ', 'O  ', 'C  ', 'N  ', 
     -              'C  ', 'N  ', '   ', 'C  ', 'N  ', '   ',
     -              'C  ', 'O  ', '   ',
     -              'CA ', 'C  ', '   ', 'CA ', 'C  ', '   ',
     -              'CA ', 'CB ', '   ', 'CA ', 'CB ', '   ',
     -              'CA ', 'CB ', '   ',
     -              'N  ', 'CA ', '   ', 'N  ', 'CA ', '   ',
     -              'N  ', 'CA ', '   ',
     -              'C  ', 'N  ', 'CA ', 'C  ', 'N  ', 'CA ',
     -              'C  ', 'N  ', 'CA ',
     -              'CA ', 'C  ', 'O  ', 'CA ', 'C  ', 'O  ',
     -              'CB ', 'CA ', 'C  ', 'CB ', 'CA ', 'C  ',
     -              'CB ', 'CA ', 'C  ',
     -              'N  ', 'CA ', 'C  ', 'N  ', 'CA ', 'C  ',
     -              'N  ', 'CA ', 'C  ',
     -              'N  ', 'CA ', 'CB ', 'N  ', 'CA ', 'CB ',
     -              'N  ', 'CA ', 'CB ', 'N  ', 'CA ', 'CB ' /

C---- Engh & Huber names corresponding to each of the triplet in the
C     BONDAT array above
      DATA ENGNAM  / 
     -              'CH1E-C-NH1   ', 'CH2G*-C-NH1  ', 'CH1E-C-N     ',
     -              'O-C-NH1      ', 'O-C-N        ',
     -              'C-NH1        ', 'C-N          ',
     -              'C-O          ',
     -              'CH1E-C       ', 'CH2G*-C      ',
     -              'CH1E-CH3E    ', 'CH1E-CH1E    ', 'CH1E-CH2E    ',
     -              'NH1-CH1E     ', 'NH1-CH2G*    ', 'N-CH1E       ',
     -              'C-NH1-CH1E   ', 'C-NH1-CH2G*  ', 'C-N-CH1E     ',
     -              'CH1E-C-O     ', 'CH2G*-C-O    ',
     -              'CH3E-CH1E-C  ', 'CH1E-CH1E-C  ', 'CH2E-CH1E-C  ',
     -              'NH1-CH1E-C   ', 'NH1-CH2G*-C  ', 'N-CH1E-C     ',
     -              'NH1-CH1E-CH3E', 'NH1-CH1E-CH1E', 'N-CH1E-CH2E  ',
     -              'NH1-CH1E-CH2E'
     -            /

C---- Mean values for each of the 31 bond lengths and angles defined
C     above.
      DATA ENGMEA / 
     -              116.2, 116.4, 116.9,
     -              123.0, 122.0,
     -              1.329, 1.341,
     -              1.231,
     -              1.525, 1.516,
     -              1.521, 1.540, 1.530,
     -              1.458, 1.451, 1.466,
     -              121.7, 120.6, 122.6,
     -              120.8, 120.8,
     -              110.5, 109.1, 110.1,
     -              111.2, 112.5, 111.8,
     -              110.4, 111.5, 103.0, 110.5
     -            /

C---- Standard deviations of the ENGMEA values above.
      DATA ENGSTD /
     -              2.0, 2.1, 1.5,
     -              1.6, 1.4,
     -              0.014, 0.016,
     -              0.020,
     -              0.021, 0.018,
     -              0.033, 0.027, 0.020,
     -              0.019, 0.016, 0.015,
     -              1.8, 1.7, 5.0,
     -              1.7, 2.1,
     -              1.5, 2.2, 1.9,
     -              2.8, 2.9, 2.5,
     -              1.5, 1.7, 1.1, 1.7
     -            /

C---- Exception categories for each of the above set of triplets.
C     The exceptions themselves (numbered 1-10) are defined in the
C     EXTYPE table below.
      DATA EXCEPT / 
     -              8, 4, 2,
     -              3, 2,
     -              3, 2,
     -              1,
     -              5, 4,
     -              6, 7, 9,
     -              8, 4, 2,
     -              8, 4, 2,
     -              5, 4,
     -              6, 7, 9,
     -              8, 4, 2,
     -              6, 7, 2,10 /

C---- The 10 exception definitions, showing which residues are
C     included/excluded when computing a particular bond length/angle.
C     The EXCEPT table above defines which of these categories applies
C     to each of the 31 main-chain bond length/angle types.
      DATA EXTYPE / 'ALL', '   ', '   ', '   ', '   ', '   ',
     -              '   ', 'PRO', '   ', '   ', '   ', '   ',
     -              'NOT', 'PRO', '   ', '   ', '   ', '   ',
     -              '   ', 'GLY', '   ', '   ', '   ', '   ',
     -              'NOT', 'GLY', '   ', '   ', '   ', '   ',
     -              '   ', 'ALA', '   ', '   ', '   ', '   ',
     -              '   ', 'ILE', 'THR', 'VAL', '   ', '   ',
     -              'NOT', 'GLY', 'PRO', '   ', '   ', '   ',
     -              'NOT', 'ILE', 'THR', 'VAL', 'ALA', '   ',
     -              'NOT', 'ILE', 'THR', 'VAL', 'ALA', 'PRO' /

C---- Descriptions corresponding to each of the exception categories
C     defined in the EXTYPE table above.
      DATA EXDESC / '                ',
     -              '(Pro)           ',
     -              '(except Pro)    ',
     -              '(Gly)           ',
     -              '(except Gly)    ',
     -              '(Ala)           ',
     -              '(Ile,Thr,Val)   ',
     -              '(except Gly,Pro)',
     -              '(the rest)      ',
     -              '(the rest)      ' /

C---- Definitions of the main-chain atoms.
      DATA  MATOM / 'C  ', 'N  ', 'O  ', 'CA ', 'CB ' /

C---- Definitions of the 20 amino-acid 3-letter codes, with space for
C     any unusual ones encountered in the structure.
      DATA RTYPE  / 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
     -              'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
     -              'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 80*'   ' /
      END

C--------------------------------------------------------------------------
