C-----------------------------------------------------------------------------
C
C SECSTR.F - Amended version of sstruc.f for use with structure-checking suite
C            of programs, PROCHECK
C
C            Amendments by R A Laskowski, April 1992, identified by
C            CHECK v.m.n--> and CHECK v.m.n<-- comment lines, where m.n is
C            the version number.
C
C v.2.0    - Rewrite to allow entry of any filename for the coordinates file
C            rather than just the 4-letter Brookhaven code
C                                               Roman Laskowski (Jul/Aug 1992)
C v.2.1    - Bug-fix so that input file is the .new file, rather than the .pdb
C            file
C                                               Roman Laskowski ( 9 Dec 1992)
C v.2.1.3  - Addition of asterisks to error prints so that can be plucked
C            out of the log files and displayed by the script file.
C                                               Roman Laskowski (26 Mar 1993)
C            Bug-fix to prevent crashing in Subroutine Mktb in with divide
C            by zero (don't know why). Fix merely checks for zero and skips
C            the appropriate part of the code
C                                               Roman Laskowski (14 Apr 1993)
C v.2.1.4  - Amendment to change name of residue information files from .res
C            to .rin so as not to clash with the .rin files produced by SHELX.
C            Amendment to FORMAT statements to change x to 1x throughout as
C            the IBM RS/6000 compiler objects to these.
C            Crude change to limit number of 'unknown amino acid' error
C            messages displayed (particularly irritating for lots of waters).
C                                               Roman Laskowski (12 May 1993)
C
C v.3.0   - Tiny change so that VAX .COM routines don't show error if no
C           asterisks found in log file.
C           Amendment so that Ooi numbers are written out to the .rin file.
C           Amendment so that Chi-3, 4, and 5 dihedrals are written out to
C           .rin file.
C                                            Roman Laskowski (14-26 Nov 1993)
C
C v.3.0.1 - Changed order of DATA statements in some routines so that they
C           follow the variable declarations (was causing problems when
C           compiling with the f2c compiler).
C           Commented out unreferenced labels (producing compiler warnings
C           when compiling on Convex).
C                                            Roman Laskowski (29/30 Mar 1994)
C
C v.3.1   - Amendments for PROCHECK-NMR to deal with each MODEL in an
C           ensemble separately, writing each to a separate .rin file.
C                                         Roman Laskowski (18/19 Apr 1994)
C
C v.3.2   - Calculation of main-chain and side-chain B-value st. devs.
C                                            Roman Laskowski (30 Apr 1994)
C           Minor amendments to various statements to make them
C           acceptable to f2c, and to deal with various uninitialised
C           variables. All tabs replaced by spaces. SAVE statements
C           for variables moved prior to declaration of those variables
C           in the routine. (Amendments supplied by Dave Love at
C           Daresbury).
C                                 David Love/Roman Laskowski (13 Oct 1994)
C           Amendment for case where one one model in the PDB file and
C           hence no MODEL and ENDMDL records
C                                            Roman Laskowski (25 Oct 1994)
C           Reinstatement of facility to process a list of PDB files in
C           a single run. (For use with PROCOMP when comparing a set of
C           similar structures).
C                                             Roman Laskowski (7 Nov 1994)
C
C v.3.3.1 - Amendment to write .rin files for CA-only files.
C           Increased number of chain-breaks allowed, and checked for
C           overflow. (Can get a large number of apparent chain-breaks
C           in structures where one chain is normal, and the other has
C           C-alpha coordinates only).
C                                             Roman Laskowski (5 Jun 1995)
C v.3.4   - Bug-fix to run PROCHECK correctly on the 1st model only of a
C           multiple-model file.
C           Addition of model-number to output NMR .rin files.
C                                            Roman Laskowski (20 Mar 1996)
C           Bug-fixes suggested by Dave Love.
C                                            Roman Laskowski (26 Apr 1996)
C
C v.3.4.4 - Increase in MAXRES to 8000 residues.
C                                            Roman Laskowski (15 Oct 1996)
C
C v.3.5   - Amendment to allow zero-occupancy atoms to be included in the
C           analysis
C                                            Roman Laskowski ( 3 Nov 1996)
C           Amendment so that residue types not hard-coded into the
C           code are not completely ignored but are treated as UNK residues.
C           Replacement of initialisation DATA statements for the large
C           by initialising loops, as for clean.f
C                                            Roman Laskowski (20 Nov 1996)
C           Amendment to cope with the first model being empty.
C                                            Roman Laskowski ( 7 May 1997)
C           Reinstatement of loop variables commented out in v.3.4.4 which
C           resulted in compiler errors under DEC alpha (as pointed out by
C           Mitch Miller).
C                                            Roman Laskowski (30 May 1997)
C v.3.5.1 - Amendment to deal with crap PDB files which don't properly
C           distinguish different models using MODEL and ENDMDL records
C           (or even TER records!)
C                                            Roman Laskowski ( 3 Sep 1998)
C v.3.5.2 - Bug-fix. Uninitialised variable clslen in routine RDBRK causing
C           core-dumps under linux. Fixed by setting to zero at start of
C           the routine. (Error found by Kai Zhi Yue).
C                                            Roman Laskowski (31 Mar 2000)
C v.3.5.3 - Bug-fix. Uninitialised variables soulen and cmplen in routine
C           RDBRK causing core-dumps under linux when PDB file has no
C           header records. Fixed by setting to zero at start of the
C           routine. (Error found by Andrew Martin).
C                                            Roman Laskowski (23 Nov 2000)
C v.3.6   - Amendment to ignore waters when picking up residues in the case
C           where the PDB file has no SEQRES records.
C                                            Roman Laskowski ( 5 Jul 2001)
C v.3.6.1 - Amendment to residue-number writing routine to .rin file where 
C           number of residues exceeds 99999 (eg as in 1o18).
C                                            Roman Laskowski (30 May 2003)
C v.3.6.2 - Amendment for PDBsum generation - write-out of disulphide
C           bonds into a separate file
C                                            Roman Laskowski (22 Jan 2004)
C
C Files:
C
C  2 mplot.in       - List of NMR ensemble .rin files, for use by PROCHECK-NMR
C                     program mplot.f
C  8 <filename>.rin - Output file holding all the residue-by-residue
C                     information required by program PPLOT.F
C 11 <filename>.new - Input .new file, being the cleaned-up version of the
C                     input PDB file.
C
C-----------------------------------------------------------------------------


      PROGRAM SECSTR

CHECK v.1.0-->
C      program sstruc
CHECK v.1.0<--

c     copyright by David Keith Smith, 1989
c****************************************************************************
c     filename:      sstruc.f
c     last modified: 11th September 1991
c****************************************************************************
c     global variables
c
      INCLUDE 'sstruc.par'

      character aacod1*(abetsz), aacod3*(3*abetsz), atname*(4*mnchna-4),
     +          ssbond(nstruc,maxres)*1, sssumm(maxres)*1,
     +          summpl(maxres)*1, seqbcd(maxres)*6, 
     +          autsum(maxres)*1, ssord*21, brksym(maxres)*1,
     +          headln*132, altcod(maxres)*1, aastd(maxres)*1,
     +          atpres(3,maxres)*1, thmprs(maxres)*1, hetcod*63,
     +          outsty*1, dsdef(maxres)*1, dudcod*105, modcod*(nmod*3),
     +          stdmod*(nmod+1), nstdmd*(nmod+1), modpl*(nmod),
     +          batnam*(4*mnchna+4)
      integer   seqcod(maxres), seqlen, hbond(maxbnd,maxres),
     +          chnsz(maxchn), chnct, bridpt(2,maxres), hedlen,
     +          dspart(maxres), chnend(maxchn), oois(2,maxres)
      real      mcaaat(ncoord,mnchna,maxres), hbonde(maxbnd,maxres),
     +          scaaat(ncoord,sdchna,maxres), mcangs(nmnang,maxres),
     +          scangs(nsdang,maxres), cadsts(11,maxres), 
     +          dsdsts(3,maxres)
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres), caonly, 
     +          readok, astruc

CHECK v.3.0.1--> (Order of statements changed)
c
c
c     local variables
c
      integer   finish, opnera, opnerb, opnerc, fiend, ok, endnfi
      character indirn*1
      parameter (finish = 2, opnera = 1, opnerb = 3, fiend = 4, ok = 0, 
     +           opnerc = 5, endnfi = 6,
     +           indirn='@') 
      logical   finuse
CHECK v.3.2-->
C      integer   fistat, namlen, iocode, i, iostat,
CHECK v.2.1-->
C     +          corepl, corend, logpl, dirpl, nodepl, outlen
C     +          corepl, corend, outlen
CHECK v.2.1<--
      integer   fistat, namlen, iocode, i, iostat
CHECK v.3.2<--
CHECK v.3.2-->
C      character*(fnamln) finam, brknam, outnam, outdir
      character*(fnamln) brknam, finam
CHECK v.3.2<--
CHECK v.2.1-->
C      integer   ifail, outln, lastsl, brklen 
CHECK v.3.2-->
C      integer   outln
CHECK v.3.2<--
CHECK v.2.1<--

CHECK v.3.5-->
      CHARACTER*3   SAVNAM(MAXRES)
      INTEGER       J, K
CHECK v.3.5<--

CHECK v.1.0-->
      CHARACTER*4   BRCODE
      REAL          BVALUE(MAXRES)
CHECK v.1.0<--

CHECK v.2.0-->
CHECK v.3.6.2-->
C      CHARACTER*78  FILRES, PDBFIL
      CHARACTER*78  FILDIS, FILRES, PDBFIL
CHECK v.3.6.2<--
      INTEGER       IEND, ILEN, ISTART, NMCATM(MAXRES), NSCATM(MAXRES)
      LOGICAL       IERROR
      REAL          MCBVAL(MAXRES), SCBVAL(MAXRES)
CHECK v.2.0<--
CHECK v.3.0.1<--

CHECK v.3.1-->
      CHARACTER*3   CMODEL
CHECK v.3.4-->
C      INTEGER       IMODEL, NMODEL
      INTEGER       IMODEL, MODLNO, NMODEL
CHECK v.3.4<--
CHECK v.3.2-->
C      LOGICAL       ENDFIL, FIRST, NMR
      LOGICAL       ENDFIL, ENSEMB, FIRST, NMR
CHECK v.3.2<--
CHECK v.3.1<--

CHECK v.3.2-->
      REAL          MCBSTD(MAXRES), SCBSTD(MAXRES)
CHECK v.3.2<--

CHECK v.3.5-->
      INTEGER       AAPLN, SQRSCT
      PARAMETER    (AAPLN = 13)
      CHARACTER*1   GAPIN(AAPLN)
      CHARACTER*3   AAIN(AAPLN)
      CHARACTER inseqr*(3*MAXTYP)
CHECK v.3.5<--

      data aacod1(1:29)  /'ABCDEFGHIXKLMNXPQRSTXVWXYZX--'/
      data aacod3( 1:39) /'ALAASXCYSASPGLUPHEGLYHISILEXXXLYSLEUMET'/,
     +     aacod3(40:78) /'ASNXXXPROGLNARGSERTHRXXXVALTRPXXXTYRGLX'/,
     +     aacod3(79:87) /'UNKPCAINI'/
      data hetcod(1:39)  /'AIBPHLSECALMMPRFRDPYRLYMGLMPPHPGLOLETPH'/,
     +     hetcod(40:63) /'ABANLEB2VB2IB1FBNOB2AB2F'/
      data modcod(1:39)  /'ACEANIBZOCBZCLTFORNH2PHORHATFATOS NHMYR'/,
     +     modcod(40:45) /'BOCMSU'/
      data dudcod( 1:39) /'  A  C  G  T  U1MA5MCOMC1MG2MGM2G7MGOMG'/,
     +     dudcod(40:78) /' YG  I +UH2U5MUPSU +C +AGALAGLGCUASGMAN'/
     +     dudcod(79:105) /'GLCCEGG4SAGSNAGGLSNGSNAMEXC'/
      data stdmod /'ABCDEFGHIJKLMNZY'/
      data nstdmd /'abcdefghijklmnzy'/
      data modpl  /'+0++++0++++0+++'/
      data ssord /'HEheBbGgIiJKLMNOPTtS '/
      data atname /' N   CA  C   O  '/
      data batnam /' N   CA  B       O1  O2 '/
CHECK v.3.5-->
C      data atomin /allatm * .false./, mcaaat /allcrd * 0/
CHECK v.3.5<--
      data fistat/0/, finuse/.false./
CHECK v.3.2-->
      data bvalue /maxres*0.0/
CHECK v.3.2<--


CHECK v.3.5-->
C---- Initialise all the variables previously initialised in DATA
C     statements
      DO 50, k = 1, maxres
          DO 40, j = 1, mnchna
              atomin(j,k) = .FALSE.
              DO 30, i = 1, ncoord
                  mcaaat(i,j,k) = 0.0
 30           CONTINUE
 40       CONTINUE
 50   CONTINUE
      inseqr = ' '
      sqrsct = 0
CHECK v.3.5<--

c
c     input file control loop
c
CHECK v.3.1-->
CHECK v.3.4-->
      IMODEL = 0
CHECK v.3.4<--
      NMODEL = 0
CHECK v.3.1<--

CHECK v.3.5-->
      KEEPAL = .FALSE.
CHECK v.3.5<--

CHECK v.1.0-->
      WRITE(*,*)
CHECK v.1.0<--

      write(*,*)'Secondary structure calculation program - ',
     +          'copyright by David Keith Smith, 1989'
      write(*,*)

CHECK v.1.0-->
      WRITE(*,*) 'Amended by R A Laskowski, 1992'
      WRITE(*,*)

CHECK v.2.0-->
COUT---- Accept 4-letter Brookhaven code
COUT      WRITE (*, 1010)
COUT 1010 FORMAT ('Enter 4-letter Brookhaven code of structure:')
COUT      READ (*,1020) BRCODE
COUT 1020 FORMAT(A)
COUT      BRKNAM = 'p' // BRCODE // '.new'
COUT      WRITE (*,1022) BRKNAM(1:NAMLEN+5)
COUT 1022 FORMAT(//,
COUT     -    'Processing file: ',A,/)

COUT---- Form the names of the output files
COUT      OUTNAM = 'p' // BRCODE // '.sst'
COUT      OUTLEN = 9

COUT---- Open the file that will hold the dihedral angles
COUT      OPEN (UNIT=8,FILE='p'//BRCODE//'.res',STATUS='UNKNOWN',ERR=900,
COUTVAX
COUT     -    CARRIAGECONTROL='LIST'
COUTVAX
COUT     -    )
CHECK v.1.0<--

CHECK v.2.0-->
C---- Accept name of original .pdb file holding the structure
      PRINT*, 'Enter filename containing coordinates of structure'
CHECK v.3.2-->
      PRINT*, '  (for file containing ensemble of NMR structures en',
     -    'ter @filename;'
      PRINT*, '   for set of separate PDB files to be processed, en',
     -    'ter %filelist,'
      PRINT*, '   where filelist contains a list of PDB files to be',
     -    ' cleaned up)'
CHECK v.3.2<--
      READ(*,10) PDBFIL
10    FORMAT(A)

CHECK v.3.1-->
C---- If the very first character of the file-name is an @, then the
C     PDB file is an NMR file containing several models of the structure
CHECK v.3.2-->
C      IF (PDBFIL(1:1).EQ.'@') THEN
C          NMR = .TRUE.
C          PDBFIL = PDBFIL(2:)
C      ELSE
C          NMR = .FALSE.
C      ENDIF
C      FIRST = .TRUE.
      IF (PDBFIL(1:1).EQ.'@') THEN
          NMR = .TRUE.
          PDBFIL = PDBFIL(2:)
          FIRST = .TRUE.
          ENSEMB = .FALSE.
      ELSE IF (PDBFIL(1:1).EQ.'%') THEN
          ENSEMB = .TRUE.
          FIRST = .TRUE.
          NMR = .FALSE.
      ELSE
          ENSEMB = .FALSE.
          FIRST = .TRUE.
          NMR = .FALSE.
      ENDIF
CHECK v.3.2<--

CHECK v.3.5-->
C---- If the next character of the PDB file is a +, then want to include
C     all the zero-occupancy atoms in the analysis
      IF (PDBFIL(1:1).EQ.'+') THEN
          PDBFIL = PDBFIL(2:)
          KEEPAL = .TRUE.
      ENDIF
CHECK v.3.5<--

C---- If this an NMR file, then open mplot.in file which will contain
C     the name of all the individual .new file output (one per model)
CHECK v.3.2-->
C      IF (NMR) THEN
      IF (ENSEMB .OR. NMR) THEN
CHECK v.3.2<--
          OPEN(UNIT=MPLTFI, FILE='mplot.in', STATUS='UNKNOWN',
     -        FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -    CARRIAGECONTROL='LIST',RECL=132,
     -    ERR=902)
      ENDIF
CHECK v.3.1<--

CHECK v.3.2-->
C---- If processing a set of PDB files listed in the %filelist specified
C     by the user, open the file-list
      fistat = ok
      IF (ENSEMB) THEN
          finam = PDBFIL
          finuse = .true.
          namlen = INDEX(finam,space)
          do 180 i = 1, namlen - 1
              finam(i:i) = finam(i+1:i+1)
  180     continue
          finam(namlen:namlen) = space
          open (namfi, file=finam, status='old',
CVAX     -        READONLY,
     -        iostat=iocode) 
          if (iocode .ne. 0) fistat = opnera
          if (fistat.eq.opnera) GO TO 902
      ENDIF
CHECK v.3.2<--

CHECK v.3.2-->
C---- If processing a list of structures, get the name of the next one
C     from the file list
 200  CONTINUE
      IF (ENSEMB .AND. fistat .eq. ok) then
          read (namfi, 10, iostat=iocode) brknam
          if (iocode .lt. 0) then
             close (namfi) 
             finuse = .false.
             fistat = endnfi
          else
             namlen = index(brknam,space) - 1
             write (*, 1030) brknam(1:namlen)
 1030        format (' Working on Brookhaven file ',A,' . . . ')
          end if
          PDBFIL = brknam
          NMR = .FALSE.
      ENDIF
CHECK v.3.2<--

C---- Peel off directory path and extension
      CALL GETNAM(PDBFIL,ISTART,IEND,IERROR)
      IF (IERROR) GO TO 999

C---- Form names of other files that will be required in default directory
      ILEN = IEND - ISTART + 1
      BRCODE = ' '
      BRKNAM = PDBFIL
      NAMLEN = 78
CHECK v.3.2-->
C      OUTLEN = 78
CHECK v.3.2<--
CHECK v.2.1.4-->
C      FILRES = PDBFIL(ISTART:IEND) // '.res'
      FILRES = PDBFIL(ISTART:IEND) // '.rin'
CHECK v.2.1.4<--
CHECK v.2.1-->
      BRKNAM = PDBFIL(ISTART:IEND) // '.new'
CHECK v.2.1<--
CHECK v.3.6.2-->
      FILDIS = PDBFIL(ISTART:IEND) // '.dis'
CHECK v.3.6.2<--


CHECK v.3.1-->
CC---- Open the file that will hold the dihedral angles
C      OPEN (UNIT=8,FILE=FILRES,STATUS='UNKNOWN',ERR=900
CCVAX     -    ,CARRIAGECONTROL='LIST'
C     -    )
CHECK v.3.1<--

C---- Initialise B-values
      DO 20, I = 1, MAXRES
          MCBVAL(I) = 0.0
          SCBVAL(I) = 0.0
CHECK v.3.2-->
          MCBSTD(I) = 0.0
          SCBSTD(I) = 0.0
CHECK v.3.2<--
          NMCATM(I) = 0
          NSCATM(I) = 0
 20   CONTINUE
CHECK v.2.0<--

c     write(*,*)'Do you require Standard or Full output (<S> or F)'
c     read(*,'(A1)')outsty
c     if (outsty .eq. 'F' .or. outsty .eq. 'f') then
c        Outsty = 'F'
c     else
c        outsty = 'S'
c     end if
      Outsty = 'F'

CHECK v.1.0-->
C      call getstr('Enter output directory', 
C     +            '/home/bsm/naylor/data/', outdir, outln, ifail)
C      if (ifail .ne. 0)  then
C         write(6,*) 'Error in GETSTR: Ifail = ', ifail
C         stop
C      end if
C  100 continue
CHECK v.1.0<--

CHECK v.3.2-->
C      fistat = ok
CHECK v.3.2<--
c
c     get the next brk file or file of names from the user
c     finish is set when they've had enough
c

CHECK v.1.0-->
C      write(*,*)'Enter the name of the next file to be processed.'
C      write(*,*)'  Just a file name for a Brookhaven format file'
C      write(*,'(a,a1,a,a)')'   ',Indirn,'file name for a file holding ',
C     +                     'names of Brookhaven format files.'
C      write(*,*)'  If you''ve had enough a blank line let''s',
C     +          ' you escape!'
C      read(*,'(A)')finam
C      namlen = index(finam,space) - 1
C      if (namlen .le. 0) then
C         fistat = finish
C      else
C         if (finam(1:1) .eq. indirn) then
C            finuse = .true.
C            do 180 I = 1,namlen-1
C            finam(i:i) = finam(I+1:I+1)
C  180       continue
C            finam(namlen:namlen) = space
C            open(unit=namfi,file=finam,status='old',iostat=iocode)
C            if (iocode .ne. 0) fistat = opnera
C         end if
C      end if
CHECK v.1.0<--

c
c     get the brookhaven file name into the appropriate place
c     either from the file of names or the input name
c

CHECK v.1.0-->
C  200 continue
C      if (fistat .eq. ok) then
C         if (finuse) then
C            read(namfi,'(A)',iostat=iocode) brknam
C            if (iocode .lt. 0) then
C               close(unit=namfi)
C               finuse = .false.
C               fistat = endnfi
C            else
C               namlen = index(brknam,space) - 1
C               write(*,*)brknam(1:namlen)
C            end if
C         else
C            brknam = finam
C         end if
C      end if
CHECK v.1.0<--

c
c     open the brookhaven file if everything ok
c
      If (fistat .eq. ok) then
CHECK v.1.0-->
C         open (Unit=brkfi, file=brknam, status='old', iostat=iocode)
CHECK v.1.0<--
         open (Unit=brkfi, file=brknam, status='old', iostat=iocode
CVAX     -         , READONLY
     -        )
         if (iocode .ne. 0) fistat = opnerb
      end if

CHECK v.1.0-->
C      if (fistat .eq. ok) then
C         nodepl = index(brknam,'::')
C         if (nodepl .ne. 0) then
C            logpl = nodepl + 1 + index(brknam(nodepl+2:),':')
C         else
C            logpl = index(brknam,':')
C         end if
C         dirpl = index(brknam,']')
C         corepl = max(logpl,dirpl) + 1
C         corepl = lastsl(brknam) + 1
C         corend = corepl - 1 + index(brknam(corepl:),'.') - 1
C         if (corend .eq. corepl - 2) corend = namlen
C         if (corend .lt. corepl) then
C            write (6,*) brknam(1:namlen)
C            fistat = opnerc
C         end if
C         brklen = corend - corepl + 1
C         outnam = outdir(1:outln)//brknam(corepl:corend)// '.sst'
C         outlen = outln + brklen + 4
C      end if
CHECK v.1.0<--

CHECK v.3.1-->
C---- If this is an NMR structure, then will loop back to here until
C     all models have been individually processed
 100  CONTINUE

CHECK v.3.2-->
C      IF (NMR) THEN
      if (fistat.eq.ok) then
CHECK v.3.2<--

C----     Loop until first ATOM record encountered, or a MODEL record
C         is found (indicating an NMR structure) or end-of-file encountered
CHECK v.3.2-->
C          CALL GETMOD(BRKFI,IMODEL,MXMODL,ENDFIL,IERROR)
CHECK v.3.4-->
C          CALL GETMOD(BRKFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR)
CHECK v.3.4.3-->
C          CALL GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,IERROR,NMR)
CHECK v.3.5-->
C          CALL GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,IERROR)
          CALL GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,AAPLN,AAIN,
     -        GAPIN,INSEQR,MAXTYP,SQRSCT,IERROR)
CHECK v.3.5<--
CHECK v.3.4.3<--
CHECK v.3.4<--
CHECK v.3.2<--

C----     Check whether end of file or data error encountered
          IF (ENDFIL) GO TO 500
          IF (IERROR) GO TO 999

C----     Set up filename for .new output file for this model
CHECK v.3.2-->
          IF (NMR) THEN
CHECK v.3.2<--
CHECK v.3.4-->
C              WRITE(CMODEL,120) IMODEL
              WRITE(CMODEL,120) MODLNO
CHECK v.3.4<--
 120          FORMAT(I3)
              IF (CMODEL(1:1).EQ.' ') CMODEL(1:1) = '0'
              IF (CMODEL(2:2).EQ.' ') CMODEL(2:2) = '0'
              FILRES = PDBFIL(ISTART:IEND) // '_' // CMODEL // '.rin'
              NMODEL = NMODEL + 1
          ENDIF
CHECK v.3.2<--

CHECK v.3.2-->
C----     Write this name out to the mplot.in file
C          WRITE(MPLTFI,140) FILRES(1:ILEN + 8)
C 140      FORMAT(A)
C----     Write this name out to the mplot.in file
          IF (ENSEMB .OR. NMR) THEN
              WRITE(MPLTFI,140) FILRES(1:ILEN + 8)
 140          FORMAT(A)
          ENDIF
CHECK v.3.2<--

C----     Close previous .new file, if it exists
          IF (.NOT.FIRST) CLOSE(8)
          FIRST = .FALSE.
CHECK v.3.2-->
          IF (NMR) THEN
CHECK v.3.2<--
              PRINT*
CHECK v.3.4-->
C              PRINT*, '   Processing NMR model', IMODEL
              PRINT*, '   Processing NMR model', MODLNO
CHECK v.3.4<--
CHECK v.3.2-->
          ENDIF
CHECK v.3.2<--
      ENDIF
CHECK v.3.1<--

CHECK v.3.1-->
C---- Open the file that will hold the dihedral angles
      OPEN (UNIT=8,FILE=FILRES,STATUS='UNKNOWN',ERR=900
CVAX     -    ,CARRIAGECONTROL='LIST'
     -    )
CHECK v.3.1<--

c
c     if survived this far have a file to process otherwise
c     provide error messages
c
      if (fistat .eq. ok) then
c       **Process file
         open(unit = errfi, status = 'SCRATCH')
         readok = .true.
         caonly = .true.
         astruc = .false.
         call rdbrk(mcaaat, atname, atomin, aacod3, hetcod, seqcod,
     +              seqbcd, scaaat, scatin, altcod, aastd,  atpres,
     +              thmprs, caonly, astruc, headln, hedlen, readok,
     +              dudcod, modcod, stdmod, nstdmd, modpl,  seqlen,

CHECK v.1.0-->
C     +              batnam)
CHECK v.2.0-->
C     +              batnam, BVALUE)
CHECK v.3.1-->
C     +              batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM)
CHECK v.3.2-->
C     +              batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM,NMR)
     +              batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM,NMR,
CHECK v.3.5-->
C     +              MCBSTD,SCBSTD)
     +              MCBSTD,SCBSTD,SAVNAM,inseqr,sqrsct)
CHECK v.3.5<--
CHECK v.3.2<--
CHECK v.3.1<--
CHECK v.2.0<--
CHECK v.1.0<--

CHECK v.3.3.1-->
C         if (readok .and. .not. caonly) then
         if (readok) then
CHECK v.3.3.1<--

CHECK v.1.0-->
C            open(unit=outfi,file=outnam,status='unknown',iostat=iocode)
CHECK v.1.0<--

CHECK v.3.2-->
C            if (iocode .ne. 0) fistat = opnerc
CHECK v.3.2<--
c           ** look for breaks across the peptide bonds
            call chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, chnct,
     +                  chnend, caonly, seqlen)
c           ** determine the calpha distances i+2 to i+11
            if (outsty .eq. 'F')
CHECK v.3.2-->
C     +         call mkcadt(mcaaat, atomin, cadsts, seqbcd, chnsz, chnct,
C     +                     seqlen)
     +         call mkcadt(mcaaat, atomin, cadsts, seqbcd,
     +                     seqlen)
CHECK v.3.2<--
c           ** create hydrogen atoms on the main chain nitrogens
CHECK v.3.2-->
C            call creatH(mcaaat, atomin, chnsz, chnct, seqlen)
            call creatH(mcaaat, atomin, chnsz, chnct)
CHECK v.3.2<--
c           ** determine all hydrogen bonds
            call mkhbnd(mcaaat, atomin, hbond, hbonde, seqcod, chnend, 
     +                  seqlen)
c           ** calculate OOI numbers
            call ooi(atomin, mcaaat, seqbcd, oois, seqlen)
c           ** find main chain torsion angles

CHECK v.1.0-->
C            call mkangl(mcaaat, mcangs, atomin, chnsz, chnct, caonly,
C     +                  seqlen)
            call mkangl(mcaaat, mcangs, atomin, chnsz, chnct, caonly,
     +                  seqlen, scaaat, scatin)
CHECK v.1.0<--

c           ** find sidechain chi angles
            if (outsty .eq. 'F')
     +         call mksang(mcaaat, atomin, scaaat, scangs, scatin, 
     +                     seqcod, aacod3, hetcod, seqlen)
c           ** determine turns, bridges and sheets
            call mktb(hbond, ssbond, mcangs, bridpt, chnend, seqlen)
c           ** determine bend residues
            call mkbend(mcangs, ssbond, seqlen)
c           ** make Kabsch and Sander type 2ndry structure summary
            call mksumm(ssbond, sssumm, summpl, seqlen)
c           ** make the authors atructure
            call mkauth(seqbcd, autsum, ssord, astruc, seqlen)
c           ** find all disulphide bonds
            if (outsty .eq. 'F')
     +         call disulf(seqbcd, dsdsts, scangs, scaaat, scatin, 
     +                     atomin, mcaaat, seqcod, dsdef, dspart, 
     +                     astruc, seqlen)
c           ** produce the output file
CHECK v.3.2-->
C            call result(hbonde, hbond,  aacod3, aacod1, seqcod, seqbcd, 
C     +                  brksym, ssbond, bridpt, sssumm, summpl, autsum, 
C     +                  mcangs, altcod, aastd,  atpres, thmprs, scangs,
C     +                  cadsts, brknam, headln, hedlen, outsty, dsdsts,
            call result(hbonde, aacod3, aacod1, seqcod, seqbcd, 
     +                  brksym, ssbond, sssumm, summpl, autsum, 
     +                  mcangs, scangs,
     +                  outsty, dsdsts,
CHECK v.3.2<--

CHECK v.1.0-->
C     +                  dsdef,  dspart, hetcod, oois,   seqlen)
CHECK v.2.0-->
C     +                  dsdef,  dspart, hetcod, oois,   seqlen, BVALUE)
CHECK v.3.2-->
C     +                  dsdef,  dspart, hetcod, oois,   seqlen, BVALUE,
C     -                  MCBVAL,SCBVAL,NMCATM,NSCATM)
     +                  hetcod, oois,   seqlen, BVALUE,
CHECK v.3.5-->
C     -                  MCBVAL,SCBVAL,NMCATM,NSCATM,MCBSTD,SCBSTD)
     -                  MCBVAL,SCBVAL,NMCATM,NSCATM,MCBSTD,SCBSTD,
     -                  SAVNAM)
CHECK v.3.5<--
CHECK v.3.2<--
CHECK v.2.0<--
CHECK v.1.0<--

CHECK v.3.6.2-->
C----       Write out any dislphide bridge partners
            CALL WRISSB(FILDIS,dspart,seqlen,seqcod,seqbcd,aacod3,
     -          hetcod)
CHECK v.3.6.2<--

CHECK v.1.0-->
C            close(unit=outfi)
CHECK v.1.0<--

         else
c           close(unit=outfi)
            if (readok .and. caonly .and. outsty .eq. 'F') then
CHECK v.3.2-->
C               outnam = outdir(1:outln)//brknam(corepl:corend) // '.sca'
               iostat = ok
CHECK v.3.2<--

CHECK v.1.0-->
C               open(unit=outfi,file=outnam,status='unknown',
C     +              iostat=iocode)
CHECK v.1.0<--

               if (iostat .ne. ok) then
CHECK v.3.2-->
C                  write(*,*)'Unable to open output file ',
C     +                      outnam(1:outlen)
CHECK v.3.2<--
                  close(unit=errfi)
               else
                  call chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, 
     +                        chnct,  chnend, caonly, seqlen)
CHECK v.3.2-->
C                  call mkcadt(mcaaat, atomin, cadsts, seqbcd, chnsz,
C     +                        chnct,  seqlen)
                  call mkcadt(mcaaat, atomin, cadsts, seqbcd,
     +                        seqlen)
CHECK v.3.2<--
                  call ooi(atomin, mcaaat, seqbcd, oois, seqlen)

CHECK v.1.0-->
C                  call mkangl(mcaaat, mcangs, atomin, chnsz, chnct, 
C     +                        caonly, seqlen)
                  call mkangl(mcaaat, mcangs, atomin, chnsz, chnct,
     +                        caonly, seqlen, scaaat, scatin)
CHECK v.1.0<--

                  call mkauth(seqbcd, autsum, ssord, astruc, seqlen)
CHECK v.3.2-->
C                  call carslt(aacod3, aacod1, seqcod, seqbcd, mcangs, 
C     +                        cadsts, brksym, altcod, aastd,  brknam, 
C     +                        headln, hedlen, autsum, hetcod, oois, 
C     +                        seqlen)
                  call carslt(aacod3, aacod1, seqcod, hetcod, seqlen)
CHECK v.3.2<--
               end if
            else
               close(unit=errfi)
            end if
            if (astruc) close(unit=authfi)
CHECK v.2.1.3-->
C            if (caonly .and. readok) write(*,*)'C-alpha only file'
C            if (.not. readok) write(*,*)'Error while reading file'
            if (caonly .and. readok) write(*,*)'**** C-alpha only file'
            if (.not. readok) write(*,*)'**** Error while reading file'
CHECK v.2.1.3<--
         end if
      else
         if (fistat .eq. opnera) then
CHECK v.3.2-->
C            write(*,*)'Unable to open the file of file names ',
C     +                 finam(1:namlen)
CHECK v.3.2<--
            finuse = .false.
         else          
            if (fistat .eq. opnerb) then
               write(*,*)'Unable to open the Brookhaven format file ',
     +                   brknam(1:namlen)
            else
               if (fistat .eq. opnerc) then
CHECK v.3.2-->
C                  write(*,*)'Unable to open output file ',
C     +                      outnam(1:outlen)
CHECK v.3.2<--
               end if
            end if
         end if
      end if
c
c     if they want more we'll cooperate
c

CHECK v.1.0-->
C      if (fistat .ne. finish) then
C         fistat = ok
C         if (finuse) then
C            go to 200
C         else
C            go to 100
C         end if
C      end if
CHECK v.1.0<--

CHECK v.3.1-->
C---- If this is an NMR structure, then loop back for next model
CHECK v.3.2-->
C      IF (NMR) GO TO 100
      IF (NMR .AND. fistat.eq.ok) GO TO 100
CHECK v.3.2<--

C---- If this is an NMR structure, show how many models have been processed
 500  CONTINUE
      IF (NMR) THEN
          PRINT*
          PRINT*, '* NMR ensemble comprises', NMODEL, ' model ',
     -        'structures'
      ENDIF
CHECK v.3.1<--

CHECK v.3.2-->
      if (fistat .ne. finish) then
         fistat = ok
         if (finuse) then
            goto 200
         end if
      end if
CHECK v.3.2<--

CHECK v.1.0-->
      GO TO 999

C---- Fatal errors
 900  CONTINUE
CHECK v.2.1.4-->
C      PRINT*, '**** ERROR. Unable to open .res file'
      PRINT*, '**** ERROR. Unable to open .rin file'
CHECK v.2.1.4<--
      GO TO 999

CHECK v.3.1-->
 902  CONTINUE
      PRINT*, '*** ERROR opening mplot.in file. Program aborted'
      GO TO 999
CHECK v.3.1<--

 999  CONTINUE
CHECK v.1.0<--
CHECK v.3.0-->
      PRINT*, '* Program completed'
CHECK v.3.0<--

      END
c****************************************************************************
      SUBROUTINE RDBRK( mcaaat, atname, atomin, aacod3, hetcod, seqcod,
     +                  seqbcd, scaaat, scatin, altcod, aastd, atpres,
     +                  thmprs, caonly, astruc, headln, hedlen, readok,
     +                  dudcod, modcod, stdmod, nstdmd, modpl, seqlen,

CHECK v.1.0-->
C     +                  batnam)
CHECK v.2.0-->
C     +                  batnam, BVALUE)
CHECK v.3.1-->
C     +                  batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM)
CHECK v.3.2-->
C     +                  batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM,NMR)
     +                  batnam, BVALUE,MCBVAL,SCBVAL,NMCATM,NSCATM,NMR,
CHECK v.3.5-->
C     +                  MCBSTD,SCBSTD)
     +                  MCBSTD,SCBSTD,SAVNAM,inseqr,sqrsct)
CHECK v.3.5<--
CHECK v.3.2<--
CHECK v.3.1<--
CHECK v.2.0<--
CHECK v.1.0<--

      INCLUDE 'sstruc.par'

      Integer seqlen, seqcod(maxres), hedlen
      Character aacod3*(abetsz*3), atname*(mnchna*4-4),seqbcd(maxres)*6,
     +          headln*132, altcod(maxres)*1, aastd(maxres)*1,
     +          atpres(3,maxres)*1, thmprs(maxres)*1, hetcod*63,
     +          modcod*(nmod*3), stdmod*(nmod+1), nstdmd*(nmod+1),
     +          modpl*(nmod), dudcod*105,batnam*(mnchna*4+4)
      Real mcaaat(ncoord,mnchna,maxres), scaaat(ncoord,sdchna,maxres)
      Logical Atomin(mnchna,maxres), scatin(sdchna,maxres), Caonly,
     +        astruc, readok
      Integer aaconv

CHECK v.3.5-->
      CHARACTER*3  UNKNAM, SAVNAM(MAXRES)
CHECK v.3.5<--

CHECK v.1.0-->
      REAL      BVALUE(MAXRES)
CHECK v.1.0<--

CHECK v.2.0-->
      INTEGER   NMCATM(MAXRES), NSCATM(MAXRES)
      REAL      MCBVAL(MAXRES), SCBVAL(MAXRES)
CHECK v.2.0<--

CHECK v.3.1-->
      LOGICAL   NMR
CHECK v.3.1<--

CHECK v.3.2-->
      REAL      MCBSTD(MAXRES), SCBSTD(MAXRES)
CHECK v.3.2<--

c
c     local variables
c
CHECK v.3.2-->
      save scname, altnam, altpos, nscats, cbet
CHECK v.3.2<--

      integer   aapln
      character allspc*16
      parameter (aapln = 13, allspc = '                ')
      real      acoord(ncoord), oldocc, occup, therm, oldthm
      logical   altset, xney
      integer   atcnt,attype, I, aacnt, wkpl,
     +          clslen, cmplen, soulen, sctype, altpos(sideaa), respl,
     +          scatct, thmcnt, aacode, nscats(sideaa), hedend, sqrsct
      character aain(maxres)*3, atmnam*4, resnam*3, altloc*1,
     +          brkrec*80, oldchn*1, seqnum*6,
     +          innum*6, keywd*(keylen), oldatm*4, setloc*1,
     +          class*40, compnd*50, source*50, brcode*4, cbet*4,
     +          scname(sideaa)*16, altnam(3)*16, contch*1, incode*4,
CHECK v.3.5-->
C     +          wkstr*50, inseqr*150, setatm*4
     +          wkstr*50, setatm*4
      CHARACTER*(*) inseqr
CHECK v.3.5<--
      integer   strlen, resfnd
      data scname(1)  /allspc/, scname(2) /allspc/
      data scname(3)  /' SG             '/
      data scname(4)  /' CG  OD1 OD2    '/
      data scname(5)  /' CG  CD  OE1 OE2'/
      data scname(6)  /' CG  CD1 CD2    '/, scname(7) /allspc/
      data scname(8)  /' CG  ND1 CD2    '/
      data scname(9)  /' CG1 CD1        '/, scname(10) /allspc/
      data scname(11) /' CG  CD  CE  NZ '/
      data scname(12) /' CG  CD1 CD2    '/
      data scname(13) /' CG  SD  CE     '/
      data scname(14) /' CG  OD1 ND2    '/, scname(15) /allspc/
      data scname(16) /' CG             '/
      data scname(17) /' CG  CD  OE1 NE2'/
      data scname(18) /' CG  CD  NE  CZ '/
      data scname(19) /' OG             '/
      data scname(20) /' OG1 CG2        '/, scname(21) /allspc/
      data scname(22) /' CG1 CG2        '/
      data scname(23) /' CG  CD1 CD2    '/, scname(24) /allspc/
      data scname(25) /' CG  CD1 CD2    '/
      data scname(26) /allspc/, scname(27) /allspc/
      data scname(28) /' CG  CD  OE     '/
      data scname(29) /' CG  CD  CE  NZ '/
      data scname(30) /' CB1 CB2        '/
      data scname(31) /' CG  CD1 CD2    '/
      data scname(32) /'SEG  OD1 OD2    '/
      data scname(33) /' CM             '/
      data scname(34) /' SG             '/
      data scname(35) /' CG  CD1 CD2    '/
      data scname(36) /'                '/
      data scname(37) /' CG  CD  CE  NZ '/
      data scname(38) /' CM             '/
      data scname(39) /' CG  CD1 CD2    '/
      data scname(40) /'                '/
      data scname(41) /' CG  CD1 CD2    '/
      data scname(42) /' CG  CD1 CD2    '/
      data scname(43) /' CG             '/
      data scname(44) /' CG  CD  CE     '/
      data scname(45) /' CG1 CG2        '/
      data scname(46) /' CG1 CD1        '/
      data scname(47) /' CG  CD1 CD2    '/
      data scname(48) /' CG  CD  CE     '/
      data scname(49) /allspc/, scname(50) /' CG  CD1 CD2    '/
      data altpos / 7*0, 1, 5*0, 2, 2*0, 3, 12*0, 13*0, 8*0/
      data altnam(1)  /' CG  AD1 AD2    '/
      data altnam(2)  /' CG  AD1 AD2    '/
      data altnam(3)  /' CG  CD  AE1 AE2'/
      data nscats / 1, 4, 2, 4, 5, 7, 0, 6, 4, -1, 5, 4, 4, 4, -1,
     +              3, 5, 7, 2, 3, -1, 3, 10, -1, 8, 5, -1, 4, 12,
     +              2, 7, 4, 2, 2, 7, 8, 6, 1, 10, 3, 5, 9, 2,  4,
     +              3, 4, 7, 4, 1, 7/
      data cbet /' CB '/
CHECK v.3.2-->
C      save scname, altnam, altpos, nscats, cbet
CHECK v.3.2<--

c
c     Routine to read the brookhaven data file and extract the sequence
c     and the atom names and coordinates. Checking is also done for
c     residue compatibility between sequence and atom. C-alpha only
c     files are flagged. The sequence residues are expected together,
c     then the secondary structure records (sheets together), and then
c     the atom records which are expected in residue order.
c
CHECK v.3.5.2-->
      clslen = 0
CHECK v.3.5.2<--
CHECK v.3.5.3-->
      soulen = 0
      cmplen = 0
CHECK v.3.5.3<--
      seqlen = 0
      hedlen = 0
CHECK v.3.5-->
C      inseqr = space
C      sqrsct = 0
CHECK v.3.5<--
      brcode = '    '
      astruc = .false.
  100 continue
      read(brkfi,'(A)',end=3000) brkrec
      keywd = brkrec(1:keylen)
CHECK v.3.5-->
      IF (KEEPAL .AND. KEYWD.EQ.'ATZERO') KEYWD = 'ATOM  '
      IF (KEEPAL .AND. KEYWD.EQ.'HEZERO') KEYWD = 'HETATM'
CHECK v.3.5<--
      if (keywd .eq. seqkey) go to 200
      if (keywd .eq. hedkey) go to 500
      if (keywd .eq. cmpkey) go to 600
      if (keywd .eq. soukey) go to 700
      if (keywd .eq. atmkey) go to 900
      if (keywd .eq. htakey) go to 900
CHECK v.3.1-->
CHECK v.3.5-->
C      IF (KEYWD .EQ. ENDMDL) GO TO 3000
      IF (KEYWD .EQ. ENDMDL .AND. aacnt.GT.0) GO TO 3000
CHECK v.3.5<--
CHECK v.3.1<--
CHECK v.3.5-->
C      if (keywd .eq. modkey) then
      if (keywd .eq. modkey .AND. aacnt.GT.0) then
CHECK v.3.5<--
         readok = .false.
         seqlen = 0
CHECK v.3.1-->
         PRINT*, '*** ERROR: MODEL record encountered unexpectedly ',
     -       'in PDB file'
CHECK v.3.1<--
         go to 3000
      end if
      if (keywd .eq. hlxkey .or. keywd .eq. shtkey .or.
     +    keywd .eq. trnkey .or. keywd .eq. dsfkey) then
         if (.not. astruc) then
            open(unit = authfi, status = 'SCRATCH')
            astruc = .true.
         end if
         write(authfi,'(a)') brkrec
      end if
      go to 100

  200 continue
      read (brkrec,220) (aain(i),i = 1, aapln)
CHECK v.3.2-->
C  220 format(18x,13(x,a3))
  220 format(18x,13(1x,a3))
CHECK v.3.2<--
      do 250, i = 1, aapln
      if (aain(i) .ne. '   ') then
         if (resfnd(inseqr,aain(i)) .eq. 0) then
            sqrsct = sqrsct + 1
            inseqr((sqrsct-1)*3+1:sqrsct*3) = aain(i)
         end if
      end if
  250 continue
      go to 100
c
c     set up header records
c
  500 continue
      read(brkrec,520) contch, wkstr, incode
  520 format(9x,a1,a40,12x,a4)
      if (contch .eq. space) then
         brcode = incode
         class = wkstr
         clslen = strlen(class)
      end if
      go to 100
  600 continue
      read(brkrec,620) contch, wkstr
  620 format(9x,a1,a50)
      if (contch .eq. space) then
         compnd = wkstr
         cmplen = strlen(compnd)
      end if
      go to 100
  700 continue
      read(brkrec,720) contch, wkstr
  720 format(9x,a1,a50)
      if (contch .eq. space) then
         source = wkstr
         soulen = strlen(source)
      end if
      go to 100
c
c     set up header if any exists
c
  900 continue
      headln(1:6) = brcode // '  '
      hedlen = 6
      if (clslen .gt. 0) then
         headln(7:8+clslen) = class(1:clslen) // '  '
         hedlen = hedlen + clslen + 2
      end if
      if (cmplen .gt. 0) then
         hedend = hedlen + cmplen + 2
         if (hedend .gt. 132) hedend = 132
         headln(hedlen+1:hedend) = compnd(1:cmplen) // '  '
         hedlen = hedend
      end if
      if (soulen .gt. 0) then
         hedend = hedlen + soulen+ 2
         if (hedend .gt. 132) hedend = 132
         headln(hedlen+1:hedend)  = source(1:soulen)
         hedlen = hedend
      end if
      if (hedlen .gt. 132) hedlen = 132
c
c     look for next atom records to process
c
 1000 continue
      if (brkrec(1:keylen) .ne. atmkey .and.
     +   brkrec(1:keylen) .ne. htakey) then
         read(brkfi,'(A)',end=3000) brkrec
CHECK v.3.1-->
CHECK v.3.5-->
C         IF (BRKREC(1:KEYLEN) .EQ. ENDMDL) GO TO 3000
         IF (BRKREC(1:KEYLEN) .EQ. ENDMDL .AND. aacnt.GT.0)
     -       GO TO 3000
CHECK v.3.5<--
CHECK v.3.1<--
CHECK v.3.5-->
      IF (KEEPAL .AND. BRKREC(1:KEYLEN).EQ.'ATZERO')
     -    BRKREC(1:KEYLEN) = 'ATOM  '
      IF (KEEPAL .AND. BRKREC(1:KEYLEN).EQ.'HEZERO')
     -    BRKREC(1:KEYLEN) = 'HETATM'
CHECK v.3.5<--
         go to 1000
      end if
c
c     have first atom record for an aa
c     for each aa get the seqnum and resnam and check against the
c     read in sequence. If ok initialise to no atoms present and
c     get the coords
c
      aacnt = 0
      do 1234 i = 1, maxres
      aastd(i) = space
 1234 continue
      read(brkrec,1080)oldchn
 1080 format(21x,a1)
 
 1100 continue
      keywd = brkrec(1:keylen)
      aacode = 0
      altset = .false.
      setatm = '    '
      read(brkrec,1120) oldatm, setloc, resnam, seqnum, oldocc, oldthm
 1120 format(12x,a4,a1,a3,1x,a6,27x,2f6.2)
      if (sqrsct .ne. 0 .and. resfnd(inseqr,resnam) .eq. 0) go to 1400
      if (resfnd(dudcod,resnam) .ne. 0) go to 1400
      respl = resfnd(modcod,resnam)
      if (respl .ne. 0) then
         respl = (respl - 1) / 3 + 1
         if (Modpl(respl:respl) .eq. '+') then
            aastd(aacnt+1) = stdmod(respl:respl)
         else
            if (aastd(aacnt) .eq. 'S') then
               aastd(aacnt) = stdmod(respl:respl)
            else
               if (aastd(aacnt) .eq. 'O') then
                  aastd(aacnt) = nstdmd(respl:respl)
               else
                  if (llt(aastd(aacnt),'a')) then
                     aastd(aacnt) = stdmod(nmod+1:nmod+1)
                  else
                     aastd(aacnt) = nstdmd(nmod+1:nmod+1)
                  end if
               end if
            end if
         end if
         go to 1400
      end if

CHECK v.3.6-->
      IF (RESNAM.EQ.'HOH') GO TO 1400
CHECK v.3.6<--
      aacode = AAconv(resnam, aacod3, hetcod)

CHECK v.3.5-->
C---- If not one of the residue-types hard-coded in the program, then
C     process as though UNKnown
      unknam = 'UNK'
      IF (AACODE.LE.0) AACODE = AAconv(unknam, aacod3, hetcod)
CHECK v.3.5<--

      if (aacode .le. 0) go to 1400
      aacnt = aacnt + 1
      if (aastd(aacnt) .eq. space) then
         aastd(aacnt) = 'S'
         if (aacode .gt. stdaa) aastd(aacnt) = 'O'
      else
         if (aacode .gt. stdaa) then
            wkpl = index(stdmod,aastd(aacnt))
            aastd(aacnt) = nstdmd(wkpl:wkpl)
         end if
      end if
      seqcod(aacnt) = aacode
CHECK v.3.5-->
      savnam(aacnt) = resnam
CHECK v.3.5<--
      seqbcd(aacnt) = seqnum
      atcnt = 0
      scatct = 0
      thmcnt = 0
      thmprs(aacnt) = space
      do 1200 i = 1, mnchna
      atomin(i,aacnt) = .false.
 1200 continue
      do 1220 i = 1, sdchna
      scatin(i,aacnt) = .false.
 1220 continue
      do 1240 i = 1, 3
      atpres(i,aacnt) = space
 1240 continue
c
c     process coordinates find atom type and if main chain set in coords
c
 1300 continue
      read(brkrec,1320)atmnam,altloc,(acoord(i),i=1,3), occup, therm
 1320 format(12x,a4,a1,13x,3f8.3,2f6.2)
      if (atmnam .eq. ' OXT' .or. atmnam .eq. ' NXT' .or.
     +    atmnam(2:2) .eq. 'H' .or. atmnam(2:2) .eq. 'D')
     +    go to 1400
      if(aacode.lt.44) then
         attype = index(atname,atmnam)
      else
         if (aacode.gt.44.and.aacode.le.50) then
            attype = index(batnam,atmnam)
         endif
      endif
      if (altloc .ne. space) then
         if (setatm .eq. '    ') then
            setatm = atmnam
            setloc = altloc
            oldocc = occup
            oldthm = therm
         end if
         if (.not. altset) then
            if (setatm .eq. atmnam) then
               if (occup .gt. oldocc) then
                  setloc = altloc
                  oldocc = occup
                  if (xney(oldthm,0.0)) thmcnt = thmcnt - 1
                  oldthm = therm
                  if (attype .ne. 0) then
                     atcnt = atcnt - 1
                  else
                     scatct = scatct - 1
                  end if
               end if
            else
               altset = .true.
            end if
         end if
      end if

      if (xney(therm,0.0).and. (altloc .eq. setloc .or.
     +          altloc .eq. space)) thmcnt = thmcnt + 1
      if (attype .ne. 0 .and. (altloc .eq. setloc .or.
     +   altloc .eq. space)) then
         attype = attype / 4 + 1
         if (attype.le.4) then

CHECK v.2.0-->
            MCBVAL(AACNT) = MCBVAL(AACNT) + THERM
            NMCATM(AACNT) = NMCATM(AACNT) + 1
CHECK v.2.0<--
CHECK v.3.2-->
            MCBSTD(AACNT) = MCBSTD(AACNT) + THERM * THERM
CHECK v.3.2<--

            do 1350 i = 1, ncoord
            mcaaat(i,attype,aacnt) = acoord(i)
 1350       continue
            atomin(attype,aacnt) = .true.
         endif
         if (attype .ne. calph) caonly = .false.
         atcnt = atcnt + 1
      end if
      if (attype .eq. 0 .and. (altloc .eq. setloc .or.
     +   altloc .eq. space)) then

CHECK v.2.0-->
         SCBVAL(AACNT) = SCBVAL(AACNT) + THERM
         NSCATM(AACNT) = NSCATM(AACNT) + 1
CHECK v.2.0<--
CHECK v.3.2-->
         SCBSTD(AACNT) = SCBSTD(AACNT) + THERM * THERM
CHECK v.3.2<--

         scatct = scatct + 1
         if (aacode .le. sideaa) then
            if (atmnam .eq. cbet) then
               sctype = 1
            else
               sctype = index(scname(aacode),atmnam)
               if (sctype .eq. 0 .and. altpos(aacode) .ne. 0)
     +            sctype = index(altnam(altpos(aacode)),atmnam)
               if (sctype .ne. 0) sctype = sctype / 4 + 2
            end if
            if (sctype .ne. 0) then     
               scatin(sctype,aacnt) = .true.

CHECK v.1.0-->
               IF (ATMNAM(3:4).EQ.'G ' .OR. ATMNAM(3:4).EQ.'G1') THEN
                   BVALUE(AACNT) = THERM
               ENDIF
CHECK v.1.0<--

               do 1375 i = 1, ncoord
               scaaat(i,sctype,aacnt) = acoord(i)
 1375          continue
            end if
         end if
      end if
c
c     search for the next atom record and treat as appropriate
c     process coords only if same aa 
c     if new aa start aa process again after setting the count fields
c
 1400 continue
      read(brkfi,'(a)',end=2000) brkrec
CHECK v.3.1-->
CHECK v.3.5-->
C      IF (BRKREC(1:KEYLEN) .EQ. ENDMDL) GO TO 2000
      IF (BRKREC(1:KEYLEN) .EQ. ENDMDL .AND. aacnt.GT.0)
     -    GO TO 2000
CHECK v.3.5<--
CHECK v.3.1<--
CHECK v.3.5-->
      IF (KEEPAL .AND. BRKREC(1:KEYLEN).EQ.'ATZERO')
     -    BRKREC(1:KEYLEN) = 'ATOM  '
      IF (KEEPAL .AND. BRKREC(1:KEYLEN).EQ.'HEZERO')
     -    BRKREC(1:KEYLEN) = 'HETATM'
CHECK v.3.5<--
      If (brkrec(1:keylen) .ne. endkey) then
         If (brkrec(1:keylen) .eq. atmkey .or.
     +      brkrec(1:keylen) .eq. htakey) then
            read(brkrec,'(21x,a6)') innum
            if (innum.eq.seqnum.and.keywd.eq.brkrec(1:keylen)) then
               if (aacode .le. 0) go to 1400
               go to 1300
            else
               if (aacode .gt. 0)
     +            call cntset(atpres, thmprs, altcod, nscats, aacode, 
     +                        aacnt, atcnt, scatct, thmcnt, setloc)
               If (aacnt .ge. maxres) then
CHECK v.2.1.3-->
C                  Write(errfi,*)'Too many amino acids',
C     +                          '  - input terminated'
                  Write(errfi,*)'**** Too many amino acids',
     +                          '  - input terminated'
CHECK v.2.1.3<--
                  goto 2100
               else
                  if (aacode .gt. 0) aastd(aacnt + 1) = space
                  go to 1100
               end if
            end if
         end if
         go to 1400
      end if
 2000 continue
      if (aacode .gt. 0)
     +   call cntset(atpres, thmprs, altcod, nscats, aacode, 
     +               aacnt, atcnt, scatct, thmcnt, setloc)
 2100 continue
      Seqlen = aacnt
 3000 continue
      if (seqlen .eq. 0) readok = .false.
CHECK v.3.1-->
C      close(unit=brkfi)
      IF (.NOT.NMR) THEN
          close(unit=brkfi)
      ENDIF
CHECK v.3.1<--
      END
c*****************************************************************************
      Subroutine cntset(atpres, thmprs, altcod, nscats, aacode, 
     +                  aacnt,  atcnt,  scatct, thmcnt, setloc)
c
c     perform end of amino acid flag setting
c
      include 'sstruc.par'
      character atpres(3,maxres)*1, thmprs(maxres)*1, altcod(maxres)*1,
     +          setloc*1
      integer   nscats(sideaa), aacode, aacnt, scatct, thmcnt, atcnt
      integer   totats
      character atmchk*1

      altcod(aacnt) = setloc
      atpres(1,aacnt) = space
      atpres(2,aacnt) = space
      atpres(3,aacnt) = space
      thmprs(aacnt) = space
      totats = atcnt + scatct
      if (aacode .le. sideaa) then
         if (nscats(aacode) .ge. 0) then
            atpres(1,aacnt) = atmchk(mnchna+nscats(aacode)-1,totats)
            atpres(2,aacnt) = atmchk(mnchna-1,atcnt)
            atpres(3,aacnt) = atmchk(nscats(aacode),scatct)
            thmprs(aacnt) = atmchk(mnchna+nscats(aacode)-1,thmcnt)
         end if
      end if
      end
c*****************************************************************************
      character*1 function atmchk(Total,count)
c
c     set 'present' fields to all none or some as appropriate
c
      include 'sstruc.par'
      integer total, count
      if (total .eq. 0) then
         atmchk = space
      else
         if (count .eq. 0) then
            atmchk = 'N'
         else
            if (count .ge. total) then
               atmchk = 'A'
            else
               atmchk = 'S'
            end if
         end if
      end if
      END
c*****************************************************************************
      Integer function Strlen(string)
c
c     find out the real length of the string
c
      include 'sstruc.par'
      character*(*) string

      integer nxch
      do 200, nxch = len(string), 1, -1
        if (string(nxch:nxch) .ne. space) goto 300
  200 continue
      nxch = 0
  300 continue
      strlen = nxch
      end
c****************************************************************************
      Integer function AAconv(aacode, aacod3, hetcod)
c
c     find the numerical position of an amino acid code
c
      include 'sstruc.par'
CHECK v.3.2-->
      save lastcd
CHECK v.3.2<--
      character aacode*3, aacod3*(3*abetsz), hetcod*63
CHECK v.3.0.1--> (Order of statements changed)
      integer   resfnd
      integer   aaval
CHECK v.3.0.1<--
CHECK v.2.1.4-->
      character lastcd*3
      data lastcd / '   ' /

CHECK v.3.2-->
C      save lastcd
CHECK v.3.2<--
CHECK v.2.1.4<--

      aaval = resfnd(aacod3,aacode)
      if (aaval .eq. 0) then
         aaval = resfnd(hetcod,aacode)
         if (aaval .eq. 0) then
            aaval = unkcod
CHECK v.2.1.4-->
C            write(errfi,530) aacode, ' unknown amino acid code'
            if (aacode.ne.lastcd) then
                write(errfi,530) aacode, ' unknown amino acid code'
                lastcd = aacode
            endif
CHECK v.2.1.4<--
  530       format(a3,a)
         else
            aaval = abetsz + aaval / 3 + 1
         end if
      else
         aaval = aaval / 3 + 1
      end if
      AAconv = aaval
      END
c*****************************************************************************
      Integer function resfnd(soustr,questr)
c
c     locate query in source on 3 character boundaries
c
      character soustr*(*), questr*3
      integer   respl, nxpl

      respl = Index(soustr,questr)
  100 continue
      If (mod(respl,3) .ne. 1 .and. respl .ne. 0) then
         nxpl = index(soustr(respl+1:),questr)
         if (nxpl .ne. 0) then
            respl = respl + nxpl
         else
            respl = nxpl
         end if
         go to 100
      end if
      resfnd = respl
      End
c*****************************************************************************
      Subroutine Chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, chnct,
     +                  chnend, caonly, seqlen)
c
c     search through looking for distant peptide bonds. If any over
c     the limit are found or if a c or n atom does not exist then a
c     chain break is deemed to occur
c
      include 'sstruc.par'
      real      mcaaat(ncoord,mnchna,maxres)
      logical   atomin(mnchna,maxres), caonly
      Integer   chnsz(maxchn), chnct, seqlen, chnend(maxchn)
      character seqbcd(maxres)*6, brksym(maxres)*1
      Real      atdist

      Integer   nxres, ninchn

      do 200 nxres = 1, seqlen
      brksym(nxres) = space
  200 continue
      ninchn = 1
      chnct = 1
      do 1000, nxres = 2, seqlen
      if (caonly) then
         if (atdist(mcaaat(1,calph,nxres-1),mcaaat(1,calph,nxres))
     +      .gt. cadbnd) then
            call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                  nxres, seqlen)
         else
            ninchn = ninchn + 1
            if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) 
     +         Write(errfi,500) nxres-1, seqbcd(nxres-1),
     +                          nxres, seqbcd(nxres)
  500       format('chain letter change but no chain break between ',
     +             i4,2x,a6,' and',i4,2x,a6)
         end if
      else
         If (atomin(carb,nxres-1) .and. atomin(nmain,nxres)) then
            if (atdist(mcaaat(1,carb,nxres-1),mcaaat(1,nmain,nxres))
     +          .gt. pepbnd) then
               call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                     nxres, seqlen)
            else
               ninchn = ninchn + 1
               if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) 
     +            Write(errfi,500) nxres-1, seqbcd(nxres-1),
     +                  nxres, seqbcd(nxres)
            end if
         else
            call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                  nxres, seqlen)
         end if
      end if
 1000 continue
      call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd,
     +            seqlen+1, seqlen)
      END
c*****************************************************************************
      SUBROUTINE Chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                  nxres, seqlen)
      include 'sstruc.par'
      integer   chnsz(maxchn), chnct, ninchn, seqlen, nxres,
     +          chnend(maxchn)
      character seqbcd(maxres)*6, brksym(maxres)*1

      chnsz(chnct) = ninchn
      if (chnct .eq. 1) then
         chnend(chnct) = ninchn
      else
         chnend(chnct) = chnend(chnct-1) + ninchn
      end if
      ninchn = 1
      if (nxres .le. seqlen) then
         brksym(nxres-1) = brkch
         brksym(nxres) = brkch
         write(errfi,100) nxres-1, seqbcd(nxres-1), nxres, seqbcd(nxres)
  100    format('chain break between ',i4,'(',a6,') and ',i4,'(',a6,')')
         chnct = chnct + 1
CHECK v.3.3.1-->
         if (chnct.gt.maxchn) then
             chnct = maxchn
             print*, '*** Maximum number of chain-breaks exceeded',
     -           maxchn
         endif
CHECK v.3.3.1<--
      end if
      end
c*****************************************************************************
CHECK v.3.2-->
C      Subroutine mkcadt(mcaaat, atomin, cadsts, seqbcd, chnsz, chnct,
C     +                  seqlen)
      Subroutine mkcadt(mcaaat, atomin, cadsts, seqbcd,
     +                  seqlen)
CHECK v.3.2<--
c
c     generate the calpha distances for I+1 to I+11. take note of nulls
c     accross brookhaven chain breaks but not accross disconections at 
c     this point
c
      include 'sstruc.par'
CHECK v.3.2-->
C      integer   seqlen, chnct, chnsz(maxchn)
      integer   seqlen
CHECK v.3.2<--
      real      mcaaat(ncoord,mnchna,maxres), cadsts(11,maxres)
      logical   atomin(mnchna,maxres)
      character seqbcd(maxres)*6

      integer   nxres, nxdist
      real      mxdist, adist
      parameter (mxdist = 99.9)
      real      atdist

      do 1000, nxres = 1, seqlen
      do 500, nxdist = 1, 11
      if (nxres+nxdist .gt. seqlen .or.
     +   .not. atomin(calph,nxres)) then
         cadsts(nxdist,nxres) = -1
      else
         if (seqbcd(nxres)(1:1) .ne. seqbcd(nxres+nxdist)(1:1)
     +      .or. .not. atomin(calph,nxres+nxdist)) then
            cadsts(nxdist,nxres) = -1
         else
            adist = atdist(mcaaat(1,calph,nxres),
     +                     mcaaat(1,calph,nxres+nxdist))
            if (adist .gt. mxdist) adist = mxdist
            cadsts(nxdist,nxres) = adist
         end if
      end if
  500 continue
 1000 continue
      END
c*****************************************************************************
CHECK v.3.2-->
C      Subroutine CreatH(mcaaat, atomin, chnsz, chnct, seqlen)
      Subroutine CreatH(mcaaat, atomin, chnsz, chnct)
CHECK v.3.2<--
c
c     Procedure to create hydrogen atoms. These are deemed to be
c     one angstrom from the nitrogen atom and angled so that N-H
c     is parallel to C=O before. By taking the vector from O to C
c     the N-H bond is generated on the other side of C-N to the O.
c
      INCLUDE 'sstruc.par'

CHECK v.3.2-->
C      integer chnsz(maxchn), chnct, seqlen
      integer chnsz(maxchn), chnct
CHECK v.3.2<--
      real    mcaaat(ncoord,mnchna,maxres)
      logical atomin(mnchna,maxres)
c
c     local variables
c
      Integer nxres, i, chnst, nxchn
      Real    atvec(ncoord), colen
c
      chnst = 0
      do 2000, nxchn = 1, chnct
      do 1900, nxres = chnst + 2, chnst + chnsz(nxchn)
      If (atomin(nmain,nxres) .and. atomin(carb,nxres-1) .and.
     +   atomin(oxyg,nxres-1)) then
         Colen = 0.0
         do 300 i = 1,ncoord
         atvec(i) = mcaaat(i,carb,nxres-1) - mcaaat(i,oxyg,nxres-1)
         colen = colen + atvec(i) ** 2
  300    continue
         colen = sqrt(colen)
         do 500 i = 1, ncoord
         atvec(i) = atvec(i) / colen
  500    continue
         do 700 i = 1, ncoord
         mcaaat(i,nhyd,nxres) = mcaaat(i,nmain,nxres) + atvec(i)
  700    continue
         atomin(nhyd,nxres) = .true.
      else
         write(errfi,1100) nxres
 1100    format('hydrogen atom not generated for residue ',i4)
         atomin(nhyd,nxres) = .false.
      end if
 1900 continue
      chnst = chnst + chnsz(nxchn)
 2000 continue
      END
c*****************************************************************************
      Subroutine MkHBnd(Mcaaat, atomin, Hbond, Hbonde, seqcod, chnend,
     +                  seqlen)
c
c     Procedure to discover the possible hydrogen bonds between two
c     amino acids. Each amino acid is allowed 4 bonds - 2 from the 
c     C=O and 2 from the N-H. Should it be that 3 or more bonds are
c     possible the most favourable are kept. A hydrogen
c     bond is deemed to exist if the Kabsch and Sander energy formula
c     gives a value less than the cutoff. The angle and distance
c     approach is implicitly built in to this method. The formula
c     corresponds to a maximum distance of 5.2 A for in line atoms
c     and a maximum angular deviation of 63 degrees. A check on the
c     maximum distance is used to abandon unnecessary calculations.
c     the formula is
c     E = q1*q2*(1/d(ON)+1/d(CH)-1/d(OH)-1/d(CN))*f
c     Residues must be separated by at least one residue to bond
c     Proline cannot be a donor of a hydrogen bond
c
      INCLUDE 'sstruc.par'

      Integer Seqlen
      Real    Mcaaat(ncoord,mnchna,maxres), Hbonde(maxbnd,maxres)
      Integer Hbond(maxbnd,maxres), seqcod(maxres), chnend(maxchn)
      logical Atomin(mnchna,maxres)
      Real    Atdist
c
c     local variables
c
CHECK v.3.2-->
C      integer    ksf, mndist
C      real       mxdist, maxeng, ksq1, ksq2, maxcad, englim
      integer    ksf
      real       mxdist, maxeng, ksq1, ksq2, maxcad, englim, mndist
CHECK v.3.2<--
      parameter (mxdist = 5.2, maxeng = -0.5, mndist = 0.5,
     +           englim = -9.9, ksq1 = 0.42, ksq2 = 0.2, ksf = 332, 
     +           maxcad = 8.0)
      Integer    Nxres, Othres, i, nbonds, fstchn, othchn
      logical    nhpres, xeqy
      real       don, doh, dch, dcn, energy, cadist

      do 200 nxres = 1, seqlen
      do 200 i = 1,maxbnd
      hbond(i,nxres) = 0
      hbonde(i,nxres) = 0.0
  200 continue
      nbonds = 0
      fstchn = 1
      do 2000 nxres = 1, seqlen
      if (nxres .gt. chnend(fstchn)) fstchn = fstchn + 1
      nhpres = .false.
      if (atomin(nmain,nxres) .and. atomin(nhyd,nxres) .and.
     +   seqcod(nxres) .ne. procod) nhpres = .true.
      if (nhpres) then
         othchn = 1
         do 1900 othres = 1, seqlen
         if (othres .gt. chnend(othchn)) othchn = othchn + 1
         If ((iabs(nxres - othres) .eq. 1 .and. fstchn .ne. othchn)
     +      .or. iabs(nxres - othres) .ge. 2) then
            If (atomin(calph,nxres) .and. atomin(calph,othres)) then
               cadist = atdist(mcaaat(1,calph,othres),
     +                  mcaaat(1,calph,nxres))
               If (cadist .lt. maxcad) then
                  If (atomin(carb,othres) .and. 
     +               atomin(oxyg,othres))then
                     don = atdist(mcaaat(1,oxyg,othres),
     +                     mcaaat(1,nmain,nxres))
                     doh = atdist(mcaaat(1,oxyg,othres),
     +                     mcaaat(1,nhyd,nxres))
                     dch = atdist(mcaaat(1,carb,othres),
     +                     mcaaat(1,nhyd,nxres))
                     dcn = atdist(mcaaat(1,carb,othres),
     +                     mcaaat(1,nmain,nxres))
                     If (xeqy(don,0.0).or.xeqy(doh,0.0).or.
     +                   xeqy(dch,0.0).or.xeqy(dcn,0.0)) then
                        write(errfi,250)nxres, othres
CHECK v.3.4-->
C  250                   format('coincident atoms in hydrogen bonding'
  250                   format('coincident atoms in hydrogen bonding',
CHECK v.3.4<--
     +                         ' donor ',I4,' acceptor ',i4)
                     else
                        energy = ksq1 * ksq2 * ksf * 
     +                          (1.0/don + 1.0/dch - 1.0/doh - 1.0/dcn)
                        If (energy .lt. englim) then
                           write(errfi,300) nxres, othres, don, dch,
     +                                doh, dcn
  300                      format('atoms too close ',I4,2x,i4,
     +                           2x,4(2x,f8.3))
                           energy = englim
                        end if
                        If (energy .lt. maxeng) call engset(hbonde,
     +                     hbond, energy, nxres, othres, nbonds)
                     end if
                  end if
               end if
            end if
         end if
 1900    continue
      end if
 2000 continue
      write(errfi,2100) 'number of hydrogen bonds is ',nbonds
 2100 format(a,i5)
      END
c*****************************************************************************
      Subroutine Engset(hbonde, hbond, energy, donres, accres, bndct)
c
c     find a place to put this h-bond. discard third bonds if necessary
c
      Include 'sstruc.par'
      real    hbonde(maxbnd,maxres), energy
      integer hbond(maxbnd,maxres), donres, accres, bndct

      integer acpl, donpl, rejptr
      integer engpl

      donpl = engpl(energy,hbonde(nst,donres))
      acpl = engpl(energy,hbonde(cst,accres))
      If (acpl .eq. 0 .or. donpl .eq. 0) then
         write(errfi,400) donres, accres, energy
  400    format('third (+) Hbond (N-C) ',I4,2x,i4,
     +                          ' energy ',f6.2,' abandoned')
      else
         if (hbond(nst+1,donres) .ne. 0) then
            write(errfi,400) donres, hbond(nst+1,donres),
     +                    hbonde(nst+1,donres)
            rejptr = hbond(nst+1,donres) 
            call rejbnd(donres,hbond(cst,rejptr),hbonde(cst,rejptr))
            bndct = bndct - 1
         end if
         if (hbond(cst+1,accres) .ne. 0) then
            write(errfi,400) hbond(cst+1,accres), accres, 
     +                    hbonde(cst+1,accres)
            rejptr = hbond(cst+1,accres) 
            call rejbnd(accres,hbond(nst,rejptr),hbonde(nst,rejptr))
            bndct = bndct - 1
         end if
         If (acpl .eq. 1) then
            hbond(cst+1,accres) = hbond(cst,accres)
            hbonde(cst+1,accres) = hbonde(cst,accres)
         end if
         If (donpl .eq. 1) then
            hbond(nst+1,donres) = hbond(nst,donres)
            hbonde(nst+1,donres) = hbonde(nst,donres)
         end if
         hbond(cst+acpl-1,accres) = donres
         hbonde(cst+acpl-1,accres) = energy
         hbond(nst+donpl-1,donres) = accres
         hbonde(nst+donpl-1,donres) = energy
         bndct = bndct + 1
      end if
      END
c*****************************************************************************
      SUBROUTINE rejbnd(rejptr,bond,bonde)
      include 'sstruc.par'
      integer rejptr, bond(2)
      real    bonde(2)
      integer rejpl

      rejpl = 1
      if (bond(1).ne. rejptr) rejpl = rejpl + 1
      if (rejpl .eq. 1) then
         bond(1) = bond(2)
         bonde(1) = bonde(2)
      end if
      bond(2) = 0
      bonde(2) = 0.0
      end
c*****************************************************************************
      INTEGER FUNCTION Engpl(energy, bonde)
      include 'sstruc.par'
      real energy, bonde(2)
      engpl = 0
      If (energy .lt. bonde(1)) then
         engpl = 1
      else
         If (energy .lt. bonde(2)) engpl = 2
      end if
      end
c****************************************************************************   
      REAL Function Atdist(coorda,coordb)
c
c     routine to determine the distance in space between the two
c     atoms represented by the coordinate arrays.
c
      INCLUDE 'sstruc.par'

      Real    Coorda(ncoord),Coordb(ncoord)
c
      Integer i
      Real    Dist
      dist = 0.0
      do 200 i = 1,ncoord
      dist = dist + (coorda(i) - coordb(i)) ** 2
  200 continue
      atdist = sqrt(dist)
      END
c*****************************************************************************
      INTEGER Function Aplace(Bondpl)
c
c     find free place in bond table or return 0
c
      include 'sstruc.par'
      integer bondpl(2), i, wkpl
      do 200 i = 1,2
      If (Bondpl(i) .eq. 0) then
         wkpl = i
         go to 250
      end if
  200 continue
      wkpl = 0
  250 continue
      aplace = wkpl
      END
c*****************************************************************************
      Subroutine Ooi(atomin, mcaaat, seqbcd, oois, seqlen)
c
c     routine to calculate Ooi numbers. these are done for n8 and n14.
c     all atoms are included as there is noway to handle i+/-2 over
c     missing densities. Numbers are for within brookhaven chain only.
c
      include 'sstruc.par'
      logical   atomin(mnchna,maxres)
      real      mcaaat(ncoord,mnchna,maxres)
      character seqbcd(maxres)*6
      integer   oois(2,maxres), seqlen

      integer   nxaa, stchn, wkchn(maxchn), stpl, nxchn, nxpart,
     +          numchn
      character chnch*1
      real      wkdist
      real      atdist

      stchn = 1
      chnch = seqbcd(1)(1:1)
      numchn = 0
      do 200, nxaa = 1, seqlen
      if (atomin(calph,nxaa)) then
         oois(1,nxaa) = 0
         oois(2,nxaa) = 0
      else
         oois(1,nxaa) = -1
         oois(2,nxaa) = -1
      end if
      if (seqbcd(nxaa)(1:1) .ne. chnch) then
         numchn = numchn + 1
         wkchn(numchn) = nxaa - stchn
         chnch = seqbcd(nxaa)(1:1)
         stchn = nxaa
      end if
  200 continue
      numchn = numchn + 1
      wkchn(numchn) = seqlen + 1 - stchn
      stpl = 0
      do 3000, nxchn = 1, numchn
      do 2000, nxaa = stpl + 1, stpl + wkchn(nxchn) - 1
      do 1000, nxpart = nxaa + 1, stpl + wkchn(nxchn)
      if (atomin(calph,nxaa) .and. atomin(calph,nxpart)) then
         wkdist = atdist(mcaaat(1,calph,nxaa),
     +            mcaaat(1,calph,nxpart))
         if (wkdist .lt. ooi2rd) then
            oois(ooi2,nxaa) = oois(ooi2,nxaa) + 1
            oois(ooi2,nxpart) = oois(ooi2,nxpart) + 1
            if (wkdist .lt. ooi1rd) then
               oois(ooi1,nxaa) = oois(ooi1,nxaa) + 1
               oois(ooi1,nxpart) = oois(ooi1,nxpart) + 1
            end if
         end if
      end if
 1000 continue
 2000 continue
      stpl = stpl + wkchn(nxchn)
 3000 continue
      end      
c*****************************************************************************

CHECK v.1.0-->
C      Subroutine MkAngl(McAAat, McAngs, Atomin, chnsz, chnct, caonly,
C     +                  Seqlen)
      Subroutine MkAngl(McAAat, McAngs, Atomin, chnsz, chnct, caonly,
     +                  Seqlen, scaaat, scatin)
CHECK v.1.0<--

c
c     Procedure to calculate the main chain dihedral angles
c     if all the atoms are present.
c
      INCLUDE 'sstruc.par'

      Real    McAAat(ncoord,mnchna,maxres), McAngs(nmnang,maxres)
      Logical Atomin(mnchna,maxres), caonly
      Integer Seqlen, chnsz(maxchn), chnct
      Real    Dihed 

CHECK v.1.0-->
      real      scaaat(ncoord,sdchna,maxres)
      logical   scatin(sdchna,maxres)
CHECK v.1.0<--

c
c     local variables
c
CHECK v.3.2-->
      SAVE angrng, angatm, angoff, angchn
CHECK v.3.2<--
      Integer angrng(nmnang,3), angatm(nmnang,4), angoff(nmnang,4),
CHECK v.3.4.4-->
C     +        angchn(nmnang), I, j, Nxres, Nxset, chnst, nxchn, stang,
CHECK v.3.5-->
C     +        angchn(nmnang), Nxres, Nxset, chnst, nxchn, stang,
     +        angchn(nmnang), I, j, Nxres, Nxset, chnst, nxchn, stang,
CHECK v.3.5<--
CHECK v.3.4.4<--
     +        finang

CHECK v.1.0-->
C      data ((angrng(I,j),j=1,3), i=1,nmnang)
C     +      /2, 0, phi, 1, -1, psi, 1, -1, omega, 2, -2, chiral,
C     +       3, -2, kappa, 2, 0, tco/
C      data ((angatm(I,j),j=1,4), i=1,nmnang)
C     +     /carb, nmain, calph, carb, nmain, calph, carb, nmain,
C     +      calph, carb, nmain, calph, calph, calph, calph, calph,
C     +      calph, calph, calph, calph, carb, oxyg, carb, oxyg/
C      data ((angoff(I,j),j=1,4), i=1,nmnang)
C     +     /-1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, -1, 0, 1, 2,
C     +      -2, 0, 0, 2, 0, 0, -1, -1/
C      data angchn / 2, 2, 2, 4, 5, 2/
CHECK v.1.0<--

CHECK v.1.0-->
      data ((angrng(I,j),j=1,3), i=1,nmnang)
     +      /2, 0, phi, 1, -1, psi, 1, -1, omega, 2, -2, chiral,
     +       1, 0, IMPLD, 3, -2, kappa, 2, 0, tco /
      data ((angatm(I,j),j=1,4), i=1,nmnang)
     +     /carb, nmain, calph, carb, nmain, calph, carb, nmain,
     +      calph, carb, nmain, calph, calph, calph, calph, calph,
     -      CALPH, NMAIN, CARB, MCBETA,
     +      calph, calph, calph, calph, carb, oxyg, carb, oxyg /
      data ((angoff(I,j),j=1,4), i=1,nmnang)
     +     /-1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, -1, 0, 1, 2,
     +       0, 0, 0, 0,-2, 0, 0, 2, 0, 0, -1, -1 /
      data angchn / 2, 2, 2, 4, 1, 5, 2 /
CHECK v.1.0<--

CHECK v.3.2-->
C      SAVE angrng, angatm, angoff, angchn
CHECK v.3.2<--
c

CHECK v.1.0-->
      DO 100, NXRES = 1, SEQLEN
          ATOMIN(MCBETA,NXRES) = SCATIN(CBETA,NXRES)
          MCAAAT(1,MCBETA,NXRES) = SCAAAT(1,CBETA,NXRES)
          MCAAAT(2,MCBETA,NXRES) = SCAAAT(2,CBETA,NXRES)
          MCAAAT(3,MCBETA,NXRES) = SCAAAT(3,CBETA,NXRES)
 100  CONTINUE
CHECK v.1.0<--

      stang = 1
      finang = nmnang
      if (caonly) then
         stang = chiral
         finang = kappa
      end if
      do 4000, Nxset = stang, finang
      chnst = 0
      do 3000, nxchn = 1, chnct
      if (chnsz(nxchn) .lt. angchn(nxset)) then
         do 1000, nxres = chnst + 1, chnst + chnsz(nxchn)
         McAngs(angrng(nxset,3),Nxres) = Null
 1000    continue
      else
         do 2000, Nxres = chnst + angrng(nxset,1),
     +                    chnst + chnsz(nxchn) + angrng(nxset,2)
         If (atomin(angatm(nxset,1), nxres+angoff(nxset,1)) .and.
     +      atomin(angatm(nxset,2), nxres+angoff(nxset,2)) .and.
     +      atomin(angatm(nxset,3), nxres+angoff(nxset,3)) .and.
     +      atomin(angatm(nxset,4), nxres+angoff(nxset,4))) then
            McAngs(angrng(nxset,3),Nxres) = Dihed(nxset,
     +          Mcaaat(1,angatm(nxset,1), nxres+angoff(nxset,1)),
     +          Mcaaat(1,angatm(nxset,2), nxres+angoff(nxset,2)),
     +          Mcaaat(1,angatm(nxset,3), nxres+angoff(nxset,3)),
     +          Mcaaat(1,angatm(nxset,4), nxres+angoff(nxset,4)))
         else
            McAngs(angrng(nxset,3),Nxres) = Null
         end if
 2000    continue
         do 2400, nxres = chnst + 1, chnst + angrng(nxset,1) - 1
         McAngs(angrng(nxset,3),nxres) = Null
 2400    continue
         do 2600, nxres = chnst + chnsz(nxchn) + angrng(nxset,2) + 1,
     +                    chnst + chnsz(nxchn)
         McAngs(angrng(nxset,3),nxres) = Null
 2600    continue
      end if
      chnst = chnst + chnsz(nxchn)
 3000 continue
 4000 continue
      END
c*****************************************************************************
      subroutine Mksang(mcaaat, atomin, scaaat, scangs, scatin, 
     +                  seqcod, aacod3, hetcod, seqlen)
c
c     calculate chi angles from atom data. check for that alternates
c     are closest to c alpha
c
      include 'sstruc.par'
CHECK v.3.2-->
      save natsdc, extprc
CHECK v.3.2<--
      real      mcaaat(ncoord,mnchna,maxres), 
     +          scaaat(ncoord,sdchna,maxres),
     +          scangs(nsdang,maxres)
      character aacod3*(abetsz*3), hetcod*63
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres)
      integer   seqcod(maxres), seqlen

      integer   natsdc(sideaa), nxang, i, j, nxchi, nxres, nchi, aacod,
     +          extprc(sideaa), swatm, extatm, chipl, nswap, wkswap
      logical   athere(maxchi+3), xney
      real      relats(ncoord,maxchi+3), distmn, distal, ang1, ang2
      character outerr*80, wkres*3
      real      atdist, dihed
      data natsdc / 0, 0, 1, 2, 3, 2, 0, 2, 2, 0, 4, 2, 3,
     +              2, 0, 1, 3, 4, 1, 1, 0, 1, 2, 0, 2, 0, 0, 3, 4,
     +              0, 2, 1, 1, 1, 2, 0, 4, 0, 2, 0, 2, 2, 1, 3, 1,
     +              2, 2, 3, 0, 2/
      data extprc / 0, 0, 0, 4, 4, 4, 0, 0, 0, 0, 0, 2, 0,
     +              0, 0, 0, 0, 0, 0, 3, 0, 2, 0, 0, 4, 0, 0, 0, 0,
     +              0, 4, 0, 0, 0, 4, 0, 0, 0, 4, 0, 2, 4, 0, 0, 2,
     +                    0, 4, 0, 0, 4/
CHECK v.3.2-->
C      save natsdc, extprc
CHECK v.3.2<--

      nswap = 0
      outerr = space
      do 5000 nxres = 1, seqlen
      aacod = seqcod(nxres)
      nchi = 0
      nchi = natsdc(aacod)
      do 100, nxang = nchi + 1, nsdang
      scangs(nxang,nxres) = null
  100 continue
      if (nchi .gt. 0) then
         athere(1) = atomin(nmain,nxres)
         athere(2) = atomin(calph,nxres)
         do 200, i = 1, ncoord
         relats(i,1) = mcaaat(i,nmain,nxres)
         relats(i,2) = mcaaat(i,calph,nxres)
  200    continue
         extatm = 0
         if (extprc(aacod) .gt. 1) extatm = 1
         do 400 i = 1, nchi + 1 + extatm
         athere(2+i) = scatin(i,nxres)
         do 300, j = 1, ncoord
         relats(j,2+i) = scaaat(j,i,nxres)
  300    continue
  400    continue
         swatm = 0
         if (extprc(aacod) .eq. 1) then
            if (atomin(calph,nxres) .and. scatin(nchi+1,nxres)
     +         .and. scatin(nchi+2,nxres)) then
               distmn = atdist(mcaaat(1,calph,nxres),
     +                  scaaat(1,nchi+1,nxres))
               distal = atdist(mcaaat(1,calph,nxres),
     +                  scaaat(1,nchi+2,nxres))
               if (distal .lt. distmn) then 
                  swatm = 1
                  athere(nchi+3) = scatin(nchi+2,nxres)
                  do 500, i = 1, ncoord
                  relats(i,nchi+3) = scaaat(i,nchi+2,nxres)
  500             continue
               end if
            end if
         end if
         do 1000, nxchi = 1, nchi
         If (athere(nxchi) .and. athere(nxchi+1) .and.
     +      athere(nxchi+2) .and. athere(nxchi+3)) then
            scangs(nxchi,nxres) = dihed(nmnang+nxchi,
     +            relats(1,nxchi), relats(1,nxchi+1),
     +            relats(1,nxchi+2), relats(1,nxchi+3))
         else
            scangs(nxchi,nxres) = null
         end if
 1000    continue
         if (extprc(aacod) .gt. 1 .and. 
     +      xney(scangs(nchi,nxres), null)) then
            If (athere(nchi) .and. athere(nchi+1) .and.
     +         athere(nchi+2) .and. athere(nchi+3+extatm)) then
               scangs(nchi+1,nxres) = dihed(nmnang+nchi+extatm,
     +            relats(1,nchi), relats(1,nchi+1),
     +            relats(1,nchi+2), relats(1,nchi+3+extatm))
            else
               scangs(nchi+1,nxres) = null
            end if
            if (xney(scangs(nchi+1,nxres), null)) then
               If (extprc(aacod) .eq. 2 .or. extprc(aacod) .eq. 3) then
                  ang1 = scangs(nchi,nxres)
                  ang2 = scangs(nchi+1,nxres)
                  if (ang1 .lt. 0.0) ang1 = 360.0 + ang1
                  if (ang2 .lt. 0.0) ang2 = 360.0 + ang2
                  chipl = nchi
                  if (abs(ang1 - ang2) .gt. 180.0) then
                     if (ang2 .gt. ang1 .and. extprc(aacod) .eq. 2)
     +                                                chipl = nchi + 1
                     if (ang1 .gt. ang2 .and. extprc(aacod) .eq. 3)
     +                                                chipl = nchi + 1
                  else
                     if (ang2 .lt. ang1 .and. extprc(aacod) .eq. 2)
     +                                                chipl = nchi + 1
                     if (ang1 .lt. ang2 .and. extprc(aacod) .eq. 3)
     +                                                chipl = nchi + 1
                  end if
                  if (chipl .ne. nchi) then
                     swatm = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               if (extprc(aacod) .eq. 4) then
                  if (abs(scangs(nchi+1,nxres)) .lt. 
     +                         abs(scangs(nchi,nxres))) then
                     swatm = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               scangs(nchi+1,nxres) = null
            end if
         end if
         if (swatm .eq. 1) then
            if (nswap .eq. 0) 
     +         write(errfi,'(A)') 'side chain atoms swapped for '
            wkswap = mod(nswap,8) + 1
            if (seqcod(nxres) .le. abetsz) then
               wkres = aacod3(seqcod(nxres)*3-2:seqcod(nxres)*3)
            else
               wkres = hetcod((seqcod(nxres)-abetsz)*3-2:
     +                       (seqcod(nxres)-abetsz)*3)
            end if
CHECK v.2.1.4-->
C            write(outerr((wkswap-1)*10+1:wkswap*10),'(a3,x,i4,2x)')
            write(outerr((wkswap-1)*10+1:wkswap*10),'(a3,1x,i4,2x)')
CHECK v.2.1.4<--
     +        wkres, nxres
            if (wkswap .eq. 8) then
               write(errfi,'(A)') outerr
               outerr = space
            end if
            nswap = nswap + 1
         end if
      end if
 5000 continue
      if (mod(nswap,8) .ne. 0) write(errfi,'(A)') outerr
      END
c**************************************************************************
      Real Function Dihed(angnum, AtomA, AtomB, AtomC, AtomD)
c
c     routine to calculate the dihedral angle between the 4 atoms
c     given
c
      INCLUDE 'sstruc.par'

      Real    Atoma(ncoord), Atomb(ncoord), Atomc(ncoord), atomd(ncoord)
      integer angnum
      Real    Atdist
      Logical xeqy
      Real    dihatm(ncoord,dihdat), codist(ncoord,dihdat-1),
     +        atmdst(dihdat-1), dotprd(dihdat-1,dihdat-1),
     +        ang, cosang, detant
      Integer i, j, k

      do 100, I = 1, ncoord
      Dihatm(i,1) = atoma(i)
      Dihatm(i,2) = atomb(i)
      Dihatm(i,3) = atomc(i)
      Dihatm(i,4) = atomd(i)
  100 continue

      do 200, i = 1, dihdat-1
      atmdst(i) = atdist(Dihatm(1,i),Dihatm(1,I+1))
      do 150, j = 1,ncoord
      codist(j,i) = dihatm(j,I+1) - dihatm(j,i)
  150 continue
  200 continue

      do 300, i = 1, dihdat-1
      do 300, j = 1, dihdat-1
      dotprd(i,j) = 0.0
  300 continue
      do 400, i = 1, dihdat - 2
      do 400, j = i+1, dihdat - 1
      do 400, k = 1, ncoord
      dotprd(i,j) = dotprd(i,j) + codist(k,i) * codist(k,j)
  400 continue
      do 500, i = 1, dihdat-2
      do 500, j = i+1, dihdat - 1
      if (xeqy(atmdst(i), 0.0) .or. xeqy(atmdst(j), 0.0)) then
         dotprd(i,j) = 1.0
      else
         dotprd(i,j) = dotprd(i,j) / atmdst(i) / atmdst(j)
      end if
  500 continue

      if (angnum .gt. ndihed .and. angnum .le. nmnang) then
         cosang =  dotprd(1,3)
      else
         if (xeqy(dotprd(1,2), 1.0) .or. xeqy(dotprd(2,3), 1.0)) then
            cosang = 1.0
         else
            cosang = (dotprd(1,2) * dotprd(2,3) - dotprd(1,3)) /
     +               (sqrt(1.0-dotprd(1,2)**2) * 
     +                sqrt(1.0-dotprd(2,3)**2))
         end if
      end if
      if (abs(cosang) .gt. 1.0) cosang = anint(cosang/abs(cosang))
      ang = acos(cosang) * radian

      if (angnum .le. ndihed .or. angnum .gt. nmnang) then
         detant = codist(1,1) * codist(2,2) * codist(3,3) -
     +            codist(1,1) * codist(2,3) * codist(3,2) +
     +            codist(1,2) * codist(2,3) * codist(3,1) -
     +            codist(1,2) * codist(2,1) * codist(3,3) +
     +            codist(1,3) * codist(2,1) * codist(3,2) -
     +            codist(1,3) * codist(2,2) * codist(3,1)
         If (detant .lt. 0.0) ang = -ang
      end if

      dihed = ang
      END
c****************************************************************************
      Subroutine MkTB(Hbond, SSbond, mcangs, bridpt, chnend, seqlen)
c
c     routine to establish the turns and bridges
c
      INCLUDE 'sstruc.par'

CHECK v.3.2-->
      SAVE turnch, turnsz, turnpl, bsofst, bsrng, pbrdg, stchap, stchpa
CHECK v.3.2<--

      Integer   Hbond(maxbnd,maxres), seqlen, bridpt(2,maxres),
     +          chnend(maxchn)
      Real      Mcangs(nmnang,maxres)
      Character SSbond(nstruc,maxres)*1
      Integer   strdpl

      Integer   nrules, parlel
      Parameter (nrules = 3, parlel = 1)
      Integer   nxres, nxbond, nxturn, turnsz(nturns), turnpl(nturns),
CHECK v.3.4.4-->
C     +          bsofst(nrules,4), i, j, nxrule, bndi, bndj, bndptr, 
CHECK v.3.5-->
C     +          bsofst(nrules,4), i, nxrule, bndi, bndj, bndptr, 
     +          bsofst(nrules,4), i, j, nxrule, bndi, bndj, bndptr, 
CHECK v.3.5<--
CHECK v.3.4.4<--
     +          brpl, brdgpt, bsrng(nrules,2), pbrdg(nrules), brpla,
     +          brdg(4,maxres), brptr, brdir, newptr, nxptr, nxpl,
     +          exrule(2,nrules), brptpl, nwptpl, nxstrn, 
     +          lststd, stsht, finsht, nxstr, chcode, stchpl, thstpl,
     +          nxstpl, nxblg, nxbrdg, strpl, brdpt, nxbrpl,
     +          strcd(2,maxres), bstval, shtcod(maxres), fstchn

      Character turnch(nturns)*1, wkch*1, strdch*1, brdch*1,
     +          stchap*(nstchs), stchpa*(nstchs),
     +          nstchp*(nstchs*2+10)
      logical   xney
c
c     for turn calculations
c
      data turnch / '3', '4', '5'/
      data turnsz / 3, 4, 5/
      data turnpl / thrtnh, alphah, pih/
c
c     for bridge calculations
c
      data ((bsofst(i,j),j = 1,4),i = 1,nrules)
     +     / -1, 0, 0, 1, 0, 0, 0, 0, -1, -1, -2, 1/
      data ((bsrng(i,j),j = 1,2), i = 1,nrules)
     +     / 2, -1, 1, 0, 2, -1/
      data pbrdg / 1, -1, -1/
      data exrule /-1, 1, 1, 1, -1, -1/
c
c     for strand and bridge labelling
c
      data stchap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data stchpa /'abcdefghijklmnopqrstuvwxyz'/
      data nstchp(1:26) /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/,
     +     nstchp(27:52) /'abcdefghijklmnopqrstuvwxyz'/,
     +     nstchp(53:62) /'1234567890'/

CHECK v.3.2-->
C      SAVE turnch, turnsz, turnpl, bsofst, bsrng, pbrdg, stchap, stchpa
CHECK v.3.2<--

      do 100, nxres = 1, seqlen
      do 100, i = 1, nstruc
      ssbond(i,nxres) = space
  100 continue

      fstchn = 1
      do 1000, nxres = 1, seqlen
      if (nxres .gt. chnend(fstchn)) fstchn = fstchn + 1
      do 800, nxbond = 1, 2
      do 600, nxturn = 1, nturns
      if (Hbond(nxbond,nxres) .ne. 0) then
         If (Hbond(nxbond,nxres) - nxres .eq. turnsz(nxturn) .and. 
     +      Hbond(nxbond,nxres) .le. chnend(fstchn)) then
            If (SSbond(turnpl(nxturn),nxres) .eq. bndend) then
               SSbond(turnpl(nxturn),nxres) = bndbth
            else
               SSbond(turnpl(nxturn),nxres) = bndbeg
            end if
            SSbond(turnpl(nxturn),Hbond(nxbond,nxres)) = bndend
            do 400, i = nxres + 1, nxres + turnsz(nxturn) - 1
            If (SSbond(turnpl(nxturn),i) .eq. space)
     +         SSbond(turnpl(nxturn),i) = turnch(nxturn)
  400       continue
            do 500 i = nxres, nxres + turnsz(nxturn)
            SSbond(turn,i) = trnch
  500       continue
         end if
      end if
  600 continue
  800 continue
 1000 continue

      do 2000, nxres = 1,seqlen
      do 2000, nxbond = 1, 4
      brdg(nxbond,nxres) = 0
 2000 continue

      do 4000, nxrule = 1, nrules
      do 3000, Nxres = bsrng(nxrule,1), seqlen + bsrng(nxrule,2)
      do 2400, bndi = 1, 2
      bndptr = hbond(bndi,nxres+bsofst(nxrule,1))
      if (nxrule .le. parlel .or. bndptr .gt. nxres) then
         brdgpt = bndptr + bsofst(nxrule,2)
         If (bndptr .gt. - bsofst(nxrule,3) .and.
     +      iabs(nxres - brdgpt) .gt. bdgsep)  then
            do 2200, bndj = 1, 2
            If (hbond(bndj,bndptr+bsofst(nxrule,3)) .eq. 
     +         nxres+bsofst(nxrule,4)) then
               brpl = 1
               If (brdg(brpl,nxres) .ne. 0) brpl = 2
               Brdg(brpl,nxres) = pbrdg(nxrule) * brdgpt
               brdg(brpl+2,nxres) = exrule(1,nxrule)
               brpla = 1
               If (brdg(brpla,brdgpt) .ne. 0) brpla = 2
               Brdg(brpla,brdgpt) = pbrdg(nxrule) * nxres
               brdg(brpla+2,brdgpt) = exrule(2,nxrule)
            end if
 2200       continue
         end if
      end if
 2400 continue
 3000 continue
 4000 continue

      do 4500, nxres = 1, seqlen
      if (xney(mcangs(chiral,nxres), null)) then
         ssbond(chisgn,nxres) = '+'
         if (mcangs(chiral,nxres) .lt. 0.0) ssbond(chisgn,nxres) = '-'
      end if
 4500 continue

CHECK v.3.0.1-->
C 5000 continue
CHECK v.3.0.1<--
      do 5100, nxres = 1, seqlen
      shtcod(nxres) = 0
      strcd(1,nxres) = 0
      strcd(2,nxres) = 0
      bridpt(1,nxres) = 0
      bridpt(2,nxres) = 0
 5100 continue
      nxstrn = 0

      do 7000, nxres = 1, seqlen
      do 6900, brpl = 1, 2
      If (brdg(brpl,nxres) .ne. 0 .and. 
     +   iabs(brdg(brpl,nxres)) .gt. nxres) then
         brptr = iabs(brdg(brpl,nxres))
         brdir = brdg(brpl,nxres) / brptr
         brptpl = 2
         if (iabs(brdg(1,brptr)) .eq. nxres) brptpl = 1
         do 5600, nxptr = nxres + 1, min(seqlen,nxres+blgszl)
         do 5500, nxpl = 1,2
         If (brdg(nxpl,nxptr) * brdir .gt. 0) then
            newptr = iabs(brdg(nxpl,nxptr))
CHECK v.2.1.3-->
            IF (iabs(newptr-brptr).NE.0) THEN
CHECK v.2.1.3-->
            if ((newptr-brptr)/iabs(newptr-brptr) .eq. brdir) then
               if ((nxptr - nxres .gt. blgszs .and.
     +            iabs(newptr - brptr) .le. blgszs) .OR.
     +            (nxptr - nxres .le. blgszs .and.
     +            iabs(newptr - brptr) .le. blgszl)) then
                  nwptpl = 2
                  if (iabs(brdg(1,newptr)) .eq. nxptr) nwptpl = 1
                  do 5300, i = nxres, nxptr
                  if (ssbond(sheet,i) .eq. space)
     +               ssbond(sheet,i) = shtsym
 5300             continue
                  if (brdg(brpl+2,nxres) .lt. 0)
     +               ssbond(sheet,nxres) = shtsma
                  if (brdg(nxpl+2,nxptr) .lt. 0)
     +               ssbond(sheet,nxptr) = shtsma
                  do 5400, i = brptr, newptr, brdir
                  if (ssbond(sheet,i) .eq. space)
     +               ssbond(sheet,i) = shtsym
 5400             continue
                  if (brdg(brptpl+2,brptr) .lt. 0)
     +               ssbond(sheet,brptr) = shtsma
                  if (brdg(nwptpl+2,newptr) .lt. 0)
     +               ssbond(sheet,newptr) = shtsma
                  if (strcd(brpl,nxres) .eq. 0) then
                     nxstrn = nxstrn + 1
                     strcd(brpl,nxres) = nxstrn * brdir
                     strcd(brptpl,brptr) = nxstrn * brdir
                  end if
                  strcd(nxpl,nxptr) = strcd(brpl,nxres)
                  strcd(nwptpl,newptr) = strcd(brptpl,brptr)
                  go to 5700
               end if
            end if
CHECK v.2.1.3-->
            ENDIF
CHECK v.2.1.3-->
         end if
 5500    continue
 5600    continue
 5700    continue
         wkch = pbch
         if (brdir .lt. 0) wkch = apbch
         ssbond(bridge,nxres) = wkch
         ssbond(bridge,brptr) = wkch
      end if
 6900 continue
 7000 continue

      nxres = 1
 7500 continue
      if (ssbond(sheet,nxres) .eq. space) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 10000
         goto 7500
      end if

      stsht = nxres
      nxres = stsht + 1
 7700 continue
      if (ssbond(sheet,nxres) .ne. space) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 7800
         goto 7700
      end if
 7800 continue
      finsht = nxres - 1
      lststd = 0

 8000 continue
      call whstrd(strcd, stsht, finsht, nxstr, nxpl, lststd, bstval)
      if (bstval .ne. 0) then
         lststd = abs(strcd(nxpl,nxstr))
         chcode = mod(lststd,nstchs)
         if (chcode .eq. 0) chcode = nstchs
         if (strcd(nxpl,nxstr) .lt. 0) then
            strdch = stchap(chcode:chcode)
         else
            strdch = stchpa(chcode:chcode)
         end if
         stchpl = bridg1
         if (ssbond(bridg1,nxstr) .ne. space) stchpl = bridg2
         if (stchpl .ne. bridg2) then
            thstpl = nxstr
 8200       continue
            nxstpl = strdpl(thstpl,finsht,strcd(nxpl,nxstr),strcd)
            if (nxstpl .le. 0) goto 8300
            if (ssbond(bridg1,nxstpl) .ne. space) stchpl = bridg2
            thstpl = nxstpl
            if (stchpl .ne. bridg2) goto 8200
         end if
 8300    continue
         ssbond(stchpl,nxstr) = strdch
         bridpt(stchpl-bridg1+1,nxstr) = abs(brdg(nxpl,nxstr))
         thstpl = nxstr
 8500    continue
         nxstpl = strdpl(thstpl,finsht,strcd(nxpl,nxstr),strcd)
         if (nxstpl .le. 0) goto 8700
         do 8600, nxblg = thstpl + 1, nxstpl - 1
         ssbond(stchpl,nxblg) = bulgch
 8600    continue
         ssbond(stchpl,nxstpl) = strdch
         strpl = 1
         if (strcd(1,nxstpl) .ne. strcd(nxpl,nxstr)) strpl = 2
         bridpt(stchpl-bridg1+1,nxstpl) = abs(brdg(strpl,nxstpl))
         thstpl = nxstpl
         go to 8500
 8700    continue
         go to 8000
      end if
      nxres = finsht + 1
      if (nxres .le. seqlen) goto 7500
10000 continue
      If (lststd .gt. nstchs)
     +   Write(errfi,10100) lststd
10100 format('there are ',I3,' strands. Labels have restarted')

      Call shtlab(ssbond, bridpt, shtcod, seqlen)
      do 10500 nxres = 1, seqlen
      if (shtcod(nxres) .ne. 0) then
c        chcode = mod(shtcod(nxres),nstchs)
         chcode = shtcod(nxres)
c        if (chcode .eq. 0) chcode = nstchs
         ssbond(sheet1,nxres) = nstchp(chcode:chcode)
      end if
10500 continue

      nxbrdg = 0
      do 11000, nxres = 1, seqlen
      do 10900, nxpl = 1, 2
      If (brdg(nxpl,nxres) .ne. 0 .and. strcd(nxpl,nxres) .eq. 0
     +   .and. abs(brdg(nxpl,nxres)) .gt. nxres) then
         nxbrdg = nxbrdg + 1
         brdpt = abs(brdg(nxpl,nxres))
         chcode = mod(nxbrdg,nstchs)
         if (chcode .eq. 0) chcode = nstchs
         if (brdg(nxpl,nxres) .lt. 0) then
            brdch = stchap(chcode:chcode)
         else
            brdch = stchpa(chcode:chcode)
         end if
         nxbrpl = bridg1
         if (ssbond(bridg1,nxres) .ne. space) nxbrpl = bridg2
         ssbond(nxbrpl,nxres) = brdch
         bridpt(nxbrpl-bridg1+1,nxres) = brdpt
         nxbrpl = bridg1
         if (ssbond(bridg1,brdpt) .ne. space) nxbrpl = bridg2
         ssbond(nxbrpl,brdpt) = brdch
         bridpt(nxbrpl-bridg1+1,brdpt) = nxres
      end if
10900 continue
11000 continue
      If (nxbrdg .gt. nstchs)
     +  Write(errfi,11100) nxbrdg
11100 format('there are ',I3,' bridges. Labels have restarted')

      END
c***************************************************************************
      Integer Function strdpl(stpl, finpl, thcode, strand)
c
c     find the next occurrence of the code if any in the array
c
      include 'sstruc.par'
      integer stpl, finpl, thcode, strand(2,maxres)
      integer nxpl

      do 1000, nxpl = stpl + 1, finpl
      If (strand(1,nxpl) .eq. thcode) go to 1100
      If (strand(2,nxpl) .eq. thcode) go to 1100
 1000 continue
      nxpl = 0
 1100 continue
      strdpl = nxpl
      END
c****************************************************************************
      Subroutine Whstrd(strcd, stsht, finsht, stres, stpl, lststd,
     +                  bstval)
c
c     find the starting position of the lowest numbered strand that
c     is greater than the last
c
      Include 'sstruc.par'
      integer strcd(2,maxres), stsht, finsht, stres, stpl, lststd,
     +        bstval

      integer  nxstr, nxpl

      bstval = 0
      do 1000, nxstr = stsht, finsht
      do 1000, nxpl = 1, 2
      If (abs(strcd(nxpl,nxstr)) .gt. lststd) then
         if (bstval .eq. 0) then
            bstval = abs(strcd(nxpl,nxstr))
            stpl = nxpl
            stres = nxstr
         else
            if (abs(strcd(nxpl,nxstr)) .lt. bstval) then
               bstval = abs(strcd(nxpl,nxstr))
               stpl = nxpl
               stres = nxstr
            end if
         end if
      end if
 1000 continue
      END
c****************************************************************************
      Subroutine Shtlab(ssbond, bridpt, shtcod, seqlen)
c
c     label the sheets using numbers which are turned into
c     alphabetic symbols later
c
      include 'sstruc.par'
      character ssbond(nstruc,maxres)*1
      integer   bridpt(2,maxres), shtcod(maxres), seqlen

      integer   nxres, stsht, finsht, shtnum, nxisht, nxpl
      logical   found

      do 100, nxres = 1, seqlen
      shtcod(maxres) = 0
  100 continue

      shtnum = 0
      nxres = 1
 1000 continue
      call fdshus(nxres, ssbond, shtcod, stsht, finsht, seqlen)
      if (stsht .ne. 0) then
         shtnum = shtnum + 1
         do 1500, nxisht = stsht, finsht
         shtcod(nxisht) = shtnum
         do 1400, nxpl = 1,2
         If (bridpt(nxpl,nxisht) .ne. 0) then
            if (shtcod(bridpt(nxpl,nxisht)) .eq. 0)
     +         call setsht(bridpt(nxpl,nxisht), shtnum, shtcod, 
     +                     ssbond, seqlen)
         end if
 1400    continue
 1500    continue

 2000    continue
         call fdshps(finsht+1, shtnum, ssbond, shtcod, bridpt, found,
     +               seqlen)
         if (found) go to 2000

         nxres = finsht + 1
         if (nxres .le. seqlen) go to 1000
      end if
      If (shtnum .gt. nstchs)
     +   Write(errfi,2100) shtnum
 2100 format('there are ',I3,' sheets. Labels have restarted')
      END
c*****************************************************************************
      Subroutine fdshus(stpt, ssbond, shtcod, stsht, finsht, seqlen)
c
c     find the first unset sheet starting at stpt
c
      include 'sstruc.par'
      integer   stpt, stsht, finsht, shtcod(maxres), seqlen
      character ssbond(nstruc,maxres)*1

      integer   nxres

      stsht = 0
      nxres = stpt
 1000 continue
      if (ssbond(sheet,nxres) .eq. space .or. shtcod(nxres) .ne. 0) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) go to 2000
         go to 1000
      end if
      stsht = nxres
 1500 continue
      nxres = nxres + 1
      if (nxres .le. seqlen) then
         if (ssbond(sheet,nxres) .ne. space) goto 1500
      end if
      finsht = nxres - 1
 2000 continue
      END
c****************************************************************************
      subroutine fdshps(stpt, shtnum, ssbond, shtcod, bridpt, found,
     +                  seqlen)
c
c     find the next partially set sheet, if any (neg shtnum), set it
c     to fully set and then set its partners to partly set if not
c     already set
c
      include 'sstruc.par'
      integer   stpt, shtnum, shtcod(maxres), bridpt(2,maxres), seqlen
      character ssbond(nstruc,maxres)
      logical   found

      integer   nxres, stsht, finsht, nxpl

      found = .false.
      nxres = stpt
 1000 continue
      if (shtcod(nxres) .ne. - shtnum) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 5000
         goto 1000
      end if
      stsht = nxres
      found = .true.
 1500 continue
      if (shtcod(nxres) .eq. -shtnum) then
         nxres = nxres + 1
         if (nxres .le. seqlen) goto 1500
      end if
      finsht = nxres - 1

      do 2000, nxres = stsht, finsht
      shtcod(nxres) = shtnum
      do 1900, nxpl = 1,2
      if (bridpt(nxpl,nxres) .ne. 0) then
         if (shtcod(bridpt(nxpl,nxres)) .eq. 0)
     +      call setsht(bridpt(nxpl,nxres), shtnum, shtcod, ssbond, 
     +                  seqlen)
      end if
 1900 continue
 2000 continue
 5000 continue
      END
c*****************************************************************************
      Subroutine Setsht(refpl, shtnum, shtcod, ssbond, seqlen)
c
c     finds the limits of the super sheet containing refpl and
c     sets this to partially set
c
      include 'sstruc.par'
      integer   refpl, shtnum, shtcod(maxres), seqlen
      character ssbond(nstruc,maxres)*1

      integer   nxres, stsht, finsht

      nxres = refpl
  200 continue
      if (ssbond(sheet,nxres) .ne. space) then
         nxres = nxres + 1
         if (nxres .le. seqlen) go to 200
      end if
      finsht = nxres - 1
      nxres = refpl
  400 continue
      if (ssbond(sheet,nxres) .ne. space) then
         nxres = nxres - 1
         if (nxres .ge. 1) go to 400
      end if
      stsht = nxres + 1

      do 800, nxres = stsht, finsht
      shtcod(nxres) = -shtnum
  800 continue
      END
c******************************************************************************
      Subroutine MkBend(mcangs, ssbond, seqlen)
c
c     routine to determine which main chain residues are bent at over
c     the requisite angle
c
      INCLUDE 'sstruc.par'

      Real      mcangs(nmnang,maxres)
      Character Ssbond(nstruc,maxres)*1
      Integer   seqlen
c
      real      bendsz
      parameter (bendsz = 70.0)
      Integer   nxres
      logical   xney

      do 1000, nxres = 1, seqlen
      if (mcangs(kappa,nxres) .gt. bendsz .and.
     +   xney(mcangs(kappa,nxres), null)) ssbond(bend,nxres) = bendch
 1000 continue
      END
c******************************************************************************
      Subroutine MkSumm(ssbond, sssumm, summpl, seqlen)
c
c     generate a structure summary using the priority rules
c
      include 'sstruc.par'
CHECK v.3.2-->
      SAVE stsym, altsym
CHECK v.3.2-->
      character*1 ssbond(nstruc,maxres), sssumm(maxres), summpl(maxres)
      integer     seqlen
CHECK v.3.0.1--> (Order of statements changed)
      Integer     nxres
CHECK v.3.0.1<--

      character   stsym(nstruc)*1, altsym(nstruc)*1
      data stsym /'H', ' ', 'E', 'B', ' ', ' ', 'G', 'I', 'T', 'S', '+'/
      data altsym/'h', ' ', 'e', 'b', ' ', ' ', 'g', 'i', 't', 's', '+'/
CHECK v.3.2-->
C      SAVE stsym, altsym
CHECK v.3.2-->


      do 200, nxres = 1, seqlen
      SSsumm(nxres) = space
      summpl(nxres) = space
  200 continue
      call hlxset(ssbond, sssumm, summpl, alphah, ahsz, stsym(alphah),
     +            altsym(alphah), seqlen)
      do 400, nxres = 1, seqlen
      If (ssbond(sheet,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) = stsym(sheet)
         if (summpl(nxres) .eq. space .or. summpl(nxres) .eq.
     +      altsym(alphah) .or. summpl(nxres) .eq. altsym(sheet))
     +      summpl(nxres) = stsym(sheet)
         if (ssbond(sheet,nxres) .eq. shtsma) then
            if (summpl(nxres-1) .eq. space .or. summpl(nxres-1) .eq.
     +         stsym(bridge)) summpl(nxres-1) = altsym(sheet)
            if (summpl(nxres+1) .eq. space .or. summpl(nxres-1) .eq.
     +         stsym(bridge)) summpl(nxres+1) = altsym(sheet)
         end if
      end if
      If (ssbond(bridge,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) = stsym(bridge)
         if (summpl(nxres) .eq. space) summpl(nxres) = stsym(bridge)
      end if
  400 continue
      call hlxset(ssbond, sssumm, summpl, thrtnh, thrtsz, stsym(thrtnh),
     +            altsym(thrtnh), seqlen)
      call hlxset(ssbond, sssumm, summpl, pih, pihsz, stsym(pih),
     +            altsym(pih), seqlen)
      call hlxres(sssumm, stsym(thrtnh), altsym(thrtnh), stsym(turn),
     +            altsym(turn), thrtsz, seqlen)
      call hlxres(summpl, stsym(thrtnh), altsym(thrtnh), stsym(turn),
     +            altsym(turn), thrtsz, seqlen)
      call hlxres(sssumm, stsym(pih), altsym(pih), stsym(turn),
     +            altsym(turn), pihsz, seqlen)
      call hlxres(summpl, stsym(pih), altsym(pih), stsym(turn),
     +            altsym(turn), pihsz, seqlen)
      call turnst(ssbond, sssumm, summpl, stsym(turn), altsym(turn),
     +            seqlen)
      call bendst(ssbond, sssumm, summpl, seqlen)

      END
c*****************************************************************************
      Subroutine Hlxset(ssbond, sssumm, summpl, hlxpos, hlxsz, hlxch,
     +                  altch,  seqlen)
c
c     go through the bond tables and note any helices into the
c     summary line
c
      include 'sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), hlxch, altch,
     +            summpl(maxres)
      integer     hlxpos, hlxsz, seqlen

      integer     nxres, nxsym

      do 1000, nxres = 2, seqlen - hlxsz
      if (ssbond(hlxpos,nxres) .eq. bndbeg .or.
     +   ssbond(hlxpos,nxres) .eq. bndbth) then
         if (ssbond(hlxpos,nxres-1) .eq. bndbeg .or.
     +      ssbond(hlxpos,nxres-1) .eq. bndbth) then
            do 500, nxsym = 0, hlxsz - 1
            if (sssumm(nxres+nxsym) .eq. space)
     +         sssumm(nxres+nxsym) = hlxch
            if (summpl(nxres+nxsym) .eq. space .OR. 
     +         summpl(nxres+nxsym) .eq. altch)
     +         summpl(nxres+nxsym) = hlxch
  500       continue
            if (summpl(nxres-1) .eq. space) summpl(nxres-1) = altch
            if (summpl(nxres+hlxsz) .eq. space)
     +             summpl(nxres+hlxsz) = altch
         end if
      end if
 1000 continue
      END
c*****************************************************************************
      Subroutine hlxres(summ, hlxch, altch, turnch, alturn, hlxsz,
     +                  seqlen)
c
c     to remove helix character runs of less than minimal helix size
c
      include 'sstruc.par'
      character*1 summ(maxres), hlxch, turnch, altch, alturn
      integer     hlxsz, seqlen

      integer     nxres, stpos, nxsym

      stpos = 0
      do 1000, nxres = 1, seqlen
      if (summ(nxres) .eq. hlxch .or. summ(nxres) .eq. altch) then
         If (stpos .eq. 0) stpos = nxres
      else
         if (stpos .ne. 0) then
            if (nxres - stpos .lt. hlxsz) then
               do 500, nxsym = stpos, nxres - 1
               if (summ(nxsym) .eq. hlxch) then
                  summ(nxsym) = turnch
               else
                  summ(nxsym) = alturn
               end if
  500          continue
            end if
            stpos = 0
         end if
      end if
 1000 continue
      if (stpos .ne. 0) then
         if (seqlen - stpos + 1 .lt. hlxsz) then
            do 1500, nxsym = stpos, seqlen
            if (summ(nxres) .eq. hlxch) then
               summ(nxsym) = turnch
            else
               summ(nxsym) = alturn
            end if
 1500       continue
         end if
      end if
      END
c**************************************************************************
      Subroutine turnst(ssbond, sssumm, summpl, turnch, alturn, seqlen)
c
c     search for solitary turn type bonds and set in the symbol
c
      include 'sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), turnch,
     +            summpl(maxres), alturn
      integer     seqlen

      integer     nxres, nxsym, nxhlx, hlxpri(nhelix), hlxsz(nhelix)

      hlxpri(1) = alphah
      hlxpri(2) = thrtnh
      hlxpri(3) = pih
      hlxsz(1) = ahsz
      hlxsz(2) = thrtsz
      hlxsz(3) = pihsz

      do 2000, nxhlx = 1, nhelix
      If (ssbond(hlxpri(nxhlx),1) .eq. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),2) .ne. bndbeg) then
         do 200, nxsym = 2, hlxsz(nxhlx)
         If (sssumm(nxsym) .eq. space) sssumm(nxsym) = turnch
         If (summpl(nxsym) .eq. space .or. summpl(nxsym) .eq.
     +      alturn) summpl(nxsym) = turnch
  200    continue
         If (summpl(1) .eq. space) summpl(1) = alturn
         If (summpl(hlxsz(nxhlx)+1) .eq. space)
     +      summpl(hlxsz(nxhlx)+1) = alturn
      end if
      do 1000, nxres = 2, seqlen - hlxsz(nxhlx)
      If ((ssbond(hlxpri(nxhlx),nxres) .eq. bndbeg .or.
     +   ssbond(hlxpri(nxhlx),nxres) .eq. bndbth) .AND.
     +   (ssbond(hlxpri(nxhlx),nxres-1) .ne. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),nxres-1) .ne. bndbth) .AND.
     +   (ssbond(hlxpri(nxhlx),nxres+1) .ne. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),nxres+1) .ne. bndbth)) then
         do 500, nxsym = nxres + 1, Nxres + hlxsz(nxhlx) - 1
         If (sssumm(nxsym) .eq. space) sssumm(nxsym) = turnch
         If (summpl(nxsym) .eq. space .or. summpl(nxsym) .eq.
     +      alturn) summpl(nxsym) = turnch
  500    continue
         If (summpl(nxres) .eq. space) summpl(nxres) = alturn
         If (summpl(nxres+hlxsz(nxhlx)) .eq. space)
     +      summpl(nxres+hlxsz(nxhlx)) = alturn
      end if
 1000 continue
 2000 continue
      END
c*****************************************************************************
      Subroutine Bendst(ssbond, sssumm, summpl, seqlen)
c
c     enter bend symbols into summary line
c
      Include 'sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), summpl(maxres)
      integer     seqlen
      Integer     nxres
c
      do 1000, nxres = 1, seqlen
      If (ssbond(bend,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) =
     +      ssbond(bend,nxres)
         If (summpl(nxres) .eq. space) summpl(nxres) =
     +      ssbond(bend,nxres)
      end if
 1000 continue
      END
c****************************************************************************
      Subroutine mkauth(seqbcd, autsum, ssord, astruc, seqlen)
c
c     reread the author structure records and create the structure
c
      include 'sstruc.par'
CHECK v.3.2-->
      save hlxch
CHECK v.3.2<--
      character seqbcd(maxres)*6, autsum(maxres)*1, ssord*21
      logical   astruc
      integer   seqlen
      integer   respl

      integer   nxsym, stpos, finpos, lstpos, hlxno, hlxcls, nxres,
     +          strno, sense, strpri, hlxcds
      parameter (hlxcds = 10)
      character stres*6, finres*6, autrec*80, strkch*1, shtid*3,
     +          hlxch(hlxcds)*1
      data hlxch/'H','J','I','K','G','L','M','N','O','P'/
CHECK v.3.2-->
C      save hlxch
CHECK v.3.2<--

      do 100, nxres = 1, seqlen
      autsum(nxres) = space
  100 continue

      if (astruc) then
         Rewind(unit=authfi)
         lstpos = 0
  200    continue
         read(authfi,'(a)',end=2000) autrec
         if (autrec(1:keylen) .eq. hlxkey) then
            read(autrec,400) hlxno, stres(1:1), stres(2:6), finres(1:1),
     +          finres(2:6), hlxcls
CHECK v.2.1.4-->
C  400       format(7x,i3,9x,a1,x,a5,5x,a1,x,a5,i2)
  400       format(7x,i3,9x,a1,1x,a5,5x,a1,1x,a5,i2)
CHECK v.2.1.4<--
            strkch = hlxch(1)
            if (hlxcls .le. hlxcds .and. hlxcls .gt. 0) 
     +         strkch = hlxch(hlxcls)
         else
            if (autrec(1:keylen) .eq. shtkey) then
               read(autrec,600) strno, shtid, stres, finres, sense
CHECK v.2.1.4-->
C  600          format(7x,i3,x,a3,7x,a6,5x,a6,i2)
  600          format(7x,i3,1x,a3,7x,a6,5x,a6,i2)
CHECK v.2.1.4<--
               strkch = shtsym
            else
               if (autrec(1:keylen) .eq. trnkey) then
                  read(autrec,800) stres, finres
  800             format(19x,a6,5x,a6)
                  strkch = trnch
               else
                  go to 200
               end if
            end if
         end if
         stpos = respl(seqbcd, stres, lstpos, seqlen)
         finpos = respl(seqbcd, finres, lstpos, seqlen)
         If (finpos .lt. 0 .or. stpos .lt. 0 .or.
     +      finpos .lt. stpos .or. (strkch .eq. hlxch(1) .and. 
     +      (hlxcls .le. 0 .or. hlxcls .gt. hlxcds))) then
            Write(errfi,*)'structure record in error'
            Write(errfi,*)autrec
         else
            strpri = index(ssord,strkch)
            do 1000, nxsym = stpos, finpos
            if (index(ssord,autsum(nxsym)) .gt. strpri) 
     +         autsum(nxsym) = strkch
 1000       continue
         end if
         go to 200
 2000    continue
      end if
      end
c*****************************************************************************
      Subroutine disulf(seqbcd, dsdsts, scangs, scaaat, scatin, atomin,
     +                  mcaaat, seqcod, dsdef, dspart,  astruc, seqlen)
c
c     read the author file and extract any author definitions of
c     disulphide bridges. check for possible bridges on a distance
c     criterion as well.
c
      include 'sstruc.par'
      character seqbcd(maxres)*6, dsdef(maxres)*1
      real      scangs(nsdang,maxres), dsdsts(3,maxres),
     +          scaaat(ncoord,sdchna,maxres), 
     +          mcaaat(ncoord,mnchna,maxres)
      integer   dspart(maxres), seqcod(maxres), seqlen
      logical   astruc, scatin(sdchna,maxres), atomin(mnchna,maxres)
      real      atdist
      integer   respl

      character autdef*1, disdef*1
      parameter (autdef = 'A', disdef = 'D')
      integer   nxdist, nxres, pos1, pos2, errset, partpl, nxcys
      real      partds, ssdist
      logical   cystyp, xney
      character autrec*80, code1*6, code2*6, typefl*1

      do 100, nxres = 1, seqlen
      dsdef(nxres) = space
      dspart(nxres) = 0
      do 100, nxdist = 1, 3
      dsdsts(nxdist,nxres) = 0.0
  100 continue

      If (astruc) then
         rewind(unit=authfi)
         pos1 = 0
  200    continue
         if (pos1 .lt. 0) pos1 = 0
         read (authfi,'(A)', end = 2000) autrec
         if (autrec(1:keylen) .ne. dsfkey) go to 200
         read(autrec,400) code1(1:1), code1(2:6), code2(1:1), code2(2:6)
CHECK v.2.1.4-->
C  400    format(15x,a1,x,a5,7x,a1,x,a5)
  400    format(15x,a1,1x,a5,7x,a1,1x,a5)
CHECK v.2.1.4<--
         pos1 = respl(seqbcd, code1, pos1, seqlen)
         errset = 0
         if (pos1 .gt. 0) pos2 = respl(seqbcd, code2, pos1, seqlen)
         If (pos1 .gt. 0 .and. pos2 .gt. 0) then
            if (Cystyp(seqcod(pos1)) .and. cystyp(seqcod(pos2))) then
               typefl = autdef
               call dsset(pos1, pos2, dsdsts, scangs, scaaat, scatin,
CHECK v.3.2-->
C     +                    atomin, mcaaat, dsdef, dspart, typefl, seqlen)
     +                    atomin, mcaaat, dsdef, dspart, typefl)
CHECK v.3.2<--
            else
               errset = 1
            end if
         else
            errset = 1
         end if
         If (errset .eq. 1) then
            write(errfi,'(A)') 'error in SSBOND record'
            write(errfi,'(A)') autrec
         end if
         go to 200

 2000    continue
         close(unit = authfi)
      end if

      nxres = 1
 2500 continue
      if (cystyp(seqcod(nxres))) then
         if (scatin(sgamma,nxres)) then
            partpl = 0
            partds = 0.0
            do 3000, nxcys = nxres + 1, seqlen
            if (cystyp(seqcod(nxcys))) then
               if (scatin(sgamma,nxcys)) then
                  ssdist = atdist(scaaat(1,sgamma,nxres),
     +                     scaaat(1,sgamma,nxcys))
                  if (ssdist .lt. sssep) then
                     if (xney(partds, 0.0)) then
                        if (partds .lt. ssdist) then
                           Write(errfi,2700) nxres, nxcys, ssdist
                        else
                           Write(errfi,2700) nxres, partpl, partds
 2700                      format('multifurcation of disulphides ',
CHECK v.3.2-->
C     +                            I4, x, i4, x, f3.1)
     +                            I4, 1x, i4, 1x, f3.1)
CHECK v.3.2<--
                           partpl = nxcys
                           partds = ssdist
                        end if
                     else
                        partpl = nxcys
                        partds = ssdist
                     end if
                  end if
               end if
            end if
 3000       continue
            if (partpl .ne. 0) then
               typefl = disdef
               call dsset(nxres, partpl, dsdsts, scangs, scaaat, scatin,
CHECK v.3.2-->
C     +                    atomin, mcaaat, dsdef, dspart, typefl, seqlen)
     +                    atomin, mcaaat, dsdef, dspart, typefl)
CHECK v.3.2<--
            end if
         end if
      end if
      nxres = nxres + 1
      if (nxres .le. seqlen) go to 2500
      END
c****************************************************************************
      logical function cystyp(seqcod)
c
c     check that we do have a potential disulfide former
c
      include 'sstruc.par'
      integer  seqcod

      If (seqcod .eq. cyscod .or. seqcod .eq. mprcod) then
         cystyp = .true.
      else
         cystyp = .false.
      end if
      end
c*****************************************************************************
      Subroutine dsset(part1,  part2,  dsdsts, scangs, scaaat, scatin,
CHECK v.3.2-->
C     +                 atomin, mcaaat, dsdef,  dspart, typefl, seqlen)
     +                 atomin, mcaaat, dsdef,  dspart, typefl)
CHECK v.3.2<--
c
c     calculate all the features for a disulphide bridge. that is 3 chi
c     angles, calpha, cbeta, and ss distances. set in the partner and
c     definition fields. Check for conflict with author defined records
c     or previously defined bridges
c
      include 'sstruc.par'
CHECK v.3.2-->
C      integer   part1, part2, dspart(maxres), seqlen
      integer   part1, part2, dspart(maxres)
CHECK v.3.2<--
      real      dsdsts(3,maxres), scangs(nsdang,maxres), 
     +          scaaat(ncoord,sdchna,maxres), 
     +          mcaaat(ncoord,mnchna,maxres)
      logical   scatin(sdchna,maxres), atomin(mnchna,maxres)
      character dsdef(maxres)*1, typefl*1
      real      atdist, dihed
  
      character*1 botdef, autdef
      real        mxdis
      parameter  (botdef = 'B', autdef = 'A', mxdis = 9.9)
      integer     nxdis

      If ((dspart(part1) .ne. 0 .and. dspart(part1) .ne. part2) .or.
     +   (dspart(part1) .eq. 0 .and. dspart(part2) .ne. 0)) then
         if (dsdef(part1) .ne. space .or. dsdef(part2) .ne. space) then
            write(errfi,100) part1, part2
  100       format('conflict with author defined disulphide -',
CHECK v.2.1.4-->
C     +             ' abandonning ',I4,x,I4)
     +             ' abandoning ',I4,1x,I4)
CHECK v.2.1.4<--
         else
            write(errfi,200) part1, part2
CHECK v.3.4-->
C  200       format('conflict with previously defined disulphide -'
  200       format('conflict with previously defined disulphide -',
CHECK v.3.4<--
CHECK v.2.1.4-->
C     +             ' abandoning ',I4,x,I4)
     +             ' abandoning ',I4,1x,I4)
CHECK v.2.1.4<--
         end if
      else
         if (dsdef(part1) .eq. autdef) then
            dsdef(part1) = botdef
            dsdef(part2) = botdef
         else
            dspart(part1) = part2
            dspart(part2) = part1
            dsdef(part1) = typefl
            dsdef(part2) = typefl
            if (atomin(calph,part1) .and. scatin(cbeta,part1) .and. 
     +         scatin(sgamma,part1) .and. scatin(sgamma,part2))
     +         scangs(2,part1) = dihed(nmnang+2,
     +         mcaaat(1,calph,part1), scaaat(1,cbeta,part1),
     +         scaaat(1,sgamma,part1), scaaat(1,sgamma,part2))
            if (scatin(cbeta,part1) .and. scatin(sgamma,part1) .and. 
     +         scatin(sgamma,part2) .and. scatin(cbeta,part2)) then
               scangs(3,part1) = dihed(nmnang+3,
     +         scaaat(1,cbeta,part1), scaaat(1,sgamma,part1), 
     +         scaaat(1,sgamma,part2), scaaat(1,cbeta,part2))
               scangs(3,part2) = scangs(3,part1)
            end if
            if (atomin(calph,part2) .and. scatin(cbeta,part2) .and. 
     +         scatin(sgamma,part2) .and. scatin(sgamma,part1))
     +         scangs(2,part2) = dihed(nmnang+2,
     +         mcaaat(1,calph,part2), scaaat(1,cbeta,part2),
     +         scaaat(1,sgamma,part2), scaaat(1,sgamma,part1))
            do 1000, nxdis = 1, 3
            dsdsts(nxdis,part1) = -1.0
 1000       continue
            If (atomin(calph,part1) .and. atomin(calph,part2))
     +         dsdsts(1,part1) = atdist(mcaaat(1,calph,part1),
     +         mcaaat(1,calph,part2))
            If (scatin(cbeta,part1) .and. atomin(cbeta,part2))
     +         dsdsts(2,part1) = atdist(scaaat(1,cbeta,part1),
     +         scaaat(1,cbeta,part2))
            If (scatin(sgamma,part1) .and. scatin(sgamma,part2))
     +         dsdsts(3,part1) = atdist(scaaat(1,sgamma,part1),
     +         scaaat(1,sgamma,part2))
            do 1200, nxdis = 1, 3
            if (dsdsts(nxdis,part1) .gt. mxdis) then
               write(errfi,1100) part1, part2, dsdsts(nxdis,part1)
 1100          format('disulphide distance too large for',
CHECK v.2.1.4-->
C     +                I4,x,I4,x,f5.1)
     +                I4,1x,I4,1x,f5.1)
CHECK v.2.1.4<--
               dsdsts(nxdis,part1) = mxdis
            end if
 1200       continue
            do 1500, nxdis = 1,3
            dsdsts(nxdis,part2) = dsdsts(nxdis,part1)
 1500       continue
         end if
      end if
      END
c*****************************************************************************
      Integer function respl(seqbcd, theres, lstpos, seqlen)
c
c     search through the residue identifiers to find the position
c     of the relevant code. Codes will normally be in order so the
c     search starts from the last position.
c
      include 'sstruc.par'
      character seqbcd(maxres)*6, theres*6
      integer   lstpos, seqlen

      integer   nxres

      do 100, nxres = lstpos + 1, seqlen
      if (seqbcd(nxres) .eq. theres) go to 300
  100 continue
      do 200, nxres = 1, lstpos
      if (seqbcd(nxres) .eq. theres) go to 300
  200 continue
      nxres = -1
  300 continue
      respl = nxres
      end
c*****************************************************************************
CHECK v.3.2-->
C      Subroutine result(hbonde, hbond,  aacod3, aacod1, seqcod, seqbcd, 
C     +                  brksym, ssbond, bridpt, sssumm, summpl, autsum, 
C     +                  mcangs, altcod, aastd,  atpres, thmprs, scangs,
C     +                  cadsts, brknam, headln, hedlen, outsty, dsdsts,
      Subroutine result(hbonde, aacod3, aacod1, seqcod, seqbcd, 
     +                  brksym, ssbond, sssumm, summpl, autsum, 
     +                  mcangs, scangs,
     +                  outsty, dsdsts,
CHECK v.3.2<--

CHECK v.1.0-->
C     +                  dsdef,  dspart, hetcod, oois,   seqlen)
CHECK v.2.0-->
C     +                  dsdef,  dspart, hetcod, oois, seqlen, BVALUE)
CHECK v.3.2-->
C     +                  dsdef,  dspart, hetcod, oois, seqlen, BVALUE,
C     -                  MCBVAL,SCBVAL,NMCATM,NSCATM)
     +                  hetcod, oois, seqlen, BVALUE,
CHECK v.3.5-->
C     -                  MCBVAL,SCBVAL,NMCATM,NSCATM,MCBSTD,SCBSTD)
     -                  MCBVAL,SCBVAL,NMCATM,NSCATM,MCBSTD,SCBSTD,
     -                  SAVNAM)
CHECK v.3.5<--
CHECK v.3.2<--
CHECK v.2.0<--
CHECK v.1.0<--

c
c     Results are in 4 parts. The standard output set of details
c     one line to an amino acid. The second set makes the full set
c     again one line per amino acid. The third set is the structure
c     patterns and summaries at 100 amino acids to the line
c     with finally error and information messages
c
      include 'sstruc.par'

CHECK v.3.2-->
      save bpl, strpts, lefthd, rthd
CHECK v.3.2<--

CHECK v.3.5-->
      CHARACTER*3   SAVNAM(MAXRES)
CHECK v.3.5<--

CHECK v.1.0-->
      REAL          BVALUE(MAXRES)
CHECK v.1.0<--

CHECK v.2.0-->
      INTEGER       NMCATM(MAXRES), NSCATM(MAXRES)
      REAL          MCBVAL(MAXRES), SCBVAL(MAXRES)
CHECK v.2.0<--

CHECK v.3.2-->
      REAL          MCBSTD(MAXRES), SCBSTD(MAXRES)
CHECK v.3.2<--

CHECK v.3.2-->
C      integer   hbond(maxbnd,maxres), seqcod(maxres), seqlen,
C     +          bridpt(2,maxres), hedlen, dspart(maxres), oois(2,maxres)
      integer   seqcod(maxres), seqlen,
     +          oois(2,maxres)
CHECK v.3.2<--
CHECK v.3.6-->
      INTEGER   RESNO
CHECK v.3.6<--
      real      hbonde(maxbnd,maxres), mcangs(nmnang,maxres),
CHECK v.3.2-->
C     +          scangs(nsdang,maxres), cadsts(11,maxres), 
     +          scangs(nsdang,maxres),
CHECK v.3.2<--
     +          dsdsts(3,maxres)
      character aacod3*(abetsz*3), ssbond(nstruc,maxres)*1,
CHECK v.3.2-->
C     +          sssumm(maxres)*1, brknam*(fnamln), summpl(maxres),
C     +          seqbcd(maxres)*6, aacod1*(abetsz), autsum(maxres)*1,
C     +          brksym(maxres)*1, headln*132, altcod(maxres)*1,
C     +          aastd(maxres)*1, atpres(3,maxres), thmprs(maxres),
C     +          outsty*1, dsdef(maxres)*1, hetcod*63
     +          sssumm(maxres)*1, summpl(maxres),
     +          seqbcd(maxres)*6, aacod1*(abetsz), autsum(maxres)*1,
     +          brksym(maxres)*1,
     +          outsty*1, hetcod*63
CHECK v.3.2<--

      integer   nxres, i, bpl(maxbnd), staa, nxstr, nxang, nxst, nxrow,
CHECK v.2.1-->
C     +          aathln, strpts(9), j, hednum
     +          aathln, strpts(9), hednum
CHECK v.2.1<--
      character aa3*3, aa1*1, lefthd(13)*9, rthd(13)*9,
     +          wkln1*(aaonln), wkln2*122, sumln(3)*(aaonln)
      data bpl /3,1,4,2/
      data strpts/chisgn, bend, turn, pih, thrtnh, bridg2, bridg1,
     +            sheet1, alphah/
      data lefthd /'chirality', '    bends', '    turns', '  5-turns',
     +             '  3-turns', ' bridge-2', ' bridge-1', '   sheets',
     +             '  4-turns', '   author', 'Kabs/Sand', '  summary',
     +             ' sequence'/
      data rthd   /'chirality', 'bends    ', 'turns    ', '5-turns  ',
     +             '3-turns  ', 'bridge-2 ', 'bridge-1 ', 'sheets   ',
     +             '4-turns  ', 'author   ', 'Kabs/Sand', 'summary  ',
     +             'sequence '/
CHECK v.3.2-->
C      save bpl, strpts, lefthd, rthd
CHECK v.3.2<--

      hednum = 1
CHECK v.3.2-->
C      call hedset(brknam, headln, hedlen, hednum, seqlen)
CHECK v.3.2<--
      do 500, nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         staa = (seqcod(nxres) - 1) * 3 + 1
         aa3 = aacod3(staa:staa+2)
         aa1 = aacod1(seqcod(nxres):seqcod(nxres))
      else
         staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
         aa3 = hetcod(staa:staa+2)
         aa1 = '-'
      end if
CHECK v.3.5-->
         aa3 = SAVNAM(NXRES)
CHECK v.3.5<--

CHECK v.1.0-->
C      write(outfi,300) nxres, brksym(nxres), seqbcd(nxres),
C     +     altcod(nxres), aa3, aa1, autsum(nxres), sssumm(nxres),
C     +     summpl(nxres), (ssbond(nxstr,nxres), nxstr = alphah, sheet1),
C     +     (ssbond(nxstr,nxres), nxstr = bridg1, chisgn),
C     +     (bridpt(nxstr,nxres), nxstr = 1, 2),
C     +     (mcangs(nxang,nxres), nxang = 1, nmnang),
C     +     (hbond(bpl(i),nxres), hbonde(bpl(i),nxres), i=1,4),
C     +     oois(ooi1,nxres), oois(ooi2,nxres)
C  300 format(1x,i4,a1,a6,x,a1,x,a3,x,a1,3(x,a1),x,9a1,2(x,i4),
C     +       4(x,f6.1),2(x,f5.1),4(x,i4,x,f4.1),2(x,i2))

C---- Write out to output .rin file
CHECK v.2.0-->
C      WRITE(8,300) NXRES, AA3, BRKSYM(NXRES), SEQBCD(NXRES),
C     -    SUMMPL(NXRES), (MCANGS(NXANG,NXRES), NXANG = 1, 3),
C     -    (SCANGS(I,NXRES), I = 1, 3), HBONDE(BPL(1),NXRES),
C     -    DSDSTS(3,NXRES), MCANGS(IMPLD,NXRES), BVALUE(NXRES)
C 300  FORMAT(I4,A3,A1,A6,A1,10F7.2)

C---- Calculate the mean and standard deviation for the main-chain
C     and side-chain B-values
      IF (NMCATM(NXRES).NE.0) THEN
          MCBVAL(NXRES) = MCBVAL(NXRES) / NMCATM(NXRES)
CHECK v.3.2-->
          MCBSTD(NXRES) = MCBSTD(NXRES) / NMCATM(NXRES)
     -        - MCBVAL(NXRES) * MCBVAL(NXRES)
          IF (NMCATM(NXRES).GT.1 .AND. MCBSTD(NXRES).GT.0.0) THEN
              MCBSTD(NXRES) = SQRT(MCBSTD(NXRES))
          ELSE
              MCBSTD(NXRES) = 0.0
          ENDIF
CHECK v.3.2<--
      ELSE
          MCBVAL(NXRES) = 0.0
CHECK v.3.2-->
          MCBSTD(NXRES) = 0.0
CHECK v.3.2<--
      ENDIF
      IF (NSCATM(NXRES).NE.0) THEN
          SCBVAL(NXRES) = SCBVAL(NXRES) / NSCATM(NXRES)
CHECK v.3.2-->
          SCBSTD(NXRES) = SCBSTD(NXRES) / NSCATM(NXRES)
     -        - SCBVAL(NXRES) * SCBVAL(NXRES)
          IF (NSCATM(NXRES).GT.1 .AND. SCBSTD(NXRES).GT.0.0) THEN
              SCBSTD(NXRES) = SQRT(SCBSTD(NXRES))
          ELSE
              SCBSTD(NXRES) = 0.0
          ENDIF
CHECK v.3.2<--
      ELSE
          SCBVAL(NXRES) = 0.0
CHECK v.3.2-->
          SCBSTD(NXRES) = 0.0
CHECK v.3.2<--
      ENDIF

C---- Write out to output .rin file
CHECK v.3.6.1-->
C     WRITE(8,300) NXRES, AA3, BRKSYM(NXRES), SEQBCD(NXRES),
      RESNO = NXRES
      IF (RESNO.GT.9999) RESNO = 9999
      WRITE(8,300) RESNO, AA3, BRKSYM(NXRES), SEQBCD(NXRES),
CHECK v.3.6.1<--
     -    SUMMPL(NXRES), (MCANGS(NXANG,NXRES), NXANG = 1, 3),
CHECK v.3.0-->
C     -    (SCANGS(I,NXRES), I = 1, 3), HBONDE(BPL(1),NXRES),
C     -    DSDSTS(3,NXRES), MCANGS(IMPLD,NXRES), BVALUE(NXRES),
C     -    MCBVAL(NXRES), SCBVAL(NXRES)
C 300  FORMAT(I4,A3,A1,A6,A1,10F7.2,2F7.3)
     -    (SCANGS(I,NXRES), I = 1, 4), HBONDE(BPL(1),NXRES),
     -    DSDSTS(3,NXRES), MCANGS(IMPLD,NXRES), BVALUE(NXRES),
     -    MCBVAL(NXRES), SCBVAL(NXRES), OOIS(OOI1,NXRES),
CHECK v.3.2-->
C     -    OOIS(OOI2,NXRES)
C 300  FORMAT(I4,A3,A1,A6,A1,11F7.2,2F7.3,2(1X,I2))
     -    OOIS(OOI2,NXRES), MCBSTD(NXRES), SCBSTD(NXRES)
CHECK v.3.5.1-->
C 300  FORMAT(I4,A3,A1,A6,A1,11F7.2,2F7.3,2(1X,I2),2F7.3)
 300  FORMAT(I4,A3,A1,A6,A1,11F7.2,2F7.3,2(I3),2F7.3)
CHECK v.3.5.1<--
CHECK v.3.2<--
CHECK v.3.0<--
CHECK v.2.0<--
CHECK v.1.0<--

  500 continue

CHECK v.1.0-->
C      write(outfi,*)char(12)
CHECK v.1.0<--

      if (outsty .eq. 'F') then
         hednum = 2
CHECK v.3.2-->
C         call hedset(brknam, headln, hedlen, hednum, seqlen)
CHECK v.3.2<--
         do 1000, nxres = 1, seqlen
         if (seqcod(nxres) .le. abetsz) then
            staa = (seqcod(nxres) - 1) * 3 + 1
            aa3 = aacod3(staa:staa+2)
            aa1 = aacod1(seqcod(nxres):seqcod(nxres))
         else
            staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
            aa3 = hetcod(staa:staa+2)
            aa1 = '-'
         end if
CHECK v.3.0.1-->
C  800    continue
CHECK v.3.0.1<--

CHECK v.1.0-->
C         write(outfi,900)nxres, brksym(nxres), seqbcd(nxres), aa3, aa1,
C     +      aastd(nxres), (atpres(i,nxres),i=1,3),
C     +      thmprs(nxres), (scangs(i,nxres),i=1,4), dspart(nxres), 
C     +      dsdef(nxres), (dsdsts(i,nxres),i=1,3), 
C     +      (cadsts(i,nxres),i=1,11)
C  900    format(x,i4,a1,a6,x,a3,x,a1,5(x,a1),4(x,f6.1),x,i4,x,a1,
C     +          3(x,f4.1),11(x,f4.1))
CHECK v.1.0<--

 1000    continue

CHECK v.1.0-->
C         write(outfi,*) char(12)
CHECK v.1.0<--

      end if

CHECK v.1.0-->
C      write(outfi,1001) brknam
C 1001 format(x,a)
C      write(outfi,1001) headln(1:hedlen)
CHECK v.1.0<--

      do 2000, nxst = 1, seqlen, aaonln

CHECK v.1.0-->
C      write(outfi,*)
CHECK v.1.0<--

      aathln = aaonln
      if (nxst + aathln - 1 .gt. seqlen) aathln = seqlen - nxst + 1
      do 1100, nxstr = 1,aathln
      sumln(1)(nxstr:nxstr) = autsum(nxst+nxstr-1)
      sumln(2)(nxstr:nxstr) = sssumm(nxst+nxstr-1)
 1100 continue
      do 1300, nxrow = 1, 2
      wkln2 = lefthd(9+nxrow) // '  ' // sumln(nxrow)(1:aathln) //
     +        '  ' // rthd(9+nxrow)

CHECK v.1.0-->
C      write(outfi,1200)wkln2(1:aathln+22)
CHECK v.1.0<--

CHECK v.3.0.1-->
C 1200 format(x,A)
CHECK v.3.0.1<--
 1300 continue
      do 1500, nxrow = 1, 9
      do 1400, nxstr = 1, aathln
      wkln1(nxstr:nxstr) = ssbond(strpts(nxrow),nxst+nxstr-1)
 1400 continue
      wkln2 = lefthd(nxrow) // '  ' // wkln1(1:aathln) // '  '
     +        // rthd(nxrow)

CHECK v.1.0-->
C      write(outfi,1200)wkln2(1:aathln+22)
CHECK v.1.0<--

 1500 continue
      do 1600, nxstr = 1,aathln
      sumln(3)(nxstr:nxstr) = summpl(nxst+nxstr-1)
      if (seqcod(nxst+nxstr-1) .le. abetsz) then
         wkln1(nxstr:nxstr) = 
     +   aacod1(seqcod(nxst+nxstr-1):seqcod(nxst+nxstr-1))
      else
         wkln1(nxstr:nxstr) = '-'
      end if
 1600 continue
      wkln2 = lefthd(12) // '  ' // sumln(3)(1:aathln) // '  '
     +        // rthd(12)

CHECK v.1.0-->
C      write(outfi,1200) wkln2(1:aathln+22)
CHECK v.1.0<--

      wkln2 = lefthd(13) // '  ' // wkln1(1:aathln) // '  ' //rthd(13)

CHECK v.1.0-->
C      write(outfi,1200) wkln2(1:aathln+22)
CHECK v.1.0<--


CHECK v.1.0-->
C      write(outfi,1700) (nxst-1+j*10,j=1,aathln/10)
CHECK v.1.0<--

CHECK v.3.0.1-->
C 1700 format(12x,10(i10))
CHECK v.3.0.1<--
 2000 continue

      call msgprt

CHECK v.1.0-->
C      write(outfi,*) char(12)
CHECK v.1.0<--

      end
CHECK v.3.6.2-->
C*****************************************************************************
C
C  SUBROUTINE WRISSB  -  If there are any disulphide bonds in the structure
C                        write out to a disulphide bonds file
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE WRISSB(FILDIS,dspart,seqlen,seqcod,seqbcd,aacod3,
     -    hetcod)

      include 'sstruc.par'

      CHARACTER*3   aa3(2)
      CHARACTER*78  FILDIS
      character     aacod3*(abetsz*3), hetcod*63, seqbcd(maxres)*6
      INTEGER       dspart(maxres), ires, NDISUL, nxres, seqcod(maxres),
     -              seqlen, staa

C---- Initialise variables
      NDISUL = 0

C---- Loop through all the residues to count the total number of
C     disulphide bonds
      DO 100, nxres = 1, seqlen

C----     If this disulphide pointer is non-zero, have a bond
          IF (dspart(nxres).GT.0) THEN

C----         Increment count only if it points to a higher-numbered
C             residue in the sequence
              IF (dspart(nxres).GT.nxres) NDISUL = NDISUL + 1
          ENDIF
 100  CONTINUE

C---- If have any disulphide bonds, open the output file
      IF (NDISUL.GT.0) THEN

C----     Open the disulphides output file
          OPEN (UNIT=15,FILE=FILDIS,STATUS='UNKNOWN',ERR=900
CVAX     -    ,CARRIAGECONTROL='LIST'
     -        )

C----     Loop through all the residues to write out the disulphide
C         bonds
          DO 200, nxres = 1, seqlen

C----         If this disulphide pointer is non-zero, have a bond
              IF (dspart(nxres).GT.0) THEN

C----             Write out only if it points to a higher-numbered
C                 residue in the sequence
                  IF (dspart(nxres).GT.nxres) THEN

C----                 Get the amino acid name of this residue
                      if (seqcod(nxres) .le. abetsz) then
                          staa = (seqcod(nxres) - 1) * 3 + 1
                          aa3(1) = aacod3(staa:staa+2)
                      else
                          staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
                          aa3(1) = hetcod(staa:staa+2)
                      end if

C----                 Get the amino acid name of other residue
                      ires = dspart(nxres)
                      if (seqcod(ires) .le. abetsz) then
                          staa = (seqcod(ires) - 1) * 3 + 1
                          aa3(2) = aacod3(staa:staa+2)
                      else
                          staa = (seqcod(ires) - abetsz - 1) * 3 + 1
                          aa3(2) = hetcod(staa:staa+2)
                      end if

C----                 Write this bond out
                      WRITE(15,120) seqbcd(nxres), aa3(1),
     -                    seqbcd(ires), aa3(2)
 120                  FORMAT(A6,1X,A3,' -> ',A6,1X,A3)

                  ENDIF
              ENDIF
 200      CONTINUE


C----     Close the output file
          CLOSE(15)
      ENDIF

      GO TO 999

C---- Error in opening file
900   CONTINUE

999   CONTINUE
      RETURN
      END

C---------------------------------------------------------------------------
CHECK v.3.6.2<--
c*****************************************************************************
CHECK v.3.2-->
C      subroutine carslt(aacod3, aacod1, seqcod, seqbcd, mcangs, cadsts,
C     +                  brksym, altcod, aastd,  brknam, headln, hedlen,
C     +                  autsum, hetcod, oois,   seqlen)
      subroutine carslt(aacod3, aacod1, seqcod, hetcod, seqlen)
CHECK v.3.2<--
c
c     results for c-alpha only files. these are placed in a differently
c     named version of the output file
c
      include 'sstruc.par'
CHECK v.3.2-->
C      character aacod3*(abetsz*3), aacod1*(abetsz), seqbcd(maxres)*6,
C     +          altcod(maxres)*1, aastd(maxres)*1, headln*132,
C     +          brknam*(fnamln), brksym(maxres)*1, autsum(maxres)*1,
C     +          hetcod*63
C      integer   seqcod(maxres), hedlen, seqlen, oois(2,maxres)
C      real      mcangs(nmnang,maxres), cadsts(11,maxres)
      character aacod3*(abetsz*3), aacod1*(abetsz),
     +          hetcod*63
      integer   seqcod(maxres), seqlen
CHECK v.3.2<--

CHECK v.2.1-->
C      integer   hednum, nxres, staa, nxang, nxdis
      integer   hednum, nxres, staa
CHECK v.2.1<--
      character aa3*3, aa1*1

      hednum = 3
CHECK v.3.2-->
C      call hedset(brknam, headln, hedlen, hednum, seqlen)
CHECK v.3.2<--
      do 1000, nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         staa = (seqcod(nxres) - 1) * 3 + 1
         aa3 = aacod3(staa:staa+2)
         aa1 = aacod1(seqcod(nxres):seqcod(nxres))
      else
         staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
         aa3 = hetcod(staa:staa+2)
         aa1 = '-'
      end if

CHECK v.1.0-->
C      write(outfi,300) nxres, brksym(nxres), seqbcd(nxres),
C     +      altcod(nxres), aa3, aa1, aastd(nxres), autsum(nxres),
C     +      (mcangs(nxang,nxres), nxang = chiral,kappa),
C     +      (cadsts(nxdis,nxres), nxdis = 1,11),
C     +      oois(ooi1,nxres), oois(ooi2,nxres)
C  300 format(x,i4,a1,a6,x,a1,x,a3,3(x,a1),x,f6.1,x,f5.1,
C     +       11(x,f4.1),2(x,i2))
CHECK v.1.0<--

 1000 continue
      call msgprt

CHECK v.1.0-->
C      write(outfi,*) char(12)
CHECK v.1.0<--

      end
c*****************************************************************************
      SUBROUTINE msgprt
c
c     write the error messages to the output file
c
      include 'sstruc.par'
      character wkrec*80
      integer   errflg

      rewind(unit=errfi)
      errflg = 0
  200 continue
      read(errfi,300,end=1000) wkrec
  300 format(A)
      if (errflg .eq. 0) then

CHECK v.1.0-->
C         write(outfi,500)'Messages'
CHECK v.1.0<--

CHECK v.2.1.4-->
C  500    format(////,x,A)
CHECK v.3.0.1-->
C  500    format(////,1x,A)
CHECK v.3.0.1<--
CHECK v.2.1.4<--
      end if
      errflg = 1

CHECK v.1.0-->
C      write(outfi,*) wkrec
CHECK v.1.0<--

      write(*,*) wkrec
      go to 200
 1000 continue
      close(unit=errfi)
      END
c****************************************************************************
CHECK v.3.2-->
C      Subroutine hedset(brknam, headln, hedlen, hednum, seqlen)
      Subroutine hedset
CHECK v.3.2<--
c
c     produce headings from results parts 1 and 2 and carslt
c
      include 'sstruc.par'
CHECK v.3.2-->
C      character brknam*(fnamln), headln*132
C      Integer   hedlen, hednum, seqlen
CHECK v.3.2<--


CHECK v.1.0-->
C      write(outfi,101)
C 101  format(' Secondary structure calculation program - ',
C     +              'copyright by David Keith Smith, 1989')
C      write(outfi,100) brknam
C  100 format(x,a)
C      write(outfi,102) headln(1:hedlen)
C 102  format(x,a)
C      write(outfi,120) seqlen
C  120 format(x,'Sequence length -',I5)
C 
C      if (hednum .eq. 1) then
C         write(outfi,1)
C 1       format
C     +  ('             A       A K K                                  ',
C     +   '                                      hydrogen bonding',
C     +   '             Ooi''s')
C         write(outfi,2)
C 2       format
C     + (' strk chain/ l amino u & S structure   bridge        dihedral',
C     +  ' angles     ',
C     +  '                donor    acceptor   donor    acceptor  N  N')
C         write(outfi,3)
C 3       format
C     +  (' num  seq.no t acids t S + patterns   partners    phi    psi',
C     +  '  omega  alpha ',
C     +  'kappa   tco to/energy fr/energy to/energy fr/energy  8 14')
C      else if (hednum .eq. 2) then
C         Write(outfi,4)
C 4       format('                   S Atoms T')
C         Write(outfi,5)
C 5       format
C     +  (' strk chain/ amino t prsnt h  side chain dihedral angles ',
C     +  '  disulphide bridges                  C alpha distances for')
C         Write(outfi,6)
C 6       format
C     +  (' num  seq.no acids d W M S m  Chi1   Chi2   Chi3   Chi4  ',
C     +  ' ptr T CACA CBCB  S/S  i+1  i+2  i+3  i+4  i+5  i+6  i+7',
C     +  '  i+8  i+9 i+10 i+11')
C      else if (hednum .eq. 3) then
C         write (outfi,7)
C 7       format
C     +  ('             A       S A    dihedral    ',
C     +  '                                                     Ooi''s')
C         write(outfi,8)
C 8       format
C     +  (' strk chain/ l amino t u     angles                   ',
C     +  'C alpha distances for                   N  N')
C         write(outfi,9)
C 9       format
C     +  (' num  seq.no t acids d t  alpha kappa  i+1  i+2  i+3  i+4  ',
C     +  'i+5  i+6  i+7  i+8  i+9 i+10 i+11  8 14')
C      end if
CHECK v.1.0<--

      end
c
c*******************************************************************************
c
      logical  function xeqy(x, y)
c     implicit none
c     implicit complex (a-z)
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xeqy = (abs(x-y) .lt. acurcy) 
      return
      end
c
c*******************************************************************************
c
      logical  function xney(x, y)
c     implicit none
c     implicit complex (a-z)
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are not "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xney = (abs(x-y) .gt. acurcy) 
      return
      end
c
c*******************************************************************************
c
      subroutine getstr (prompt, defans, string, ilen, ifail)
c
c     implicit complex(a-z)
c     implicit none
c===============================================================================
c     General purpose routine to prompt the user to type in a character
c     string. The prompt that will be written out to the user should be passed
c     to the routine in variable PROMPT, along with a default response in
c     DEFANS. The variable STRING will return the user's input to the calling
c     program, and ILEN will hold its length. The flag IFAIL will have one of
c     the following values:
c        IFAIL = 0 : everything OK;
c        IFAIL = 1 : variable PROMPT is empty;
c        IFAIL = 2 : variable DEFANS is empty;
c        IFAIL = 3 : user's response is too long to fit in STRING;
c        IFAIL = 4 : error occurred writing out prompt to terminal;
c        IFAIL = 5 : error occurred reading in reply from keyboard;
c        IFAIL = 6 : EOF detected reading in reply from keyboard.
c===============================================================================
      character*(*) prompt, defans, string
      integer       ilen, ifail
c
      character     line*255
      integer       lenp, lens, strlen, lend
c
CHECK v.3.0.1-->
C 1000 format (A, ' <', A, '> - ')
CHECK v.3.0.1<--
 1010 format (A)
c
      ilen  = 0
      ifail = 0
c
      lenp  = strlen(prompt)
      if (lenp .le. 0)  then
         ifail = 1
         return
      end if
c
      lend  = strlen(defans)
      if (lend .le. 0)  then
         ifail = 2
         return
      end if
c
      lens  = len(string)
c
CHECK v.3.0.1-->
C   10 write (6, 1000, err = 904) prompt(1:lenp), defans(1:lend)
CHECK v.3.0.1<--
      read  (5, 1010, err = 905, end = 906) line
      ilen = strlen(line)
      if (ilen .eq. 0) then
         line = defans(1:lend)
         ilen = lend
      end if
      if (ilen .gt. lens)  then
         ifail = 3
         return
      end if
      string = line(1:ilen)
      ifail = 0
      return
c
CHECK v.3.0.1-->
C  904 ifail = 4
C      return
CHECK v.3.0.1<--
  905 ifail = 5
      return
  906 ifail = 6
      return
c
      end 
c
c*******************************************************************************
c
      integer function lastsl (string)
c     implicit complex(a-z)
c     implicit none
      character *(*) string
      integer   i, lens
c
      lens = len(string)
      do 10 i = lens, 1, -1
      if (string(i:i).eq.'/')  then
         lastsl = i
         return
      end if
   10 continue
      lastsl = 0
      return
      end
c
c*****************************************************************************
c
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

C---------------------------------------------------------------------------
CHECK v.2.0<--

CHECK v.3.1-->
C*****************************************************************************
C
C  SUBROUTINE GETMOD  -  Read through the PDB file until hit the next
C                        MODEL record
C
C----------------------------------------------------------------------+--- 

CHECK v.3.2-->
C      SUBROUTINE GETMOD(BRKFI,IMODEL,MXMODL,ENDFIL,IERROR)
CHECK v.3.4-->
C      SUBROUTINE GETMOD(BRKFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR)
CHECK v.3.4.3-->
C      SUBROUTINE GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,IERROR,NMR)
CHECK v.3.5-->
C      SUBROUTINE GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,IERROR)
      SUBROUTINE GETMOD(BRKFI,IMODEL,MXMODL,MODLNO,ENDFIL,AAPLN,AAIN,
     -    GAPIN,INSEQR,MAXTYP,SQRSCT,IERROR)
CHECK v.3.5<--
CHECK v.3.4.3<--
CHECK v.3.4<--
CHECK v.3.2<--

      SAVE

      CHARACTER*80  IREC
CHECK v.3.4-->
C      INTEGER       BRKFI, IMODEL, LINE, MXMODL
      INTEGER       BRKFI, IMODEL, LINE, MODLNO, MXMODL
CHECK v.3.4<--
CHECK v.3.2-->
C      LOGICAL       ENDFIL, FOUND, IERROR
CHECK v.3.4.3-->
C      LOGICAL       ENDFIL, FOUND, IERROR, NMR
      LOGICAL       ENDFIL, FOUND, IERROR
CHECK v.3.4.3<--
CHECK v.3.2<--

CHECK v.3.5-->
      INTEGER       AAPLN

      CHARACTER*1   GAPIN(AAPLN)
      CHARACTER*3   AAIN(AAPLN)
      CHARACTER*(*) INSEQR
      INTEGER       I, MAXTYP, NBAD, RESFND, SQRSCT
      LOGICAL       BADSEQ
CHECK v.3.5<--

C---- Initialise variables
CHECK v.3.5-->
      BADSEQ = .FALSE.
CHECK v.3.5<--
      ENDFIL = .FALSE.
      FOUND = .FALSE.
      IERROR = .FALSE.
      LINE = 0

C---- Loop through the file until hit a MODEL record, the first ATOM
C     record, or the EOF
 100  CONTINUE

C----     Read in the next record
          READ(BRKFI,120,ERR=900,END=500) IREC
 120      FORMAT(A)

C----     Test if this is a MODEL record
          IF (IREC(1:5).EQ.'MODEL') THEN
              LINE = LINE + 1
CHECK v.3.4-->
C              READ(IREC,140,ERR=902) IMODEL
              READ(IREC,140,ERR=902) MODLNO
              IMODEL = IMODEL + 1
CHECK v.3.4<--
 140          FORMAT(11X,I3)
              FOUND = .TRUE.
              IF (IMODEL.GT.MXMODL) GO TO 904
CHECK v.3.2-->
CHECK v.3.4-->
C              NMR = .TRUE.
CHECK v.3.4<--
CHECK v.3.5-->
C          ELSE IF (IREC(1:6).EQ.'ATOM  ') THEN
          ELSE IF (IREC(1:6).EQ.'ATOM  ' .OR.
     -        IREC(1:6).EQ.'ATZERO') THEN
CHECK v.3.5<--
              IMODEL = 1
              FOUND = .TRUE.
              BACKSPACE(BRKFI)
CHECK v.3.2<--

CHECK v.3.5-->
C----     If this is a SEQRES record: store the details
          ELSE IF (IREC(1:6).EQ.'SEQRES') THEN
              READ(IREC,220) (GAPIN(I), AAIN(I), I = 1, AAPLN)
 220          FORMAT(18X,13(A1,A3))
              DO 250 I = 1, AAPLN
                  IF (AAIN(I).NE.'   ') THEN
                      IF (RESFND(INSEQR,AAIN(I)).EQ.0) THEN
                          SQRSCT = SQRSCT + 1
                          IF (SQRSCT.GT.MAXTYP) THEN
                              PRINT*, '*** Warning. Maximum number of ',
     -                            'different types of residue in SEQRE',
     -                            'S records, MAXTYP, exceeded', MAXTYP
                              SQRSCT = MAXTYP
                          ENDIF
                          INSEQR(((SQRSCT-1)*3)+1:SQRSCT*3) = AAIN(I)
                      ENDIF
                   ENDIF
 250          CONTINUE

C----         If any of the characters between the residues are not
C             blank, then discard all the SEQRES data read in (it's
C             probably a comment, or something!)
              NBAD = 0
              DO 300, I = 1, AAPLN
                  IF (GAPIN(I).NE.' ') THEN
                      BADSEQ = .TRUE.
                      NBAD = NBAD + 1
                  ENDIF
 300          CONTINUE
              IF (NBAD.GT.0) THEN
                  PRINT*, '*** Warning. Invalid data in SEQRES record.',
     -                ' SEQRES data will be ignored'
                  PRINT*, '*** ', IREC(1:66)
              ENDIF
CHECK v.3.5<--
          ENDIF

C---- Loop back for next record
      IF (.NOT.FOUND) GO TO 100
CHECK v.3.5-->
C      GO TO 999
      GO TO 800
CHECK v.3.5<--

C---- End of file reached
 500  CONTINUE
      ENDFIL = .TRUE.

CHECK v.3.5-->
C---- If any bad SEQRES records were encountered, then discard all the
C     SEQRES data
 800  CONTINUE
      IF (BADSEQ) THEN
          SQRSCT = 0
      ENDIF
CHECK v.3.5<--

      GO TO 999

C---- Error in file
900   CONTINUE
      PRINT*, ' *** ERROR reading NMR file while looking for MODEL',
     -    ' record'
      IF (IMODEL.EQ.0) THEN
          PRINT*, ' *** Error at line', LINE
      ELSE
          PRINT*, ' *** Error at line', LINE, '  after ENDMDL record',
     -        ' for model', IMODEL
      ENDIF
      GO TO 990

 902  CONTINUE
      PRINT*, ' *** Data error in MODEL number'
      IF (IMODEL.EQ.0) THEN
          PRINT*, ' *** Error at line', LINE
      ELSE
          PRINT*, ' *** Error at line', LINE, '  after ENDMDL record',
     -        ' for model', IMODEL
      ENDIF
      GO TO 990

 904  CONTINUE
      PRINT*, ' *** Maximum number of models exceeded'
      PRINT*, ' *** Model found:', IMODEL, '     Maxmimum allowed:',
     -    MXMODL
      GO TO 990


 990  CONTINUE
      IERROR = .TRUE.

999   CONTINUE
      RETURN
      END

C---------------------------------------------------------------------------
CHECK v.3.1<--
