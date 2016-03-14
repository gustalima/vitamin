C**************************************************************************
C
C  TPLOT.F - Program to plot torsion angle distributions, including
C            Ramachandran plots, and to calculate torsion angle G-factors
C
C            Written by Roman Laskowski, University College, London,
C            November 1993.
C
C            Original version was part of v.3.0 of the PROCHECK suite
C            of programs.
C
C            Subsequent amendments will be labelled by CHECK v.m.n-->
C            and CHECK v.m.n<-- where m.n is the version number
C            corresponding to the change
C
C  v.3.0.1 - Minor bug-fix: array out of bounds causing program to hang
C            on SG machines.
C            Correction of integer variables incorrectly defined as reals.
C            Correction of another array-bound problem.
C            Commented out unreferenced labels (producing compiler warnings
C            when compiling on Convex).
C                                             Roman Laskowski (29 Mar 1994)
C
C  v.3.1   - Transfer of FINKEY, GETCOL and GETNAM routines to ps.f. Unit
C            number for procheck.prm file changed to 10 for consistency.
C            Addition of plot description to display of plot file name.
C                                           Roman Laskowski (8-19 Apr 1994)
C
C  v.3.2   - Addition of identifying plot handle to each plot filename.
C            Addition of descriptive plot title to header of PostScript
C            file.
C            Addition of .sum file holding a summary page.
C            Addition of option to print plot filename on the plot itself.
C            Chain ID made case-insensitive.
C                                          Roman Laskowski (22-30 Apr 1994)
C            Addition of option allowing plots of same type to be combined
C            in a single PostScript file.
C            Amendments for PROCHECK-NMR:
C               Pick up of parameters from the procheck_nmr.prm file
C               Read in the individual .rin files listed in mplot.in.
C               Print in title of number of models in the ensemble.
C                                          Roman Laskowski (11-18 Oct 1994)
C            Minor amendments to various statements to make them
C            acceptable to f2c, and to deal with various uninitialised
C            variables.(Amendments supplied by Dave Love at Daresbury).
C                                  David Love/Roman Laskowski (18 Oct 1994)
C            Print of model-numbers for ensembles.
C                                              Roman Laskowski (3 Nov 1994)
C            Addition of PROCHECK-COMP heading for plots when the
C            ensemble version is run, and use of procheck_comp.prm
C            parameter file.
C                                             Roman Laskowski (10 Nov 1994)
C
C  v.3.3   - Additon of option to split the Ramachandrans onto a separate
C            page for each member of the ensemble.
C                                             Roman Laskowski (13 Feb 1995)
C
C  v.3.3.2 - Output of summary statistics in html format
C                                              Roman Laskowski (21 Aug 1995)
C
C  v.3.4   - Additional options for Ramachandran plot based on suggestions
C            by Ton Rullmann, including smaller boxes, filled boxes and
C            removal of model numbers.
C                      Ton Rullmann (Jun 1995) Roman Laskowski (12 Feb 1996)
C            Implementation of option to split the Ramachandrans onto a
C            separate page for each member of the ensemble.
C            Addition of range-selection allowing outputs to be limited
C            for certain models and certain ranges of residue-numbers.
C                                              Roman Laskowski (22 Mar 1996)
C            Amendments of plot headings where model- or residue-ranges
C            apply.
C                                              Roman Laskowski (25 Mar 1996)
C            Increase of number of colours that can be user-defined.
C            Bug-fix to datapoint colours on Ramachandran plot when printed
C            in colour.
C                                              Roman Laskowski (24 Apr 1996)
C
C  v.3.4.3 - Transfer of GETDAT routine to ps.f.
C            Colouring of points on chi1-chi2 plots and individual residue
C            Ramachandran plots according to G-factor value. Smearing of
C            these distributions.
C                                           Roman Laskowski (22-23 May 1996)
C            New version of routine PROCPT to plot points in order of
C            worsening G-factor so that the model-numbers of the outliers
C            will be more clearly visible.
C                                              Roman Laskowski (28 May 1996)
C            New option in parameter file to allow "publication version" of
C            Ramachandran plot to be generated (ie without the headings,
C            borders and statistics). [Suggested by Gerard Kleywegt].
C                                              Roman Laskowski (22 Jul 1996)
C            Output of HTML file only if being run for WWW pages.
C                                              Roman Laskowski ( 1 Oct 1996)
C
C  v.3.4.4 - Bug-fix to Ramachandran plot for Black-and-white option.
C                                              Roman Laskowski (30 Oct 1996)
c
C  v.3.5   - Bug-fix to Ramachandran plot - residue labels being printed
C            in white text.
C                                              Roman Laskowski (21 Nov 1996)
C
C  v.3.5.2 - Extra parameter to call to GETDAT (required elsewhere).
C                                              Roman Laskowski (20 May 1999)
C            Generalisation of arrays holding the torsion angle distributions
C            to allow addition of other data.
C                                           Roman Laskowski (22-25 May 1999)
C            Addition of CA-angle/CA-torsion angle distribution plots.
C                                              Roman Laskowski (28 May 1999)
C
C--------------------------------------------------------------------------
C
C Compilation and linking (on unix)
C -----------------------
C
C f77 -u -c tplot.f
C f77 -u -c ps.f
C f77 -o tplot tplot.o ps.o
C
C Compilation and linking (on VAX VMS)
C -----------------------
C
C FORT TPLOT
C FORT PS
C LINK TPLOT, PS
C
C--------------------------------------------------------------------------
C
C     Files
C     -----
C
C  1  mplot.in       - File holding list of input .rin file making up
C                      an ensemble of NMR structures
C  2  prodata        - File holding the torsion angle distributions for
C                      phi-psi and chi1-chi2 combinations. Data generated
C                      by program gentors from a given database of high
C                      resolution protein structures. (Latest version is
C                      163 non-homologous structures from the Apr 93
C                      release of the Brookhaven databank)
C  3  <filename>.rin - File generated by program SECSTR, holding all the
C                      residue-by-residue torsion angle information
C  4  <filename>.lan - File generated by program ANGLEN, holding all the
C                      mainchain bond lengths and bond angles
C  7  <filename>.sdh - Output file holding the residue-by-residue log-odds
C                      scores calculated here
C  10 procheck.prm   - Input parameter file containing user-defined
C                      options governing the plots produced. Name of
C                      file might be:
C                         procheck_nmr.prm - when NMR option being run
C                         procheck_comp.prm - when COMP option being run
C                                            on an ensemble of structures
C  11 <filename>_nn.ps - Output PostScript files, numbered nn = 01, 02, ...
C  12  ps.number     - Input file holding last-used number for the
C                      PostScript files
C  14 <filename>.sum - Output summary page showing summarised results from
C                      all the plots
C  15 <filename>.html- Output html file giving summary statistics in html
C                      format
C  16 <user-defined> - Optional input file holding the model- and residue-
C                      number ranges to be included in the plots.
C
C--------------------------------------------------------------------------
C
C     Subroutine calling tree
C     -----------------------
C
C     MAIN    --> INITS
C
C                 Get filename and resolution of structure, and read in
C                 plot parameters
C             --> GETCOD  --> GETNAM
C                         --> GETFIL
C             --> PARAMS  --> FINKEY
C                         --> GETCOL
C             --> GETRNG  --> GETOKN
C                         --> INTOKN  --> LENSTR
C                         --> DELTOK
C                         --> STOTOK
C                         --> PRNRNG
C             --> GETWNT  --> INMODL
C             --> NORMGS  --> CALCOV
C             --> OPESUM
C             --> OPEHTM
C
C                 Read in and uncompress the torsion angle distributions
C             --> GETDAT
C             --> SMEAR
C             --> CALCLO
C
C                 Read in the data points from the structure
C             --> GETPTS  --> RAMOPE   --> PSNAME
C                                      --> PSOPEN
C                         --> RAMACH   --> PSPAGE
C                                      --> PSTEXT
C                                      --> PSTRAN
C                                      --> PSCTXT
C                                      --> RAMSHD  --> PSLWID
C                                                  --> PSLINE
C                                                  --> RAMFIL  --> PSLWID
C                                                              --> RAMREG
C                                                              --> PSHADE
C                                                              --> PSUBOX
C                                                              --> PSCOLR
C                                                  --> RAMLIN  --> PSLWID
C                                                              --> RAMREG
C                                                              --> PSLINE
C                                                  --> INCLET  --> PSHADE
C                                                              --> PSCTXT
C                                      --> RAMAXE  --> AXES    --> PSLWID
C                                                              --> PSLINE
C                                                              --> PSCTXT
C                                                  --> PSLWID
C                                                  --> PSDASH
C                                                  --> PSLINE
C                                                  --> PSCTXT
C                                                  --> PSRTXT
C                          --> INRANG  --> LENSTR
C                          --> RAMREG
C                          --> RAMPLT  --> PSMARK  --> PSBBOX
C                                                  --> PSTRIA
C                                      --> PSCTXT
C                          --> FILKEY  --> PSLWID
C                                      --> PSHADE
C                                      --> PSBBOX
C                                      --> PSCTXT
C                                      --> PSTEXT
C                          --> MSTATS
C                          --> RAMEND  --> RMSTAT  --> RAMPUT  --> PSCTXT
C                                                              --> PSTEXT
C                                                  --> PSCTXT
C                                      --> PSENDP
C                          --> PSCLOS
C
C                 Plot the torsion angle distributions
C             --> PLOTS   --> PSNAME
C                         --> PSOPEN
C                         --> PSPAGE
C                         --> PLTGRD  --> PLTDIS --> PSLWID
C                                                --> PSHADE
C                                                --> PSUBOX
C                                                --> PSCALE
C                                                --> PSLWID
C                                                --> PSDASH
C                                                --> PSLINE
C                                     --> AXES   --> PSLWID
C                                                --> PSLINE
C                                                --> PSCTXT
C                                     --> PSCTXT
C                                     --> PSRTXT
C                         --> PROCPT  --> CALC1D
C                                     --> CALC2D
C                                     --> PSLWID
C                                     --> PSLINE
C                                     --> PSCALE
C                                     --> PSBBOX
C                                     --> PSCTXT
C                                     --> PSCOLB
C                         --> FILKEY  --> PSLWID
C                                     --> PSHADE
C                                     --> PSBBOX
C                                     --> PSCTXT
C                                     --> PSTEXT
C                         --> PLTEND  --> PSTEXT
C                                     --> PSCTXT
C                         --> PSENDP
C                         --> PSCLOS
C
C                 Calculate covalent quality scores from main-chain bond
C                 lengths and angles
C             --> GETLAN
C             --> COVSCO  --> CALCOV
C                         --> RESCOR
C
C                 Write out the calculated, residue-by-residue log-odds
C                 scores
C             --> WRISDH
C             --> PUTPSN
C
C----------------------------------------------------------------------+---


      PROGRAM TPLOT

      INCLUDE 'tplot.inc'

CHECK v.3.2-->
C      INTEGER       DISTRB
CHECK v.3.5.2-->
C      INTEGER       DISTRB, LOOPS
      INTEGER       DISTRB, LOOPS, NCELL1, NCELL2
CHECK v.3.5.2<--
      LOGICAL       BOTHND
CHECK v.3.2<--

CHECK v.3.5.2-->
      REAL          ENERGY(MXCELL*(NAMINO+1)), NOBSER(MXCELL*(NAMINO+1))
CHECK v.3.5.2<--

C---- Initialise variables
      CALL INITS
      IF (IFAIL) GO TO 999

CHECK v.3.2-->
C---- Read in the code of the Brookhaven file and the structure's
C     resolution
      CALL GETCOD
      IF (IFAIL) GO TO 999
CHECK v.3.2<--

C---- Read in the program parameters
      CALL PARAMS
      IF (IFAIL) GO TO 999

CHECK v.3.4-->
C---- Read in the model- and residue-ranges in the supplied ranges file,
C     if there is one
      IF (HAVRAN) THEN
          CALL GETRNG(FILRNG,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,MRANGE,
     -        NRANGE,BOTHND,IFAIL)
          IFAIL = .FALSE.

C----     Check whether the user has selected a range of residues
          IF (RESFRM(1).NE.'*ALL  ') RSELEC = .TRUE.
      ENDIF

C---- If this is an NMR structure, then find which of the files
C     contain the wanted models
      IF (NMR) THEN
          CALL GETWNT(MXFILE,NFILE,FILRIN,FFILE,MAXRNG,MODFRM,MODTO,
     -        MRANGE,RSELEC,MWANT,NMODEL,MODNUM,ACTNUM,TOPMOD,
     -        BRCODE,NAMLEN,BLEN,IFAIL)
          IF (IFAIL) GO TO 999

C---- Model ranges only apply to NMR ensembles, so set to default
C     if this is not an NMR structure
      ELSE
          MRANGE = 1
          MODFRM(MRANGE) = -99999
          MODTO(MRANGE) = 99999
      ENDIF
CHECK v.3.4<--

C---- Calculate normalization factors for the Gaussian distribs used for
C     bond lengths and bond angles
      CALL NORMGS

CHECK v.3.2-->
CC---- Read in the code of the Brookhaven file and the structure's
CC     resolution
C      CALL GETCOD
C      IF (IFAIL) GO TO 999
CHECK v.3.2<--

CHECK v.3.2-->
C---- Open the summary file
      CALL OPESUM(FILSUM,.TRUE.,IFAIL)
      IF (IFAIL) GO TO 999

C---- For NMR or ensemble, only need to loop over the distributions to
C     be plotted
      IF (NMR .OR. ENSEMB) THEN
          LOOPS = 2
      ELSE
          LOOPS = NDISTR
      ENDIF
CHECK v.3.2<--

CHECK v.3.3.2-->
C---- Open the html file
CHECK v.3.4.3-->
      IF (WWWOUT) THEN
CHECK v.3.4.3<--
          CALL OPEHTM(FILHTM,.TRUE.,IFAIL)
          IF (IFAIL) GO TO 999
CHECK v.3.4.3-->
      ENDIF
CHECK v.3.4.3<--
CHECK v.3.3.2<--

C---- Loop for the different types of distributions to be plotted
CHECK v.3.2-->
C      DO 100, DISTRB = 1, NDISTR
      DO 100, DISTRB = 1, LOOPS
CHECK v.3.2<--

C----     Read in the torsion angle distributions and uncompress the
C         data
CHECK v.3.4.3-->
C          CALL GETDAT(DISTRB)
CHECK v.3.5.2-->
C          CALL GETDAT(DISTRB,TWODEE,NOBSER,NCELL,NCELL1,NAMINO,VALBEG,
C     -        VALEND,STEP,NCOUNT,NRMEAN,NRMSTD,2,IFAIL)
          CALL GETDAT(DISTRB,TWODEE,NOBSER,MXCELL,NCELL1,NCELL2,NAMINO,
     -        VALBEG,VALEND,STEP,NCOUNT,NRMEAN,NRMSTD,2,.FALSE.,IFAIL)
CEXTRACT
          IF (DISTRB.EQ.1) THEN
              CALL WRIDST(NOBSER,NCELL1,NCELL2)
          ENDIF
CEXTRACT
CHECK v.3.5.2<--
CHECK v.3.4.3<--
          IF (IFAIL) GO TO 999

CHECK v.3.4.3-->
C----     Smear the observations over the 2D distributions
          IF (TWODEE) THEN
CHECK v.3.5.2-->
C              CALL SMEAR(NOBSER,ENERGY,NCELL,NCOUNT,NAMINO)
              CALL SMEAR(NOBSER,ENERGY,NCELL1,NCELL2,NCOUNT,NAMINO)
CHECK v.3.5.2<--
          ENDIF

C----     Calculate log-odds scores
CHECK v.3.5.2-->
C          CALL CALCLO(NOBSER,ENERGY,NCELL,NAMINO,NCOUNT)
          CALL CALCLO(NOBSER,ENERGY,NCELL1,NCELL2,NAMINO,NCOUNT)
CHECK v.3.5.2<--
CHECK v.3.4.3<--

C----     Read in the data points from the appropriate file

CHECK v.3.5.2-->
C----     If looking at the CA-angle/CA torsion angle distribution
C         read in the data from the .ca file
          IF (DISTRB.EQ.7) THEN

C----         Read in the CA angle data
              CALL GETCAD(DISTRB)

C----     Otherwise have main-chain and side-chain torsion angles
          ELSE
CHECK v.3.5.2<--

C----         Read in the data points from the .rin file
              CALL GETPTS(DISTRB)
CHECK v.3.5.2-->
          ENDIF
CHECK v.3.5.2<--
          IF (IFAIL) GO TO 999

C----     Plot the data superimposed on the distributions
CHECK v.3.5.2-->
C          CALL PLOTS(DISTRB)
          CALL PLOTS(DISTRB,NOBSER,ENERGY,NCELL1,NCELL2)
CHECK v.3.5.2<--

 100  CONTINUE

C---- Read through the main-chain bond lengths and bond angles file to
C     calculate mean values
CHECK v.3.2-->
      IF (.NOT.NMR .AND. .NOT.ENSEMB) THEN
CHECK v.3.2<--
          CALL GETLAN
          IF (IFAIL) GO TO 999

C----     Calculate covalent quality scores from main-chain bond lengths
C         and angles
          CALL COVSCO

C----     Write out the calculated residue-by-residue log-odds scores
          CALL WRISDH
CHECK v.3.2-->
      ENDIF
CHECK v.3.2<--

C---- Write out the current PostScript plot number
      CALL PUTPSN(IPLOT)

 999  CONTINUE
      IF (IFAIL) THEN
          PRINT*, '*** Program tplot terminated with error'
CHECK v.3.2-->
          WRITE(14,*)
          WRITE(14,*) '*** Program tplot terminated with error. See lo',
     -        'g file:  tplot.log'
          WRITE(14,*)
CHECK v.3.2<--
      ELSE
         PRINT*, '* Program complete'
      ENDIF
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE INITS  -  Initialise variables
C
C----------------------------------------------------------------------+---

      SUBROUTINE INITS

      INCLUDE 'tplot.inc'

      CHARACTER*3   AMNAME(NAMINO)
CHECK v.3.4-->
C      INTEGER       IAMINO, IDISTR, IRESID, ITYPE
      INTEGER       IAMINO, IDISTR, IFILE, IRESID, ITYPE
CHECK v.3.4<--

      DATA AMNAME /'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
     -             'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
     -             'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' /

C---- Initialise variables
      ALLRAM = .FALSE.
      DO 50, IAMINO = 1, NAMINO
          AMINO(IAMINO) = AMNAME(IAMINO)
 50   CONTINUE
      CHAIN = ' '
CHECK v.3.2-->
      COMBPS = .FALSE.
CHECK v.3.2<--
      CRONLY = .FALSE.
      DO 60, ITYPE = 1, 3
          COVAVE(ITYPE) = 0.0
          NUMAVE(ITYPE) = 0
 60   CONTINUE
      DO 80, IDISTR = 1, NDISTR
CHECK v.3.5.2-->
C          IF (IDISTR.GT.2) THEN
          IF (IDISTR.GT.2 .AND. IDISTR.LT.7) THEN
CHECK v.3.5.2<--
              DOPLOT(IDISTR) = .FALSE.
          ELSE
              DOPLOT(IDISTR) = .TRUE.
          ENDIF
 80   CONTINUE
CHECK v.3.4-->
      FFILE = 1
      HAVRAN = .FALSE.
CHECK v.3.4<--
      IFAIL = .FALSE.
      INCOLR = .FALSE.
      INCOLG = .FALSE.
      INCOLC = .FALSE.
      IPLOT = 0
      LABRES = 1
      LIMCHI = -2.0
      LIMGP = -2.0
      INCLET = .TRUE.
CHECK v.3.2-->
      NFILE = 1
      ENSEMB = .FALSE.
      NMR = .FALSE.
CHECK v.3.2<--
      NOLINE = .FALSE.
      NOSHAD = .FALSE.
CHECK v.3.4-->
      MRANGE = 1
      MODFRM(MRANGE) = -99999
      MODTO(MRANGE) = 99999
      DO 100, IFILE = 1, MXFILE
          MWANT(IFILE) = .TRUE.
 100  CONTINUE
      NRANGE = 1
      NMODEL = 1
      RESFRM(NRANGE) = '*ALL  '
      RESTO(NRANGE) = 'XXXXXX'
CHECK v.3.4<--
CHECK v.3.2-->
      PLABEL = .TRUE.
CHECK v.3.2<--
      PLOTRM = .TRUE.
CHECK v.3.4-->
      RBXFIL = .FALSE.
      RBXMOD = .TRUE.
      RBXSIZ = 1.0
      RSELEC = .FALSE.
CHECK v.3.4<--
      DO 300, IRESID = 1, MXRES
CHECK v.3.2-->
          MODEL(IRESID) = 0
CHECK v.3.2<--
          SAVRES(IRESID) = 0
          DO 150, ITYPE = 1, 3
              SCOVAL(ITYPE,IRESID) = 0.0
 150      CONTINUE
          DO 200, IDISTR = 1, NDISTR
              SCORE(IDISTR,IRESID) = 999.99
 200      CONTINUE
 300  CONTINUE
      USENGH = .TRUE.
CHECK v.3.2-->
      WITHAN = .FALSE.
CHECK v.3.2<--
CHECK v.3.4.3-->
      WWWOUT = .FALSE.
CHECK v.3.4.3<--

C---- Open data file, prodata
      OPEN(UNIT=2, FILE='prodata', STATUS='OLD',
     -     FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST', READONLY,
     -     ERR=900)

      GO TO 999

900   CONTINUE
      PRINT*, '*** ERROR. Unable to open parameter file, prodata'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETCOD  -  Read in the name of the .rin file to be processed
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCOD

      INCLUDE 'tplot.inc'

CHECK v.3.2-->
C      CHARACTER*80  PDBFIL
C      INTEGER       LINE
      CHARACTER*1   LOWERA, YESNO
      CHARACTER*3   NUMBER
      CHARACTER*26  UPPER
CHECK v.3.4-->
C      INTEGER       IPOS, LINE, N
      INTEGER       IPOS, LENSTR, LINE, N
CHECK v.3.4<--
CHECK v.3.2<--
      LOGICAL       IERROR

CHECK v.3.2-->
      DATA  LOWERA / 'a' /
      DATA  UPPER  / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
CHECK v.3.2<--

C---- Initialise variables
      LINE = 0

C---- Accept name of original .pdb file holding the structure
CHECK v.3.4-->
C      PRINT*, 'Enter filename of coordinates'
      PRINT*, 'Enter filename of coordinates (prefix with a % for',
     -    ' a non-NMR ensemble'
CHECK v.3.4<--
      READ(*,110,ERR=900) PDBFIL
CHECK v.3.2-->
      IF (PDBFIL(1:1).EQ.'%') THEN
          PDBFIL = PDBFIL(2:)
          ENSEMB = .TRUE.
          NMR = .FALSE.
      ENDIF
CHECK v.3.2<--
 110  FORMAT(A)

C---- Accept chain-ID
      PRINT*, 'Enter required chain-ID, or leave blank for all'
      READ(*,110,ERR=906) CHAIN

CHECK v.3.2-->
C---- Convert chain ID to upper-case if necessary
      N = ICHAR(CHAIN) - ICHAR(LOWERA) + 1
      IF (N.GE.1 .AND. N.LE.26) CHAIN = UPPER(N:N)
CHECK v.3.2<--

C---- Accept structure's resolution
      PRINT*, 'Enter resolution'
      READ(*,*,ERR=904) RESOL

CHECK v.3.2-->
C---- Accept whether processing an ensemble of structures
      PRINT*, 'Is this an ensemble of structures (Y/N)?'
      READ(*,110) YESNO
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          IF (.NOT.ENSEMB) NMR = .TRUE.
      ENDIF
CHECK v.3.2<--

C---- Peel off directory path and extension
      CALL GETNAM(PDBFIL,ISTART,IEND,IERROR)
      IF (IERROR) GO TO 990

C---- Form names of other files that will be required in default directory
      BRCODE = PDBFIL(ISTART:IEND)
      FILLAN = PDBFIL(ISTART:IEND) // '.lan'
CHECK v.3.2-->
C      FILRIN = PDBFIL(ISTART:IEND) // '.rin'
      FILRIN(1) = PDBFIL(ISTART:IEND) // '.rin'
CHECK v.3.2<--
      FILSDH = PDBFIL(ISTART:IEND) // '.sdh'
CHECK v.3.2-->
      FILSUM = PDBFIL(ISTART:IEND) // '.sum'
CHECK v.3.2<--
CHECK v.3.3.2-->
      FILHTM = PDBFIL(ISTART:IEND) // '.html'
CHECK v.3.3.2<--
      FILPS = BRCODE
      PSLEN = IEND - ISTART + 1
      BLEN = PSLEN
      IF (CHAIN.NE.' ') THEN
          BRCODE = BRCODE(1:BLEN) // ' - Chain ' // CHAIN
          BLEN = BLEN + 10
          IF (BLEN.GT.78) BLEN = 78
      ENDIF

CHECK v.3.4-->
C---- If separate Ramachandran for each model, then will label
C     each one with the appropriate model number
      RAMHED = BRCODE(1:BLEN) // ' (Model '
      RLEN = BLEN + 14

C---- Save the length of the PDB code plus chain ID (if relevant)
      NAMLEN = BLEN
CHECK v.3.4<--

CHECK v.3.2-->
C---- If this is an ensemble, read in the list of .rin filenames
      IF (NMR .OR. ENSEMB) THEN
          CALL GETFIL
          IF (IFAIL) GO TO 999
          
C----     Add number of models to plot title
          IF (BLEN.LT.65) THEN
              WRITE(NUMBER,120) NFILE
 120          FORMAT(I3)
              IPOS = 1
              IF (NUMBER(1:1).EQ.' ') IPOS = 2
              IF (NUMBER(2:2).EQ.' ') IPOS = 3
              BRCODE = BRCODE(1:BLEN) // ' (' // NUMBER(IPOS:) //
     -            ' models)'
              BLEN = BLEN + 14 - IPOS
          ENDIF

C---- Otherwise, if only a single .rin file, then open it
      ELSE
CHECK v.3.2-->
C      OPEN(UNIT=3,FILE=FILRIN,STATUS='OLD',ACCESS='SEQUENTIAL',
CCVAX     -     READONLY,
C     -    FORM='FORMATTED',ERR=902)
          OPEN(UNIT=3,FILE=FILRIN(1),STATUS='OLD',ACCESS='SEQUENTIAL',
CVAX     -     READONLY,
     -    FORM='FORMATTED',ERR=902)
      ENDIF
CHECK v.3.2<--

CHECK v.3.4-->
C---- Read in (optional) name of ranges file
      FILRNG = ' '
      IF (NMR) THEN
          PRINT*, 'Enter name of ranges file (blank if none)'
          READ(*,110,ERR=900) FILRNG
      ENDIF
      IF (FILRNG.EQ.' ') THEN
          HAVRAN = .FALSE.
      ELSE
          HAVRAN = .TRUE.
      ENDIF
CHECK v.3.4<--

      GO TO 999

C---- Fatal errors
 900  CONTINUE
      PRINT*, '*** ERROR. Data error in entered filename'
      GO TO 990

 902  CONTINUE
      PRINT*, '*** Error opening .rin file.'
CHECK v.3.2-->
CHECK v.3.4-->
C      PRINT*, FILRIN(1)
      PRINT*, '* [', FILRIN(1)(1:LENSTR(FILRIN(1))), ']'
CHECK v.3.4<--
CHECK v.3.2<--
      GO TO 990

 904  CONTINUE
      PRINT*, '*** Error in entered resolution'
      GO TO 990

 906  CONTINUE
      PRINT*, '*** ERROR. Data error in entered chain-ID'
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2-->
C**************************************************************************
C     
C  SUBROUTINE GETFIL  -  Read in the list of input .rin files in mplot.in
C
C----------------------------------------------------------------------+---
      
      SUBROUTINE GETFIL
      
      INCLUDE 'tplot.inc'
      
      CHARACTER*80  FNAME
      INTEGER       IPOS
      LOGICAL       FOUND
      
C---- Initialise variables
      NFILE = 0

C---- Open input file, mplot.in
      OPEN(UNIT=1, FILE='mplot.in', STATUS='OLD', FORM='FORMATTED',
     -     ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST', READONLY,
     -     ERR=900)
      
C---- Read through the file
 100  CONTINUE
      
C----     Read in the next filename
          READ(1,120,END=500,ERR=902) FNAME
 120      FORMAT(A)
          IF (FNAME(1:1).EQ.' ') GO TO 100

C----     Initialise variables for extraction of filename
          IEND = 80
          IPOS = 0
          FOUND = .FALSE.

C----     Search the line for end of filename
 200      CONTINUE
              IPOS = IPOS + 1
              IF (FNAME(IPOS:IPOS).EQ.' ') THEN
                  FOUND = .TRUE.
                  IEND = IPOS - 1
              ENDIF
          IF (.NOT.FOUND .AND. IPOS.LT.IEND) GO TO 200

C----     If have located the filename, add extension and store
          IF (FOUND) THEN
              PDBFIL = FNAME(1:IEND)
              CALL GETNAM(PDBFIL,ISTART,IEND,IFAIL)
              IF (IFAIL) GO TO 999
      
C----         Form full filename
              NFILE = NFILE + 1
              IF (NFILE.GT.MXFILE) GO TO 908
              FILRIN(NFILE) = PDBFIL(ISTART:IEND) // '.rin'

C----         Store identifying code
              LENID(NFILE) = IEND - ISTART + 1
              IF (LENID(NFILE).GT.MXLEND) THEN
                  LENID(NFILE) = MXLEND
                  FILID(NFILE) = PDBFIL(IEND - MXLEND + 1:IEND)
              ELSE
                  FILID(NFILE) = PDBFIL(ISTART:IEND)
              ENDIF
          ENDIF
      GO TO 100

C---- Close the mplot.in file
 500  CONTINUE
      CLOSE(1)
      IF (NFILE.EQ.0) GO TO 910
CHECK v.3.4-->
      NMODEL = NFILE
      TOPMOD = NFILE
CHECK v.3.4<--

      GO TO 999

C---- Fatal errors
 900  CONTINUE
      PRINT*, '**** Unable to open input file, mplot.in'
      GO TO 990
 
 902  CONTINUE
      PRINT*, '**** File error reading file mplot.in at line',
     -    NFILE + 1
      GO TO 990
 
 908  CONTINUE
      PRINT*, '**** Maximum number of files exceeded. Amend parameter ',
     -    'MXFILE =', MXFILE
      GO TO 990
 
 910  CONTINUE
      PRINT*, 'No files found in mplot.in list'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.
 
 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2<--
C**************************************************************************
C
C  SUBROUTINE PARAMS  -  Read in program parameters from parameter file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PARAMS

      INCLUDE 'tplot.inc'

      CHARACTER*1   YESNO
CHECK v.3.2-->
CHECK v.3.4-->
C      CHARACTER*80  FNAME
      CHARACTER*80  FNAME, IREC
CHECK v.3.4<--
CHECK v.3.2<--
      INTEGER       I, ICOL, INUMB, LINE
CHECK v.3.4-->
C      LOGICAL       ALLCOL, FINERR
      LOGICAL       ALLCOL, ENDCOL, FINERR
CHECK v.3.4<--
      REAL          RNUMB

C---- Open parameter file
CHECK v.3.1-->
C      OPEN(UNIT=1, FILE='procheck.prm', STATUS='OLD',
CHECK v.3.2-->
C      OPEN(UNIT=10, FILE='procheck.prm', STATUS='OLD',
      IF (NMR) THEN
          FNAME = 'procheck_nmr.prm'
      ELSE IF (ENSEMB) THEN
          FNAME = 'procheck_comp.prm'
      ELSE
          FNAME = 'procheck.prm'
      ENDIF
      OPEN(UNIT=10, FILE=FNAME, STATUS='OLD',
CHECK v.3.2<--
CHECK v.3.1<--
     -     FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST', READONLY,
     -     ERR=900)
      FINERR = .FALSE.
      LINE = 0

C---- Read in the parameters

C---- Check that have the right version number in the parameter file
CHECK v.3.2-->
C      CALL FINKEY('PROCHECK v.3.0',14,LINE,FINERR)
      IF (NMR) THEN
CHECK v.3.4-->
C          CALL FINKEY('PROCHECK-NMR. PROCHECK v.3.3',28,LINE,FINERR)
C      ELSE IF (ENSEMB) THEN
C          CALL FINKEY('PROCHECK-COMP. PROCHECK v.3.3',29,LINE,FINERR)
          CALL FINKEY('PROCHECK-NMR. PROCHECK v.3.4',28,LINE,FINERR)
      ELSE IF (ENSEMB) THEN
          CALL FINKEY('PROCHECK-COMP. PROCHECK v.3.4',29,LINE,FINERR)
CHECK v.3.4<--
      ELSE
          CALL FINKEY('PROCHECK v.3.3',14,LINE,FINERR)
          IF (FINERR) THEN
              CALL FINKEY('PROCHECK v.3',12,LINE,FINERR)
              IF (FINERR) THEN
                  GO TO 901
              ELSE
                  PRINT*, '* Warning. Parameter file not up-to-date'
                  PRINT*, '*          Defaults will be used for',
     -                ' missing items'
              ENDIF
          ENDIF
      ENDIF
CHECK v.3.2<--
      IF (FINERR) GO TO 901

C---- Find the colours key-word
      CALL FINKEY('Colours',7,LINE,FINERR)
      IF (FINERR) GO TO 990

C---- Read in all the RGB colours and corresponding colour names
CHECK v.3.4-->
C      DO 100, ICOL = 1, MXCOLR
C          LINE = LINE + 1
CCHECK v.3.1-->
CC          READ(1,*,END=902,ERR=904) (RGB(I, ICOL), I = 1, 3),
C          READ(10,*,END=902,ERR=904) (RGB(I, ICOL), I = 1, 3),
CCHECK v.3.1<--
C     -        COLNAM(ICOL)
C 100  CONTINUE
      ENDCOL = .FALSE.
      NCOLOR = MXCOLR
      DO 100, ICOL = 1, MXCOLR
          IF (.NOT.ENDCOL) THEN
              LINE = LINE + 1
              READ(10,20,END=902,ERR=904) IREC
 20           FORMAT(A)
              IF (IREC.EQ.' ') THEN
                  ENDCOL = .TRUE.
              ELSE
                  READ(IREC,*,ERR=904) (RGB(I,ICOL), I = 1, 3),
     -                COLNAM(ICOL)
                  NCOLOR = ICOL
              ENDIF
          ENDIF

C----     If have reached end of colours then insert default
          IF (ENDCOL) THEN
              RGB(1,ICOL) = 0.0
              RGB(2,ICOL) = 0.0
              RGB(3,ICOL) = 0.0
              COLNAM(ICOL) = 'WHITE'
          ENDIF
 100  CONTINUE
CHECK v.3.4<--

C---- Determine whether all plots are to be in colour

C---- Find the colour-all-plots keywords
      CALL FINKEY('Colour all plots?',17,LINE,FINERR)
      IF (FINERR) GO TO 990

C---- Determine whether all plots are to be in colour
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          ALLCOL = .TRUE.
CHECK v.3.4.3-->
      ELSE IF (YESNO.EQ.'W' .OR. YESNO.EQ.'w') THEN
          ALLCOL = .TRUE.
          WWWOUT = .TRUE.
CHECK v.3.4.3<--
      ELSE
          ALLCOL = .FALSE.
      ENDIF

C---- Determine which plots are to be produced

C---- Find the "Which plots" keyword
CHECK v.3.4-->
C      CALL FINKEY('Which plots',11,LINE,FINERR)
      IF (NMR) THEN
          CALL FINKEY('a. Geometry plots',17,LINE,FINERR)
      ELSE
          CALL FINKEY('Which plots',11,LINE,FINERR)
      ENDIF
CHECK v.3.4<--
      IF (FINERR) GO TO 990

C---- 1. Produce the main Ramachandran plot (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
 120  FORMAT(A)
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          PLOTRM = .TRUE.
      ELSE
          PLOTRM = .FALSE.
      ENDIF

C---- 2. Plot the Gly & Pro Ramachandran plots? (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          DOPLOT(1) = .TRUE.
      ELSE
          DOPLOT(1) = .FALSE.
      ENDIF

C---- 3. Plot the chi1-chi2 plots? (Y/N)
      LINE = LINE + 1
CHECK v.3.1<--
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1-->
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          DOPLOT(2) = .TRUE.
      ELSE
          DOPLOT(2) = .FALSE.
      ENDIF

C---- Read in the plot parameters

C---- Find the Ramachandran plot key-word
      CALL FINKEY('1. Ramachandran plot',20,LINE,FINERR)
      IF (FINERR) GO TO 990

C---- Ramachandran plot options
CHECK v.3.4-->
      RMSPLT = .FALSE.
CHECK v.3.4<--

CHECK v.3.3-->
C---- Show the Ramachandran plots on separate pages
      IF (NMR .OR. ENSEMB) THEN
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              RMSPLT = .TRUE.
          ELSE
              RMSPLT = .FALSE.
          ENDIF
      ENDIF

C---- If not splitting Ramachandrans, then just use standard subheading
C     for the Ramachandran plot
      IF (.NOT.RMSPLT) THEN
          RAMHED = BRCODE
          RLEN = BLEN
      ENDIF
CHECK v.3.3<--

C---- Shade in the regions of the Ramachandran plot (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          NOSHAD = .FALSE.
      ELSE
          NOSHAD = .TRUE.
      ENDIF

C---- Print letter-codes in the regions of the Ramachandran plot (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          INCLET = .TRUE.
      ELSE
          INCLET = .FALSE.
      ENDIF

C---- Draw borders round the regions of the Ramachandran plot (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          NOLINE = .FALSE.
      ELSE
          NOLINE = .TRUE.
      ENDIF

C---- Show core region of Ramachandran plot only (Y/N)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          CRONLY = .TRUE.
      ELSE
          CRONLY = .FALSE.
      ENDIF

C---- Which points are to be labelled on Ramachandran plot: (0 = disallowed,
C     1 = generous, 2 = allowed, 3 = core)
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,*,END=902,ERR=904) INUMB
      READ(10,*,END=902,ERR=904) INUMB
CHECK v.3.1<--
      LABRES = INUMB

CHECK v.3.4-->
C---- Read in new parameters specific to ensemble plots
      IF (NMR .OR. ENSEMB) THEN

C----     Get size of data points (0.0-2.0)
          LINE = LINE + 1
          READ(10,*,END=902,ERR=904) RNUMB
          RBXSIZ = RNUMB
          IF (RBXSIZ.LT.0.0) RBXSIZ = 0.0
          IF (RBXSIZ.GT.2.0) RBXSIZ = 2.0

C----     Fill in the data points (Y/N)
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              RBXFIL = .TRUE.
          ELSE
              RBXFIL = .FALSE.
          ENDIF

C----     Show model numbers inside the data points (Y/N)
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              RBXMOD = .TRUE.
          ELSE
              RBXMOD = .FALSE.
          ENDIF
      ENDIF
CHECK v.3.4<--

C---- Produce a black-and-white or colour PostScript file
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          INCOLR = .TRUE.
      ELSE
          INCOLR = .FALSE.
      ENDIF

C---- Get each of the user-defined colours for this plot
      DO 200, ICOL = 1, 7
          CALL GETCOL(COLRAM(ICOL),LINE,COLNAM,MXCOLR)
 200  CONTINUE

CHECK v.3.4.3-->
C---- Generate "publication version" of the Ramachandran plot (Y/N)
      LINE = LINE + 1
      READ(10,120,END=902,ERR=904) YESNO
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          RAMPUB = .TRUE.
      ELSE
          RAMPUB = .FALSE.
      ENDIF
CHECK v.3.4.3<--

C---- Gly & Pro Ramachandran plots

C---- Find the Gly & Pro Ramachandran plots key-word
      CALL FINKEY('2. Gly & Pro',12,LINE,FINERR)
      IF (FINERR) GO TO 990

C---- Cutoff value for labelling of residues on the plots
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,*,END=902,ERR=904) RNUMB
      READ(10,*,END=902,ERR=904) RNUMB
CHECK v.3.1<--
      LIMGP = RNUMB
      IF (LIMGP.LT.-9.9) LIMGP = -9.9
      IF (LIMGP.GT.99.9) LIMGP = 99.9

C---- Plot all 20 Ramachandran plots or just the Gly & Pro
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          ALLRAM = .TRUE.
      ELSE
          ALLRAM = .FALSE.
      ENDIF

C---- Produce a black-and-white or colour PostScript file
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          INCOLG = .TRUE.
      ELSE
          INCOLG = .FALSE.
      ENDIF

C---- Get each of the user-defined colours for this plot
      DO 300, ICOL = 1, 5
          CALL GETCOL(COLRGP(ICOL),LINE,COLNAM,MXCOLR)
 300  CONTINUE

C---- Chi1-Chi2 plots

C---- Find the Chi1-Chi2 plots key-word
      CALL FINKEY('3. Chi1-Chi2 plots',18,LINE,FINERR)
      IF (FINERR) GO TO 990

C---- Cutoff value for labelling of residues on the plots
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,*,END=902,ERR=904) RNUMB
      READ(10,*,END=902,ERR=904) RNUMB
CHECK v.3.1<--
      LIMCHI = RNUMB
      IF (LIMCHI.LT.-9.9) LIMCHI = -9.9
      IF (LIMCHI.GT.99.9) LIMCHI = 99.9

C---- Produce a black-and-white or colour PostScript file
      LINE = LINE + 1
CHECK v.3.1-->
C      READ(1,120,END=902,ERR=904) YESNO
      READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
      IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
          INCOLC = .TRUE.
      ELSE
          INCOLC = .FALSE.
      ENDIF

C---- Get each of the user-defined colours for this plot
      DO 400, ICOL = 1, 5
          CALL GETCOL(COLCHI(ICOL),LINE,COLNAM,MXCOLR)
 400  CONTINUE

C---- Find G-factor parameters
CHECK v.3.2-->
      IF (.NOT.NMR .AND. .NOT.ENSEMB) THEN
CHECK v.3.2<--
          CALL FINKEY('G-factors',9,LINE,FINERR)
          IF (FINERR) GO TO 990

C----     See whether Engh & Huber means to be used for G-factor
C         calculations
          LINE = LINE + 1
CHECK v.3.1-->
C          READ(1,120,END=902,ERR=904) YESNO
          READ(10,120,END=902,ERR=904) YESNO
CHECK v.3.1<--
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              USENGH = .TRUE.
          ELSE
              USENGH = .FALSE.
          ENDIF
CHECK v.3.2-->
      ENDIF
CHECK v.3.2<--

CHECK v.3.2-->
C---- Find whether file-handles are required or not
      CALL FINKEY('File-handles',12,LINE,FINERR)
      IF (.NOT.FINERR) THEN

C----     See whether the file-handle is required
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              WITHAN = .TRUE.
          ELSE
              WITHAN = .FALSE.
          ENDIF

C----     See whether plot filename to be printed on the plot itself
          LINE = LINE + 1
          READ(10,*,END=902,ERR=904)
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              PLABEL = .TRUE.
          ELSE
              PLABEL = .FALSE.
          ENDIF

C----     See whether all pages of same plot to be combined in one
C         paginated PostScript file
          LINE = LINE + 1
          READ(10,120,END=902,ERR=904) YESNO
          IF (YESNO.EQ.'Y' .OR. YESNO.EQ.'y') THEN
              COMBPS = .TRUE.
          ELSE
              COMBPS = .FALSE.
          ENDIF
      ENDIF
CHECK v.3.2<--

C---- If all plots are required to be in colour, then set all the
C     appropriate flags to TRUE
      IF (ALLCOL) THEN
          INCOLR = .TRUE.
          INCOLC = .TRUE.
          INCOLG = .TRUE.
      ENDIF

      GO TO 999


C---- Errors reading parameter file
900   CONTINUE
CHECK v.3.2-->
C      PRINT*, '*** ERROR Parameters file (procheck.prm) not found.'
      IF (NMR) THEN
          PRINT*, '*** ERROR Parameters file (procheck_nmr.prm) not',
     -        ' found.'
      ELSE IF (ENSEMB) THEN
          PRINT*, '*** ERROR Parameters file (procheck_comp.prm) not',
     -        ' found.'
      ELSE
          PRINT*, '*** ERROR Parameters file (procheck.prm) not found.'
      ENDIF
CHECK v.3.2<--
      GO TO 990

 901  CONTINUE
CHECK v.3.2-->
C       PRINT*, '*** Failed to find correct version number in ',
C     -        'parameters file, procheck.prm'
      IF (NMR) THEN
CHECK v.3.4-->
C          PRINT*, '*** Failed to find correct version number in ',
C     -            'parameters file, procheck_nmr.prm'
C      ELSE IF (ENSEMB) THEN
C          PRINT*, '*** Failed to find correct version number in ',
C     -            'parameters file, procheck_comp.prm'
C      ELSE
C          PRINT*, '*** Failed to find correct version number in ',
C     -            'parameters file, procheck.prm'
          PRINT*, '*** Old version of parameter file, ',
     -            'procheck_nmr.prm, found'
      ELSE IF (ENSEMB) THEN
          PRINT*, '*** Old version of parameter file, ',
     -            'procheck_comp.prm, found'
      ELSE
          PRINT*, '*** Old version of parameter file, ',
     -            'procheck.prm, found'
      ENDIF
      PRINT*, '*** Please delete the file and re-run the program'
      PRINT*, '*** ---------------------------------------------'
CHECK v.3.4<--
CHECK v.3.2<--
      GO TO 990

902   CONTINUE
CHECK v.3.2-->
C      PRINT*, '*** Premature end of parameters file, procheck.prm, ',
C     -        'encountered at line', LINE
      IF (NMR) THEN
          PRINT*, '*** Premature end of parameters file, procheck_nmr',
     -            '.prm, encountered at line', LINE
      ELSE IF (ENSEMB) THEN
          PRINT*, '*** Premature end of parameters file, procheck_comp',
     -            '.prm, encountered at line', LINE
      ELSE
          PRINT*, '*** Premature end of parameters file, procheck',
     -            '.prm, encountered at line', LINE
      ENDIF
CHECK v.3.2<--
CHECK v.3.4-->
      PRINT*, '*** Please delete the file and re-run the program'
      PRINT*, '*** ---------------------------------------------'
CHECK v.3.4<--
      GO TO 990

904   CONTINUE
CHECK v.3.2-->
C      PRINT*, '*** Error reading parameter file, procheck.prm, at',
C     -    ' line', LINE
      IF (NMR) THEN
          PRINT*, '*** Error reading parameter file, procheck_nmr.prm',
     -        ', at line', LINE
      ELSE IF (ENSEMB) THEN
          PRINT*, '*** Error reading parameter file, procheck_comp.prm',
     -        ', at line', LINE
      ELSE
          PRINT*, '*** Error reading parameter file, procheck.prm, at',
     -        ' line', LINE
      ENDIF
CHECK v.3.2<--
CHECK v.3.4-->
      PRINT*, '*** Please delete the file and re-run the program'
      PRINT*, '*** ---------------------------------------------'
CHECK v.3.4<--
      GO TO 990

C---- Close the parameter file
 990  CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE NORMGS  -  Calculate the normalization factors for the
C                        Gaussian distributions used for the bond lengths
C                        and bond angles
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE NORMGS

      INCLUDE 'tplot.inc'

      REAL          MAXSTD
      PARAMETER    (MAXSTD = 6.0)
 
      INTEGER       ISTEP, NSTEPS
      REAL          AVSCOR, CALCOV, COUNT, PROBAB, RSCORE, STDEV,
     -              TCOUNT, TOTAV2

C---- Calculate number of steps required
      AVSCOR = 0.0
      NSTEPS = 2.0 * MAXSTD / GSTEP + 1
      STDEV = -MAXSTD + GSTEP / 2.0
      TCOUNT = 0.0
      TOTAV2 = 0.0

C---- Loop for the required number of steps
      DO 500, ISTEP = 1, NSTEPS
          RSCORE = CALCOV(STDEV,GSTEP)
          PROBAB = EXP(RSCORE)
          COUNT = 1000.0 * PROBAB
          AVSCOR = AVSCOR + RSCORE * COUNT
          TCOUNT = TCOUNT + COUNT
          TOTAV2 = TOTAV2 + RSCORE * RSCORE * COUNT

C----     Increment value of standard deviation
          STDEV = STDEV + GSTEP
 500  CONTINUE

C---- Calculate average score and standard deviation
      AVSCOR = AVSCOR / TCOUNT
      TOTAV2 = TOTAV2 / TCOUNT
      STDEV = SQRT(TOTAV2 - AVSCOR * AVSCOR)
      GAUMEA = AVSCOR
      GAUSTD = STDEV

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4.3-->
C Routine GETDAT transferred to ps.f
CHECK v.3.4.3<--
C**************************************************************************
C
C  SUBROUTINE GETPTS  -  Read in the data from the .rin file for the
C                        appropriate distribution
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETPTS(DISTRB)

      INCLUDE 'tplot.inc'

      CHARACTER*1   INCHN, SECSTR
      CHARACTER*2   REGION
      CHARACTER*3   RESDUE
CHECK v.3.2-->
C      CHARACTER*5   SEQNO
CHECK v.3.4-->
C      CHARACTER*5   CNUMB, SEQNO
C      CHARACTER*79  BLNKLN, SUMLIN
      CHARACTER*5   CNUMB, NUMBER, SEQNO
      CHARACTER*79  BLNKLN, IREC, SUMLIN
CHECK v.3.4<--
CHECK v.3.2<--
CHECK v.3.2-->
C      INTEGER       IRESID, DISTRB, REGNO, REGTYP, SERNO
CHECK v.3.4-->
C      INTEGER       IFILE, IPOS, IRESID, DISTRB, MXIPOS, RLABEL, REGNO,
C     -              REGTYP, SERNO
      INTEGER       DISTRB, IERR, IFILE, IMODEL, IPAGE, IPOS, IRESID,
     -              LENSTR, MALLOW, MCORE, MENDCH, MGENER, MGLY, MONPG,
     -              MOUTSI, MPRO, MRAMPL, MRESID, MXIPOS, RLABEL, REGNO,
     -              REGTYP, SERNO
      LOGICAL       INRANG, WANTED
CHECK v.3.4<--
CHECK v.3.2<--
      REAL          CGAMB, CHI1, CHI2, CHI3, CHI4, DISULF, ENHB, MCBVAL,
     -              OMEGA, PHI, PSI, SCBVAL, ZETA
CHECK v.3.2-->
      DATA BLNKLN( 1:40) / ' |                                      ' /
      DATA BLNKLN(41:79) / '                                      |'  /
CHECK v.3.2<--

C---- Initialise variables
      IRESID = 0
CHECK v.3.4-->
      IPAGE = 0
      MALLOW = 0
      MCORE = 0
      MENDCH = 0
      MGENER = 0
      MGLY = 0
      MONPG = 0
      MOUTSI = 0
      MPRO = 0
      MRAMPL = 0
      MRESID = 0
CHECK v.3.4<--
      NALLOW = 0
      NCORE = 0
      NENDCH = 0
      NGENER = 0
      NGLY = 0
      NONPG = 0
      NOUTSI = 0
      NPRO = 0
      NRAMPL = 0
      NRESID = 0
CHECK v.3.2-->
      RLABEL = 0
CHECK v.3.2<--

CHECK v.3.2-->
C---- If this is an ensemble, open the first of the .rin files
      IF (NMR .OR. ENSEMB) THEN
CHECK v.3.4-->
C          IFILE = 1
          IFILE = FFILE
CHECK v.3.4<--
          OPEN(UNIT=3,FILE=FILRIN(IFILE),STATUS='OLD',
     -        ACCESS='SEQUENTIAL',
CVAX     -     READONLY,
     -        FORM='FORMATTED',ERR=900)

C---- Otherwise rewind the single .rin file
      ELSE
CHECK v.3.2<--
          REWIND(3)
CHECK v.3.2-->
      ENDIF
CHECK v.3.2<--

CHECK v.3.4-->
C---- Check whether the first record contains the actual model number
      READ(3,120,END=900) IREC
 120  FORMAT(A)
      IF (IREC(1:5).EQ.'MODEL') THEN
          READ(IREC,160,IOSTAT=IERR) IMODEL
 160      FORMAT(9X,I5)
          IF (IERR.NE.0) THEN
              IMODEL = IFILE
          ENDIF
      ELSE
          IMODEL = IFILE
          BACKSPACE(3)
      ENDIF
CHECK v.3.4<--

C---- Initialise Ramachandran plot
      IF (PLOTRM) THEN
CHECK v.3.4-->
C----     If plotting a separate Ramachandran plot for each model,
C         prepare header information for the next .rin file
          IF (RMSPLT) THEN
              WRITE(NUMBER,520) IFILE
 520          FORMAT(I5)
              IPOS = 1
              IF (NUMBER(1:1).EQ.' ') IPOS = 2
              IF (NUMBER(2:2).EQ.' ') IPOS = 3
              IF (NUMBER(3:3).EQ.' ') IPOS = 4
              IF (NUMBER(4:4).EQ.' ') IPOS = 5
              RAMHED = RAMHED(1:RLEN - 6) // NUMBER(IPOS:) // ')'
              IF (HAVRAN .AND. RESFRM(1).NE.'*ALL  ') THEN
                  RAMHED = RAMHED(1:RLEN - IPOS + 1) // '**'
              ENDIF
          ELSE
              RAMHED = BRCODE
              RLEN = BLEN
          ENDIF
CHECK v.3.4.3-->
          IF (RAMPUB) THEN
              RAMHED = ' '
              RLEN = 1
          ENDIF
CHECK v.3.4.3<--

C----     Open the Ramachandran file
          CALL RAMOPE
CHECK v.3.4<--

C----     Plot the Ramachandran file headings
CHECK v.3.4-->
C          CALL RAMACH
          CALL RAMACH(IPAGE)
CHECK v.3.4<--
      ENDIF

C---- Read through the residues file
 100  CONTINUE
          READ(3,110,END=500,ERR=902) SERNO, RESDUE, INCHN, SEQNO,
     -        SECSTR, PHI, PSI, OMEGA, CHI1, CHI2, CHI3, CHI4, ENHB,
     -        DISULF, ZETA, CGAMB, MCBVAL, SCBVAL
 110      FORMAT(I4,A3,1X,A1,A5,A1,11F7.2,2F7.3)

CHECK v.3.4-->
C----     Test whether the current residue is within one of the
C         user-defined ranges
          WANTED = INRANG(SEQNO,INCHN,RESFRM,RESTO,MAXRNG,NRANGE)
CHECK v.3.4<--

C----     Only process this residue if it belongs to the required chain
C         or residue is within one of the selected ranges
CHECK v.3.4-->
C          IF (CHAIN.EQ.' ' .OR. INCHN.EQ.CHAIN) THEN
          IF (WANTED .AND. (CHAIN.EQ.' ' .OR. INCHN.EQ.CHAIN)) THEN
CHECK v.3.4<--

C----         Increment residue-count
              IRESID = IRESID + 1
              IF (IRESID.GT.MXRES) GO TO 904
CHECK v.3.4-->
              MRESID = MRESID + 1
CHECK v.3.4<--

C----         Save the serial number of this residue
              SAVRES(IRESID) = SERNO

CHECK v.3.2-->
C----         For ensemble structures save the model number
              MODEL(IRESID) = IFILE
CHECK v.3.2<--

C----         Chi-1 torsion angle
              IF (RESDUE.EQ.'PRO') CHI1 = 999.9
              IF (CHI1.LT.0.0) CHI1 = CHI1 + 360.0
              IF (CHI1.EQ.360.0) CHI1 = 0.0

C----         Chi-2 torsion angle
              IF (CHI2.LE.0.0) CHI2 = CHI2 + 360.0
              IF (CHI2.EQ.360.0) CHI2 = 0.0

C----         Chi-3 torsion angle
              IF (CHI3.LE.0.0) CHI3 = CHI3 + 360.0
              IF (CHI3.EQ.360.0) CHI3 = 0.0

C----         Chi-4 torsion angle
              IF (CHI4.LE.0.0) CHI4 = CHI4 + 360.0
              IF (CHI4.EQ.360.0) CHI4 = 0.0

C----         Omega torsion angle
              IF (OMEGA.LE.0.0) OMEGA = OMEGA + 360.0
              IF (OMEGA.EQ.360.0) OMEGA = 0.0

C----         Zeta torsion angle
              IF (ZETA.LE.0.0) ZETA = ZETA + 360.0
              IF (ZETA.EQ.360.0) ZETA = 0.0

C----         Store required value(s)
              IF (DISTRB.EQ.1) THEN
                  VALUE(1,IRESID) = PHI
                  VALUE(2,IRESID) = PSI
              ELSE IF (DISTRB.EQ.2) THEN
                  VALUE(1,IRESID) = CHI1
                  VALUE(2,IRESID) = CHI2
              ELSE IF (DISTRB.EQ.3) THEN
                  VALUE(1,IRESID) = CHI1
              ELSE IF (DISTRB.EQ.4) THEN
                  VALUE(1,IRESID) = CHI3
              ELSE IF (DISTRB.EQ.5) THEN
                  VALUE(1,IRESID) = CHI4
              ELSE IF (DISTRB.EQ.6) THEN
                  VALUE(1,IRESID) = OMEGA
              ENDIF

C----         Save the residue name
              RESID(IRESID) = INCHN // ' ' // RESDUE // SEQNO

C----         Determine which region of Ramachandran plot residue lies in
              CALL RAMREG(PHI,PSI,REGION,REGNO,REGTYP)

C----         Increment count of different types of residues
              IF (RESDUE.EQ.'GLY') NGLY = NGLY + 1
              IF (RESDUE.EQ.'PRO') NPRO = NPRO + 1
CHECK v.3.4-->
              IF (RESDUE.EQ.'GLY') MGLY = MGLY + 1
              IF (RESDUE.EQ.'PRO') MPRO = MPRO + 1
CHECK v.3.4<--

C----         If residue is neither a glycine nor a proline, then determine
C             in which region of the Ramachandran plot its phi-psi values
C             fall
              IF (RESDUE.NE.'GLY' .AND. RESDUE.NE.'PRO') THEN

C----             Check whether this is an end-residue at a chain break
                  IF (PHI.EQ.999.9 .OR. PSI.EQ.999.9) THEN
                      NENDCH = NENDCH + 1
CHECK v.3.4-->
                      MENDCH = MENDCH + 1
CHECK v.3.4<--

C----             If valid phi-psi, then determine region type
                  ELSE
                      NONPG = NONPG + 1
CHECK v.3.4-->
                      MONPG = MONPG + 1
CHECK v.3.4<--

C----                 Increment counts of residues in each region type
                      IF (REGION.EQ.'XX') THEN
                          NOUTSI = NOUTSI + 1
CHECK v.3.4-->
                          MOUTSI = MOUTSI + 1
CHECK v.3.4<--
                          DISALL = .TRUE.
                      ELSE IF (REGION.EQ.'B ' .OR. REGION.EQ.'A ' .OR.
     -                    REGION.EQ.'L ') THEN
                          NCORE = NCORE + 1
CHECK v.3.4-->
                          MCORE = MCORE + 1
CHECK v.3.4<--
                      ELSE IF (REGION.EQ.'b ' .OR. REGION.EQ.'a ' .OR.
     -                    REGION.EQ.'l ' .OR. REGION.EQ.'p ') THEN
                          NALLOW = NALLOW + 1
CHECK v.3.4-->
                          MALLOW = MALLOW + 1
CHECK v.3.4<--
                      ELSE IF (REGION.EQ.'~b' .OR. REGION.EQ.'~a' .OR.
     -                    REGION.EQ.'~l' .OR. REGION.EQ.'~p') THEN
                          NGENER = NGENER + 1
CHECK v.3.4-->
                          MGENER = MGENER + 1
CHECK v.3.4<--
                      ENDIF
                  ENDIF
              ENDIF

C----         Plot this residue on the Ramachandran plot
              IF (PLOTRM .AND. PHI.NE.999.9 .AND. PSI.NE.999.9) THEN
                  NRAMPL = NRAMPL + 1
CHECK v.3.4-->
                  MRAMPL = MRAMPL + 1
CHECK v.3.4<--
CHECK v.3.2-->
C                  CALL RAMPLT(PHI,PSI,RESDUE,REGTYP,SEQNO)
                  CALL RAMPLT(PHI,PSI,RESDUE,REGTYP,SEQNO,IFILE,RLABEL,
     -                INCHN)
CHECK v.3.2<--
              ENDIF

          ENDIF

C---- Loop back for next record in file
      GO TO 100

C---- End of file reached
 500  CONTINUE
      NRESID = IRESID

C---- Check that file wasn't empty
      IF (NRESID.EQ.0) GO TO 906
CHECK v.3.2-->
      IF (.NOT.NMR .AND. .NOT.ENSEMB) THEN
CHECK v.3.2<--
          PRINT*, NRESID, ' residues read in'
CHECK v.3.2-->
      ENDIF
CHECK v.3.2<--

CHECK v.3.2-->
C---- If this is an ensemble, compute statistics for model just processed
C     and then get the next .rin file
      IF (NMR .OR. ENSEMB) THEN

C----     Close the .rin file
          CLOSE(3)

CHECK v.3.4-->
C----     If plotting a separate Ramachandran plot for each model,
C         then close off this one and open the next
          IF (PLOTRM .AND. RMSPLT) THEN

C----         Calculate Ramachandran plot statistics for last model
CHECK v.3.4.3-->
C              CALL MSTATS(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
C     -            MPRO,MRAMPL,MRESID)
              CALL MSTATS(MALLOW,MCORE,MGENER,MONPG,MOUTSI)
CHECK v.3.4.3<--

              IF (ENSEMB) THEN
                  CALL FILKEY(XRAMC2 + 4.0 * HWIDX,YRAMC2 - 2.0 * HWIDY,
     -                2.0 * HWIDX,2.0 * HWIDY,8.0,5.0)
              ENDIF

C----         Print statistics and close file
CHECK v.3.4.3-->
C              CALL RAMEND(BBOXX1,BBOXX2,BBOXY1,BBOXY2,MALLOW,MCORE,
C     -            MENDCH,MGENER,MGLY,MONPG,MOUTSI,MPRO,MRAMPL,MRESID)
              CALL RAMEND(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,
     -            MOUTSI,MPRO,RAMPUB)
CHECK v.3.4.3<--

C----         Reinitialise the counts
              MALLOW = 0
              MCORE = 0
              MENDCH = 0
              MGENER = 0
              MGLY = 0
              MONPG = 0
              MOUTSI = 0
              MPRO = 0
              MRAMPL = 0
              MRESID = 0
          ENDIF
CHECK v.3.4<--

C----     If there's another .rin file to go, open it
CHECK v.3.4-->
 550      CONTINUE
CHECK v.3.4<--
          IF (IFILE.LT.NFILE) THEN
              IFILE = IFILE + 1
CHECK v.3.4-->
              IF (MWANT(IFILE)) THEN
CHECK v.3.4<--
                  OPEN(UNIT=3,FILE=FILRIN(IFILE),STATUS='OLD',
     -                ACCESS='SEQUENTIAL',
CVAX     -             READONLY,
     -                FORM='FORMATTED',ERR=900)
CHECK v.3.4-->
C                  GO TO 100

C----             If plotting a separate Ramachandran plot for each model,
C                 prepare header information for the next .rin file
                  IF (PLOTRM .AND. RMSPLT) THEN
                      WRITE(NUMBER,520) IFILE
                      IPOS = 1
                      IF (NUMBER(1:1).EQ.' ') IPOS = 2
                      IF (NUMBER(2:2).EQ.' ') IPOS = 3
                      IF (NUMBER(3:3).EQ.' ') IPOS = 4
                      IF (NUMBER(4:4).EQ.' ') IPOS = 5
                      RAMHED = RAMHED(1:RLEN - 6) // NUMBER(IPOS:) //
     -                    ')'
                      IF (HAVRAN .AND. RESFRM(1).NE.'*ALL  ') THEN
                          RAMHED = RAMHED(1:RLEN - IPOS + 1) // '**'
                      ENDIF
                      CALL RAMACH(IPAGE)
                  ENDIF
                  GO TO 100
              ELSE
                  GO TO 550
              ENDIF
CHECK v.3.4<--
          ENDIF
          PRINT*, NRESID, ' residues read in for entire ensemble'
      ENDIF
CHECK v.3.2<--

C---- Calculate percentages for residues in the different regions of the
C     Ramachandran plot
      IF (NONPG.GT.0) THEN
          ALLOWP = 100.0 * REAL(NALLOW) / REAL(NONPG)
          COREPC = 100.0 * REAL(NCORE) / REAL(NONPG)
          GENERP = 100.0 * REAL(NGENER) / REAL(NONPG)
          OUTSIP = 100.0 * REAL(NOUTSI) / REAL(NONPG)
      ELSE
          ALLOWP = 0.0
          COREPC = 0.0
          GENERP = 0.0
          OUTSIP = 0.0
      ENDIF

C---- Print statistics for Ramachandran plot
      IF (PLOTRM) THEN
CHECK v.3.2-->
C----     If this is an ensemble of structures from separate files
C         (ie rather than a set of NMR models in a single file), then
C         print key
          IF (ENSEMB) THEN
              CALL FILKEY(XRAMC2 + 4.0 * HWIDX,YRAMC2 - 2.0 * HWIDY,
     -            2.0 * HWIDX,2.0 * HWIDY,8.0,5.0)
          ENDIF
CHECK v.3.2<--

C----     Print statistics and close file
CHECK v.3.4-->
C          CALL RAMEND(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
          IF (.NOT.RMSPLT) THEN
CHECK v.3.4.3-->
C              CALL RAMEND(BBOXX1,BBOXX2,BBOXY1,BBOXY2,NALLOW,NCORE,
C     -            NENDCH,NGENER,NGLY,NONPG,NOUTSI,NPRO,NRAMPL,NRESID)
              CALL RAMEND(NALLOW,NCORE,NENDCH,NGENER,NGLY,NONPG,
     -            NOUTSI,NPRO,RAMPUB)
CHECK v.3.4.3<--
          ENDIF

C----     Close the Ramachandran plot file
          CALL PSCLOS(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
CHECK v.3.4<--
      ENDIF

C---- Set Ramachandran plot flag off so that don't plot this again
      PLOTRM = .FALSE.

CHECK v.3.2-->
C---- If this is the first pass through this routine, write out the
C     Ramachandran statistics to the .sum summary file
      IF (DISTRB.EQ.1) THEN

C----     File name, resolution, chain-ID and number of residues

C----     Blank out the output line
          SUMLIN = BLNKLN

C----     Determine length of filename of original PDB file
          IPOS = 79
 600      CONTINUE
              IPOS = IPOS - 1
          IF (IPOS.GT.1 .AND. PDBFIL(IPOS:IPOS).EQ.' ') GO TO 600

C----     Determine the maximum length of filename that we can fit in
          IF (CHAIN.EQ.' ') THEN
              MXIPOS = 54
          ELSE
              MXIPOS = 51
          ENDIF
          IF (IPOS.GT.MXIPOS) IPOS = MXIPOS

C----     Fill the output line with the relevant data
          WRITE(CNUMB,620) RESOL
 620      FORMAT(F5.1)
          SUMLIN(4:4 + IPOS - 1) = PDBFIL(1:IPOS)
          SUMLIN(4 + IPOS + 1:4 + IPOS + 5) = CNUMB
          IF (CHAIN.NE.' ') SUMLIN(4 + IPOS + 8:4 + IPOS + 8) = CHAIN
          WRITE(CNUMB,640) NRESID
 640      FORMAT(I5)
          SUMLIN(64:68) = CNUMB
          SUMLIN(70:77) = 'residues'

C----     Write the line out to the summary file
          WRITE(14,660) SUMLIN
 660      FORMAT(A)
          WRITE(14,660) BLNKLN

C----     Ramachandran plot statisticsm

C----     Blank out the output line
          SUMLIN = BLNKLN

C----     Fill the line with the percentages in the various regions of the
C         Ramachandran plot
          SUMLIN(4:77) = 'Ramachandran plot:       % core       % al' //
     -        'low       % gener       % disall'
          WRITE(CNUMB,680) COREPC
 680      FORMAT(F5.1)
          SUMLIN(24:28) = CNUMB
          WRITE(CNUMB,680) ALLOWP
          SUMLIN(37:41) = CNUMB
          WRITE(CNUMB,680) GENERP
          SUMLIN(51:55) = CNUMB
          WRITE(CNUMB,680) OUTSIP
          SUMLIN(65:69) = CNUMB

C----     Determine whether to highlight the line with a + or a *
          IF (GENERP.GT.0.0 .OR. COREPC.LT.90.0) THEN
              SUMLIN(1:1) = '+'
          ENDIF
          IF (OUTSIP.GT.0.0 .OR. COREPC.LT.80.0) THEN
              SUMLIN(1:1) = '*'
          ENDIF

C----     Write the line out to the summary file
          WRITE(14,660) SUMLIN
          WRITE(14,660) BLNKLN
      ENDIF
CHECK v.3.2<--

      GO TO 999

C---- Fatal errors
CHECK v.3.0.1-->
C900   CONTINUE
C      PRINT*, '*** ERROR. Unable to open data file'
C      GO TO 990
CHECK v.3.0.1<--

CHECK v.3.2-->
 900  CONTINUE
CHECK v.3.4-->
C      PRINT*, '*** Error opening .rin file.', FILRIN(IFILE)
      PRINT*, '*** Error opening .rin file.'
      PRINT*, '* [', FILRIN(IFILE)(1:LENSTR(FILRIN(IFILE))), ']'
CHECK v.3.4<--
      GO TO 990
CHECK v.3.2<--

902   CONTINUE
      PRINT*, '*** ERROR. Data error in .rin file at line:', IRESID + 1
      GO TO 990

 904  CONTINUE
      PRINT*, '*** ERROR. Maximum number of residues exceeded:', MXRES
      GO TO 990

 906  CONTINUE
      PRINT*, '*** Warning. No data points point in .rin file.'
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5.2-->
C**************************************************************************
C
C  SUBROUTINE GETCAD  -  Read in the CA angle data from the .ca file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCAD(DISTRB)

      INCLUDE 'tplot.inc'

      INTEGER       DISTRB

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5.2<--
CHECK v.3.4-->
C**************************************************************************
C
C  SUBROUTINE RAMOPE  -  Open the Ramachandran plots file
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMOPE

      INCLUDE 'tplot.inc'

      CHARACTER*30  PLDESC
      CHARACTER*60  PTITLE

C---- Open PostScript file
      PLDESC = 'Main Ramachandran plot'
      CALL PSNAME(FILPS,PSLEN,IPLOT,PLDESC,'ramachand',WITHAN)
      PTITLE = 'Ramachandran plot'
CHECK v.3.4-->
C      CALL PSOPEN(FILPS,MXCOLR,RGB,INCOLR,PTITLE,1)
      CALL PSOPEN(FILPS,MXCOLR,RGB,NCOLOR,INCOLR,PTITLE,1)
CHECK v.3.4<--

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4<--
C**************************************************************************
C
C  SUBROUTINE RAMACH  -  Print headings, etc for Ramachandran plot
C
C----------------------------------------------------------------------+---

CHECK v.3.4-->
C      SUBROUTINE RAMACH
      SUBROUTINE RAMACH(IPAGE)
CHECK v.3.4<--

      INCLUDE 'tplot.inc'

CHECK v.3.4-->
      CHARACTER*7   PAGE
CHECK v.3.4<--
CHECK v.3.1-->
CHECK v.3.4-->
C      CHARACTER*30  PLDESC
CHECK v.3.4<--
CHECK v.3.1<--
CHECK v.3.2-->
CHECK v.3.4-->
C      CHARACTER*60  PTITLE
CHECK v.3.4<--
CHECK v.3.2<--
CHECK v.3.4-->
CHECK v.3.4.3-->
C      INTEGER       IPAGE
      INTEGER       BAKCOL, IPAGE
CHECK v.3.4.3<--
CHECK v.3.4<--
CHECK v.3.4.3-->
      LOGICAL       HLABEL
CHECK v.3.4.3<--

C---- Open PostScript file
CHECK v.3.4-->
CCHECK v.3.1-->
CC      CALL PSNAME(FILPS,PSLEN,IPLOT)
C      PLDESC = 'Main Ramachandran plot'
CCHECK v.3.2-->
CC      CALL PSNAME(FILPS,PSLEN,IPLOT,PLDESC)
C      CALL PSNAME(FILPS,PSLEN,IPLOT,PLDESC,'ramachand',WITHAN)
CCHECK v.3.2<--
CCHECK v.3.1<--
CCHECK v.3.2-->
CC      CALL PSOPEN(FILPS,BBOXX1,BBOXX2,BBOXY1,BBOXY2,MXCOLR,RGB,
CC     -    INCOLR,COLRAM(1))
C      PTITLE = 'Ramachandran plot'
C      CALL PSOPEN(FILPS,MXCOLR,RGB,INCOLR,PTITLE,1)
CHECK v.3.4<--
CHECK v.3.4.3-->
      BAKCOL = COLRAM(1)
      HLABEL = PLABEL
      IF (RAMPUB) THEN
          BAKCOL = -COLRAM(1)
          HLABEL = .FALSE.
      ENDIF
CHECK v.3.4.3<--
      CALL PSPAGE(FILPS,BBOXX1,BBOXX2,BBOXY1,BBOXY2,MXCOLR,RGB,
CHECK v.3.4-->
C     -    INCOLR,COLRAM(1),1,PLABEL)
CHECK v.3.4.3-->
C     -    INCOLR,COLRAM(1),1,PLABEL,RSELEC)
     -    INCOLR,BAKCOL,IPAGE + 1,HLABEL,RSELEC)
CHECK v.3.4.3<--
CHECK v.3.4<--

C---- Print program name in top left-hand corner
CHECK v.3.4.3-->
      IF (.NOT.RAMPUB) THEN
CHECK v.3.4.3<--
          IF (NMR) THEN
              CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -            'PROCHECK-NMR')
          ELSE IF (ENSEMB) THEN
              CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -            'PROCHECK-COMP')
          ELSE
              CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -            'PROCHECK')
          ENDIF
CHECK v.3.4.3-->
      ENDIF
CHECK v.3.4.3<--
CHECK v.3.2<--

CHECK v.3.4-->
C---- If splitting Ramachandran plots, print page number
      IF (RMSPLT) THEN
          IPAGE = IPAGE + 1
          WRITE(PAGE,120) IPAGE
 120      FORMAT('Page ',I2)
          CALL PSCTXT(BBOXX2 - 40.0,BBOXY2 - 20.0,12.0,PAGE)
      ENDIF
CHECK v.3.4<--

C---- Shade in the different regions of the Ramachandran plot
CHECK v.3.4.3-->
      IF (RAMPUB) THEN
          CALL PSTRAN(0.0,-80.0)
      ENDIF
CHECK v.3.4.3<--
      CALL RAMSHD

C---- Generate axes and labels for Ramachandran plot
      CALL RAMAXE

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4-->
C**************************************************************************
C
C  SUBROUTINE MSTATS  -  Calculate Ramachandran plot statistics for current
C                        model only
C
C----------------------------------------------------------------------+---

CHECK v.3.4.3-->
C      SUBROUTINE MSTATS(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
C     -    MPRO,MRAMPL,MRESID)
      SUBROUTINE MSTATS(MALLOW,MCORE,MGENER,MONPG,MOUTSI)
CHECK v.3.4.3<--

      INCLUDE 'tplot.inc'

CHECK v.3.4.3-->
C      INTEGER       MALLOW, MCORE, MENDCH, MGENER, MGLY, MONPG, MOUTSI,
C     -              MPRO, MRAMPL, MRESID
      INTEGER       MALLOW, MCORE, MGENER, MONPG, MOUTSI
CHECK v.3.4.3<--

C---- Calculate percentages for residues in the different regions of the
C     Ramachandran plot
      IF (NONPG.GT.0) THEN
          ALLOWP = 100.0 * REAL(MALLOW) / REAL(MONPG)
          COREPC = 100.0 * REAL(MCORE) / REAL(MONPG)
          GENERP = 100.0 * REAL(MGENER) / REAL(MONPG)
          OUTSIP = 100.0 * REAL(MOUTSI) / REAL(MONPG)
      ELSE
          ALLOWP = 0.0
          COREPC = 0.0
          GENERP = 0.0
          OUTSIP = 0.0
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4<--
C**************************************************************************
C
C  SUBROUTINE RAMSHD  -  Mark the different Ramachandran plot regions, using
C                        lines/borders and lettering, as required
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMSHD

      INCLUDE 'tplot.inc'

      INTEGER        MAPSIZ
      PARAMETER     (MAPSIZ = 36)

      REAL           PGAP, XGAP, YGAP

C---- Write out background shade for the plot
      CALL PSLWID(0.1)
      IF (.NOT.NOSHAD) THEN
          CALL PSHADE(0.975,COLRAM(2),RGB,MXCOLR,INCOLR)
          CALL PSUBOX(XRAMC1,YRAMC1,XRAMC1,YRAMC2,XRAMC2,YRAMC2,
     -        XRAMC2,YRAMC1)
      ENDIF

C---- Initialise variables for shading
      PGAP = 360.0 / MAPSIZ
      XGAP = (XRAMC2 - XRAMC1) / MAPSIZ
      YGAP = (YRAMC2 - YRAMC1) / MAPSIZ

C---- Shade in the different regions, if required
      CALL RAMFIL(PGAP,XGAP,YGAP)

C---- If borders round regions required, then draw them in
      IF (.NOT.NOLINE) THEN
          CALL RAMLIN(PGAP,XGAP,YGAP)
      ENDIF

C---- Print labels for each of the major regions, if required
      IF (INCLET) THEN
          CALL RAMLET(XGAP,YGAP)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMFIL  -  Shade in the different Ramachandran plot regions
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMFIL(PGAP,XGAP,YGAP)

      INCLUDE 'tplot.inc'

      INTEGER        MAPSIZ
      PARAMETER     (MAPSIZ = 36)

      CHARACTER*2    REGION
      INTEGER        I, J, IREG, LSTREG, RCOL(12), REGTYP, RTYPE
      REAL           COLOUR(3,12), FACTOR(12), PGAP, PHI, PSI,
     -               SHADE(12), X, XGAP, XLEFT, Y, YGAP, YLEFT

      DATA FACTOR  / 1.0, 0.95, 1.0, 0.95,0.95, 1.0, 1.0, 1.0,
     -               1.0,0.95,0.95,0.95 /
      DATA RCOL    /   2,   5,   4,    5,   4,   5,   4,   4,
     -                 3,   3,   3,   3 /
      DATA SHADE   / 1.0, 0.45, 0.7, 0.45, 0.6, 0.5, 0.7, 0.7,
     -               0.9, 0.8, 0.8, 0.8 /


C---- If no shading required, then blank out all the different
C     shade-codes
      DO 100, I = 1, 12
          IF (NOSHAD) THEN
              COLOUR(1,I) = 1.0
              COLOUR(2,I) = 1.0
              COLOUR(3,I) = 1.0
              SHADE(I) = 1.0
          ELSE
              COLOUR(1,I) = FACTOR(I) * RGB(1,COLRAM(RCOL(I)))
              COLOUR(2,I) = FACTOR(I) * RGB(2,COLRAM(RCOL(I)))
              COLOUR(3,I) = FACTOR(I) * RGB(3,COLRAM(RCOL(I)))
          ENDIF
 100  CONTINUE

C---- Loop through all the map-points
      CALL PSLWID(0.0)
      RTYPE = 0
      Y = YRAMC1
      PSI = -180.0 + PGAP / 2.0
      DO 300, I = 1, MAPSIZ
          X = XRAMC1
          XLEFT = X
          YLEFT = Y
          LSTREG = 0
          PHI = -180.0 + PGAP / 2.0
          DO 200, J = 1, MAPSIZ

C----         Determine which region this point is in
              CALL RAMREG(PHI,PSI,REGION,IREG,REGTYP)

C----         If region is a new one, write out strip corresponding to
C             previous region
              IF (LSTREG.NE.0 .AND. IREG.NE.LSTREG) THEN

C----             Write out the appropriate shade
                  IF (LSTREG.GT.1) THEN
                      IF (.NOT.CRONLY .OR. RTYPE.EQ.4) THEN
                          IF (INCOLR) THEN
                              CALL PSCOLR(COLOUR(1,LSTREG),
     -                            COLOUR(2,LSTREG),COLOUR(3,LSTREG))
                          ELSE
                              CALL PSHADE(SHADE(LSTREG),0,RGB,MXCOLR,
     -                            INCOLR)
                          ENDIF
                          CALL PSUBOX(XLEFT, Y, X, Y, X,
     -                        Y + YGAP, XLEFT, Y + YGAP)
                      ENDIF
                  ENDIF
                  XLEFT = X
              ENDIF
              LSTREG = IREG
              RTYPE = REGTYP
              X = X + XGAP
              PHI = PHI + PGAP
 200      CONTINUE

C----     Write out strip corresponding to last region, if required
          IF (LSTREG.GT.1) THEN
              IF (.NOT.CRONLY .OR. RTYPE.EQ.4) THEN
                  IF (INCOLR) THEN
                      CALL PSCOLR(COLOUR(1,LSTREG),
     -                    COLOUR(2,LSTREG),COLOUR(3,LSTREG))
                  ELSE
                      CALL PSHADE(SHADE(LSTREG),0,RGB,MXCOLR,
     -                    INCOLR)
                  ENDIF
                  CALL PSUBOX(XLEFT, Y, X, Y, X,
     -                Y + YGAP, XLEFT, Y + YGAP)
              ENDIF
          ENDIF

C----     Increment y-value
          Y = Y + YGAP
          PSI = PSI + PGAP
 300  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMLIN  -  Draw in the lines separating the different
C                        Ramachandran plot regions
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMLIN(PGAP,XGAP,YGAP)

      INCLUDE 'tplot.inc'

      INTEGER        MAPSIZ, NHATCH
      PARAMETER     (MAPSIZ = 36, NHATCH=50)

      CHARACTER*2    REGION
      INTEGER        I, J, IREG, ITYPE, LSTREG, LTYPE(MAPSIZ), REGTYP
      REAL           PGAP, PHI, PSI, X, XGAP, Y, YGAP

C---- Initialise region-types (ie core, allowed, generous, disallowed)
      DO 100, I = 1, MAPSIZ
          LTYPE(I) = 0
 100  CONTINUE

C---- Loop through all the map-points
      IF (NOSHAD) THEN
          CALL PSLWID(0.2)
      ELSE
          CALL PSLWID(0.8)
      ENDIF
      Y = YRAMC1
      PSI = -180.0 + PGAP / 2.0
      DO 300, I = 1, MAPSIZ
          X = XRAMC1
          LSTREG = 0
          PHI = -180.0 + PGAP / 2.0
          DO 200, J = 1, MAPSIZ

C----         Determine which region this point is in
              CALL RAMREG(PHI,PSI,REGION,IREG,REGTYP)

C----         Determine region-type
              ITYPE = REGTYP

C----         If region is a new one, then may need to draw a line between
C             it and the last
              IF (J.GT.1 .AND. ITYPE.NE.LTYPE(J - 1)) THEN

C----             Write out the appropriate line
                  IF (NOSHAD .OR. ITYPE.LT.LTYPE(J - 1)) THEN
                      IF (.NOT.CRONLY .OR. (ITYPE.EQ.4 .OR.
     -                    LTYPE(J - 1).EQ.4)) THEN
                          CALL PSLINE(X, Y, X, Y + YGAP)
                      ENDIF
                  ENDIF
              ENDIF

C----         Check whether to draw line between region and one above it
              IF (I.GT.1) THEN
                  IF (ITYPE.NE.LTYPE(J)) THEN
                      IF (NOSHAD .OR. ITYPE.GT.LTYPE(J)) THEN
                          IF (.NOT.CRONLY .OR. (ITYPE.EQ.4 .OR.
     -                        LTYPE(J).EQ.4)) THEN
                              CALL PSLINE(X, Y, X + XGAP, Y)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF

C----         Store current details and prepare for next square
              LSTREG = IREG
              LTYPE(J) = ITYPE
              X = X + XGAP
              PHI = PHI + PGAP
 200      CONTINUE

C----     Increment y-value
          Y = Y + YGAP
          PSI = PSI + PGAP
 300  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMLET  -  Write in the letter-codes identifying the
C                        different Ramachandran plot regions
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMLET(XGAP,YGAP)

      INCLUDE 'tplot.inc'

      REAL          XGAP, YGAP

C---- Print labels for each of the major regions

C---- Core regions
      CALL PSCTXT(XRAMC1 +  2.0 * XGAP,YRAMC2 -  1.0 * YGAP,15.0,'B')
      CALL PSCTXT(XRAMC1 +  7.0 * XGAP,YRAMC2 - 17.0 * YGAP,15.0,'A')
      CALL PSCTXT(XRAMC1 + 23.5 * XGAP,YRAMC2 - 13.5 * YGAP,10.0,'L')

C---- Other regions, of required
      IF (CRONLY) GO TO 999
      CALL PSCTXT(XRAMC1 +  1.0 * XGAP,YRAMC2 -  5.0 * YGAP,15.0,'b')
      CALL PSCTXT(XRAMC1 +  2.0 * XGAP,YRAMC2 - 15.0 * YGAP,15.0,'a')
      CALL PSCTXT(XRAMC1 + 21.5 * XGAP,YRAMC2 - 10.0 * YGAP,10.0,'l')
      CALL PSCTXT(XRAMC1 + 23.5 * XGAP,YRAMC2 - 33.0 * YGAP,10.0,'p')
      CALL PSCTXT(XRAMC1 + 21.5 * XGAP,YRAMC2 - 31.0 * YGAP,10.0,'~p')
      CALL PSCTXT(XRAMC1 + 16.5 * XGAP,YRAMC2 -  5.0 * YGAP,10.0,'~b')
      CALL PSCTXT(XRAMC1 + 16.5 * XGAP,YRAMC2 - 19.5 * YGAP,10.0,'~a')
      CALL PSCTXT(XRAMC1 + 22.5 * XGAP,YRAMC2 -  7.0 * YGAP,10.0,'~l')
      CALL PSCTXT(XRAMC1 +  2.0 * XGAP,YRAMC2 - 33.0 * YGAP,10.0,'b')
      CALL PSCTXT(XRAMC1 +  0.5 * XGAP,YRAMC2 - 31.0 * YGAP,10.0,'~b')
      CALL PSCTXT(XRAMC1 + 35.0 * XGAP,YRAMC2 -  2.0 * YGAP,10.0,'b')
      CALL PSCTXT(XRAMC1 + 31.5 * XGAP,YRAMC2 -  1.0 * YGAP,10.0,'~b')
      CALL PSCTXT(XRAMC1 + 33.5 * XGAP,YRAMC2 - 34.0 * YGAP,10.0,'~b')

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMAXE  -  Draw axes and axis-labels of Ramachandran plot,
C                        writing out to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMAXE

      INCLUDE 'tplot.inc'

      REAL           XCENTR, YCENTR

C---- Initialise variables
      XCENTR = (XRAMC1 + XRAMC2) / 2.0
      YCENTR = (YRAMC1 + YRAMC2) / 2.0

C---- Draw box round graph and label axes
      CALL AXES(XRAMC1,XRAMC2,YRAMC1,YRAMC2,8,8,-180.0,180.0,
     -    -180.0,180.0,15.0,0,0,22.0,.TRUE.,.TRUE.,.FALSE.,.FALSE.,
     -    .FALSE.)

C---- Draw axes
      CALL PSLWID(0.2)
      CALL PSDASH(1)
      CALL PSLINE(XCENTR,YRAMC1,XCENTR,YRAMC2)
      CALL PSLINE(XRAMC1,YCENTR,XRAMC2,YCENTR)
      CALL PSDASH(0)

C---- Graph heading and x-axis title
CHECK v.3.2-->
C      CALL PSCTXT(XCENTR,YRAMC2 + 45.0,30.0,'Ramachandran Plot')
CHECK v.3.4.3-->
      IF (.NOT.RAMPUB) THEN
CHECK v.3.4.3<--
          IF (ENSEMB) THEN
              CALL PSCTXT(XCENTR,YRAMC2 + 45.0,30.0,
     -            'Ensemble Ramachandran Plot')
          ELSE
              CALL PSCTXT(XCENTR,YRAMC2 + 45.0,30.0,'Ramachandran Plot')
          ENDIF
CHECK v.3.4.3-->
      ENDIF
CHECK v.3.4.3<--
CHECK v.3.2<--
CHECK v.3.4-->
C      CALL PSCTXT(XCENTR,YRAMC2 + 16.0,25.0,BRCODE(1:BLEN))
      CALL PSCTXT(XCENTR,YRAMC2 + 16.0,25.0,RAMHED(1:RLEN + 2))
CHECK v.3.4<--
      CALL PSCTXT(XCENTR,YRAMC1 - 30.0,15.0,'Phi (degrees)')

C---- y-axis title
CHECK v.3.4.3-->
C      CALL PSRTXT(XRAMC1 - 48.0,YCENTR,15.0,'Psi (degrees)')
      CALL PSRCTX(XRAMC1 - 48.0,YCENTR,15.0,'Psi (degrees)')
CHECK v.3.4.3<--

C---- Set line-thickness for point-markers
      CALL PSLWID(0.2)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMPLT  -  Print point on Ramachandran plot
C
C----------------------------------------------------------------------+---

CHECK v.3.2-->
C      SUBROUTINE RAMPLT(PHI,PSI,RESDUE,RTYPE,SEQNO)
      SUBROUTINE RAMPLT(PHI,PSI,RESDUE,RTYPE,SEQNO,IMODEL,RLABEL,INCHN)
CHECK v.3.2<--

      INCLUDE 'tplot.inc'

CHECK v.3.2-->
C      CHARACTER*3   RESDUE
      CHARACTER*1   INCHN
      CHARACTER*3   NUMBER, RESDUE
CHECK v.3.2<--
      CHARACTER*5   SEQNO
CHECK v.3.2-->
C      CHARACTER*5   SEQNO
C      INTEGER       RTYPE
C      REAL          PHI, PSI, SHADE, X, Y
      CHARACTER*13  STRING
CHECK v.3.5-->
C      INTEGER       ICOLR, IMODEL, IPOS, LEN, RLABEL, RTYPE
      INTEGER       ICOLR, IMODEL, IPOS, LEN, RLABEL, RTYPE, SCOLR
CHECK v.3.5<--
      REAL          DX, DY, PHI, PSI, SHADE, TSIZE, X, Y
CHECK v.3.2<--

C---- Calculate appropriate position on graph
      X = XRAMC1 + (XRAMC2 - XRAMC1) * (PHI + 180.0) / 360.0
      Y = YRAMC1 + (YRAMC2 - YRAMC1) * (PSI + 180.0) / 360.0
CHECK v.3.4.4-->
      ICOLR = 1
      CALL PSCOLB(0.0,0.0,0.0)
CHECK v.3.4.4<--

C---- Adjust colour of marker, as appropriate
      IF (INCOLR) THEN
          IF (RTYPE.LE.LABRES + 1 .AND. RESDUE.NE.'GLY') THEN
              CALL PSCOLB(RGB(1,COLRAM(7)),RGB(2,COLRAM(7)),
     -            RGB(3,COLRAM(7)))
CHECK v.3.4-->
              ICOLR = COLRAM(7)
CHECK v.3.4<--
          ELSE
              CALL PSCOLB(RGB(1,COLRAM(6)),RGB(2,COLRAM(6)),
     -            RGB(3,COLRAM(6)))
CHECK v.3.4-->
              ICOLR = COLRAM(6)
CHECK v.3.4<--
          ENDIF
      ENDIF
CHECK v.3.5-->
      SCOLR = ICOLR
CHECK v.3.5<--

C---- Plot appropriate marker depending on whether glycine or otherwise
CHECK v.3.2-->
C      IF (RESDUE.EQ.'GLY') THEN
C          CALL PSMARK(X,Y,2,HWIDX,HWIDY)
C      ELSE
C          IF (RTYPE.EQ.1) THEN
C              CALL PSMARK(X,Y,3,HWIDX,HWIDY)
C          ELSE
C              CALL PSMARK(X,Y,1,HWIDX,HWIDY)
C          ENDIF
C      ENDIF

C---- For ensembles, print empty squares and boxes
CHECK v.3.4-->
C      IF ((ENSEMB .OR. NMR) .AND. NFILE.GT.1) THEN
C          SHADE = 1.0
C          ICOLR = COLRAM(1)
C          DX = 2.0 * HWIDX
C          DY = 2.0 * HWIDY
      IF ((ENSEMB .OR. NMR) .AND. NMODEL.GT.1) THEN
          DX = 2.0 * HWIDX * RBXSIZ
          DY = 2.0 * HWIDY * RBXSIZ
          IF (RBXFIL) THEN
              SHADE = 0.0
CHECK v.3.5-->
C              ICOLR = COLRAM(7)
              SCOLR = COLRAM(7)
CHECK v.3.5<--
          ELSE
              SHADE = 1.0
CHECK v.3.5-->
C              ICOLR = COLRAM(1)
              SCOLR = COLRAM(1)
CHECK v.3.5<--
          ENDIF
CHECK v.3.4<--

C---- For single structure, print filled-in squares and triangles
      ELSE
          SHADE = 0.0
CHECK v.3.4-->
C          ICOLR = COLRAM(7)
C          DX = HWIDX
C          DY = HWIDY
          DX = HWIDX * RBXSIZ
          DY = HWIDY * RBXSIZ
CHECK v.3.4<--
CHECK v.3.4.3-->
          IF (INCOLR) CALL PSCOLB(0.0,0.0,0.0)
CHECK v.3.4.3<--
      ENDIF

C---- Plot the marker
      CALL PSLWID(0.2)
CHECK v.3.5-->
C      CALL PSHADE(SHADE,ICOLR,RGB,MXCOLR,INCOLR)
      CALL PSHADE(SHADE,SCOLR,RGB,MXCOLR,INCOLR)
CHECK v.3.5<--
      IF (RESDUE.EQ.'GLY') THEN
          CALL PSTRIA(X - DX,Y - 0.6 * DY,
     -        X + DX,Y - 0.6 * DY,X,Y + 1.2 * DY)
      ELSE
          CALL PSBBOX(X - DX,Y - DY,X - DX,Y + DY,
     -                X + DX,Y + DY,X + DX,Y - DY)
      ENDIF

C---- For ensembles, print the model number inside the marker
CHECK v.3.4-->
C      IF ((ENSEMB .OR. NMR) .AND. NFILE.GT.1) THEN
      IF ((ENSEMB .OR. NMR) .AND. RBXMOD .AND. NMODEL.GT.1) THEN
CHECK v.3.4<--
          WRITE(NUMBER,110) IMODEL
 110      FORMAT(I3)
          IF (NUMBER(2:2).EQ.' ') THEN
              IPOS = 3
          ELSE IF (NUMBER(1:1).EQ.' ') THEN
              IPOS = 2
          ELSE
              IPOS = 1
          ENDIF
CHECK v.3.4-->
C          IF (NFILE.LT.10) THEN
          IF (TOPMOD.LT.10) THEN
CHECK v.3.4<--
              TSIZE = 8.0
          ELSE
              TSIZE = 5.0
          ENDIF
          IF (RESDUE.EQ.'GLY') TSIZE = TSIZE * 0.75
CHECK v.3.4-->
          TSIZE = TSIZE * RBXSIZ
CHECK v.3.4<--
          CALL PSCTXT(X,Y,TSIZE,NUMBER(IPOS:))
      ENDIF
CHECK v.3.2<--

C---- Label point, if necessary
      IF (RTYPE.LE.LABRES + 1 .AND. RESDUE.NE.'GLY') THEN
CHECK v.3.2-->
C          CALL PSCTXT(X,Y + 3.5 * DY,10.0,RESDUE // SEQNO)
          IF (SEQNO(3:3).EQ.' ') THEN
              IPOS = 4
          ELSE IF (SEQNO(2:2).EQ.' ') THEN
              IPOS = 3
          ELSE IF (SEQNO(1:1).EQ.' ') THEN
              IPOS = 2
          ELSE
              IPOS = 1
          ENDIF
          STRING = RESDUE // ' ' // SEQNO(IPOS:)
          LEN = 4 + (6 - IPOS)
          IF (INCHN.NE.' ') THEN
              STRING(LEN+1:LEN+3) = '(' // INCHN // ')'
              LEN = LEN + 3
          ENDIF
CHECK v.3.4.3-->
CHECK v.3.4.4-->
          IF (INCOLR) THEN
CHECK v.3.4.4<--
              CALL PSCOLB(RGB(1,ICOLR),RGB(2,ICOLR),RGB(3,ICOLR))
CHECK v.3.4.4-->
          ENDIF
CHECK v.3.4.4<--
CHECK v.3.4.3<--
          CALL PSCTXT(X,Y + 3.5 * DY,10.0,STRING(1:LEN))
          RLABEL = RLABEL + 1
CHECK v.3.2<--
      ENDIF

C---- Reset background colour to black
      IF (INCOLR) CALL PSCOLB(0.0,0.0,0.0)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMEND  -  Print statistics for Ramachandran plot
C
C----------------------------------------------------------------------+---

CHECK v.3.4-->
C      SUBROUTINE RAMEND(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
CHECK v.3.4.3-->
C      SUBROUTINE RAMEND(BBOXX1,BBOXX2,BBOXY1,BBOXY2,MALLOW,MCORE,MENDCH,
C     -    MGENER,MGLY,MONPG,MOUTSI,MPRO,MRAMPL,MRESID)
      SUBROUTINE RAMEND(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
     -    MPRO,RAMPUB)
CHECK v.3.4.3<--
CHECK v.3.4<--

CHECK v.3.0.1-->
C      REAL           BBOXX1, BBOXX2, BBOXY1, BBOXY2
CHECK v.3.4.3-->
C      INTEGER        BBOXX1, BBOXX2, BBOXY1, BBOXY2
CHECK v.3.4.3<--
CHECK v.3.4-->
      INTEGER       MALLOW, MCORE, MENDCH, MGENER, MGLY, MONPG, MOUTSI,
CHECK v.3.4.3-->
C     -              MPRO, MRAMPL, MRESID
     -              MPRO
      LOGICAL       RAMPUB
CHECK v.3.4.3<--
CHECK v.3.4<--
CHECK v.3.0.1<--

C---- Print the overall statistics on the Ramachandran plot
CHECK v.3.4-->
C      CALL RMSTAT
CHECK v.3.4.3-->
C      CALL RMSTAT(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
C     -    MPRO,MRAMPL,MRESID)
      IF (.NOT.RAMPUB) THEN
          CALL RMSTAT(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
     -        MPRO)
      ENDIF
CHECK v.3.4.3<--
CHECK v.3.4<--

C---- End of PostScript file
CHECK v.3.2-->
      CALL PSENDP
CHECK v.3.2<--
CHECK v.3.4-->
C      CALL PSCLOS(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
CHECK v.3.4<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RMSTAT  -  Print the overall statistics on the Ramachandran
C                        plot
C
C----------------------------------------------------------------------+---

CHECK v.3.4-->
C      SUBROUTINE RMSTAT
CHECK v.3.4.3-->
C      SUBROUTINE RMSTAT(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
C     -    MPRO,MRAMPL,MRESID)
      SUBROUTINE RMSTAT(MALLOW,MCORE,MENDCH,MGENER,MGLY,MONPG,MOUTSI,
     -    MPRO)
CHECK v.3.4.3<--
CHECK v.3.4<--

      INCLUDE 'tplot.inc'

      CHARACTER*100 IREC
CHECK v.3.4-->
C      INTEGER       NTOTAL
      INTEGER       MALLOW, MCORE, MENDCH, MGENER, MGLY, MONPG, MOUTSI,
CHECK v.3.4.3-->
C     -              MPRO, MRAMPL, MRESID, NTOTAL
     -              MPRO, NTOTAL
CHECK v.3.4.3<--
CHECK v.3.4<--
CHECK v.3.4-->
C      REAL          X, XCENTR, XPOS1, XPOS2, XSTART, Y, YCENTR, YGAP,
C     -              YSTART
      REAL          X, XCENTR, XFPOS1, XFPOS2, XPOS1, XPOS2, XSTART, Y,
     -              YCENTR, YGAP, YSTART
CHECK v.3.4<--

CHECK v.3.4-->
C      DATA XSTART, XPOS1, XPOS2, YGAP / 50.0, 300.0, 340.0, 10.0 /
      DATA XSTART, XFPOS1, XFPOS2, YGAP / 50.0, 300.0, 340.0, 10.0 /
CHECK v.3.4<--

C---- Initialise variables
CHECK v.3.4-->
C      NTOTAL = NCORE + NALLOW + NGENER + NOUTSI + NGLY + NPRO + NENDCH
      NTOTAL = MCORE + MALLOW + MGENER + MOUTSI + MGLY + MPRO + MENDCH
CHECK v.3.4<--
      XCENTR = (XRAMC1 + XRAMC2) / 2.0
      YCENTR = (YRAMC1 + YRAMC2) / 2.0
      X = XRAMC1 + XSTART
CHECK v.3.4-->
C      XPOS1 = XPOS1 + XRAMC1
C      XPOS2 = XPOS2 + XRAMC1
      XPOS1 = XFPOS1 + XRAMC1
      XPOS2 = XFPOS2 + XRAMC1
CHECK v.3.4<--
      YSTART = YRAMC1 - 52.0
      Y = YSTART

C---- Print statistics
      CALL PSCTXT(XCENTR,Y,12.0,'Plot statistics')
      Y = Y - 2 * YGAP

C---- Numbers in each of the four different region types
CHECK v.3.4-->
C      WRITE(IREC,100) NCORE, COREPC
      WRITE(IREC,100) MCORE, COREPC
CHECK v.3.4<--
 100  FORMAT('Residues in most favoured regions  [A,B,L]           ',
     -    I6,6X,F5.1,'%')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
CHECK v.3.4-->
C      WRITE(IREC,120) NALLOW, ALLOWP
      WRITE(IREC,120) MALLOW, ALLOWP
CHECK v.3.4<--
 120  FORMAT('Residues in additional allowed regions  [a,b,l,p]    ',
     -    I6,6X,F5.1,'%')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
CHECK v.3.4-->
C      WRITE(IREC,140) NGENER, GENERP
      WRITE(IREC,140) MGENER, GENERP
CHECK v.3.4<--
 140  FORMAT('Residues in generously allowed regions  [~a,~b,~l,~p]',
     -    I6,6X,F5.1,'%')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
CHECK v.3.4-->
C      WRITE(IREC,160) NOUTSI, OUTSIP
      WRITE(IREC,160) MOUTSI, OUTSIP
CHECK v.3.4<--
 160  FORMAT('Residues in disallowed regions                       ',
     -    I6,6X,F5.1,'%')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      WRITE(IREC,180)
 180  FORMAT(55X,'----',6X,'------')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
CHECK v.3.4-->
C      WRITE(IREC,200) NONPG, 100.0
      WRITE(IREC,200) MONPG, 100.0
CHECK v.3.4<--
 200  FORMAT('Number of non-glycine and non-proline residues       ',
     -    I6,6X,F5.1,'%')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      Y = Y - YGAP / 2.0

C---- Numbers of end-residues
CHECK v.3.4-->
C      WRITE(IREC,300) NENDCH
      WRITE(IREC,300) MENDCH
CHECK v.3.4<--
 300  FORMAT('Number of end-residues (excl. Gly and Pro)           ',
     -    I6)
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      Y = Y - YGAP / 2.0

C---- Numbers oF glycines and prolines
CHECK v.3.4-->
C      WRITE(IREC,320) NGLY
      WRITE(IREC,320) MGLY
CHECK v.3.4<--
 320  FORMAT('Number of glycine residues (shown as triangles)      ',
     -    I6)
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
CHECK v.3.4-->
C      WRITE(IREC,340) NPRO
      WRITE(IREC,340) MPRO
CHECK v.3.4<--
 340  FORMAT('Number of proline residues                           ',
     -    I6)
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      WRITE(IREC,360)
 360  FORMAT(55X,'----')
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      WRITE(IREC,380) NTOTAL
 380  FORMAT('Total number of residues                             ',
     -    I6)
      CALL RAMPUT(X,XPOS1,XPOS2,Y,YGAP,IREC)
      Y = Y - 2.0 * YGAP

C---- Print final explanatory text
      IREC = 'Based on an analysis of 118 structures'//
     -    ' of resolution of at least 2.0 Angstroms'
      CALL PSCTXT(XCENTR,Y,8.0,IREC(1:78))
      Y = Y - YGAP
      IREC = 'and R-factor no greater than 20%,' //
     -    ' a good quality model would be expected'
      CALL PSCTXT(XCENTR,Y,8.0,IREC(1:73))
      Y = Y - YGAP
      IREC = 'to have over 90% in the most favoured regions.'
      CALL PSCTXT(XCENTR,Y,8.0,IREC(1:46))

CHECK v.3.2-->
CHECK v.3.4-->
C      IF ((ENSEMB .OR. NMR) .AND. NFILE.GT.1) THEN
      IF ((ENSEMB .OR. NMR) .AND. NMODEL.GT.1) THEN
CHECK v.3.4<--
          Y = Y - YGAP
          IREC = 'Model numbers shown inside each data point.'
          CALL PSCTXT(XCENTR,Y,8.0,IREC(1:43))
      ENDIF
CHECK v.3.2<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RAMPUT  -  Print statistics line for Ramachandran plot
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMPUT(X,X1,X2,Y,YGAP,IREC)

      CHARACTER*(*) IREC
      REAL          X, X1, X2, Y, YGAP

C---- Print the text and numbers
      CALL PSTEXT(X,Y,9.0,IREC(1:53))
      CALL PSCTXT(X1,Y,9.0,IREC(54:59))
      CALL PSCTXT(X2,Y,9.0,IREC(66:71))

C---- Increment y-value
      Y = Y - YGAP

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2-->
C**************************************************************************
C     
C  SUBROUTINE FILKEY  -  Print out a key cross-referencing the model
C                        numbers to the corresponding filenames
C
C----------------------------------------------------------------------+---
      
      SUBROUTINE FILKEY(XSTART,YSTART,DX,DY,TBIG,TSMALL)
      
      INCLUDE 'tplot.inc'
      
      CHARACTER*3   NUMBER
      INTEGER       IFILE, IPOS
      REAL          DX, DY, TBIG, TSIZE, TSMALL, X, XSEP, XSTART, Y,
     -              YSEP, YSTART

C---- Initialise coordinate positions and text size
      X = XSTART
      Y = YSTART
      XSEP = 2.0 * DX
      IF (NFILE.LT.30) THEN
          YSEP = DY
      ELSE
          YSEP = 0.5 * DY
      ENDIF
CHECK v.3.4-->
C      IF (NFILE.LT.10) THEN
      IF (TOPMOD.LT.10) THEN
CHECK v.3.4<--
          TSIZE = TBIG
      ELSE
          TSIZE = TSMALL
      ENDIF

C---- Print heading
      CALL PSTEXT(X - DX,Y + 2.0 * DY + YSEP,TSIZE,'Files:')

C---- Loop through all the files in the ensemble
      DO 500, IFILE = 1, NFILE

C----     Draw marker box
          CALL PSLWID(0.2)
          CALL PSHADE(1.0,0,RGB,MXCOLR,.FALSE.)
          CALL PSBBOX(X - DX,Y - DY,X - DX,Y + DY,
     -                X + DX,Y + DY,X + DX,Y - DY)

C----     Form and print the model number inside the marker
          WRITE(NUMBER,110) IFILE
 110      FORMAT(I3)
          IF (NUMBER(2:2).EQ.' ') THEN
              IPOS = 3
          ELSE IF (NUMBER(1:1).EQ.' ') THEN
              IPOS = 2
          ELSE
              IPOS = 1
          ENDIF
          CALL PSCTXT(X,Y,TSIZE,NUMBER(IPOS:))

C----     Print file name
          CALL PSTEXT(X + XSEP,Y,TSIZE,
     -        FILID(IFILE)(1:LENID(IFILE)))

C----     Move y-position down for next file
          Y = Y - 2.0 * DY - YSEP
 500  CONTINUE
      
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2<--
C**************************************************************************
C
C  SUBROUTINE WRIDST  -  Write out the data making up the given
C                        distribution
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRIDST(NOBSER,NCELL1,NCELL2)

      INCLUDE 'tplot.inc'

      INTEGER        NCELL1, NCELL2

      CHARACTER*3    AMNAME(NAMINO+1)
      INTEGER        IAMINO, ICELL, JCELL
      REAL           NOBSER(NCELL1,NCELL2,NAMINO+1)

      DATA AMNAME /'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu',
     -             'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe',
     -             'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'ALL' /

C---- Loop over all the amino acids
      DO 800, IAMINO = 1, NAMINO

C----     Write out this amino acid code
          WRITE(27,110) IAMINO, AMNAME(IAMINO)
 110      FORMAT('Amino acid:',I3,'. ',A3)

C----     Loop through all the lines
          DO 200, JCELL = 1, NCELL2

C----         Write out the counts
              WRITE(27,120)
     -            (NINT(NOBSER(ICELL,JCELL,IAMINO)), ICELL = 1, NCELL1)
 120          FORMAT(45I3)
 200      CONTINUE
          WRITE(27,*)
 800  CONTINUE
      
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PLOTS  -  Plot the distributions and data values
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      SUBROUTINE PLOTS(DISTRB)
      SUBROUTINE PLOTS(DISTRB,NOBSER,ENERGY,NCELL1,NCELL2)
CHECK v.3.5.2<--

      INCLUDE 'tplot.inc'

CHECK v.3.5.2-->
      INTEGER        NCELL1, NCELL2
CHECK v.3.5.2<--

CHECK v.3.2-->
      CHARACTER*4    CNUMB
CHECK v.3.2<--
      CHARACTER*7    PAGE
CHECK v.3.2-->
      CHARACTER*9    PLHAND
CHECK v.3.2<--
CHECK v.3.1-->
      CHARACTER*30   PLDESC
CHECK v.3.1<--
      CHARACTER*40   HEADIN
CHECK v.3.2-->
      CHARACTER*60   PTITLE
      CHARACTER*79   BLNKLN, SUMLIN
CHECK v.3.2<--
      INTEGER        BAKCOL, HLEN, IAMINO, ICELL, ICOL, IPAGE, IROW,
     -               JCELL, NCOL, NROW, DISTRB
CHECK v.3.2-->
C      LOGICAL        GENPLT, INCOL
C      REAL           PERCEN, XCENTR, XCORN1, XCORN2, YCORN1, YCORN2
      LOGICAL        FIRST, GENPLT, INCOL, LAST
CHECK v.3.5.2-->
C      REAL           KEYPSX, KEYPSY, PERCEN, XCENTR, XCORN1, XCORN2,
C     -               YCORN1, YCORN2
      REAL           ENERGY(NCELL1,NCELL2,NAMINO+1), KEYPSX, KEYPSY,
     -               NOBSER(NCELL1,NCELL2,NAMINO+1), PERCEN, XCENTR,
     -               XCORN1, XCORN2, YCORN1, YCORN2
CHECK v.3.5.2<--
CHECK v.3.2<--
CHECK v.3.2-->
      DATA BLNKLN( 1:40) / ' |                                      ' /
      DATA BLNKLN(41:79) / '                                      |'  /
CHECK v.3.2<--

C---- Initialise variables
      IF (DISTRB.EQ.1) THEN
          BAKCOL = COLRGP(1)
          IF (ALLRAM) THEN
CHECK v.3.4-->
C              HEADIN = 'Ramachandran plots for all residues'
              HEADIN = 'Ramachandran plots for all residue types'
CHECK v.3.4<--
              HLEN = 40
CHECK v.3.1-->
              PLDESC = 'All-residue Ramachandran plots'
CHECK v.3.1<--
CHECK v.3.2-->
              PLHAND = 'allramach'
              PTITLE = PLDESC
CHECK v.3.2<--
          ELSE
              HEADIN = 'Ramachandran plots for Gly & Pro'
              HLEN = 32
CHECK v.3.1-->
              PLDESC = 'Gly & Pro Ramachandran plots'
CHECK v.3.1<--
CHECK v.3.2-->
              PLHAND = 'ramglypro'
              PTITLE = PLDESC
CHECK v.3.2<--
          ENDIF
          INCOL = INCOLG
      ELSE IF (DISTRB.EQ.2) THEN
          BAKCOL = COLCHI(1)
          HEADIN = 'Chi1-Chi2 plots'
          HLEN = 15
          INCOL = INCOLC
CHECK v.3.1-->
          PLDESC = 'All-residue chi1-chi2 plots'
CHECK v.3.1<--
CHECK v.3.2-->
          PLHAND = 'chi1_chi2'
          PTITLE = PLDESC
CHECK v.3.2<--

CHECK v.3.5.2-->
C---- CA-angle/CA torsion angle distributions
      ELSE IF (DISTRB.EQ.7) THEN
          BAKCOL = COLCHI(1)
          HEADIN = 'CA-angle vs CA torsion angle'
          HLEN = 28
          INCOL = INCOLC
          PLDESC = 'All-residue CA angle plots'
          PLHAND = 'calph_dih'
          PTITLE = PLDESC
CHECK v.3.5.2<--
      ENDIF
      IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
          NCOL = 1
          NROW = 2
      ELSE
          NCOL = 3
          NROW = 3
      ENDIF
      ICOL = 0
      IPAGE = 0
      IROW = 1
      XCENTR = (BBOXX1 + BBOXX2) / 2.0

CHECK v.3.2-->
C---- Initialise counts of labelled, and total, points on the graphs
      NLABEL = 0
      NPOINT = 0
      FIRST = .TRUE.
      KEYPSX = 0.0
      KEYPSY = 0.0
      LAST = .FALSE.
CHECK v.3.2<--

C---- Determine the maximum and minimum values across all graphs
      PERMAX = -999999.99
      PERMIN = 999999.99
      IF (DOPLOT(DISTRB)) THEN
          DO 300, IAMINO = 1, NAMINO
              IF (NCOUNT(IAMINO).NE.0) THEN
CHECK v.3.5.2-->
C                  DO 200, JCELL = 1, NCELL
C                      DO 100, ICELL = 1, NCELL
                  DO 200, JCELL = 1, NCELL2
                      DO 100, ICELL = 1, NCELL1
CHECK v.3.5.2<--
CHECK v.3.4.3-->
C                          PERCEN = ARRAY(ICELL,JCELL,IAMINO)
                          PERCEN = NOBSER(ICELL,JCELL,IAMINO)
CHECK v.3.4.3<--
     -                        * 100.0 / REAL(NCOUNT(IAMINO))
                          PERMIN = MIN(PERMIN,PERCEN)
                          PERMAX = MAX(PERMAX,PERCEN)
100                   CONTINUE
200               CONTINUE
              ENDIF
300       CONTINUE
      ENDIF

C---- Loop through all the amino acid types to be plotted
      DO 500, IAMINO = 1, NAMINO

C----     Determine whether 2D distribution is to be plotted
          GENPLT = .FALSE.
          IF (DOPLOT(DISTRB) .AND. NCOUNT(IAMINO).GT.0) THEN
              GENPLT = .TRUE.
          ENDIF

C----     For phi-psi distribution, plot only the Gly and Pro distrib
C         as have already plotted the Ramachandran plot for the other
C         residue types
          IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM .AND. IAMINO.NE.8 .AND.
     -        IAMINO.NE.15) THEN
              GENPLT = .FALSE.
          ENDIF

C----     Plot if the plot is required
          IF (GENPLT) THEN

C----         Increment column number on page (1-3) throwing new page if
C             necessary
              ICOL = ICOL + 1
              IF (ICOL.GT.NCOL) THEN
                  IROW = IROW + 1
                  ICOL = 1
              ENDIF
              IF (IPAGE.EQ.0 .OR. IROW.GT.NROW) THEN

C----             If this is not the first page, then close current
C                 PostScript file
                  IF (IPAGE.GT.0) THEN
CHECK v.3.2-->
C                      CALL PLTEND(DISTRB)

C----                 For ensemble, plot a file-key to the right of
C                     the Gly/Pro Ramachandran plot
                      IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
                          IF (ENSEMB) THEN
                              CALL FILKEY(KEYPSX,KEYPSY,MKSZ,MKSZ,
     -                            8.0,5.0)
                          ENDIF
                      ENDIF
                      CALL PLTEND(DISTRB,LAST)
CHECK v.3.2<--
                  ENDIF

C----             Open PostScript file and print page-number
CHECK v.3.1-->
C                  CALL PSNAME(FILPS,PSLEN,IPLOT)
CHECK v.3.2-->
C                  CALL PSNAME(FILPS,PSLEN,IPLOT,PLDESC)
C                  PLDESC = ' '
C                  CALL PSOPEN(FILPS,BBOXX1,BBOXX2,BBOXY1,BBOXY2,MXCOLR,
C     -                RGB,INCOL,BAKCOL)
                  IF (.NOT.COMBPS .OR. FIRST) THEN
                      CALL PSNAME(FILPS,PSLEN,IPLOT,PLDESC,PLHAND,
     -                    WITHAN)
CHECK v.3.4-->
C                      CALL PSOPEN(FILPS,MXCOLR,RGB,INCOL,PTITLE,
C     -                    IPAGE + 1)
                      CALL PSOPEN(FILPS,MXCOLR,RGB,NCOLOR,INCOL,PTITLE,
     -                    IPAGE + 1)
CHECK v.3.4<--
                      FIRST = .FALSE.
                  ENDIF
                  PLDESC = ' '
                  CALL PSPAGE(FILPS,BBOXX1,BBOXX2,BBOXY1,BBOXY2,
CHECK v.3.4-->
C     -                MXCOLR,RGB,INCOL,BAKCOL,IPAGE + 1,PLABEL)
     -                MXCOLR,RGB,INCOL,BAKCOL,IPAGE + 1,PLABEL,RSELEC)
CHECK v.3.4<--
CHECK v.3.2<--
CHECK v.3.1<--
                  IPAGE = IPAGE + 1
                  WRITE(PAGE,420) IPAGE
 420              FORMAT('Page ',I2)
                  CALL PSCTXT(BBOXX2 - 40.0,BBOXY2 - 20.0,12.0,PAGE)
                  IROW = 1

CHECK v.3.2-->
C----             Print program name in top left-hand corner
                  IF (NMR) THEN
                      CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -                    'PROCHECK-NMR')
                  ELSE IF (ENSEMB) THEN
                      CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -                    'PROCHECK-COMP')
                  ELSE
                      CALL PSTEXT(BBOXX1 + 15.0,BBOXY2 - 10.0,10.0,
     -                    'PROCHECK')
                  ENDIF
CHECK v.3.2<--

C----             Print graph headings
                  CALL PSCTXT(XCENTR,BBOXY2 - 40.0,25.0,HEADIN(1:HLEN))
CHECK v.3.4-->
C                  CALL PSCTXT(XCENTR,BBOXY2 - 70.0,18.0,BRCODE(1:BLEN))
                  CALL PSCTXT(XCENTR,BBOXY2 - 70.0,18.0,
     -                BRCODE(1:BLEN + 2))
CHECK v.3.4<--
              ENDIF

C----         Initialise variables
              IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
                  XCORN1 = XLEFT1
                  XCORN2 = XCORN1 + XWID1
                  YCORN1 = YBOTT1 + (NROW - IROW) * (YHIGH1 + YSEP1)
                  YCORN2 = YCORN1 + YHIGH1
CHECK v.3.2-->
                  KEYPSX = XCORN2 + 10.0 * MKSZ
                  KEYPSY = YBOTT1 + 2.0 * YHIGH1 + YSEP1 - MKSZ
CHECK v.3.2<--
              ELSE
                  XCORN1 = XLEFT3 + (ICOL - 1) * (XWID3 + XSEP3)
                  XCORN2 = XCORN1 + XWID3
                  YCORN1 = YBOTT3 + (NROW - IROW) * (YHIGH3 + YSEP3)
                  YCORN2 = YCORN1 + YHIGH3
              ENDIF

C----         Plot the current graph
CHECK v.3.5.2-->
C              CALL PLTGRD(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,ICOL,
C     -            DISTRB)
              CALL PLTGRD(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,ICOL,
     -            DISTRB,NOBSER,NCELL1,NCELL2)
CHECK v.3.5.2<--
CHECK v.3.2-->
          ELSE
              XCORN1 = 0.0
              XCORN2 = 0.0
              YCORN1 = 0.0
              YCORN2 = 0.0
CHECK v.3.2<--
          ENDIF

C----     Plot all the data points on the current graph, calculating
C         log-odds score for each one
          IF (NCOUNT(IAMINO).GT.0) THEN
CHECK v.3.5.2-->
C              CALL PROCPT(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,DISTRB,
C     -            GENPLT)
              CALL PROCPT(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,DISTRB,
     -            GENPLT,ENERGY,NCELL1,NCELL2)
CHECK v.3.5.2<--
          ENDIF
 500  CONTINUE

C---- Close the current PostScript file
      IF (DOPLOT(DISTRB)) THEN
CHECK v.3.2-->
C          CALL PLTEND(DISTRB)

C----     For ensemble, plot a file-key to the right of
C         the Gly/Pro Ramachandran plot
          IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
              IF (ENSEMB) THEN
                  CALL FILKEY(KEYPSX,KEYPSY,MKSZ,MKSZ,
     -                8.0,5.0)
              ENDIF
          ENDIF
          LAST = .TRUE.
          CALL PLTEND(DISTRB,LAST)
CHECK v.3.2<--
      ENDIF

CHECK v.3.2-->
C---- Write out the distribution statistics to the .sum summary file
      IF (DOPLOT(DISTRB)) THEN

C----     Blank out the output line
          SUMLIN = BLNKLN

C----     Set up line label according to graph type
          IF (DISTRB.EQ.1) THEN

C----         All 20 amino acid Ramachandrans
              IF (ALLRAM) THEN
                  SUMLIN(4:21) = 'All Ramachandrans:'

C----         Gly & Pro Ramachandran plots
              ELSE
                  SUMLIN(4:21) = 'Gly & Pro Ramach: '
              ENDIF
          ELSE IF (DISTRB.EQ.2) THEN
              SUMLIN(4:21) = 'Chi1-chi2 plots:  '
CHECK v.3.5.2-->
          ELSE IF (DISTRB.EQ.2) THEN
              SUMLIN(4:21) = 'CA angle plots:  '
CHECK v.3.5.2<--
          ENDIF

C----     Show the number of labelled residues found in the low-probability
C         regions of the plots
          WRITE(CNUMB,620) NLABEL
 620      FORMAT(I4)
          SUMLIN(23:26) = CNUMB
          SUMLIN(28:57) = 'labelled residues (out of    )'

C----     Show the total number of values plotted
          WRITE(CNUMB,620) NPOINT
          SUMLIN(53:56) = CNUMB

C----     Determine whether to highlight the line with a + or a *
          IF (NLABEL.GT.0) THEN
              SUMLIN(1:1) = '+'
          ENDIF
          IF (NPOINT.GT.0) THEN
              IF (100.0 * REAL(NLABEL) / REAL(NPOINT).GT.3.0) THEN
                  SUMLIN(1:1) = '*'
              ENDIF
          ENDIF

C----     Write the line out to the summary file
          WRITE(14,660) SUMLIN
 660      FORMAT(A)
      ENDIF

C---- Write out blank line to summary file after chi1-ch2 line
CHECK v.3.5.2-->
C      IF (DISTRB.EQ.2) THEN
      IF (DISTRB.EQ.7) THEN
CHECK v.3.5.2<--
          WRITE(14,660) BLNKLN
      ENDIF
CHECK v.3.2<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PLTGRD  -  Plot graph boxes containing the distributions as
C                        obtained from a high-resolution data set, and the
C                        individual values seen in the structure
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      SUBROUTINE PLTGRD(X1,X2,Y1,Y2,IAMINO,ICOL,DISTRB)
      SUBROUTINE PLTGRD(X1,X2,Y1,Y2,IAMINO,ICOL,DISTRB,NOBSER,NCELL1,
     -    NCELL2)
CHECK v.3.5.2<--

      INCLUDE 'tplot.inc'

CHECK v.3.5.2-->
      INTEGER        NCELL1, NCELL2
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      CHARACTER*6    XAXIS, YAXIS
      CHARACTER*16   XAXIS, YAXIS
CHECK v.3.5.2<--
      INTEGER        ICOL, IAMINO, DISTRB
CHECK v.3.5.2-->
C      REAL           LDROP, LSIZE, TWID, XCENTR, XL1, XL2, X1, X2,
C     -               YCENTR, YL1, YL2, Y1, Y2
      REAL           LDROP, LSIZE, NOBSER(MXCELL,NAMINO+1), TWID,
     -               XCENTR, XL1, XL2, X1, X2, YCENTR, YL1, YL2, Y1, Y2
CHECK v.3.5.2<--

C---- Initialise variables
      XCENTR = (X1 + X2) / 2.0
      YCENTR = (Y1 + Y2) / 2.0
      IF (DISTRB.EQ.1) THEN
          XAXIS = 'Phi'
          YAXIS = 'Psi'
      ELSE IF (DISTRB.EQ.2) THEN
          XAXIS = 'Chi-1'
          YAXIS = 'Chi-2'
CHECK v.3.5.2-->
      ELSE IF (DISTRB.EQ.7) THEN
          XAXIS = 'CA angle'
          YAXIS = 'CA torsion angle'
CHECK v.3.5.2<--
      ENDIF
      IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
          LDROP = 24.0
          LSIZE = 11.0
          TWID = 18.0
      ELSE
          LDROP = 20.0
          LSIZE = 8.0
          TWID = 12.0
      ENDIF

C---- Plot the distribution if not empty
      IF (NCOUNT(IAMINO).GT.0) THEN
CHECK v.3.5.2-->
C          CALL PLTDIS(X1,X2,Y1,Y2,IAMINO,DISTRB)
          CALL PLTDIS(X1,X2,Y1,Y2,IAMINO,DISTRB,NOBSER,NCELL1,NCELL2)
CHECK v.3.5.2<--
      ENDIF

C---- Draw box round graph and label axes
      CALL AXES(X1,X2,Y1,Y2,NTICKX,NTICKY,VALBEG(1),VALEND(1),
     -    VALBEG(2),VALEND(2),LSIZE,0,0,TWID,.TRUE.,.TRUE.,.FALSE.,
     -    .FALSE.,.FALSE.)

C---- Print x-axis title
      CALL PSCTXT(XCENTR,Y1 - LDROP,12.0,XAXIS)

C---- y-axis title
      IF (ICOL.EQ.1) CALL PSRTXT(X1 - LDROP - 5.0,YCENTR,12.0,YAXIS)

C---- For Chi1-Chi2 plots, draw in dotted lines for the favourable
C     positions
      IF (DISTRB.EQ.2) THEN
          XL1 = X1 + (X2 - X1) / 6.0
          XL2 = X1 + (X2 - X1) * (5.0 / 6.0)
          YL1 = Y1 + (Y2 - Y1) / 6.0
          YL2 = Y1 + (Y2 - Y1) * (5.0 / 6.0)
          CALL PSDASH(3)
          CALL PSLWID(0.1)
          CALL PSLINE(X1,YL1,X2,YL1)
          CALL PSLINE(X1,YCENTR,X2,YCENTR)
          CALL PSLINE(X1,YL2,X2,YL2)
          CALL PSLINE(XL1,Y1,XL1,Y2)
          CALL PSLINE(XCENTR,Y1,XCENTR,Y2)
          CALL PSLINE(XL2,Y1,XL2,Y2)
          CALL PSDASH(0)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PLTDIS  -  Plot the distribution as different-shaded regions
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      SUBROUTINE PLTDIS(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,DISTRB)
      SUBROUTINE PLTDIS(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,DISTRB,
     -    NOBSER,NCELL1,NCELL2)
CHECK v.3.5.2<--

      INCLUDE 'tplot.inc'

CHECK v.3.5.2-->
      INTEGER        NCELL1, NCELL2
CHECK v.3.5.2<--

      INTEGER        DISTRB, IAMINO, ICELL, JCELL, MAXCOL, MINCOL
      LOGICAL        INCOL
CHECK v.3.5.2-->
C      REAL           ANGLE, LSHADE, PERCEN, PI, PIBY2, SHADE, X, XCORN1,
C     -               XCORN2, XGAP, XLEFT, X1, Y, YCORN1, YCORN2, YGAP,
C     -               YLEFT
      REAL           ANGLE, LSHADE, NOBSER(NCELL1,NCELL2,NAMINO+1),
     -               PERCEN, PI, PIBY2, SHADE, X, XCORN1, XCORN2, XGAP,
     -               XLEFT, X1, Y, YCORN1, YCORN2, YGAP, YLEFT
CHECK v.3.5.2<--

      PARAMETER     (PI = 3.141592654, PIBY2 = PI / 2.0)

C---- Initialise variables for shading
      IF (DISTRB.EQ.1) THEN
          INCOL = INCOLG
          MINCOL = COLRGP(2)
          MAXCOL = COLRGP(3)
      ELSE
          INCOL = INCOLC
          MINCOL = COLCHI(2)
          MAXCOL = COLCHI(3)
      ENDIF
CHECK v.3.5.2-->
C      XGAP = (XCORN2 - XCORN1) / NCELL
C      YGAP = (YCORN2 - YCORN1) / NCELL
      XGAP = (XCORN2 - XCORN1) / NCELL1
      YGAP = (YCORN2 - YCORN1) / NCELL2
CHECK v.3.5.2<--

CHECK v.3.4.3-->
C---- Plot the graph background
      CALL PSLWID(0.0)
      IF (.NOT.NOSHAD) THEN
          CALL PSHADE(1.0,MINCOL,RGB,MXCOLR,INCOLR)
          CALL PSUBOX(XCORN1,YCORN1,XCORN1,YCORN2,XCORN2,YCORN2,
     -        XCORN2,YCORN1)
      ENDIF
CHECK v.3.4.3<--

C---- Loop through all the cells in this graph
      CALL PSLWID(0.0)
      Y = YCORN1
CHECK v.3.5.2-->
C      DO 300, JCELL = 1, NCELL
      DO 300, JCELL = 1, NCELL2
CHECK v.3.5.2<--
          LSHADE = 1.0
          XLEFT = XCORN1
          X = XLEFT + XGAP
          YLEFT = Y
          X1 = XLEFT
CHECK v.3.5.2-->
C          DO 200, ICELL = 1, NCELL
          DO 200, ICELL = 1, NCELL1
CHECK v.3.5.2<--

C----         Determine the shade to be plotted
              IF (NCOUNT(IAMINO).NE.0) THEN
CHECK v.3.4.3-->
C                  PERCEN = 100.0 * ARRAY(ICELL,JCELL,IAMINO)
                  PERCEN = 100.0 * NOBSER(ICELL,JCELL,IAMINO)
CHECK v.3.4.3<--
     -                / REAL(NCOUNT(IAMINO))
                  ANGLE = PIBY2 * SQRT(PERCEN / PERMAX)
                  SHADE = 1.0 - SIN(ANGLE)
                  IF (SHADE.LT.0.0) SHADE = 0.0
                  IF (SHADE.GT.1.0) SHADE = 1.0
              ELSE
                  SHADE = 0.0
              ENDIF

C----         Write out the appropriate shade
              IF (SHADE.NE.LSHADE) THEN
CHECK v.3.4.3-->
C----             Write out if this square isn't the background colour
                  IF (LSHADE.NE.1.0) THEN
CHECK v.3.4.3<--
                      CALL PSCALE(LSHADE,INCOL,MXCOLR,RGB,MINCOL,MAXCOL)
                      CALL PSUBOX(X1, Y, XLEFT, Y, XLEFT, Y + YGAP,
     -                    X1, Y + YGAP)
CHECK v.3.4.3-->
                  ENDIF
CHECK v.3.4.3<--
                  LSHADE = SHADE
                  X1 = XLEFT
              ENDIF

C----         If the final cell of this row, then write out
CHECK v.3.5.2-->
C              IF (ICELL.EQ.NCELL) THEN
              IF (ICELL.EQ.NCELL1) THEN
CHECK v.3.5.2<--
CHECK v.3.4.3-->
C----             Write out if this square isn't the background colour
                  IF (LSHADE.NE.1.0) THEN
CHECK v.3.4.3<--
                      CALL PSCALE(LSHADE,INCOL,MXCOLR,RGB,MINCOL,MAXCOL)
                      CALL PSUBOX(X1, Y, X, Y, X, Y + YGAP,
     -                    X1, Y + YGAP)
CHECK v.3.4.3-->
                  ENDIF
CHECK v.3.4.3<--
              ENDIF

C----         Increment x-value
              XLEFT = X
              X = X + XGAP
 200      CONTINUE

C----     Increment y-value
          Y = Y + YGAP
 300  CONTINUE
      CALL PSLWID(0.2)

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4.3-->
C**************************************************************************
C
C  SUBROUTINE PROCPT  -  Process the individual data points for this amino
C                        acid type
C
C----------------------------------------------------------------------+---

      SUBROUTINE PROCPT(XCORN1,XCORN2,YCORN1,YCORN2,IAMINO,DISTRB,
CHECK v.3.5.2-->
C     -    GENPLT)
     -    GENPLT,ENERGY,NCELL1,NCELL2)
CHECK v.3.5.2<--

      INCLUDE 'tplot.inc'

CHECK v.3.5.2-->
      INTEGER        NCELL1, NCELL2
CHECK v.3.5.2<--

      CHARACTER*3    AMNAME(NAMINO+1), NUMBER, RESNAM, RESRQD
      CHARACTER*6    CHPTS
      CHARACTER*15   ACTEXT, RIDTXT
      CHARACTER*17   HEADTX
      INTEGER        DISTRB, HEND, IAMINO, ICOLB, ICOLBK, ICOLG, IPOINT,
     -               IPOS, IRESID, I1, I2, NBLANK, NDONE, NINDEX,
     -               SINDEX(MXINDX)
      LOGICAL        GENPLT, INCOL, VALOK
CHECK v.3.5.2-->
C      REAL           CALC1D, CALC2D, CBAD(3), CGOOD(3), LIMSCR, MSIZ,
C     -               SC, SCALEX, SCALEY, SHADE, TOPSHD, TSIZ, TSIZE, TZ,
C     -               X, XANG, XCENTR, XCORN1, XCORN2, Y, YANG, YCORN1,
C     -               YCORN2
      REAL           CALC1D, CALC2D, CBAD(3), CGOOD(3),
     -               ENERGY(NCELL1,NCELL2,NAMINO+1), LIMSCR, MSIZ,
     -               SC, SCALEX, SCALEY, SHADE, TOPSHD, TSIZ, TSIZE, TZ,
     -               X, XANG, XCENTR, XCORN1, XCORN2, Y, YANG, YCORN1,
     -               YCORN2
CHECK v.3.5.2<--

      DATA AMNAME /'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu',
     -             'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe',
     -             'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'ALL' /

C---- Initialise variables
      IF (DISTRB.EQ.1) THEN
          CGOOD(1) = RGB(1,COLRGP(4))
          CGOOD(2) = RGB(2,COLRGP(4))
          CGOOD(3) = RGB(3,COLRGP(4))
          CBAD(1) = RGB(1,COLRGP(5))
          CBAD(2) = RGB(2,COLRGP(5))
          CBAD(3) = RGB(3,COLRGP(5))
          ICOLBK = COLRGP(1)
          ICOLG = COLRGP(4)
          ICOLB = COLRGP(5)
          INCOL = INCOLG
          LIMSCR = LIMGP
          MSIZ = MKSZ
          TSIZ = TXTSZ * 1.75
      ELSE
          CGOOD(1) = RGB(1,COLCHI(4))
          CGOOD(2) = RGB(2,COLCHI(4))
          CGOOD(3) = RGB(3,COLCHI(4))
          CBAD(1) = RGB(1,COLCHI(5))
          CBAD(2) = RGB(2,COLCHI(5))
          CBAD(3) = RGB(3,COLCHI(5))
          ICOLBK = COLCHI(1)
          ICOLG = COLCHI(4)
          ICOLB = COLCHI(5)
          INCOL = INCOLC
          LIMSCR = LIMCHI
          MSIZ = MKSZ * 0.75
          TSIZ = TXTSZ
      ENDIF
      NDONE = 0
      NINDEX = 0
      RESRQD = AMINO(IAMINO)
      SCALEX = (XCORN2 - XCORN1) / (VALEND(1) - VALBEG(1))
      SCALEY = (YCORN2 - YCORN1) / (VALEND(2) - VALBEG(2))
      IF (INCOL) THEN
          TOPSHD = 0.0
      ELSE
          TOPSHD = MAXSHD
      ENDIF
      IF ((ENSEMB .OR. NMR) .AND. NMODEL.GT.1) THEN
          MSIZ = MSIZ * 1.5
      ENDIF

C---- Loop through all the residues in the structure
      DO 100, IRESID = 1, NRESID

C----     Check whether this residue is one of those required
          RESNAM = RESID(IRESID)(3:5)
          IF (RESNAM.EQ.RESRQD) THEN
              XANG = VALUE(1,IRESID)
              IF (XANG.GE.VALBEG(1) .AND. XANG.LE.VALEND(1)) THEN
                  IF (TWODEE) THEN
                      YANG = VALUE(2,IRESID)
                      IF (YANG.GE.VALBEG(2) .AND.
     -                    YANG.LE.VALEND(2)) THEN
                          VALOK = .TRUE.
                      ELSE
                          VALOK = .FALSE.
                      ENDIF
                  ELSE
                      VALOK = .TRUE.
                  ENDIF
              ELSE
                  VALOK = .FALSE.
              ENDIF

C----         If falls within plot area, then calculate its G-factor
C             and store
              IF (VALOK) THEN

C----             Calculate the log-odds score for this value
                  IF (TWODEE) THEN
CHECK v.3.5.2-->
C                      SC = CALC2D(IAMINO,XANG,YANG,NAMINO,NCELL,
C     -                    ENERGY,STEP,VALBEG)
                      SC = CALC2D(IAMINO,XANG,YANG,NAMINO,NCELL1,
     -                    NCELL2,ENERGY,STEP,VALBEG)
CHECK v.3.5.2<--
                  ELSE
CHECK v.3.5.2-->
C                      SC = CALC1D(IAMINO,XANG,NAMINO,NCELL,NCELL1,
C     -                    ENERGY,STEP,VALBEG)
                      SC = CALC1D(IAMINO,XANG,NAMINO,NCELL1,ENERGY,
     -                    STEP,VALBEG)
CHECK v.3.5.2<--

C----                 If this is a Chi-1 distribution and already have
C                     a score for this residue from the Chi1-Chi2
C                     distribution, then discard this score
                      IF (DISTRB.EQ.3 .AND. SCORE(2,IRESID).LT.999.0)
     -                    SC = 0.0
                  ENDIF

C----             Apply normalisation factor to the score
                  IF (SC.NE.0.0 .AND.NRMSTD(IAMINO).NE.0.0) THEN
                      SC = (SC - NRMEAN(IAMINO)) / NRMSTD(IAMINO)
                  ELSE
                      SC = 999.99
                  ENDIF
                  SCORE(DISTRB,IRESID) = SC

C----             Store pointer to this value for sorting by G-factor
                  NINDEX = NINDEX + 1
                  NDONE = NDONE + 1
                  IF (NINDEX.LE.MXINDX) THEN
                      SINDEX(NINDEX) = IRESID
                  ENDIF
              ENDIF
          ENDIF
 100  CONTINUE

C---- Sort all the data-points by G-factor so that the worst ones are
C     plotted last
      NINDEX = MIN(NINDEX,MXINDX)
      IF (NINDEX.GT.1) THEN
          CALL GSORT(DISTRB,NINDEX,SINDEX)
      ENDIF

C---- If data points to be plotted, then do the plot
      IF (GENPLT) THEN

C----     Loop through all the sorted data-points
          DO 300, IPOINT = 1, NINDEX

C----         Get the corresponding residue number and x- and y-values
              IRESID = SINDEX(IPOINT)
              XANG = VALUE(1,IRESID)
              IF (TWODEE) THEN
                  YANG = VALUE(2,IRESID)
              ENDIF
              SC = SCORE(DISTRB,IRESID)

C----         Calculate plot coordinates
              X = XCORN1 + SCALEX * (XANG - VALBEG(1))
              Y = YCORN1 + SCALEY * (YANG - VALBEG(2))
              IF (SC.GT.VALLOW) THEN
                  SHADE = MINSHD
              ELSE IF (SC.LE.VALUPP) THEN
                  SHADE = TOPSHD
              ELSE
                  SHADE = MINSHD - (MINSHD - TOPSHD)
     -                * (SC - VALLOW) / (VALUPP - VALLOW)
              ENDIF
              CALL PSLWID(0.05)
              CALL PSLWID(0.1)
              CALL PSCALE(SHADE,INCOL,MXCOLR,RGB,ICOLG,
     -            ICOLB)
              CALL PSBBOX(X - MSIZ,Y - MSIZ,X - MSIZ,Y + MSIZ,
     -                    X + MSIZ,Y + MSIZ,X + MSIZ,Y - MSIZ)

C----         For ensembles, print the model number inside the marker
              IF ((ENSEMB .OR. NMR) .AND. NMODEL.GT.1) THEN
                  WRITE(NUMBER,110) MODEL(IRESID)
 110              FORMAT(I3)
                  IF (NUMBER(2:2).EQ.' ') THEN
                      IPOS = 3
                  ELSE IF (NUMBER(1:1).EQ.' ') THEN
                      IPOS = 2
                  ELSE
                      IPOS = 1
                  ENDIF
                  IF (TOPMOD.LT.10) THEN
                      TSIZE = 8.0
                  ELSE
                      TSIZE = 5.0
                  ENDIF
                  CALL PSCTXT(X,Y,TSIZE,NUMBER(IPOS:))
              ENDIF

C----         If the point is in an unfavourable region, then print
C             its residue ID
              IF (SC.LT.LIMSCR) THEN
                  Y = Y + 2.5 * MSIZ
                  RIDTXT = RESID(IRESID)

C----             Check for blank spaces at the start and end of
C                 the residue ID
                  IF (RIDTXT(1:1).EQ.' ') THEN
                      I1 = 3
                  ELSE
                      I1 = 1
                  ENDIF
                  IF (RIDTXT(9:10).EQ.'  ') THEN
                      I2 = 8
                  ELSE IF (RIDTXT(10:10).EQ.' ') THEN
                      I2 = 9
                  ELSE
                      I2 = 10
                  ENDIF

C----             For the chi1-chi2 plot cut out the residue name
C                 from the ID
                  IF (DISTRB.EQ.1 .AND. .NOT.ALLRAM) THEN
                      ACTEXT = RIDTXT
                  ELSE
                      ACTEXT = RIDTXT(1:1) // RIDTXT(6:15)
                      IF (I1.EQ.3) I1 = 2
                      I2 = I2 - 4
                  ENDIF
                  TZ = TSIZ
                  IF (I2 - I1.GT.3) TZ = TSIZ * 0.75
                  CALL PSCTXT(X,Y,TZ,ACTEXT(I1:I2))

C----             Increment count of labelled residues
                  NLABEL = NLABEL + 1
              ENDIF
              IF (GENPLT .AND. INCOL) CALL PSCOLB(0.0,0.0,0.0)
 300      CONTINUE
      ENDIF

C---- Print graph heading showing amino acid name and number of data
C     points
      IF (DOPLOT(DISTRB) .AND. GENPLT) THEN
          WRITE(CHPTS,20) NDONE
 20       FORMAT(I6)
          IF (NDONE.EQ.0) CHPTS = '     -'
          NBLANK = 5
          IF (CHPTS(5:5).NE.' ') NBLANK = 4
          IF (CHPTS(4:4).NE.' ') NBLANK = 3
          IF (CHPTS(3:3).NE.' ') NBLANK = 2
          IF (CHPTS(2:2).NE.' ') NBLANK = 1
          IF (CHPTS(1:1).NE.' ') NBLANK = 0
          HEADTX = AMNAME(IAMINO) // ' (' // CHPTS(NBLANK + 1:6) // ')'
          HEND = 6 + (6 - NBLANK)
          XCENTR = (XCORN1 + XCORN2) / 2.0
          CALL PSCTXT(XCENTR,YCORN2 + YTITGP,15.0,HEADTX(1:HEND))

C----     Increment count of points plotted on the distributions
          NPOINT = NPOINT + NDONE
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GSORT  -  Sort of data-points by G-factor, using the SHELL
C                       sort method.
C
C                       [Adapted from Numerical Recipes, Press et al.]
C
C--------------------------------------------------------------------------

      SUBROUTINE GSORT(DISTRB,NINDEX,SINDEX)

      INCLUDE 'tplot.inc'

      REAL          ALN2I, TINY
      PARAMETER (ALN2I=1.4426950, TINY=1.E-5)

      INTEGER       NINDEX
      INTEGER       DISTRB, IRESID, SINDEX(NINDEX)

      INTEGER       I, ISWAP, J, K, L, LOGNB2, M, N, NN
      REAL          G1, G2

C---- Initialise variables

C---- Perform the sort
      N = NINDEX
      LOGNB2 = INT(ALOG(FLOAT(N)) * ALN2I + TINY)
      M = N
      DO 12, NN = 1, LOGNB2
          M = M / 2
          K = N - M
          DO 11, J = 1, K
              I = J
3             CONTINUE
              L = I + M

C----         Retrieve the first G-factor
              IRESID = SINDEX(L)
              G1 = SCORE(DISTRB,IRESID)

C----         Retrieve the second G-factor
              IRESID = SINDEX(I)
              G2 = SCORE(DISTRB,IRESID)

C----         If order the wrong way round, swap the indices
              IF (G1.GT.G2) THEN
                  ISWAP = SINDEX(I)
                  SINDEX(I) = SINDEX(L)
                  SINDEX(L) = ISWAP
                  I = I - M
                  IF (I.GE.1) GO TO 3
              ENDIF
11        CONTINUE
12    CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4.3<--
CHECK v.3.4.3-->
C Modified functions CALC1D and CALC2D transferred to ps.f
CHECK v.3.4.3<--
C**************************************************************************
C
C  SUBROUTINE PLTEND  -  Finish current page
C
C----------------------------------------------------------------------+---

CHECK v.3.2-->
C      SUBROUTINE PLTEND(DISTRB)
      SUBROUTINE PLTEND(DISTRB,LAST)
CHECK v.3.2<--

      INCLUDE 'tplot.inc'

      CHARACTER*120 IREC
      INTEGER       DISTRB
CHECK v.3.2-->
      LOGICAL       LAST
CHECK v.3.2<--
      REAL          LIMSCR, X1, Y1

C---- Print explanatory note about highlighted residues
      IF (DISTRB.EQ.1) THEN
          LIMSCR = LIMGP
      ELSE
          LIMSCR = LIMCHI
      ENDIF
      X1 = BBOXX1 + 30.0
      Y1 = BBOXY1 + 30.0
      WRITE(IREC,20) LIMSCR
 20   FORMAT('Numbers of residues are shown in brackets. Those in unfa',
     -    'vourable conformations (score < ',F5.2,') are labelled.')
      CALL PSTEXT(X1,Y1,10.0,IREC)

C---- Print explanatory note at bottom of page
      X1 = BBOXX1 + 30.0
      Y1 = BBOXY1 + 20.0
      CALL PSTEXT(X1,Y1,10.0,'Shading shows favourable conforma' //
     -    'tions as obtained from an analysis of 163 structures' //
     -    ' at resolution 2.0A or better.')

CHECK v.3.2-->
CHECK v.3.4-->
C      IF ((ENSEMB .OR. NMR) .AND. NFILE.GT.1) THEN
      IF ((ENSEMB .OR. NMR) .AND. NMODEL.GT.1) THEN
CHECK v.3.4<--
          X1 = BBOXX1 + 30.0
          Y1 = BBOXY1 + 10.0
          CALL PSTEXT(X1,Y1,10.0,'Model numbers shown inside ea' //
     -        'ch data point.')
      ENDIF
CHECK v.3.2<--

C---- Close the current PostScript file
CHECK v.3.2-->
C      CALL PSCLOS(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
      CALL PSENDP
      IF (.NOT.COMBPS .OR. LAST) THEN
          CALL PSCLOS(BBOXX1,BBOXX2,BBOXY1,BBOXY2)
      ENDIF
CHECK v.3.2<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETLAN  -  Read through the bond lengths and angles file
C                        to calculate mean values
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE GETLAN

      INCLUDE 'tplot.inc'
 
      CHARACTER*1   INCHN
      CHARACTER*3   RESNAM
      CHARACTER*80  IREC
      INTEGER       IMAIN, IRESID, LINEA
      REAL          LENANG

C---- Initialise variables
      LINEA = 0
      DO 20, IMAIN = 1, NMAIN
          MCHMEA(IMAIN) = 0.0
          MCHNUM(IMAIN) = 0
 20   CONTINUE

C---- Open lengths and angles file, <filename>.lan
      OPEN(UNIT=4, FILE=FILLAN, STATUS='OLD', FORM='FORMATTED',
     -    ACCESS='SEQUENTIAL',
CVAX     -    READONLY,
     -    ERR=900)

C---- Read through the header records
      DO 50, IMAIN = 1, NALP
          READ(4,*)
 50   CONTINUE

C---- Loop while reading in records
 100  CONTINUE
          READ(4,200,END=500,ERR=904) IREC
 200      FORMAT(A)
          LINEA = LINEA + 1

C----     Only process this record if it is not a planar group
          IF (IREC(29:33).NE.'PLANE') THEN
              READ(IREC,210,ERR=906) IRESID, INCHN, RESNAM, IMAIN,
     -            LENANG
 210          FORMAT(I6,A1,5X,A3,I2,F9.4)

C----         Only process this residue if it belongs to the required
C             chain
              IF (CHAIN.EQ.' ' .OR. INCHN.EQ.CHAIN) THEN

C----             Accumulate mean values
                  MCHMEA(IMAIN) = MCHMEA(IMAIN) + LENANG
                  MCHNUM(IMAIN) = MCHNUM(IMAIN) + 1
              ENDIF
          ENDIF

C---- Loop back for next record in file
      GO TO 100

C---- End of file reached
 500  CONTINUE

C---- Calculate mean values
      DO 800, IMAIN = 1, NMAIN
          IF (MCHNUM(IMAIN).GT.0) THEN
              MCHMEA(IMAIN) = MCHMEA(IMAIN) / REAL(MCHNUM(IMAIN))
          ENDIF
 800  CONTINUE

      GO TO 999
 
C---- Fatal errors
 900  CONTINUE
      PRINT*, '*** ERROR. Unable to open .lan file'
      GO TO 990
 
 904  CONTINUE
      PRINT*, '*** ERROR. Error reading .lan file at line:', LINEA + 1
      GO TO 990
 
 906  CONTINUE
      PRINT*, '*** ERROR. Data error in .lan file at line:', LINEA + 1
      PRINT*, IREC
      GO TO 990
 
990   CONTINUE
      IFAIL = .TRUE.
 
999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE COVSCO  -  Read in the main-chain bond lengths and bond
C                        angles and calculate a covalent quality score
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE COVSCO

      INCLUDE 'tplot.inc'
 
      CHARACTER*1   AL(NMAIN), INCHN, OLDCHN
      CHARACTER*3   OLDRES, RESNAM
      CHARACTER*5   OLDSEQ, SEQNO
      CHARACTER*8   INLANG, LASTLA
      CHARACTER*13  ENCODE(NMAIN)
      CHARACTER*16  INDESC
      CHARACTER*132 AREC, IREC
      INTEGER       BANUMB, BLNUMB, IMAIN, IRESID, JRESID, LINE,
     -              SRESID
      REAL          BAMEAN, BLMEAN, CALCOV, ENGMEA(NMAIN),
     -              ENGSTD(NMAIN), LENANG, NSTDEV, SVAL

C---- Initialise variables
      AREC = ' '
      BAMEAN = 0.0
      BANUMB = 0
      BLMEAN = 0.0
      BLNUMB = 0
      JRESID = 0
      LASTLA = ' '
      LINE = 0
      OLDCHN = ' '
      OLDRES = ' '
      OLDSEQ = ' '
      SRESID = 0
      IREC = ' '
      SVAL = 0.0

C---- Rewind the main-chain bond-lengths and angles file
      REWIND(4)

C---- Read through the header records and pick up the small-molecule
C     data
      DO 100, IMAIN = 1, NMAIN
          READ(4,40,ERR=904) AL(IMAIN), INLANG, ENCODE(IMAIN),
     -        INDESC, ENGMEA(IMAIN), ENGSTD(IMAIN)
 40       FORMAT(A1,1X,A8,6X,A13,4X,A16,2F8.3)
 100  CONTINUE

C---- Skip all the records giving data on the planar groups
      DO 200, IMAIN = 1, NPLANE
          READ(4,*)
 200  CONTINUE

C---- Loop while reading in records
 300  CONTINUE
          READ(4,320,END=800,ERR=904) IREC
 320      FORMAT(A)
          LINE = LINE + 1

C----     Only process this record if it is not a planar group
          IF (IREC(29:33).NE.'PLANE') THEN

C----         Read in the required data
              READ(IREC,340,ERR=906) IRESID, INCHN, SEQNO, RESNAM,
     -            IMAIN, LENANG
 340          FORMAT(I6,A1,A5,A3,I2,F9.4)

C----         Only process this record if it belongs to the required chain
              IF (CHAIN.EQ.' ' .OR. INCHN.EQ.CHAIN) THEN

C----             Calculate number of standard deviations that value differs
C                 either from the small-molecule mean, or from protein mean
                  IF (USENGH) THEN
                      NSTDEV = ABS(LENANG - ENGMEA(IMAIN))
     -                    / ENGSTD(IMAIN)
                  ELSE
                      NSTDEV = ABS(LENANG - MCHMEA(IMAIN))
     -                    / ENGSTD(IMAIN)
                  ENDIF

C----             Calculate log-odds score for this value and normalize it
                  SVAL = CALCOV(NSTDEV,GSTEP)
                  SVAL = (SVAL - GAUMEA) / GAUSTD

C----             Accumulate bond length/angle means
                  IF (AL(IMAIN).EQ.'L') THEN
                      BLMEAN = BLMEAN + SVAL
                      BLNUMB = BLNUMB + 1
                      COVAVE(1) = COVAVE(1) + SVAL
                      NUMAVE(1) = NUMAVE(1) + 1
                  ELSE
                      BAMEAN = BAMEAN + SVAL
                      BANUMB = BANUMB + 1
                      COVAVE(2) = COVAVE(2) + SVAL
                      NUMAVE(2) = NUMAVE(2) + 1
                  ENDIF
                  COVAVE(3) = COVAVE(3) + SVAL
                  NUMAVE(3) = NUMAVE(3) + 1

C----             If residue has changed, then calculate average score
                  IF (SRESID.GT.0 .AND. IRESID.NE.SRESID) THEN
                      JRESID = JRESID + 1
                      CALL RESCOR(BANUMB,BLNUMB,BAMEAN,BLMEAN,JRESID)

C----                 Reinitialise stored values
                      SVAL = 0.0
                      BAMEAN = 0.0
                      BANUMB = 0
                      BLMEAN = 0.0
                      BLNUMB = 0
                  ENDIF

C----             Store current value
                  OLDCHN = INCHN
                  OLDRES = RESNAM
                  OLDSEQ = SEQNO
                  SVAL = LENANG
                  SRESID = IRESID
              ENDIF
          ENDIF

C---- Loop back for next record in file
      GO TO 300

C---- End of file reached
 800  CONTINUE

C---- Process final residue
      IF (SRESID.GT.0) THEN
          JRESID = JRESID + 1
          CALL RESCOR(BANUMB,BLNUMB,BAMEAN,BLMEAN,JRESID)
      ENDIF

C---- Calculate overall mean scores
      IF (NUMAVE(1).NE.0) COVAVE(1) = COVAVE(1) / NUMAVE(1)
      IF (NUMAVE(2).NE.0) COVAVE(2) = COVAVE(2) / NUMAVE(2)
      IF (NUMAVE(3).NE.0) COVAVE(3) = COVAVE(3) / NUMAVE(3)

      GO TO 999
 
C---- Fatal errors
 904  CONTINUE
      PRINT*, '*** ERROR. Error reading .lan file at line:', LINE + 1
      GO TO 990
 
 906  CONTINUE
      PRINT*, '*** ERROR. Data error in .lan file at line:', LINE + 1
      PRINT*, IREC

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION CALCOV  -  Calculate log-odds score for the supplied number
C                      of st.devs from mean and supplied interval width
C
C----------------------------------------------------------------------+---

      REAL FUNCTION CALCOV(NSTDEV,GSTEP)

      REAL           STLIM, TWOPI
      PARAMETER     (
     -               STLIM  = 12.0,
     -               TWOPI  = 2.0 * 3.14159265
     -              )

      REAL           AREA, GSTEP, NSTDEV, X1, X2, Y1, Y2

C---- Check that not an extreme outlier
      IF (NSTDEV.GT.STLIM) NSTDEV = STLIM

C---- Calculate height of Normal distribution on either side of the given
C     number of standard deviations
      X1 = NSTDEV - GSTEP / 2.0
      X2 = NSTDEV + GSTEP / 2.0
      Y1 = (1.0 / SQRT(TWOPI)) * EXP (- X1 * X1 / 2.0)
      Y2 = (1.0 / SQRT(TWOPI)) * EXP (- X2 * X2 / 2.0)

C---- Calculate approximate area between these two points
      AREA = Y1 * GSTEP - (Y1 - Y2) * GSTEP / 2.0

C---- Taking area to be proportional to the probability of the value
C     being at this point on the Normal distribution, return the log-odds
C     score
      CALCOV = LOG(AREA)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE RESCOR  -  Calculate the average score for current residue
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE RESCOR(BANUMB,BLNUMB,BAMEAN,BLMEAN,IRESID)

      INCLUDE 'tplot.inc'
 
      INTEGER       BANUMB, BLNUMB, IRESID
      REAL          BAMEAN, BLMEAN

C---- Calculate mean values for this residue
      SCOVAL(3,IRESID) = BLMEAN + BAMEAN
      IF (BLNUMB + BANUMB.NE.0)
     -    SCOVAL(3,IRESID) = SCOVAL(3,IRESID) / (BLNUMB + BANUMB)
      IF (BLNUMB.NE.0) BLMEAN = BLMEAN / BLNUMB
      IF (BANUMB.NE.0) BAMEAN = BAMEAN / BANUMB
      SCOVAL(1,IRESID) = BLMEAN
      SCOVAL(2,IRESID) = BAMEAN

      RETURN
      END

C--------------------------------------------------------------------------
C*****************************************************************************
C
C  SUBROUTINE WRISDH  -  Write out the log-odds scores for each residue
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRISDH

      INCLUDE 'tplot.inc'

      CHARACTER*1   BRCLOS, BROPEN
CHECK v.3.0.1-->
C      INTEGER       ICOUNT(NDISTR + 1), IDISTR, IPOS, IRESID, RCOUNT,
      INTEGER       ICOUNT(NDISTR + 2), IDISTR, IPOS, IRESID, RCOUNT,
CHECK v.3.0.1<--
     -              TCOUNT
      REAL          AVERGE, MEAN(NDISTR + 2), TOTAV

C---- Initialise accumulators for overall means
      DO 50, IDISTR = 1, NDISTR + 2
          ICOUNT(IDISTR) = 0
          MEAN(IDISTR) = 0.0
 50   CONTINUE

C---- Open output .sdh file
      OPEN(UNIT=7,FILE=FILSDH,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL='LIST',
     -    FORM='FORMATTED',ERR=900)

C---- Loop through all the residues, writing the log-odds scores out
      DO 200, IRESID = 1, NRESID
          AVERGE = 0.0
          RCOUNT = 0
          TOTAV = 0.0
          TCOUNT = 0

C----     Calculate average for this residue
          DO 100, IDISTR = 1, NDISTR
              IF (SCORE(IDISTR,IRESID).LT.999.0) THEN
                  AVERGE = AVERGE + SCORE(IDISTR,IRESID)
                  RCOUNT = RCOUNT + 1
                  TOTAV = TOTAV + SCORE(IDISTR,IRESID)
                  TCOUNT = TCOUNT + 1
                  IF (IDISTR.LT.5) THEN
                      IPOS = IDISTR
                  ELSE
                      IPOS = IDISTR - 1
                  ENDIF
                  ICOUNT(IPOS) = ICOUNT(IPOS) + 1
                  MEAN(IPOS) = MEAN(IPOS) + SCORE(IDISTR,IRESID)
                  ICOUNT(NDISTR) = ICOUNT(NDISTR) + 1
                  MEAN(NDISTR) = MEAN(NDISTR)
     -                + SCORE(IDISTR,IRESID)
                  ICOUNT(NDISTR+1) = ICOUNT(NDISTR+1) + 1
                  MEAN(NDISTR+1) = MEAN(NDISTR+1)
     -                + SCORE(IDISTR,IRESID)

C----             Combine the chi-3 and chi-4 scores together
                  IF (IDISTR.EQ.5) THEN
                      IF (SCORE(IDISTR - 1,IRESID).GT.999.0) THEN
                          SCORE(IDISTR - 1,IRESID)
     -                        = SCORE(IDISTR,IRESID)
                      ELSE
                          SCORE(IDISTR - 1,IRESID)
     -                        = (SCORE(IDISTR - 1,IRESID)
     -                        + SCORE(IDISTR,IRESID)) / 2.0
                      ENDIF
                  ENDIF
              ENDIF

C----         For omega distribution and above, move it back one position
CHECK v.3.5.2-->
C              IF (IDISTR.EQ.6) THEN
              IF (IDISTR.GE.6) THEN
CHECK v.3.5.2<--
                  SCORE(IDISTR - 1,IRESID) = SCORE(IDISTR,IRESID)
              ENDIF
 100      CONTINUE
          IF (RCOUNT.NE.0) AVERGE = AVERGE / RCOUNT

C----     Add in the mean bond-lengths and bond angles scores
          IF (SCOVAL(1,IRESID).LT.999.0) THEN
              TOTAV = TOTAV + SCOVAL(1,IRESID)
              TCOUNT = TCOUNT + 1
              ICOUNT(NDISTR+1) = ICOUNT(NDISTR+1) + 1
              MEAN(NDISTR+1) = MEAN(NDISTR+1) + SCOVAL(1,IRESID)
          ENDIF
          IF (SCOVAL(2,IRESID).LT.999.0) THEN
              TOTAV = TOTAV + SCOVAL(2,IRESID)
              TCOUNT = TCOUNT + 1
              ICOUNT(NDISTR+1) = ICOUNT(NDISTR+1) + 1
              MEAN(NDISTR+1) = MEAN(NDISTR+1) + SCOVAL(2,IRESID)
          ENDIF
          IF (TCOUNT.NE.0) TOTAV = TOTAV / TCOUNT

C----     Write the log-odds scores and average out
          WRITE(7,120) SAVRES(IRESID), RESID(IRESID),
     -        (SCORE(IDISTR,IRESID), IDISTR = 1, NDISTR - 1),
     -        AVERGE, (SCOVAL(IDISTR,IRESID), IDISTR = 1, 3), TOTAV
CHECK v.3.5.2-->
C 120      FORMAT(I6,1X,A10,10F9.4)
 120      FORMAT(I6,1X,A10,11F9.4)
CHECK v.3.5.2<--
 200  CONTINUE

C---- Calculate overall means
      DO 300, IDISTR = 1, NDISTR + 1
          IF (ICOUNT(IDISTR).NE.0) THEN
              MEAN(IDISTR) = MEAN(IDISTR) / ICOUNT(IDISTR)
          ENDIF
 300  CONTINUE

C---- Write the overall means out to file
      BROPEN = '('
      BRCLOS = ')'
      IF (CHAIN.EQ.' ') THEN
          BROPEN = ' '
          BRCLOS = ' '
      ENDIF
      WRITE(7,320) (MEAN(IDISTR), IDISTR = 1, NDISTR),
     -    (COVAVE(IDISTR), IDISTR = 1, 3), MEAN(NDISTR + 1), RESOL,
     -    BRCODE(1:PSLEN), BROPEN, CHAIN, BRCLOS
CHECK v.3.5.2-->
C 320  FORMAT('Means:           ',10F9.4,F6.2,2X,A,3A1)
 320  FORMAT('Means:           ',10F9.4,F6.2,2X,A,3A1)
CTEMP 320  FORMAT('Means:           ',11F9.4,F6.2,2X,A,3A1)
CHECK v.3.5.2<--

C---- Close the output file
      CLOSE(7)

      GO TO 999

 900  CONTINUE
      PRINT*, '*** ERROR. Unable to open output file: '
      PRINT*, FILSDH, '*'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------+---
