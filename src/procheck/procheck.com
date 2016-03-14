$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!  PROCHECK.COM
$!  ------------
$!  Command file for running sterochemical checking programs
$!  Roman Laskowski, August 1992
$!
$!  VAX VMS version
$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!
$! Check number of parameters entered. If less than 2, then accept missing
$! parameters from user. If there are 3, then 3rd one is the name of the
$! default directory (used when PROCHECK.COM is submitted to a batch queue by
$! PROSUB.COM)
$!
$ NULL:=""
$ P6 = " "
$ IF (P1.EQS.NULL) THEN $ GOTO CODEINP
$ IF (P2.EQS.NULL) THEN $ GOTO CHAININP
$ IF (P3.EQS.NULL) THEN $ GOTO RUNCODE
$ IF (P4.EQS.NULL) THEN $ GOTO GOTCHA
$ IF (P5.EQS.NULL) THEN $ GOTO NOCHAIN
$ P6 = P2
$ P2 = P3
$ SET DEF 'P4'
$ DEFINE PRODIR$ 'P5'
$ GOTO RUNCODE
$!
$ NOCHAIN:
$ SET DEF 'P3'
$ DEFINE PRODIR$ 'P4'
$ GOTO RUNCODE
$!
$ CODEINP:
$ READ/PROMPT="Enter name of coordinates file: " SYS$COMMAND P1
$ P1 = F$EDIT(P1,"TRIM,UPCASE")
$!
$ CHAININP:
$ READ/PROMPT="Enter chain-ID, blank for all: " SYS$COMMAND P6
$!
$ RESINP:
$ READ/PROMPT="Enter resolution: " SYS$COMMAND P2
$ GOTO RUNCODE
$!
$ GOTCHA:
$ P6 = P2
$ P2 = P3
$!
$ RUNCODE:
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT " P   R   O   C   H   E   C   K  "
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Check of stereochemical quality - PROCHECK v.3.2"
$ WRITE SYS$OUTPUT "------------------------------------------------"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Coordinates file:   [''P1']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Chain:              [''P6']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Resolution:          ''P2'"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "  "
$!
$! Check if .pdb file exists
$ PDBFILE:="''P1'"
$ FILE_SPEC = F$SEARCH(PDBFILE)
$ IF (FILE_SPEC.NES."") THEN $ GOTO RUNPROGS
$ WRITE SYS$OUTPUT "File not found: ''PDBFILE'"
$ EXIT
$!
$ RUNPROGS:
$!
$! Create temporary file to hold keyboard inputs required by programs
$ OPEN/WRITE OUTPUT_FILE PROCHECK.XXX
$ WRITE      OUTPUT_FILE "''P1'"
$ WRITE      OUTPUT_FILE "''P6'"
$ WRITE      OUTPUT_FILE "''P2'"
$ WRITE      OUTPUT_FILE "N"
$ CLOSE      OUTPUT_FILE
$!
$!    +-------------+
$!    ! C L E A N   !
$!    +-------------+
$!
$ LOGFILE:="CLEAN.LOG"
$ WRITE SYS$OUTPUT "Running clean-up on file   : ''PDBFILE'"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER PRODIR$:RESDEFS.DAT RESDEFS.DAT
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:CLEAN
$ SEARCH CLEAN.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! S E C S T R !
$!    +-------------+
$!
$ LOGFILE:="SECSTR.LOG"
$ WRITE SYS$OUTPUT "Secondary structure assignment"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:SECSTR
$ SEARCH SECSTR.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! N B         !
$!    +-------------+
$!
$ LOGFILE:="NB.LOG"
$ WRITE SYS$OUTPUT "Non-bonded interactions"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:NB
$ SEARCH NB.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! A N G L E N !
$!    +-------------+
$!
$ LOGFILE:="ANGLEN.LOG"
$ WRITE SYS$OUTPUT "Calculation of bond lengths and bond angles"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:ANGLEN
$ SEARCH ANGLEN.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! T P L O T   !
$!    +-------------+
$!
$ LOGFILE:="TPLOT.LOG"
$ WRITE SYS$OUTPUT "Phi-psi and chi1-chi2 distributions"
$ WRITE SYS$OUTPUT "  "
$! Check if procheck.prm file exists
$ PFILE:="PROCHECK.PRM"
$ FILE_SPEC = F$SEARCH(PFILE)
$ IF (FILE_SPEC.EQS."") THEN $ COPY PRODIR$:PROCHECK.PRM []
$ ASSIGN/USER PRODIR$:PROCHECK.DAT PRODATA
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:TPLOT
$ SEARCH TPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! P P L O T   !
$!    +-------------+
$!
$ LOGFILE:="PPLOT.LOG"
$ WRITE SYS$OUTPUT "Preparing stereochemical quality plots"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:PPLOT
$ SEARCH PPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! B P L O T   !
$!    +-------------+
$!
$ LOGFILE:="BPLOT.LOG"
$ WRITE SYS$OUTPUT "Main-chain bond-lengths and angles, and planar groups"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:BPLOT
$ SEARCH BPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$! Delete unwanted file at end
$ DEL PROCHECK.XXX;*
$ DEL PS.NUMBER;*
