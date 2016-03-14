$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!  PROPLOT.COM
$!  -----------
$!  Command file for running just the PPLOT program of PROCHECK 
$!  Roman Laskowski, August 1992
$!
$!  VAX VMS version
$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!
$!
$! Check that 2 parameters have been entered
$!
$ NULL:=""
$ IF (P1.EQS.NULL) THEN $ GOTO CODEINP
$ IF (P2.EQS.NULL) THEN $ GOTO RESINP
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
$!
$ GOTCHA:
$ P6 = P2
$ P2 = P3
$!
$ RUNCODE:
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Stereochemical quality plots - PROCHECK v.3.2"
$ WRITE SYS$OUTPUT "---------------------------------------------"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Coordinates file: [''P1']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Chain:            [''P6']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Resolution:        ''P2'"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "  "
$!
$ RUNPROGS:
$!
$! Create temporary file to hold keyboard inputs required by programs
$ OPEN/WRITE OUTPUT_FILE PPLOT.IN
$ WRITE      OUTPUT_FILE "''P1'"
$ WRITE      OUTPUT_FILE "''P6'"
$ WRITE      OUTPUT_FILE "''P2'"
$ WRITE      OUTPUT_FILE "N"
$ CLOSE      OUTPUT_FILE
$!
$!    +-------------+
$!    ! T P L O T   !
$!    +-------------+
$!
$ LOGFILE:="TPLOT.LOG"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Phi-psi and chi1-chi2 distributions"
$ WRITE SYS$OUTPUT "  "
$! Check if procheck.prm file exists
$ PFILE:="PROCHECK.PRM"
$ FILE_SPEC = F$SEARCH(PFILE)
$ IF (FILE_SPEC.EQS."") THEN $ COPY PRODIR$:PROCHECK.PRM []
$ ASSIGN/USER PRODIR$:PROCHECK.DAT PRODATA
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PPLOT.IN SYS$INPUT
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
$ ASSIGN/USER PPLOT.IN SYS$INPUT
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
$ ASSIGN/USER PPLOT.IN SYS$INPUT
$ RUN PRODIR$:BPLOT
$ SEARCH BPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$! Delete unwanted files at end
$ DEL PPLOT.IN;*
$ DEL PS.NUMBER;*
