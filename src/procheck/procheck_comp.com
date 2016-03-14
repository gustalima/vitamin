$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!  PROCHECK_COMP.COM
$!  -----------------
$!  Command file for running sterochemical checking programs
$!  for ensembles of structures.
$!  Roman Laskowski, November 1994
$!
$!  VAX VMS version
$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!
$! Check number of parameters entered. If less than 1, then accept missing
$! parameter1 from user. If there are 3, then 3rd one is the name of the
$! default directory (used when PROCHECK_COMP.COM is submitted to a
$! batch queue by PROSUB_COMP.COM)
$!
$ NULL:=""
$ IF (P1.EQS.NULL) THEN $ GOTO CODEINP
$ IF (P2.EQS.NULL) THEN $ GOTO RUNCODE
$ SET DEF 'P2'
$ DEFINE PRODIR$ 'P3'
$ GOTO RUNCODE
$!
$!
$ CODEINP:
$ READ/PROMPT="Enter name of file containing list of PDB files: " SYS$COMMAND P1
$ P1 = F$EDIT(P1,"TRIM,UPCASE")
$!
$!
$ RUNCODE:
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT " P   R   O   C   H   E   C   K   -   C  O  M  P"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "PROCHECK-COMP"
$ WRITE SYS$OUTPUT "-------------"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "List of PDB files:   [''P1']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "  "
$!
$! Check if filelist exists
$ FILELIST:="''P1'"
$ FILE_SPEC = F$SEARCH(FILELIST)
$ IF (FILE_SPEC.NES."") THEN $ GOTO RUNPROGS
$ WRITE SYS$OUTPUT "File not found: ''FILELIST'"
$ EXIT
$!
$ RUNPROGS:
$!
$! Create temporary file to hold keyboard inputs required by programs
$ OPEN/WRITE OUTPUT_FILE PROCHECK.XXX
$ WRITE      OUTPUT_FILE "%''P1'"
$ WRITE      OUTPUT_FILE " "
$ WRITE      OUTPUT_FILE " "
$ CLOSE      OUTPUT_FILE
$!
$!    +-------------+
$!    ! C L E A N   !
$!    +-------------+
$!
$ LOGFILE:="CLEAN.LOG"
$ WRITE SYS$OUTPUT "Running clean-up on file   : ''FILELIST'"
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
$!    ! R M S D E V !
$!    +-------------+
$!
$ LOGFILE:="RMSDEV.LOG"
$ WRITE SYS$OUTPUT "RMS deviations for ensemble"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:RMSDEV
$ SEARCH RMSDEV.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! T P L O T   !
$!    +-------------+
$!
$! Create temporary file to hold keyboard inputs required by programs
$ OPEN/WRITE OUTPUT_FILE PROCHECK.XXX
$ WRITE      OUTPUT_FILE "%''P1'"
$ WRITE      OUTPUT_FILE " "
$ WRITE      OUTPUT_FILE "0.0"
$ WRITE      OUTPUT_FILE "Y"
$ CLOSE      OUTPUT_FILE
$!
$ LOGFILE:="TPLOT.LOG"
$ WRITE SYS$OUTPUT "Phi-psi and chi1-chi2 distributions"
$ WRITE SYS$OUTPUT "  "
$! Check if procheck_comp.prm file exists
$ PFILE:="PROCHECK_COMP.PRM"
$ FILE_SPEC = F$SEARCH(PFILE)
$ IF (FILE_SPEC.EQS."") THEN $ COPY PRODIR$:PROCHECK_COMP.PRM []
$ ASSIGN/USER PRODIR$:PROCHECK.DAT PRODATA
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:TPLOT
$ SEARCH TPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$!    +-------------+
$!    ! M P L O T   !
$!    +-------------+
$!
$ LOGFILE:="MPLOT.LOG"
$ WRITE SYS$OUTPUT "Dihedral angle distributions and quality plots"
$ WRITE SYS$OUTPUT "  "
$ ASSIGN/USER PRODIR$:PROCHECK.DAT PRODATA
$ ASSIGN/USER 'LOGFILE' SYS$OUTPUT
$ ASSIGN/USER PROCHECK.XXX SYS$INPUT
$ RUN PRODIR$:MPLOT
$ SEARCH MPLOT.LOG "*"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "........................................................"
$!
$! Delete unwanted file at end
$ DEL PROCHECK.XXX;*
