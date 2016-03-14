$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!  PROSUB.COM
$!  ----------
$!  Command file for submitting PROCHECK.COM to a given batch queue
$!  Roman Laskowski, August 1992
$!
$!  Usage:-
$!
$!     PROSUB filename [chain] resolution queue_name output_directory
$!
$!  VAX VMS version
$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!
$!
$! Check that exactly four parameters have been entered
$!
$ NULL:=""
$ P6 = ' '
$ IF (P1.EQS.NULL) THEN $ GOTO CODEINP
$ IF (P2.EQS.NULL) THEN $ GOTO CHAININP
$ IF (P3.EQS.NULL) THEN $ GOTO QUEUEINP
$ IF (P4.EQS.NULL) THEN $ GOTO CHAININP
$ IF (P5.EQS.NULL) THEN $ GOTO CHAININP
$ P6 = P2
$ P2 = P3
$ P3 = P4
$ P4 = P5
$ GOTO RUNCODE
$!
$ CODEINP:
$ READ/PROMPT="Enter name of coordinates file: " SYS$COMMAND P1
$ P1 = F$EDIT(P1,"TRIM,UPCASE")
$!
$ CHAININP:
$ READ/PROMPT="Enter chain-ID, blank for all:  " SYS$COMMAND P6
$!
$ RESINP:
$ READ/PROMPT="Enter resolution:               " SYS$COMMAND P2
$!
$ QUEUEINP:
$ READ/PROMPT="Enter name of batch-queue:      " SYS$COMMAND P3
$!
$ DIRINP:
$ READ/PROMPT="Enter name of output-directory: " SYS$COMMAND P4
$!
$ RUNCODE:
$ SET DEF 'P4'
$!
$! Check if .pdb file exists
$ PDBFILE:="''P1'"
$ FILE_SPEC = F$SEARCH(PDBFILE)
$ IF (FILE_SPEC.NES."") THEN $ GOTO SUBMITJOB
$ WRITE SYS$OUTPUT "File not found: ''PDBFILE'"
$ EXIT
$!
$ SUBMITJOB:
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Check of stereochemical quality"
$ WRITE SYS$OUTPUT "-------------------------------"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Coordinates file:   [''P1']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Chain-ID:           [''P6']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Resolution:          ''P2'"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Submitting to queue: ''P3'"
$ WRITE SYS$OUTPUT "  "
$ P5=F$LOGICAL(""PRODIR$"")
$ SUBMIT PRODIR$:PROCHECK.COM /PARAM=('P1','P6','P2','P4','P5') /QUEUE='P3'
