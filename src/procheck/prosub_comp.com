$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!  PROSUB_COMP.COM
$!  ---------------
$!  Command file for submitting PROCHECK_COMP.COM to a given
$!  batch queue
$!  Roman Laskowski, April 1994
$!
$!  Usage:-        
$!      PROSUB_COMP filelist queue_name output_directory
$!
$!  VAX VMS version
$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$!
$!
$! Check that exactly three parameters have been entered
$!
$ NULL:=""
$ P6 = ' '
$ IF (P1.EQS.NULL) THEN $ GOTO CODEINP
$ IF (P2.EQS.NULL) THEN $ GOTO QUEUEINP
$ IF (P3.EQS.NULL) THEN $ GOTO QUEUEINP
$ GOTO RUNCODE
$!
$ CODEINP:
$ READ/PROMPT="Enter name of file containing list of PDB files: " SYS$COMMAND P1
$ P1 = F$EDIT(P1,"TRIM,UPCASE")
$!
$ QUEUEINP:
$ READ/PROMPT="Enter name of batch-queue:      " SYS$COMMAND P2
$!
$ DIRINP:
$ READ/PROMPT="Enter name of output-directory: " SYS$COMMAND P3
$!
$ RUNCODE:
$ SET DEF 'P3'
$!
$! Check if filelist exists
$ FILELIST:="''P1'"
$ FILE_SPEC = F$SEARCH(FILELIST)
$ IF (FILE_SPEC.NES."") THEN $ GOTO SUBMITJOB
$ WRITE SYS$OUTPUT "File not found: ''FILELIST'"
$ EXIT
$!
$ SUBMITJOB:
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Check of stereochemical quality"
$ WRITE SYS$OUTPUT "-------------------------------"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Coordinates file:   [''P1']"
$ WRITE SYS$OUTPUT "  "
$ WRITE SYS$OUTPUT "Submitting to queue: ''P2'"
$ WRITE SYS$OUTPUT "  "
$ P5=F$LOGICAL(""PRODIR$"")
$ SUBMIT PRODIR$:PROCHECK_COMP.COM /PARAM=('P1','P3','P5') /QUEUE='P2'
