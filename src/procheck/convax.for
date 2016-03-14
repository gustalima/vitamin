C**************************************************************************
C
C  CONVAX.FOR  -  Program for stripping out the CVAX characters from the
C                 source code of PROCHECK. This makes the VAX-specific
C                 lines of code active
C
C----------------------------------------------------------------------+---
C
C The files used by this program are:
C
C  1  Input file (.f)
C  2  Output file (.for)
C
C----------------------------------------------------------------------+---

      IMPLICIT NONE

      INTEGER       NFILES
      PARAMETER    (NFILES=12)


      CHARACTER*12  FILIN(NFILES), FILOUT(NFILES)
      CHARACTER*80  IREC
      INTEGER       IFILE, ILEN, LINE

      DATA FILIN  / 'anglen.f',  'bplot.f',  'clean.f',  'vplot.f ',
     -              'mplot.f',   'pplot.f',  'ps.f',     'rmsdev.f',
     -              'secstr.f',  'tplot.f',  'viol2pdb.f',
     -              'gfac2pdb.f'    /
      DATA FILOUT / 'anglen.for','bplot.for','clean.for','vplot.for ',
     -              'mplot.for', 'pplot.for','ps.for',   'rmsdev.for',
     -              'secstr.for','tplot.for','viol2pdb.for',
     -              'gfac2pdb.for'  /


C---- Loop through the files to be processed
      DO 800, IFILE = 1, NFILES

C----     Initialise variables
          LINE = 0

C----     Open input file, FILIN
          OPEN(UNIT=1,FILE=FILIN(IFILE),STATUS='OLD',
     -        ACCESS='SEQUENTIAL',CARRIAGECONTROL='LIST',
     -        FORM='FORMATTED',ERR=900)
 
C----     Open ouput file
          OPEN(UNIT=2,FILE=FILOUT(IFILE),STATUS='UNKNOWN',
     -        ACCESS='SEQUENTIAL',CARRIAGECONTROL='LIST',
     -        FORM='FORMATTED',ERR=902)

C----     Print message
          PRINT*, 'Converting file: ', FILIN(IFILE)

C----     Loop through all the records in the input file
 100      CONTINUE

C----         Read in the next record from the input file
              READ(1,110,ERR=904,END=500) IREC
 110          FORMAT(A)
              LINE = LINE + 1

C----         Determine the length of the record so that trailing blanks
C             can be stripped off
              ILEN = 80
 200          CONTINUE
              IF (IREC(ILEN:ILEN).EQ.' ') THEN
                  IF (ILEN.GT.1) THEN
                      ILEN = ILEN - 1
                      GO TO 200
                  ENDIF              
              ENDIF              

C----         If it is a CVAX record, write it out stripped of the CVAX
              IF (IREC(1:4).EQ.'CVAX') THEN
                  WRITE(2,110,ERR=906) IREC(5:ILEN)

C----         Otherwise, write it out as it stands
              ELSE
                  WRITE(2,110,ERR=906) IREC(1:ILEN)
              ENDIF

C----     Loop back for next record from input file
          GO TO 100

C----     End of file reached
 500      CONTINUE

C----     Close both files, deleting the old .F file
          CLOSE(1,STATUS='DELETE')
          CLOSE(2)
          PRINT 510, IFILE, FILIN(IFILE), FILOUT(IFILE)
 510      FORMAT(6X,'File',I2,' converted:  ',A,' to ',A)
 800  CONTINUE

      GO TO 999

C---- File-read errors
900   CONTINUE
      PRINT*, '*** ERROR - Unable to open input file:', FILIN(IFILE)
      GO TO 999

902   CONTINUE
      PRINT*, '*** ERROR - Unable to open output file:', FILOUT(IFILE)
      GO TO 999

904   CONTINUE
      PRINT*, '*** Read error on input file at line ', LINE + 1
      GO TO 999

 906  CONTINUE
      PRINT*, '*** Write error on output file at line ', LINE + 1
      GO TO 999

999   CONTINUE
      END

        
C----------------------------------------------------------------------+---
