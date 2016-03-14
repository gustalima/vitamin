C**************************************************************************
C
C  PS.F    - PostScript subroutines
C
C            Written by Roman Laskowski, University College, London,
C            November 1993.
C
C            Original version was part of v.3.0 of the PROCHECK suite
C            of programs. (Previously the routines were incorporated in
C            pplot.f).
C
C            Subsequent amendments will be labelled by CHECK v.m.n-->
C            and CHECK v.m.n<-- where m.n is the version number
C            corresponding to the change
C
C  v.3.0.1 - Change of plot-file names. All plot-file names now have the
C            form: <filename>_nn.ps  where <filename> is the name of the
C            input data file, and nn is a sequential number 01, 02, ...
C                                             Roman Laskowski (29 Mar 1994)
C
C  v.3.1   - Transfer of FINKEY, GETCOL and GETNAM routines used by the
C            plot programs to here. (Unit number of procheck.prm set to
C            10 to be consistent across the plot programs).
C            Addition of part-circle and open-circle routines required
C            by mplot.f.
C            Addition of plot description to display of plot file name.
C                                           Roman Laskowski (8-19 Apr 1994)
C
C  v.3.2   - Addition of identifying plot handle to each plot filename.
C            Addition of descriptive plot title to header of PostScript
C            file.
C            Addition of OPESUM routine to open the .sum summary file.
C            Addition of (optional) name of PostScript file to bottom
C            left-hand corner of the page for ease of reference.
C            Addition of Shell-sort routine, SHSORT.
C                                      Roman Laskowski (22 Apr-3 May 1994)
C            Amendments to PSOPEN and PSCLOS routines, plus addition of
C            PSPAGE and PSENDP routines, to allow several pages to be
C            written into the same PostScript file.
C                                         Roman Laskowski (11-18 Oct 1994)
C            
C  v.3.3.2 - Addition of OPEHTM routine to open the .html summary file.
C                                            Roman Laskowski (21 Aug 1995)
C
C  v.3.4   - Addition of footnote in PSPAGE for residue-range selection.
C                                            Roman Laskowski (25 Mar 1996)
C            Bug-fix for uncentered rotated text print.
C                                            Roman Laskowski (27 Mar 1996)
C            Addition of INSORT routine, being an integer version of
C            SHSORT for sorting arrays of integers. Inclusion of function
C            LENSTR and subroutine ADJLIM.
C                                            Roman Laskowski (29 Mar 1996)
C            Addition of range-file reading and processing routines,
C            GETRNG, GETOKN, INTOKN, DELTOKN, STOTOK, PRNRNG, GETWNT,
C            INMODL and INRANG
C                                            Roman Laskowski ( 2 Apr 1996)
C            Addition of other common plotting routines:- PINRNO, PINTIC,
C            and DHELIX.
C                                            Roman Laskowski (10 Apr 1996)
C            Amendment to PostScript header routines to plot only those
C            colours that have been defined.
C            Removal of unused variables.
C                                            Roman Laskowski (25 Apr 1996)
C            Bug-fixes suggested by Dave Love (routine PSMERK removed as
C            never used), and Adam.
C            Change to output PostScript filenames to deal with cases where
C            over 99 files generated (eg 1HBS and 2GLS - pointed out by
C            Gerard Kleywegt).
C                                            Roman Laskowski (26 Apr 1996)
C
C  v.3.4.3 - Amendment to print of selected model ranges such that default
C            range is whatever is supplied to the routine.
C                                              Roman Laskowski (16 May 1996)
C            Transfer of GETDAT routine from tplot.f and mplot.f.
C                                              Roman Laskowski (22 May 1996)
C            Addition of SMEAR and CALCLO routines, and transfer of
C            modified CALC1D and CALC2D routines from tplot.f and mplot.f.
C                                              Roman Laskowski (24 May 1996)
C            Bug-fix for problem of crashing on VMS systems when writing
C            .html file - output record not long enough - pointed out by
C            Clare Sansom of Leeds.
C                                               Roman Laskowski (6 Jun 1996)
C            Bug-fix for acceptance of BOTH keyword in lower-case as well
C            as upper-case in ranges file.
C                                              Roman Laskowski (11 Jul 1996)
C            Amendments to PostScript set-up lines to make border round plot
C            optional. Addition of PSTRAN routine.
C                                              Roman Laskowski (22 Jul 1996)
C
C  v.3.5     Bug-fix in secondary structure print for short helices.
C            Addition of Greek-letter print and arc plot.
C            Addition of PSDOT routine.
C                                            Roman Laskowski (1-10 Nov 1996)
C            Addition of PSUCIR routine.
C                                              Roman Laskowski (18 Sep 1996)
C
C  v.3.5.2 - Extra parameter to allow GETDAT to continue when normalization
C            data missing from the prodata file such that routine can be
C            called by the program that calculates the normalization factors.
C                                              Roman Laskowski (20 May 1999)
C            Generalisation of arrays holding the torsion angle distributions
C            to allow addition of other data.
C                                           Roman Laskowski (22-25 May 1999)
C
C  v.3.5.3   Addition of PostScript scissors for use by wirplot.f for
C            marking intron positions.
C                                               Roman Laskowski (4 Aug 1999)
C
C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSNAME  -  Determine file name for next PostScript file
C
C----------------------------------------------------------------------+---

CHECK v.3.1-->
C      SUBROUTINE PSNAME(FILPS,PSLEN,IPLOT)
CHECK v.3.2-->
C      SUBROUTINE PSNAME(FILPS,PSLEN,IPLOT,PLDESC)
      SUBROUTINE PSNAME(FILPS,PSLEN,IPLOT,PLDESC,PLHAND,WITHAN)
CHECK v.3.2<--
CHECK v.3.1<--

CHECK v.3.4-->
C      CHARACTER*2    CPLOT
      CHARACTER*3    CPLOT
CHECK v.3.4<--
CHECK v.3.2-->
      CHARACTER*9    PLHAND
CHECK v.3.2<--
CHECK v.3.1-->
      CHARACTER*30   PLDESC
CHECK v.3.1<--
      CHARACTER*(*)  FILPS
CHECK v.3.2-->
C      INTEGER        IPLOT, PSLEN
CHECK v.3.4-->
C      INTEGER        EXTLEN, IPLOT, PSLEN
      INTEGER        EXTLEN, IPLOT, NLEN, PSLEN
CHECK v.3.4<--
      LOGICAL        WITHAN
CHECK v.3.2<--

C---- Increment plot-number
      IPLOT = IPLOT + 1
CHECK v.3.4-->
C      IF (IPLOT.GT.99) IPLOT = IPLOT - 100
C      WRITE(CPLOT,100) IPLOT
C 100  FORMAT(I2)
C      IF (CPLOT(1:1).EQ.' ') CPLOT(1:1) = '0'
      IF (IPLOT.LT.100) THEN
          WRITE(CPLOT,100) IPLOT
 100      FORMAT(I2)
          IF (CPLOT(1:1).EQ.' ') CPLOT(1:1) = '0'
          NLEN = 2
      ELSE
          WRITE(CPLOT,120) IPLOT
 120      FORMAT(I3)
          NLEN = 3
      ENDIF
CHECK v.3.4<--

C---- Form plot filename
CHECK v.3.0.1-->
C      FILPS = FILPS(1:PSLEN) // '.d' // CPLOT
C      PRINT*, '* Plotfile: ', FILPS(1:PSLEN + 4)
CHECK v.3.2-->
C      FILPS = FILPS(1:PSLEN) // '_' // CPLOT // '.ps'
      IF (WITHAN) THEN
CHECK v.3.4-->
C          FILPS = FILPS(1:PSLEN) // '_' // CPLOT // '_' //
          FILPS = FILPS(1:PSLEN) // '_' // CPLOT(1:NLEN) // '_' //
CHECK v.3.4<--
     -        PLHAND // '.ps'
CHECK v.3.4-->
C          EXTLEN = 16
          EXTLEN = 14 + NLEN
CHECK v.3.4<--
      ELSE
CHECK v.3.4-->
C          FILPS = FILPS(1:PSLEN) // '_' // CPLOT // '.ps'
          FILPS = FILPS(1:PSLEN) // '_' // CPLOT(1:NLEN) // '.ps'
CHECK v.3.4<--
CHECK v.3.4-->
C          EXTLEN = 6
          EXTLEN = 4 + NLEN
CHECK v.3.4<--
      ENDIF
CHECK v.3.2<--
CHECK v.3.1-->
C      PRINT*, '* Plotfile: ', FILPS(1:PSLEN + 6)
CHECK v.3.2-->
C      PRINT*, PLDESC, '  * Plotfile: ', FILPS(1:PSLEN + 6)
      PRINT*, PLDESC, '  * File: ', FILPS(1:PSLEN + EXTLEN)
CHECK v.3.2<--
CHECK v.3.1<--
CHECK v.3.0.1<--

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSOPEN  -  Open PostScript file and write out header records
C
C----------------------------------------------------------------------+---

CHECK v.3.2-->
C      SUBROUTINE PSOPEN(FNAME,BBOXX1,BBOXX2,BBOXY1,BBOXY2,MXCOLR,RGB,
C     -    INCOL,BAKCOL)
CHECK v.3.4-->
C      SUBROUTINE PSOPEN(FNAME,MXCOLR,RGB,INCOL,PTITLE,PAGE)
      SUBROUTINE PSOPEN(FNAME,MXCOLR,RGB,NCOLOR,INCOL,PTITLE,PAGE)
CHECK v.3.4<--
CHECK v.3.2<--

      CHARACTER*(*)  FNAME
      CHARACTER*2    COLNO
CHECK v.3.2-->
      CHARACTER*60   PTITLE
CHECK v.3.2<--
CHECK v.3.2-->
C      INTEGER        BAKCOL, BBOXX1, BBOXX2, BBOXY1, BBOXY2, I, ICOL,
C     -               MXCOLR
CHECK v.3.4-->
C      INTEGER        I, ICOL, MXCOLR, PAGE, TITLEN
      INTEGER        I, ICOL, MXCOLR, NCOLOR, PAGE, TITLEN
CHECK v.3.4<--
CHECK v.3.2<--
      LOGICAL        INCOL
      REAL           RGB(3,MXCOLR)

C---- Open output file
      OPEN(UNIT=11, FILE=FNAME, STATUS='UNKNOWN',
     -     FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST',
     -     ERR=900)

CHECK v.3.2-->
C---- Determine the title length
      TITLEN = 61
 5    CONTINUE
          TITLEN = TITLEN - 1
      IF (PTITLE(TITLEN:TITLEN).EQ.' ' .AND. TITLEN.GT.1) GO TO 5
CHECK v.3.2<--

C---- Write out headings to PostScript file
CHECK v.3.2-->
C      WRITE(11,10)
      WRITE(11,10) PTITLE(1:TITLEN), PAGE
CHECK v.3.2<--
 10   FORMAT(
     -    '%!PS-Adobe-3.0',/,
     -    '%%Creator: Procheck',/,
     -    '%%DocumentNeededResources: font Times-Roman Symbol',/,
     -    '%%BoundingBox: (atend)',/,
     -    '%%Pages: 1',/,
CHECK v.3.2-->
     -    '%%Title: ',A,', page ',I2,/,
CHECK v.3.2<--
     -    '%%EndComments',/,
     -    '%%BeginProlog')
      WRITE(11,20)
 20   FORMAT(
     -    '/L { moveto lineto stroke } bind def',/,
     -    '/Col { setrgbcolor } bind def')

C---- Loop to define the colours
      IF (INCOL) THEN
CHECK v.3.4-->
C          DO 100, ICOL = 1, MXCOLR
          DO 100, ICOL = 1, NCOLOR
CHECK v.3.4<--
              WRITE(COLNO,40) ICOL
 40           FORMAT(I2.0)
              IF (COLNO(1:1).EQ.' ') COLNO(1:1) = '0'
              WRITE(11,60) COLNO, (RGB(I,ICOL), I = 1, 3)
 60           FORMAT('/Col',A2,' {gsave ',3F8.4,' setrgbcolor } def')
 100      CONTINUE
      ENDIF

C---- Define commands for boxes, circles, etc.
      WRITE(11,120)
 120  FORMAT(
     -    '/Poly3 { moveto lineto lineto fill grestore } bind ',
     -    'def',/,
     -    '/Pl3 { 6 copy Poly3 moveto moveto moveto closepath ',
     -    'stroke } bind def',/,
     -    '/Pline3 { 6 copy Poly3 moveto lineto lineto closepa',
     -    'th stroke } bind def')
      WRITE(11,250)
 250  FORMAT(
     -    '/Poly4 { moveto lineto lineto lineto fill grestore } bind ',
     -    'def',/,
     -    '/Pl4 { 8 copy Poly4 moveto moveto moveto moveto closepath ',
     -    'stroke } bind def',/,
     -    '/Pline4 { 8 copy Poly4 moveto lineto lineto lineto closepa',
     -    'th stroke } bind def',/,
CHECK v.3.1-->
C     -    '/Circle { gsave newpath 0.0 setgray 3 copy 0 360 arc fill ',
C     -    'stroke grestore }',/
C     -    'bind def')
     -    '/Circol { 1 setgray } def',/,
     -    '/Circle { gsave newpath 0 360 arc gsave Circol fill',
     -    ' grestore stroke grestore }',/
     -    ' bind def',/,
CHECK v.3.5-->
     -    '/Ucircle { gsave newpath 0 360 arc gsave Circol fill ',
     -    'grestore grestore } bind def',/,
CHECK v.3.5<--
     -    '/Ocircle { gsave newpath 0 360 arc stroke grestore } bind d',
     -    'ef',/,
     -     '/PC { moveto arc gsave Circol fill grestore } bind def',/,
     -     '/Pcircle { 7 copy PC moveto arc closepath stroke } bind ',
     -     'def')
CHECK v.3.1<--
      WRITE(11,300)
 300  FORMAT(
     -    '/Print { /Times-Roman findfont exch scalefont setfont show',
     -    ' } bind def',/,
     -    '/Gprint { /Symbol findfont exch scalefont setfont show } b',
     -    'ind def',/,
     -    '/Center {',/,
     -    '  dup /Times-Roman findfont exch scalefont setfont',/,
     -    '  exch stringwidth pop -2 div exch -3 div rmoveto',/,
CHECK v.3.4-->
C     -    ' } bind def',/,
     -    ' } bind def')
      WRITE(11,350)
 350  FORMAT(
CHECK v.3.4<--
     -    '/CenterRot90 {',/,
     -    '  dup /Times-Roman findfont exch scalefont setfont',/,
     -    '  exch stringwidth pop -2 div exch 3 div exch rmoveto',/,
     -    ' } bind def',/,
     -    '/UncenterRot90 {',/,
     -    '  dup /Times-Roman findfont exch scalefont setfont',/,
     -    '  exch stringwidth } bind def',/,
     -    '/Rot90 { gsave currentpoint translate 90 rotate } bind def',/,
CHECK v.3.5.3-->
     -    '/Rot270 { gsave currentpoint translate 270 rotate } bind de',
     -    'f',/,
     -    '/ZDF { /ZapfDingbats findfont 10 scalefont } def',/,
CHECK v.3.5.3<--
     -    '%%EndProlog')
      WRITE(11,400)
 400  FORMAT(
     -    '%%BeginSetup',/,
     -    '1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray [ ] 0 ',
     -    'setdash newpath',/,
CHECK v.3.2-->
C     -    '%%EndSetup',/,
C     -    '%%Page: 1 1 ')
     -    '%%EndSetup')
CHECK v.3.2<--

CHECK v.3.2-->
C      WRITE(11,500) REAL(BBOXX1), REAL(BBOXY1), REAL(BBOXX2),
C     -     REAL(BBOXY1), REAL(BBOXX2), REAL(BBOXY2), REAL(BBOXX1),
C     -     REAL(BBOXY2)
C 500  FORMAT(
C     -    '/ProcheckSave save def',/,
C     -    2F6.2,' moveto',2(2F7.2,' lineto'),/,
C     -    2F7.2,' lineto closepath gsave',/,
C     -    'gsave 1.0000 setgray fill grestore',/,
C     -    'clip newpath')
C
CC---- Draw in the background colour
C      IF (INCOL) THEN
C          CALL PSHADE(0.0,BAKCOL,RGB,MXCOLR,INCOL)
C          CALL PSUBOX(REAL(BBOXX1), REAL(BBOXY1), REAL(BBOXX2),
C     -        REAL(BBOXY1), REAL(BBOXX2), REAL(BBOXY2), REAL(BBOXX1),
C     -        REAL(BBOXY2))
C      ENDIF
CHECK v.3.2<--

      GO TO 999

C---- Errors on parameter file
 900  CONTINUE
      PRINT*, '**** WARNING. Unable to open Postscript file'
      GO TO 999

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2-->
C**************************************************************************
C
C  SUBROUTINE PSPAGE  -  Write out records for the start of a new
C                        PostScript page
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSPAGE(FNAME,BBOXX1,BBOXX2,BBOXY1,BBOXY2,MXCOLR,RGB,
CHECK v.3.4-->
C     -    INCOL,BAKCOL,PAGE,PLABEL)
     -    INCOL,BAKCOL,PAGE,PLABEL,RSELEC)
CHECK v.3.4<--

      CHARACTER*(*)  FNAME
      INTEGER        BAKCOL, BBOXX1, BBOXX2, BBOXY1, BBOXY2, MXCOLR,
     -               NLEN, PAGE
CHECK v.3.4-->
C      LOGICAL        INCOL, PLABEL
      LOGICAL        INCOL, PLABEL, RSELEC
CHECK v.3.4<--
      REAL           RGB(3,MXCOLR)

C---- Write out page-number
      WRITE(11,400) PAGE, PAGE
 400  FORMAT('%%Page: ',2I4)

C---- Draw the border around the plot area
CHECK v.3.4.3-->
C      WRITE(11,500) REAL(BBOXX1), REAL(BBOXY1), REAL(BBOXX2),
C     -     REAL(BBOXY1), REAL(BBOXX2), REAL(BBOXY2), REAL(BBOXX1),
C     -     REAL(BBOXY2)
C 500  FORMAT(
C     -    '/ProcheckSave save def',/,
C     -    2F6.2,' moveto',2(2F7.2,' lineto'),/,
C     -    2F7.2,' lineto closepath',/,
C     -    '  gsave 1.0000 setgray fill grestore',/,
C     -    '  stroke gsave')
      WRITE(11,500) -10.0, -10.0, 800.0, -10.0, 800.0, 900.0,
     -    -10.0, 900.0
 500  FORMAT(
     -    '/ProcheckSave save def',/,
     -    2F7.2,' moveto',2(2F8.2,' lineto'),/,
     -    2F8.2,' lineto closepath',/,
     -    '  gsave 1.0000 setgray fill grestore',/,
     -    '  stroke gsave')
CHECK v.3.4.3<--

C---- If required, write the name of this PostScript file in the bottom
C     left-hand corner of the page for reference
      IF (PLABEL) THEN
          NLEN = INDEX(FNAME,' ') - 1
          WRITE(11,450) REAL(BBOXX1), REAL(BBOXY1 - 10), FNAME(1:NLEN),
     -        10.0
 450      FORMAT(2F5.1,' moveto',/,'(',A,')',/,F4.1,' Print')
      ENDIF

C---- If plot covers only a selected subset of residues, then show footnote
C     in bottom right of plot
      IF (RSELEC) THEN
          WRITE(11,480) REAL(BBOXX2 - 108), REAL(BBOXY1 - 10),
     -        '** Selected residues only', 10.0
 480      FORMAT(2F5.1,' moveto',/,'(',A,')',/,F4.1,' Print')
      ENDIF

C---- Draw in the background colour if required
      IF (INCOL) THEN
CHECK v.3.4.3-->
C          CALL PSHADE(0.0,BAKCOL,RGB,MXCOLR,INCOL)
          IF (BAKCOL.LT.0) CALL PSCOLB(1.0,1.0,1.0)
          CALL PSHADE(0.0,ABS(BAKCOL),RGB,MXCOLR,INCOL)
CHECK v.3.4.3<--
          CALL PSBBOX(REAL(BBOXX1), REAL(BBOXY1), REAL(BBOXX2),
     -        REAL(BBOXY1), REAL(BBOXX2), REAL(BBOXY2), REAL(BBOXX1),
     -        REAL(BBOXY2))
CHECK v.3.4.3-->
      ELSE IF (BAKCOL.GT.0) THEN
          CALL PSHADE(1.0,BAKCOL,RGB,MXCOLR,INCOL)
          CALL PSBBOX(REAL(BBOXX1), REAL(BBOXY1), REAL(BBOXX2),
     -        REAL(BBOXY1), REAL(BBOXX2), REAL(BBOXY2), REAL(BBOXX1),
     -        REAL(BBOXY2))
CHECK v.3.4.3<--
      ENDIF
      CALL PSLWID(0.1)
CHECK v.3.4.3-->
      IF (BAKCOL.LT.0) CALL PSCOLB(0.0,0.0,0.0)
CHECK v.3.4.3<--

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2<--
C**************************************************************************
C
C  SUBROUTINE PSCLOS  -  Write final lines to PostScript file and close
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCLOS(BBOXX1,BBOXX2,BBOXY1,BBOXY2)

      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2

C---- Write out closing lines to PostScript file
      WRITE(11,100) BBOXX1 - 1, BBOXY1 - 1, BBOXX2 + 1, BBOXY2 + 1
 100  FORMAT(
CHECK v.3.2-->
C     -    'grestore stroke',/,
C     -    'ProcheckSave restore',/,
C     -    'showpage',/,
CHECK v.3.2<--
     -    '%%Trailer',/,
     -    '%%BoundingBox:',4I4,/,
     -    '%%EOF')


C---- Close the file
      CLOSE(11)

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSARC  -  Write out an arc of a circle
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSARC(X,Y,RADIUS,ANGLE1,ANGLE2)

      REAL          ANGLE1, ANGLE2, RADIUS, X, Y

      WRITE(11,100) X, Y, RADIUS, ANGLE1, ANGLE2
 100  FORMAT('newpath ',5F7.2,' arc stroke')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
CHECK v.3.2-->
C**************************************************************************
C
C  SUBROUTINE PSENDP  -  Write final lines for the current page to
C                        the PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSENDP

C---- Write out closing lines for this page
      WRITE(11,100)
 100  FORMAT(
     -    'ProcheckSave restore',/,
     -    'showpage')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2<--
C**************************************************************************
C
C  SUBROUTINE PSBBOX  -  Write out bounded box to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSBBOX(X1,Y1,X2,Y2,X3,Y3,X4,Y4)

      REAL          X1, X2, Y1, Y2, X3, Y3, X4, Y4

      WRITE(11,100) X1, Y1, X2, Y2, X3, Y3, X4, Y4
 100  FORMAT(8F7.2,' Pline4')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.1-->
C**************************************************************************
C
C  SUBROUTINE PSCCOL  -  Define colour for filled circle
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCCOL(R,G,B)

      REAL          R, G, B

      WRITE(11,100) R, G, B
 100  FORMAT('/Circol { ',3F7.4,' setrgbcolor } def')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSCSHD  -  Define shade for filled circle
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCSHD(SHADE)

      REAL          SHADE

      WRITE(11,100) SHADE
 100  FORMAT('/Circol { ',F7.4,' setgray } def')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.1<--
C**************************************************************************
C
C  SUBROUTINE PSCIRC  -  Write out a black circle to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCIRC(X,Y,RADIUS)

      REAL          X, Y, RADIUS

      WRITE(11,100) X, Y, RADIUS
 100  FORMAT(3F7.2,' Circle')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5.3-->
C**************************************************************************
C
C  SUBROUTINE PSCISS  -  Write out a pair of scissors of the given colour
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCISS(X,Y,SHADE,COLOUR,INCOL)

      CHARACTER*2   COLNO
      INTEGER       COLOUR
      LOGICAL       INCOL
      REAL          SHADE, X, Y

C---- Set the colour of the scissors
      IF (.NOT.INCOL) THEN
          WRITE(11,100) SHADE
 100      FORMAT('gsave',F7.4,' setgray')
      ELSE
          WRITE(COLNO,140) COLOUR
 140      FORMAT(I2.0)
          IF (COLNO(1:1).EQ.' ') COLNO(1:1) = '0'
          WRITE(11,200) COLNO
 200      FORMAT('Col',A2)
      ENDIF

C---- Fix scissor position
      WRITE(11,300) X, Y
 300  FORMAT(2F7.2,' moveto')

C---- Draw scissors
      WRITE(11,400)
 400  FORMAT('(c) 15.0 CenterRot90',/,
     -    'Rot270',/,
     -    'ZDF setfont',/,
     -    '/charstring ( ) dup 0 34 put def',/,
     -    'charstring show')

C---- Restore graphics state
      WRITE(11,500)
 500  FORMAT('grestore')
      WRITE(11,500)

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5.3<--
C**************************************************************************
C
C  SUBROUTINE PSCOLB  -  Write background colour level
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PSCOLB(R,G,B)

      REAL          R, G, B

      WRITE(11,100) R, G, B
 100  FORMAT(3F8.4,' Col')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSCOLR  -  Write colour level
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PSCOLR(R,G,B)

      REAL          R, G, B

      WRITE(11,100) R, G, B
 100  FORMAT('gsave ',3F8.4,' Col')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSCTXT  -  Write out centred text to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCTXT(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
CXXX      CHARACTER*158 IREC
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y
 100  FORMAT(2F7.2,' moveto')
CXXX      WRITE(IREC,150) TEXT
CXXX 150  FORMAT('(',A,')')
CXXX      IF (IREC(80:).EQ.' ') THEN
          WRITE(11,200) TEXT, SIZE
 200      FORMAT('(',A,')',/,F4.1,' Center')
          WRITE(11,300) TEXT, SIZE
 300      FORMAT('(',A,')',/,F4.1,' Print')
CXXX      ELSE
CXXX          WRITE(11,320) IREC(1:79)
CXXX 320      FORMAT(A)
CXXX          WRITE(11,320) IREC(80:)
CXXX          WRITE(11,350) SIZE
CXXX 350      FORMAT(F4.1,' Center')
CXXX          WRITE(IREC,150) TEXT
CXXX          WRITE(11,320) IREC(1:79)
CXXX          WRITE(11,320) IREC(80:)
CXXX          WRITE(11,400) SIZE
CXXX 400      FORMAT(F4.1,' Print')
CXXX      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSCGTX  -  Write out centred Greek text to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCGTX(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y
 100  FORMAT(2F7.2,' moveto')
      WRITE(11,200) TEXT, SIZE
 200  FORMAT('(',A,')',/,F4.1,' Center')
      WRITE(11,300) TEXT, SIZE
 300  FORMAT('(',A,')',/,F4.1,' Gprint')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
C**************************************************************************
C
C  SUBROUTINE PSLINE  -  Write line out to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSLINE(X1,Y1,X2,Y2)

      REAL          X1, X2, Y1, Y2

      WRITE(11,100) X1, Y1, X2, Y2
 100  FORMAT(4F7.2,' L')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSLWID  -  Write line-width out to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSLWID(LWIDTH)

      REAL          LWIDTH

      WRITE(11,100) LWIDTH
 100  FORMAT(F5.2,' setlinewidth')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.1-->
C**************************************************************************
C
C  SUBROUTINE PSMARK  -  Write point-marker to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSMARK(X,Y,MARKER,HWIDX,HWIDY)

      INTEGER       MARKER
      REAL          HWIDX, HWIDY, X, Y

C---- Print appropriate marker

C---- Oridnary residue
      IF (MARKER.EQ.1) THEN
           CALL PSBBOX(X - HWIDX,Y - HWIDY,X - HWIDX,Y + HWIDY,
     -                 X + HWIDX,Y + HWIDY,X + HWIDX,Y - HWIDY)

C---- Glycine
      ELSE IF (MARKER.EQ.2) THEN
          CALL PSTRIA(X - HWIDX,Y - 0.6 * HWIDY,X + HWIDX,
     -        Y - 0.6 * HWIDY,X,Y + 1.2 * HWIDY)

C---- Residue in disallowed region
      ELSE IF (MARKER.EQ.3) THEN
           CALL PSBBOX(X - HWIDX,Y - HWIDY,X - HWIDX,Y + HWIDY,
     -                 X + HWIDX,Y + HWIDY,X + HWIDX,Y - HWIDY)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSOCIR  -  Write out an open (ie uncoloured, unshaded)
C                        circle to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSOCIR(X,Y,RADIUS)

      REAL          X, Y, RADIUS

      WRITE(11,100) X, Y, RADIUS
 100  FORMAT(3F7.2,' Ocircle')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSPCIR  -  Write out a shaded/coloured part-circle to
C                        PostScript file (start and end-angles supplied
C                        in degrees). The shading goes anti-clockwise from
C                        ANGBEG round to ANGEND, where zero is at 3 o'clock
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSPCIR(X,Y,RADIUS,ANGBEG,ANGEND)

      REAL          ANGBEG, ANGEND, X, Y, RADIUS

      WRITE(11,100) X, Y, RADIUS, ANGBEG, ANGEND, X, Y
 100  FORMAT(7F7.2,' Pcircle')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.1<--
C**************************************************************************
C
C  SUBROUTINE PSRCTX  -  Write out text centred and rotated through 90
C                        degrees to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSRCTX(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y
 100  FORMAT(2F7.2,' moveto')
      WRITE(11,200) TEXT, SIZE
 200  FORMAT('(',A,') ',F4.1,' CenterRot90 Rot90')
      WRITE(11,300) TEXT, SIZE
 300  FORMAT('(',A,') ',F4.1,' Print')
      WRITE(11,400)
 400  FORMAT('grestore')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSRTXT  -  Write out text rotated through 90 degrees to
C                        PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSRTXT(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y
 100  FORMAT(2F7.2,' moveto')
CHECK v.3.4-->
C      WRITE(11,200) TEXT, SIZE
C 200  FORMAT('(',A,') ',F4.1,' UncenterRot90 Rot90')
      WRITE(11,200)
 200  FORMAT('Rot90')
CHECK v.3.4<--
      WRITE(11,300) TEXT, SIZE
 300  FORMAT('(',A,') ',F4.1,' Print')
      WRITE(11,400)
 400  FORMAT('grestore')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSDASH  -  Switch dashed lines on/off
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSDASH(ONOFF)

      INTEGER       ONOFF

      IF (ONOFF.EQ.0) THEN
          WRITE(11,100)
 100      FORMAT('[]  0 setdash')
      ELSE
          WRITE(11,120) ONOFF, ONOFF
 120      FORMAT('[ ', I2, ' ', I2, ' ]  0 setdash')
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSDOT  -  Switch dotted lines on/off
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSDOT(LDASH,LSPACE)

      INTEGER       LDASH, LSPACE

      WRITE(11,100) LDASH, LSPACE
 100  FORMAT('[ ', I2, ' ', I2, ' ] 0 setdash')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
C**************************************************************************
C
C  SUBROUTINE PSCALE  -  Write out colour or black-and-white shading
C                        level to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSCALE(SHADE,INCOL,MXCOLR,RGB,MINCOL,MAXCOL)

      INTEGER       ICOL, MAXCOL, MINCOL, MXCOLR
      LOGICAL       INCOL
      REAL          COLOUR(3), RGB(3,MXCOLR), SHADE

C---- If black-and-white PostScript file, then write out grey-level
      IF (.NOT.INCOL) THEN
          WRITE(11,100) SHADE
 100      FORMAT('gsave',F7.4,' setgray')

C---- Otherwise, if in colour, determine the right mix of colours
      ELSE
          DO 200, ICOL = 1, 3
              COLOUR(ICOL) = RGB(ICOL,MINCOL) + (1.0 - SHADE)
     -            * (RGB(ICOL,MAXCOL) - RGB(ICOL,MINCOL))
 200      CONTINUE

C----     Write out the 3 RGB values
          WRITE(11,220) (COLOUR(ICOL), ICOL = 1, 3)
 220      FORMAT('gsave ',3F8.4,' Col')
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSHADE  -  Write shading level or colour out to PostScript
C                        file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSHADE(SHADE,COLOUR,RGB,MXCOLR,INCOL)

      CHARACTER*2   COLNO
      INTEGER       COLOUR, MXCOLR
      LOGICAL       INCOL
CHECK v.3.4.3-->
C      REAL          RGB(3,MXCOLR), SHADE
      REAL          DUMMY, RGB(3,MXCOLR), SHADE
CHECK v.3.4.3<--

      IF (.NOT.INCOL) THEN
          WRITE(11,100) SHADE
 100      FORMAT('gsave',F7.4,' setgray')
      ELSE
          WRITE(COLNO,140) COLOUR
 140      FORMAT(I2.0)
          IF (COLNO(1:1).EQ.' ') COLNO(1:1) = '0'
          WRITE(11,200) COLNO
 200      FORMAT('Col',A2)
      ENDIF
CHECK v.3.4.3-->
      DUMMY = RGB(1,1)
CHECK v.3.4.3<--

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSREST  -  Write out a "grestore" to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSREST

      WRITE(11,100)
 100  FORMAT('grestore')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSSAVE  -  Write out a "gsave" to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSSAVE

      WRITE(11,100)
 100  FORMAT('gsave')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
C**************************************************************************
C
C  SUBROUTINE PSTEXT  -  Write out text to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSTEXT(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
CXXX      CHARACTER*158 IREC
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y - SIZE / 4.0
 100  FORMAT(2F7.2,' moveto')
CXXX      WRITE(IREC,250) TEXT
CXXX 250  FORMAT('(',A,')')
CXXX      IF (IREC(80:).EQ.' ') THEN
          WRITE(11,300) TEXT, SIZE
 300      FORMAT('(',A,')',/,F4.1,' Print')
CXXX      ELSE
CXXX          WRITE(11,320) IREC(1:79)
CXXX 320      FORMAT(A)
CXXX          WRITE(11,320) IREC(80:)
CXXX          WRITE(11,350) SIZE
CXXX 350      FORMAT(F4.1,' Print')
CXXX      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSTXTG  -  Write out Greek text to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSTXTG(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(11,100) X, Y - SIZE / 4.0
 100  FORMAT(2F7.2,' moveto')
      WRITE(11,300) TEXT, SIZE
 300  FORMAT('(',A,')',/,F4.1,' Gprint')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
CHECK v.3.4.3-->
C**************************************************************************
C
C  SUBROUTINE PSTRAN  -  Write out translation of entire picture in x and y
C                        to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSTRAN(X,Y)

      REAL          X, Y

      WRITE(11,100) X, Y
 100  FORMAT(2F8.2,' translate')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4.3<--
C**************************************************************************
C
C  SUBROUTINE PSTRIA  -  Write out bounded triangle to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSTRIA(X1,Y1,X2,Y2,X3,Y3)

      REAL          X1, X2, Y1, Y2, X3, Y3

      WRITE(11,100) X1, Y1, X2, Y2, X3, Y3
 100  FORMAT(6F7.2,' Pline3')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSUBOX  -  Write out unbounded box to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSUBOX(X1,Y1,X2,Y2,X3,Y3,X4,Y4)

      REAL          X1, X2, Y1, Y2, X3, Y3, X4, Y4

      WRITE(11,100) X1, Y1, X2, Y2, X3, Y3, X4, Y4
 100  FORMAT(8F7.2,' Pl4')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE PSUCIR  -  Write out a circle without a border to PostScript
C                        file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSUCIR(X,Y,RADIUS)

      REAL          X, Y, RADIUS

      WRITE(11,100) X, Y, RADIUS
 100  FORMAT(3F7.2,' Ucircle')

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
C**************************************************************************
C
C  SUBROUTINE PSUTRI  -  Write out unbounded triangle to PostScript file
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PSUTRI(X1,Y1,X2,Y2,X3,Y3)

      REAL          X1, X2, X3, Y1, Y2, Y3

      WRITE(11,100) X1, Y1, X2, Y2, X3, Y3
 100  FORMAT(6F7.2,' Pl3')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE AXES  -  Draw graph axes and axis-labels, writing out to
C                      PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE AXES(XLIM1,XLIM2,YLIM1,YLIM2,NPTSX,NPTSY,XMIN,XMAX,
     -    YMIN,YMAX,SIZLAB,INTRLX,INTRLY,TSIZE,NUMBEX,NUMBEY,ZEROS,
     -    NOEND,FULTIC)

      CHARACTER*5    LABEL
      INTEGER        IFROM, IMARK, INTRLX, INTRLY, IVALUE, NPTSX, NPTSY
      LOGICAL        FULTIC, NOEND, NUMBEX, NUMBEY, ZEROS
      REAL           SIZLAB, TSIZE, UGAP, VALUE, VGAP, X, XCENTR, XGAP,
     -               XLIM1, XLIM2, XMAX, XMIN, X1, X2, Y, YCENTR, YGAP,
     -               YLIM1, YLIM2, YMAX, YMIN, Y1, Y2

C Routine parameters:-
C ------------------
C XLIM1, XLIM2, YLIM1, YLIM2 - x- and y-extent of graph area (ie coords
C          of graph-box). In PostScript coords.
C NPTSX, NPTSY - Numbers of ticks along x- and y-axes
C XMIN, XMAX, YMIN, YMAX - Minimum and maximum x- and y-values on plotted
C          axes
C SIZLAB - Label-size for numbers on the axis (in PostScript coords).
C INTRLX, INTRLY - Integer values indicating how the numbers on the axes
C          are to be formatted. (O=format as integer,1=format as real with
C          1 number after decimal point,2=real with 2 numbers after decimal
C          point.
C TSIZE  - A measure of the distance of axis labels from the axis (in
C          PostScript coords).
C NUMBEX, NUMBEY - Flags indicating whether numbers are actually required
C          on the x- and y-axes
C ZEROS  - Flag indicating whether zero value on the y-axis is required
C NOEND  - Flag indicating whether the last numbers on the x- and y-axes
C          are required
C FULTIC - Flag indicating whether axis ticks are to be full-ticks,
C          extending a little to the other side of the axis, or half-ticks
C          only.

C---- Initialise variables
      IF (NPTSX.GT.0) THEN
          XGAP = (XLIM2 - XLIM1) / NPTSX
          UGAP = (XMAX - XMIN) / NPTSX
      ENDIF
      IF (NPTSY.GT.0) THEN
          YGAP = (YLIM2 - YLIM1) / NPTSY
          VGAP = (YMAX - YMIN) / NPTSY
      ENDIF
      XCENTR = (XLIM1 + XLIM2) / 2.0
      YCENTR = (YLIM1 + YLIM2) / 2.0

C---- Draw box round graph
      CALL PSLWID(0.6)
      CALL PSLINE(XLIM1,YLIM1,XLIM1,YLIM2)
      CALL PSLINE(XLIM1,YLIM2,XLIM2,YLIM2)
      CALL PSLINE(XLIM2,YLIM2,XLIM2,YLIM1)
      CALL PSLINE(XLIM2,YLIM1,XLIM1,YLIM1)

C---- x-axis point markers and labels
      X = XLIM1
      Y = YLIM1 - TSIZE / 2.0
      IF (NPTSX.GT.0) THEN
          DO 100, IMARK = 1, NPTSX + 1
              VALUE = XMIN + (IMARK - 1) * UGAP
              IF (VALUE.GE.0.0) THEN
                  IVALUE = VALUE + 0.5
              ELSE
                  IVALUE = VALUE - 0.5
              ENDIF
              IF (INTRLX.EQ.0) THEN
                  WRITE(LABEL,20) IVALUE
 20               FORMAT(I5)
              ELSE IF (INTRLX.EQ.1) THEN
                  WRITE(LABEL,40) VALUE
 40               FORMAT(F5.1)
              ELSE IF (INTRLX.EQ.2) THEN
                  WRITE(LABEL,60) VALUE
 60               FORMAT(F5.2)
              ENDIF
              IFROM = 1
              IF (LABEL(1:1).EQ.' ') IFROM = 2
              IF (LABEL(2:2).EQ.' ') IFROM = 3
              IF (LABEL(3:3).EQ.' ') IFROM = 4
              IF (LABEL(4:4).EQ.' ') IFROM = 5
              IF (NUMBEX) THEN
                  IF (.NOT.(NOEND .AND. IMARK.EQ.NPTSX + 1)) THEN
                      Y1 = YLIM1
                      Y2 = YLIM1 - SIZLAB / 3.0
                      IF (FULTIC) Y1 = Y1 + SIZLAB / 6.0
                      CALL PSCTXT(X,Y,SIZLAB,LABEL(IFROM:5))
                      CALL PSLINE(X,Y1,X,Y2)
                  ENDIF
              ENDIF
              X = X + XGAP
 100      CONTINUE
      ENDIF

C---- y-axis point markers and labels
      X = XLIM1 - TSIZE
      Y = YLIM1
      IF (NPTSY.GT.0) THEN
          DO 200, IMARK = 1, NPTSY + 1
              VALUE = YMIN + (IMARK - 1) * VGAP
              IF (VALUE.GE.0.0) THEN
                  IVALUE = VALUE + 0.5
              ELSE
                  IVALUE = VALUE - 0.5
              ENDIF
              IF (INTRLY.EQ.0) THEN
                  WRITE(LABEL,20) IVALUE
              ELSE IF (INTRLY.EQ.1) THEN
                  WRITE(LABEL,40) VALUE
              ELSE
                  WRITE(LABEL,60) VALUE
              ENDIF
              IFROM = 1
              IF (NUMBEY .AND. (ZEROS .OR. IMARK.GT.1)) THEN
                  IF (.NOT.(NOEND .AND. IMARK.EQ.NPTSY + 1)) THEN
                      X1 = XLIM1 - SIZLAB / 3.0
                      X2 = XLIM1
                      IF (FULTIC) X2 = X2 + SIZLAB / 6.0
                      CALL PSCTXT(X,Y,SIZLAB,LABEL(1:5))
                      CALL PSLINE(X1,Y,X2,Y)
                  ENDIF
              ENDIF
              Y = Y + YGAP
 200      CONTINUE
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.1-->
C**************************************************************************
C
C  SUBROUTINE FINKEY  -  Find the required keyword in the parameter file
C
C----------------------------------------------------------------------+---

      SUBROUTINE FINKEY(KEY,LEN,LINE,FINERR)

      CHARACTER*(*) KEY
      CHARACTER*80  IREC
      INTEGER       LEN, LINE
      LOGICAL       FINERR

C---- Search through the file until find the keyword
      REWIND(10)
      FINERR = .FALSE.
      LINE = 0
 100  CONTINUE
          LINE = LINE + 1
          READ(10,110,END=900) IREC
 110      FORMAT(A)
      IF (IREC(1:LEN).NE.KEY) GO TO 100

C---- Skip the next line
      LINE = LINE + 1
      READ(10,110,END=900) IREC

      GO TO 999

C---- Unable to find keyword
 900  CONTINUE
CHECK v.3.2-->
C      PRINT*, '*** ERROR. Unable to find keyword in procheck.prm: ',
      PRINT*, '* Unable to find keyword in parameter file: ',
CHECK v.3.2<--
     -    KEY
      FINERR = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETCOL  -  Determine the colour-number of this item
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETCOL(ICOL,LINE,COLNAM,MXCOLR)

      INTEGER       MXCOLR

      CHARACTER*12  COLIN, COLNAM(MXCOLR)
      INTEGER       ICOL, ILOOP, LINE


C---- Read in the next colour
      COLIN = ' '
      LINE = LINE + 1
      READ(10,20,END=100,ERR=100) COLIN
 20   FORMAT(A)
      ILOOP = 0

C---- Loop through all the colour names until find the required one
 100  CONTINUE
          ILOOP = ILOOP + 1
      IF (COLIN.NE.COLNAM(ILOOP) .AND. ILOOP.LT.MXCOLR) GO TO 100

C---- If have found the appropriate colour, store its number
      IF (ILOOP.NE.0) THEN
          ICOL = ILOOP
      ELSE
          PRINT*, '* Colour-type unknown: [', COLIN, '] at line',
     -        LINE, ' in procheck.prm file'
      ENDIF

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

      CHARACTER*1   BSLASH, CLOSEB, OPENB, PCHAR, SLASH
      CHARACTER*78  PDBFIL
      INTEGER       IEND, IPOS, ISTART, ISTATE
      LOGICAL       FINISH, GOTDOT, IERROR

      OPENB = '['
      CLOSEB = ']'
      BSLASH = '\\'
      SLASH = '/'

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
              IF (PCHAR.EQ.SLASH .OR. PCHAR.EQ.BSLASH .OR.
     -            PCHAR.EQ.CLOSEB) THEN
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
              ELSE IF (PCHAR.EQ.SLASH .OR. PCHAR.EQ.BSLASH
     -            .OR. PCHAR.EQ.CLOSEB) THEN
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

C----------------------------------------------------------------------+---
CHECK v.3.1<--
C**************************************************************************
C
C  SUBROUTINE RAMREG  -  Determine which part of the Ramachandran plot
C                        this residue is in
C
C----------------------------------------------------------------------+---

      SUBROUTINE RAMREG(PHI,PSI,REGION,REGNO,REGTYP)

      SAVE

      INTEGER        MAPSIZ
      PARAMETER     (MAPSIZ = 36)

      CHARACTER*1    CHAR
      CHARACTER*2    REGION
      CHARACTER*12   CODE
      CHARACTER*24   NEWCOD
      CHARACTER*36   MAPSTR(MAPSIZ)

      INTEGER        I, J, K, MAP(MAPSIZ,MAPSIZ), REGNO, REGTYP,
     -               RTYPE(12)
      LOGICAL        FIRST
      REAL           GAP, PHI, PSI

      DATA FIRST   / .TRUE. /

      DATA CODE    / 'oBbAaLlexyzw' /
      DATA NEWCOD  / 'XXB b A a L l p ~b~a~l~p' /
      DATA RTYPE   /  1,4,3,4,3,4,3,3,2,2,2,2 /

      DATA MAPSTR( 1) / 'bBBBBBBBBBBbbbxxooooowwwwwoooooxxxxb' /
      DATA MAPSTR( 2) / 'bBBBBBBBBBBBbbbxxooooooooooooooxxxbb' /
      DATA MAPSTR( 3) / 'bBBBBBBBBBBBbbbxxxxooooooooooooxxbbb' /
      DATA MAPSTR( 4) / 'bbBBBBBBBBBBBbbbbxxooooooooooooxxxbb' /
      DATA MAPSTR( 5) / 'bbBBBBBBBBBBBBbbxxxooooooooooooxxxxb' /
      DATA MAPSTR( 6) / 'bbBBBBBBBBBBBbbbxxxoooooooooooooxxxb' /
      DATA MAPSTR( 7) / 'bbbBBBBBBBBBbbbbxxxooozzzzzzoooooxxb' /
      DATA MAPSTR( 8) / 'bbbbBBBBBBBbbbbbbxxzzzzzzzzzoooooxxb' /
      DATA MAPSTR( 9) / 'bbbbbBbbBbbbbbbbbxxzzzzzllzzoooooxxb' /
      DATA MAPSTR(10) / 'bbbbbbbbbbbbbbxxxxxzzllllzzzoooooxxx' /
      DATA MAPSTR(11) / 'bbbbbbbbbbbbxxxxxxxzzllllzzzoooooxxx' /
      DATA MAPSTR(12) / 'xbbbbbbbbbbbbxxoooozzzlllzzzzooooxxx' /
      DATA MAPSTR(13) / 'xxbbbbbbbbbbbxxoooozzllllllzzooooxxx' /
      DATA MAPSTR(14) / 'yyaaaaaaaaaayyyoooozzllLllzzzooooxxx' /
      DATA MAPSTR(15) / 'yaaaaaaaaaaayyyoooozzzlLLlzzzooooxxx' /
      DATA MAPSTR(16) / 'yaaaaaaAaaaaayyyooozzzlllllzzooooxxx' /
      DATA MAPSTR(17) / 'yaaaaaAAAAaaayyyyooozzzlllzzzooooxxx' /
      DATA MAPSTR(18) / 'yaaaaaAAAAAaaayyyooozzzzlllzzooooxxx' /
      DATA MAPSTR(19) / 'yaaaaAAAAAAAaaayyyoozzzlzllzzooooxxx' /
      DATA MAPSTR(20) / 'yaaaaaAAAAAAAaayyyyozzzzzzzzzooooxxx' /
      DATA MAPSTR(21) / 'yyaaaaAAAAAAAAaayyyozzzzzzzzzooooxxx' /
      DATA MAPSTR(22) / 'yyaaaaaAAAAAAAaaayyyoooooooooooooxxx' /
      DATA MAPSTR(23) / 'yyyaaaaaAAAAAAAaayyyoooooooooooooxxx' /
      DATA MAPSTR(24) / 'yaaaaaaaaAAAAAAaaayyyooooooooooooxxx' /
      DATA MAPSTR(25) / 'aayaaaaaaaaAAAAaaayyyooooooooooooxxx' /
      DATA MAPSTR(26) / 'yyyaaaaaaaaaaaaaaaayyooooooooooooxxx' /
      DATA MAPSTR(27) / 'yyyyyaaaaaaaaaaaaayyyooooooooooooxxx' /
      DATA MAPSTR(28) / 'oyyyyaaaaaaaaaayyyyyyooooooooooooooo' /
      DATA MAPSTR(29) / 'oooyyyyyyyyayyyyyyyyoooooooooooooooo' /
      DATA MAPSTR(30) / 'oooxxxxbbxxxxxxxxooooooooooooooooooo' /
      DATA MAPSTR(31) / 'xxxxxbbxxxxxxxooooooowwwwwoooooooooo' /
      DATA MAPSTR(32) / 'xxxxxxbbbbxxxxooooooowwwwwoooooooooo' /
      DATA MAPSTR(33) / 'xbbbxbbbbbxxxxxoooooowwewwwooooooooo' /
      DATA MAPSTR(34) / 'xbbbbbbbbbbxxxxoooooowwewwwooooooxxx' /
      DATA MAPSTR(35) / 'xbbbbbbbbbbbbxxoooooowweewwooooooxxx' /
      DATA MAPSTR(36) / 'bbbbbbbbbbbbbxxoooooowwwwwwooooooxxb' /


C---- If this is the first call to this routine, initialise the MAP array
      IF (FIRST) THEN
          FIRST = .FALSE.
          DO 300, I = 1, MAPSIZ
              DO 200, J = 1, MAPSIZ
                  CHAR = MAPSTR(I)(J:J)
                  DO 100, K = 1, 12
                      IF (CHAR.EQ.CODE(K:K)) MAP(MAPSIZ - I + 1,J) = K
 100              CONTINUE
 200          CONTINUE
 300      CONTINUE
          GAP = 360.0 / MAPSIZ
      ENDIF

C---- Determine which part of the Ramachandran plot the residue is in
      IF (PHI.GT.180.0 .OR. PSI.GT.180.0) THEN
          REGNO = 1
          REGION = 'XX'
          REGTYP = 0
      ELSE
          I = (PSI + 180.0) / GAP + 1
          IF (I.LT.1) I = 1
          IF (I.GT.MAPSIZ) I = MAPSIZ
          J = (PHI + 180.0) / GAP + 1
          IF (J.LT.1) J = 1
          IF (J.GT.MAPSIZ) J = MAPSIZ
          REGNO = MAP(I,J)
          I = 2 * REGNO - 1
          REGION = NEWCOD(I:I+1)
          REGTYP = RTYPE(REGNO)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETPSN  -  Read in the next PostScript file number from
C                        ps.number file
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE GETPSN(IPLOT)

      INTEGER       IPLOT

C---- Open the ps.number file holding the last-used PostScript number
      OPEN(UNIT=12, FILE='ps.number', STATUS='OLD', FORM='FORMATTED',
     -    ACCESS='SEQUENTIAL',
CVAX     -    READONLY,
     -    ERR=900)

C---- Read in the PostScript number
      READ(12,*,END=900,ERR=900) IPLOT

      GO TO 999
 
C---- File not found, or empty
 900  CONTINUE
      IPLOT = 0

999   CONTINUE
      CLOSE(12)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PUTPSN  -  Write out the next PostScript file number from
C                        ps.number file
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PUTPSN(IPLOT)

      INTEGER       IPLOT

C---- Open the ps.number file holding the last-used PostScript number
      OPEN(UNIT=12,FILE='ps.number',STATUS='UNKNOWN',FORM='FORMATTED',
     -    ACCESS='SEQUENTIAL',
CVAX     -    CARRIAGECONTROL='LIST',
     -    ERR=999)

C---- Write out the PostScript number
      WRITE(12,*) IPLOT

C---- Close the file
      CLOSE(12)

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2-->
C**************************************************************************
C
C  SUBROUTINE OPESUM  -  Open the summary file, reading through to the
C                        end if the file already exists
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPESUM(FNAME,NEWFIL,IFAIL)

      CHARACTER*80   IREC
      CHARACTER*(*)  FNAME
      LOGICAL        IFAIL, NEWFIL

C---- Initialise variables
      IFAIL = .FALSE.

C---- Open the summary file
      OPEN(UNIT=14,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
CVAX     -    CARRIAGECONTROL='LIST',
     -    ERR=900)

C---- If this is an existing file, then read through all its records until
C     get to the end of file
      IF (.NOT.NEWFIL) THEN

C----     Loop through the file until reach the end
 100      CONTINUE
             READ(14,110,END=500) IREC
 110         FORMAT(A)
          GO TO 100

C---- If this is a new file, then write the header records to it
      ELSE
          WRITE(14,120)
 120      FORMAT(/,
     -        ' +----------<<<  P  R  O  C  H  E  C  K     S  U  M ',
     -        ' M  A  R  Y  >>>----------+',/,
     -           ' |                                                  ',
     -        '                          |')
      ENDIF

C---- End of file reached
 500  CONTINUE

      GO TO 999

 900  CONTINUE
      PRINT*, '*** ERROR opening summary file'
      PRINT*, '* ', FNAME
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.
      GO TO 999

 999  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------+---
C**************************************************************************
C
C  SUBROUTINE SHSORT  -  Sort routine using the SHELL sort method. Sorts the
C                        real array RARRAY, of size NSIZE, as indexed by the
C                        array INDEX. The latter array has pointers to only
C                        NINDEX of the elements in RARRAY. The routine
C                        returns the indices in INDEX rearranged in descending
C                        order of the corresponding values in RARRAY.
C                        [Adapted from Numerical Recipes, Press et al.]
C
C--------------------------------------------------------------------------

      SUBROUTINE SHSORT(NSIZE,RARRAY,NINDEX,INDEX)

CVAX      IMPLICIT NONE

      REAL          ALN2I, TINY
      PARAMETER (ALN2I=1.4426950, TINY=1.E-5)

      INTEGER       NSIZE, NINDEX
      INTEGER       INDEX(NINDEX)
      REAL          RARRAY(NSIZE)

      INTEGER       I, ISWAP, J, K, L, LOGNB2, M, N, NN

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
              IF (RARRAY(INDEX(L)).GT.RARRAY(INDEX(I))) THEN
                  ISWAP = INDEX(I)
                  INDEX(I) = INDEX(L)
                  INDEX(L) = ISWAP
                  I = I - M
                  IF (I.GE.1) GO TO 3
              ENDIF
11        CONTINUE
12    CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.2<--
CHECK v.3.4-->
C**************************************************************************
C
C  SUBROUTINE INSORT  -  Sort routine using the SHELL sort method. Sorts the
C                        integer array IARRAY, of size NSIZE, as indexed by the
C                        array INDEX. The latter array has pointers to only
C                        NINDEX of the elements in IARRAY. The routine
C                        returns the indices in INDEX rearranged in ascending
C                        order of the corresponding values in IARRAY.
C                        [Adapted from Numerical Recipes, Press et al.]
C
C--------------------------------------------------------------------------

      SUBROUTINE INSORT(NSIZE,IARRAY,NINDEX,INDEX)

CVAX      IMPLICIT NONE

      REAL          ALN2I, TINY
      PARAMETER (ALN2I=1.4426950, TINY=1.E-5)

      INTEGER       NSIZE, NINDEX
      INTEGER       INDEX(NINDEX)
      INTEGER       IARRAY(NSIZE)

      INTEGER       I, ISWAP, J, K, L, LOGNB2, M, N, NN

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
              IF (IARRAY(INDEX(L)).LT.IARRAY(INDEX(I))) THEN
                  ISWAP = INDEX(I)
                  INDEX(I) = INDEX(L)
                  INDEX(L) = ISWAP
                  I = I - M
                  IF (I.GE.1) GO TO 3
              ENDIF
11        CONTINUE
12    CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION LENSTR  -  Return length of string 'CHARS' excluding
C                      trailing blanks
C
C----------------------------------------------------------------------+---

      INTEGER FUNCTION LENSTR(CHARS)

      CHARACTER*(*) CHARS
      INTEGER J, MVAL

      MVAL = LEN(CHARS)
      DO 10, J = MVAL, 1, -1
          IF (CHARS(J:J).NE.' ') GO TO 20
   10 CONTINUE
      J = 0
   20 CONTINUE
      LENSTR = J

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE ADJLIM  -  Adjust the maximum and minimum limits for use
C                        on the graphs
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE ADJLIM(VALMIN,VALMAX,NGAPS,MUSTIN)

      INTEGER       MXSPAN
      PARAMETER    (MXSPAN = 40)

      INTEGER       ACTGAP, IGAPS, ISPAN, NBEST(MXSPAN), NFINT(MXSPAN),
     -              NGAPS
      LOGICAL       DONE, HAVGAP, MUSTIN
      REAL          GAP, VALTOP(MXSPAN), VALMAX, VALMIN, OTHER, SPAN

      DATA VALTOP / 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0,
     -              2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 8.0, 10.0, 15.0, 20.0,
     -              30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 125.0, 150.0,
     -              200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 750.0,
     -              800.0, 1000.0, 1200.0, 1500.0, 2000.0 /

      DATA NBEST  /   2,   3,   4,   5,   3,   4,   5,   6,   3,   4,
     -                5,   6,   4,   5,   6,   3,   4,   5,     5,    4,
     -                 6,    4,    5,    6,    4,     5,     5,     3,
     -                  4,     5,     6,     4,     5,     6,     3,
     -                  4,     5,       6,      3,      4 /

      DATA NFINT  /   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,
     -                0,  -1,   1,   1,   1,   0,   1,   1,     1,    1,
     -                 1,    1,    1,    1,    1,     1,     1,     1,
     -                  1,     1,     1,     1,     1,     1,     1,
     -                  1,     1,       1,      1,      1 /

C---- Calculate y-extents above and below the y=0 line
      SPAN = MAX(ABS(VALMIN),ABS(VALMAX))
      OTHER = MIN(ABS(VALMIN),ABS(VALMAX))
      DONE = .FALSE.
      DO 100, ISPAN = 1, MXSPAN
          IF (.NOT.DONE .AND. SPAN.LE.VALTOP(ISPAN)) THEN

C----         Get the number of gaps for this span
              IGAPS = NBEST(ISPAN)
              HAVGAP = .TRUE.

C----         Check whether have to have integral gaps on the axis
C             and whether this is possible
              IF (MUSTIN) THEN

C----             If this is gives gaps with integer labels, then OK
                  IF (NFINT(ISPAN).EQ.1) THEN
                      HAVGAP = .TRUE.

C----             Otherwise, if can force integral labels, then do so
                  ELSE IF (NFINT(ISPAN).EQ.-1) THEN
                      IGAPS = NINT(VALTOP(ISPAN))
                      HAVGAP = .TRUE.

                  ELSE
                      HAVGAP = .FALSE.
                  ENDIF
              ENDIF

C----         If have a suitable set of gaps, then store
              IF (HAVGAP) THEN
                  ACTGAP = IGAPS
                  SPAN = VALTOP(ISPAN)
                  DONE = .TRUE.
              ENDIF
          ENDIF
 100  CONTINUE

C---- Calculate gap size and maximum value for the other side of zero
      IF (OTHER.NE.0.0) THEN
          GAP = SPAN / ACTGAP
          NGAPS = OTHER / GAP + 1
          OTHER = NGAPS * GAP 
          NGAPS = NGAPS + ACTGAP
      ELSE
          NGAPS = ACTGAP
      ENDIF

C---- Replace maximum and minimum values by adjusted values
      IF (ABS(VALMIN).LT.ABS(VALMAX)) THEN
          VALMIN = - OTHER
          VALMAX = SPAN
      ELSE
          VALMIN = - SPAN
          VALMAX = OTHER
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETRNG  -  Read in the list of residue ranges to be included
C                        in the outputs
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETRNG(FNAME,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,MRANGE,
     -    NRANGE,BOTHND,IFAIL)

      INTEGER       MXTOKN, RECLEN, TOKLEN
      PARAMETER    (MXTOKN = 6, RECLEN = 80, TOKLEN = 5)

      INTEGER       MAXRNG

      CHARACTER*1   CH, TOKCHN(MXTOKN)
      CHARACTER*(TOKLEN) DUMMY, TOKEN(MXTOKN)
      CHARACTER*6   RESFRM(MAXRNG), RESTO(MAXRNG)
      CHARACTER*(RECLEN) IREC
      CHARACTER*(*) FNAME
      INTEGER       IPOS, LENSTR, LINE, MODFRM(MAXRNG), MODTO(MAXRNG),
     -              MRANGE, NRANGE, NTOKEN
      LOGICAL       BOTHND, IFAIL, MODEL

C---- Initialise variables
      BOTHND = .FALSE.
      IFAIL = .FALSE.
      LINE = 0
      MRANGE = 0
      NRANGE = 0

C---- Open input file defining the ranges of residue of interest
      OPEN(UNIT=16, FILE=FNAME, STATUS='OLD', FORM='FORMATTED',
     -     ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST', READONLY,
     -     ERR=900)

C---- Read through the file
 100  CONTINUE

C----     Read in the next residue range and interpret
          READ(16,120,END=1500,ERR=902) IREC
 120      FORMAT(A)
          LINE = LINE + 1

C----     Get first character to check that line hasn't been commented out
          CH = IREC(1:1)
          IF (CH.NE.'%' .AND. CH.NE.'#') THEN

C----         Check whether this line contains the BOTH keyword
              IPOS = INDEX(IREC,'BOTH')
CHECK v.3.4.3-->
              IF (IPOS.EQ.0) THEN
                  IPOS = INDEX(IREC,'both')
              ENDIF
CHECK v.3.4.3<--

C----         If keyword found, then set flag to include only those
C             restraints for which both residues are present in the
C             selected residue ranges
              IF (IPOS.GT.0) THEN
                  BOTHND = .TRUE.
                  CH = '%'
              ENDIF
          ENDIF

C----     Process the line to extract the model- and residue-ranges
          IF (CH.NE.'%' .AND. CH.NE.'#') THEN

C----         Extract all the tokens (ie separate character-strings) from the
C             line just read in
              CALL GETOKN(IREC,RECLEN,TOKEN,TOKCHN,MXTOKN,TOKLEN,NTOKEN)

C----         Interpret the tokens, modifying them accordingly
              CALL INTOKN(TOKEN,TOKCHN,MXTOKN,TOKLEN,NTOKEN,DUMMY,MODEL)

C----         Delete any unwanted tokens (marked for deletion with XXXXX)
              CALL DELTOK(TOKEN,TOKCHN,MXTOKN,NTOKEN)

C----         Store all the tokens representing current model- or
C             residue-range
              CALL STOTOK(TOKEN,TOKCHN,MXTOKN,NTOKEN,MODEL,FNAME,IREC,
     -            LINE,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,NRANGE,MRANGE,
     -            IFAIL)
              IF (IFAIL) GO TO 999
          ENDIF

C---- Loop back for the next record
      GO TO 100

C---- Close the residue ranges file
1500  CONTINUE
      CLOSE(16)

C---- If no model ranges defined, take default assumptionm
      IF (MRANGE.EQ.0) THEN
          MRANGE = 1
CHECK v.3.4.3-->
C          MODFRM(MRANGE) = -99999
C          MODTO(MRANGE) = 99999
CHECK v.3.4.3<--
      ENDIF

C---- If no residue ranges defined, assume that all residues wanted
      IF (NRANGE.EQ.0) THEN
          NRANGE = 1
          RESFRM(NRANGE) = '*ALL  '
          RESTO(NRANGE) = 'XXXXXX'
      ENDIF

C---- Print the model- and residue-ranges selected
      CALL PRNRNG(FNAME,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,NRANGE,
     -    MRANGE,BOTHND)

      GO TO 999

C---- Fatal errors
 900  CONTINUE
      PRINT*, '**** Unable to open input file: ', FNAME(1:LENSTR(FNAME))
      GO TO 990

 902  CONTINUE
      PRINT*, '**** File error reading file ', FNAME(1:LENSTR(FNAME)),
     -    ' at line', LINE + 1
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETOKN  -  Extract all the tokens (ie separate 
C                        character-strings) from the line just read in
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETOKN(IREC,RECLEN,TOKEN,TOKCHN,MXTOKN,TOKLEN,NTOKEN)

      INTEGER       MXTOKN, RECLEN, TOKLEN

      CHARACTER*1   CH, TOKCHN(MXTOKN)
      CHARACTER*(*) TOKEN(MXTOKN)
      CHARACTER*52  LETTER
      CHARACTER*(*) IREC
      INTEGER       IPOS, IPSTN, ITOKEN, IVALUE, NTOKEN

      DATA LETTER( 1:26) / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LETTER(27:52) / 'abcdefghijklmnopqrstuvwxyz'/

C---- Initialise variables for extraction of information from
C     the line just read in
      IPSTN = 0
      IPOS = 0
      DO 200, ITOKEN = 1, MXTOKN
          TOKEN(ITOKEN) = ' '
          TOKCHN(ITOKEN) = ' '
 200  CONTINUE
      ITOKEN = 0

C---- Search the line for each token and store
 400  CONTINUE
          IPOS = IPOS + 1

C----     If have a space, then check if have just ended a token
          IF (IREC(IPOS:IPOS).EQ.' ') THEN

C----         If currently saving a token, then have reached its end
              IPSTN = 0

C----     Otherwise, if this is a character, then store in the
C         current token
          ELSE

C----         If this is the start of a new token, increment
C             token-count
              IF (IPSTN.EQ.0) THEN
                  ITOKEN = ITOKEN + 1
              ENDIF

C----         Process providing we haven't exceeded number of tokens
C             allowed
              IF (ITOKEN.LE.MXTOKN) THEN

C----             Increment position in the current token
                  IPSTN = IPSTN + 1

C----             Convert the character into upper-case if necessary
                  CH = IREC(IPOS:IPOS)
                  IF (LGE(CH,'a') .AND. LLE(CH,'z')) THEN
                      IVALUE = INDEX(LETTER,CH) - 26
                      CH = LETTER(IVALUE:IVALUE)
                  ENDIF

C----             Store the current character in the token
                  IF (IPSTN.LE.TOKLEN) THEN
                      TOKEN(ITOKEN)(IPSTN:IPSTN) = CH
                  ENDIF
              ENDIF
          ENDIF
      IF (IPOS.LT.RECLEN) GO TO 400

C---- Save the number of tokens read in fom this line
      NTOKEN = MIN(ITOKEN,MXTOKN)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE INTOKN  -  Interpret the current line's tokens
C
C----------------------------------------------------------------------+---

      SUBROUTINE INTOKN(TOKEN,TOKCHN,MXTOKN,TOKLEN,NTOKEN,DUMMY,MODEL)

      INTEGER       MXTOKN, TOKLEN

      CHARACTER*1   CH, CHAIN, TOKCHN(MXTOKN)
      CHARACTER*(*) DUMMY, TOKEN(MXTOKN)
      INTEGER       ILEN, ITOKEN, JTOKEN, LENSTR, LTOKEN,
     -              NTOKEN
      LOGICAL       MODEL

C---- Initialise variables for interpretation of tokens
      CHAIN = ' '
      LTOKEN = 0
      MODEL = .FALSE.

C---- Check whether this is a MODEL range
      IF (TOKEN(1).EQ.'MODEL') THEN
          MODEL = .TRUE.
          TOKCHN(1) = 'X'
          TOKEN(1) = 'XXXXX'
      ENDIF

C---- If this is a residue line, then mark token for deletion
      IF (TOKEN(1)(1:3).EQ.'RES') THEN
          TOKCHN(1) = 'X'
          TOKEN(1) = 'XXXXX'
      ENDIF

C---- Process all the stored tokens
      DO 800, ITOKEN = 1, NTOKEN

C----     Check if this is a chain identifier
          CH = TOKEN(ITOKEN)(1:1)
          IF (TOKEN(ITOKEN)(2:2).EQ.' ' .AND.
     -        LGE(CH,'A') .AND. LLE(CH,'Z')) THEN
              CHAIN = CH

C----         If this is the first token of this chain, treat
C             as though an "ALL" entry
              IF (ITOKEN - LTOKEN.EQ.1) THEN
                  TOKEN(ITOKEN) = 'ALL'
                  TOKCHN(ITOKEN) = CHAIN

C----         Otherwise, update all prior tokens with this chain
              ELSE
                  DO 600, JTOKEN = LTOKEN + 1, ITOKEN - 1
                      IF (TOKCHN(JTOKEN).EQ.' ')
     -                    TOKCHN(JTOKEN) = CHAIN
 600              CONTINUE

C----             Mark the current token for deletion
                  TOKCHN(ITOKEN) = 'X'
                  TOKEN(ITOKEN) = 'XXXXX'
              ENDIF
              LTOKEN = ITOKEN

C----     If it is one of the Reserved words, then leave alone
          ELSE IF (TOKEN(ITOKEN).EQ.'ALL  ' .OR.
     -        TOKEN(ITOKEN).EQ.'FIRST' .OR.
     -        TOKEN(ITOKEN).EQ.'LAST ') THEN

C----     If this is a 'TO' string, then just mark for deletion
          ELSE IF (TOKEN(ITOKEN).EQ.'TO   ') THEN
              TOKCHN(ITOKEN) = 'X'
              TOKEN(ITOKEN) = 'XXXXX'

C----     Otherwise, assume this is a residue- or model-number, so
C         shift it into the correct column positions
          ELSE

C----         Right-justify the token string
              ILEN = LENSTR(TOKEN(ITOKEN))
              DUMMY = ' '
              DUMMY(TOKLEN - ILEN + 1:) = TOKEN(ITOKEN)
              TOKEN(ITOKEN) = DUMMY

C----         If this is not a model-number, or the last character is
C             a letter of an insertion code, take that to be an
C             insertion code, so keep whole token right-justified
              CH = TOKEN(ITOKEN)(TOKLEN:TOKLEN)
              IF (MODEL .OR. (LGE(CH,'A') .AND. LLE(CH,'Z'))) THEN

C----         Otherwise, shift back a character
              ELSE
                  DUMMY = TOKEN(ITOKEN)(2:)
                  TOKEN(ITOKEN) = DUMMY
              ENDIF
          ENDIF
 800   CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE DELTOK  -  Delete unwanted tokens
C
C----------------------------------------------------------------------+---

      SUBROUTINE DELTOK(TOKEN,TOKCHN,MXTOKN,NTOKEN)

      INTEGER       MXTOKN

      CHARACTER*1   TOKCHN(MXTOKN)
      CHARACTER*(*) TOKEN(MXTOKN)
      INTEGER       ITOKEN, JTOKEN, NTOKEN

C---- Initialise location of next free token
      JTOKEN = 0

C---- Loop through all the tokens, overwriting any deleted ones
      DO 1000, ITOKEN = 1, NTOKEN

C----     If token a valid one, then move to next available spot
          IF (TOKEN(ITOKEN).NE.'XXXXX') THEN
              JTOKEN = JTOKEN + 1
              IF (JTOKEN.LT.ITOKEN) THEN
                  TOKEN(JTOKEN) = TOKEN(ITOKEN)
                  TOKCHN(JTOKEN) = TOKCHN(ITOKEN)
              ENDIF
          ENDIF
 1000 CONTINUE

C---- Save the adjusted number of tokens stored
      NTOKEN = MIN(JTOKEN,MXTOKN)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE STOTOK  -  Store all the tokens representing current model-
C                        or residue-range
C
C----------------------------------------------------------------------+---

      SUBROUTINE STOTOK(TOKEN,TOKCHN,MXTOKN,NTOKEN,MODEL,FNAME,IREC,
     -    LINE,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,NRANGE,MRANGE,IFAIL)

      INTEGER       MAXRNG, MXTOKN

      CHARACTER*1   TOKCHN(MXTOKN)
      CHARACTER*6   RESFRM(MAXRNG), RESTO(MAXRNG)
      CHARACTER*(*) FNAME, IREC, TOKEN(MXTOKN)
      INTEGER       IVAL1, IVAL2, LENSTR, LINE, MODFRM(MAXRNG),
     -              MODTO(MAXRNG), MRANGE, NRANGE,
     -              NTOKEN
      LOGICAL       IFAIL, MODEL, VALID

C---- If the tokens are for a model range, then store first and last
C     model numbers
      IF (MODEL .AND. NTOKEN.GT.0) THEN
          VALID = .TRUE.

C----     If all models required, then set range as all-encompassing
          IF (TOKEN(1).EQ.'ALL  ') THEN
              IVAL1 = -99999
              IVAL2 = 99999
              MODFRM(1) = IVAL1
              MODTO(1) = IVAL2

C----     Otherwise, interpret model-number range
          ELSE

C----         If FIRST model defined, store start
              IF (TOKEN(1).EQ.'FIRST') THEN
                  IVAL1 = -99999
              ELSE
                  READ(TOKEN(1),1020,ERR=1060) IVAL1
 1020             FORMAT(I5)
                  GO TO 1100
 1060             CONTINUE
                  PRINT*, '*** Warning. Invalid model number: [',
     -                TOKEN(1), '] in line', LINE, '  of file',
     -                FNAME(1:LENSTR(FNAME))
                  PRINT*, '***          Line ignored: [',
     -                IREC(1:LENSTR(IREC))
                  VALID = .FALSE.
              ENDIF

C----         Interpret end model-number of range
 1100         CONTINUE

C----         If LAST model defined, store end
              IF (TOKEN(2).EQ.'LAST') THEN
                  IVAL2 = 99999

C----         If no second model number, then range consists of a
C             single model number
              ELSE IF (TOKEN(2).EQ.'XXXXX') THEN
                  IVAL2 = IVAL1

C----         Otherwise assume it is a model number
              ELSE
                  READ(TOKEN(2),1020,ERR=1160) IVAL2
                  GO TO 1200
 1160             CONTINUE
                  PRINT*, '*** Warning. Invalid model number: [',
     -                TOKEN(2), '] in line', LINE, '  of file',
     -                FNAME(1:LENSTR(FNAME))
                  PRINT*, '***          Line ignored: [',
     -                IREC(1:LENSTR(IREC))
                  VALID = .FALSE.
 1200             CONTINUE
              ENDIF
          ENDIF

C----     If have a valid range, then store
          IF (VALID) THEN
              MRANGE = MRANGE + 1
              IF (MRANGE.GT.MAXRNG) GO TO 910
              MODFRM(MRANGE) = IVAL1
              MODTO(MRANGE) = IVAL2
          ENDIF

C---- Otherwise, store the residue range
      ELSE IF (NTOKEN.GT.0) THEN
          NRANGE = NRANGE + 1
          IF (NRANGE.GT.MAXRNG) GO TO 908
          RESFRM(NRANGE) = TOKCHN(1) // TOKEN(1)
          RESTO(NRANGE) = TOKCHN(2) // TOKEN(2)
      ENDIF

      GO TO 999

C---- Fatal errors
 908  CONTINUE
      PRINT*, '**** Maximum number of residue ranges exceeded in file ',
     -    FNAME(1:LENSTR(FNAME)), '  MAXRNG =', MAXRNG
      GO TO 990

 910  CONTINUE
      PRINT*, '**** Maximum number of model ranges exceeded in file ',
     -    FNAME(1:LENSTR(FNAME)), '  MAXRNG =', MAXRNG
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PRNRNG  -  Print out the model- and residue ranges
C
C----------------------------------------------------------------------+---

      SUBROUTINE PRNRNG(FNAME,MODFRM,MODTO,RESFRM,RESTO,MAXRNG,NRANGE,
     -    MRANGE,BOTHND)

      INTEGER       MAXRNG

      CHARACTER*6   RESFRM(MAXRNG), RESTO(MAXRNG)
      CHARACTER*(*) FNAME
      CHARACTER*80  IREC
      INTEGER       ILEN, IRANGE, LENSTR, MODFRM(MAXRNG), MODTO(MAXRNG),
     -              MRANGE, NRANGE
      LOGICAL       BOTHND,VALID

C---- Print heading
      PRINT 20, FNAME(1:LENSTR(FNAME))
 20   FORMAT('*',/,'*',/,
     -    '*  Model and residue ranges selected in file: ',A,/,
     -    '*  ------------------------------------------',/,'*')

C---- Model-number range
      IF (MODFRM(1).EQ.-99999 .AND. MODTO(1).EQ.99999) THEN
          PRINT 120, '*      ALL models in the ensemble'
 120      FORMAT(A)
      ELSE
          DO 400, IRANGE = 1, MRANGE

C----         Retrieve start of model-number range
              IF (MODFRM(IRANGE).EQ.-99999) THEN
                  IREC = 'From FIRST model'
              ELSE
                  WRITE(IREC,140) MODFRM(IRANGE)
 140              FORMAT('From model ',I5)
              ENDIF

C----         Retrieve end of model-number range
              IF (MODTO(IRANGE).EQ.99999) THEN
                  IREC(17:) = ' to LAST model'
              ELSE
                  WRITE(IREC(17:),160) MODTO(IRANGE)
 160              FORMAT(' to model ',I5)
              ENDIF

C----         Show the interpreted range
              PRINT 180, '*      ', IREC(1:31)
 180          FORMAT(2A)
 400      CONTINUE
      ENDIF
      PRINT 120, '* '

C---- Residue-number range
      IF (RESFRM(1).EQ.'*ALL  ') THEN
          PRINT 120, '*      ALL residues'
      ELSE
          DO 800, IRANGE = 1, NRANGE

C----         Retrieve start of residue-number range
              IF (RESFRM(IRANGE)(2:).EQ.'ALL  ') THEN
                  IREC = 'All residues'
                  ILEN = 12
              ELSE IF (RESFRM(IRANGE)(2:).EQ.'FIRST') THEN
                  IREC = 'From FIRST residue'
                  ILEN = 18
              ELSE
                  WRITE(IREC,640) RESFRM(IRANGE)(2:)
 640              FORMAT('From residue ',A5)
                  ILEN = 18
              ENDIF

C----         Add chain identifier, if required
              IF (RESFRM(IRANGE)(1:1).NE.' ') THEN
                  IREC(ILEN + 1:) = ' in chain [' //
     -                RESFRM(IRANGE)(1:1) // ']'
                  ILEN = ILEN + 13
              ENDIF

C----         Retrieve end of residue-number range
              VALID = .FALSE.
              IF (RESTO(IRANGE)(2:).EQ.'LAST ') THEN
                  IREC(ILEN + 1:) = ' to LAST residue'
                  ILEN = ILEN + 16
                  VALID = .TRUE.
              ELSE IF (RESTO(IRANGE)(2:).NE.'XXXXX') THEN
                  WRITE(IREC(ILEN + 1:),660) RESTO(IRANGE)(2:)
 660              FORMAT(' to residue ',A5)
                  ILEN = ILEN + 17
                  VALID = .TRUE.
              ENDIF

C----         Add chain identifier, if required
              IF (VALID .AND. RESTO(IRANGE)(1:1).NE.' ') THEN
                  IREC(ILEN + 1:) = ' in chain [' //
     -                RESTO(IRANGE)(1:1) // ']'
                  ILEN = ILEN + 13
              ENDIF

C----         Show the interpreted range
              PRINT 180, '*      ', IREC(1:ILEN)
 800      CONTINUE
      ENDIF
      PRINT 120, '* '
      IF (BOTHND) THEN
          PRINT 120, '* Restraints only where both residues are sele' //
     -        'cted'
          PRINT 120, '* '
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETWNT  -  For an NMR ensemble, identify which files contain
C                        the wanted models (ie those specified by the user
C                        in the ranges file)
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETWNT(MXFILE,NFILE,FILRIN,FFILE,MAXRNG,MODFRM,MODTO,
     -    MRANGE,RSELEC,MWANT,NMODEL,MODNUM,ACTNUM,TOPMOD,TITLE,NAMLEN,
     -    TLEN,IFAIL)

      INTEGER       MAXRNG, MXFILE, NFILE

      CHARACTER*3   NUMBER
      CHARACTER*78  TITLE
      CHARACTER*80  FILRIN(MXFILE), IREC
CHECK v.3.4.3-->
C      INTEGER       ACTNUM(MXFILE), FFILE, IERR, IFAIL, IFILE, IMODEL,
      INTEGER       ACTNUM(MXFILE), FFILE, IERR, IFILE, IMODEL,
CHECK v.3.4.3<--
     -              IPOS, MODFRM(MAXRNG), MODNUM(MXFILE), MODTO(MAXRNG),
     -              MRANGE, NAMLEN, NMODEL, TLEN, TOPMOD
CHECK v.3.4.3-->
C      LOGICAL       FOUND, INMODL, MWANT(MXFILE), RSELEC, WANTED
      LOGICAL       FOUND, IFAIL, INMODL, MWANT(MXFILE), RSELEC, WANTED
CHECK v.3.4.3<--

C---- Initialise file number
      FOUND = .FALSE.
      FFILE = 0
      NMODEL = 0
      TOPMOD = -99

C---- Loop through all the .rin files
      DO 500, IFILE = 1, NFILE

C----     Open the corresponding file
          OPEN(UNIT=4,FILE=FILRIN(IFILE),STATUS='OLD',
     -        ACCESS='SEQUENTIAL',
CVAX     -        READONLY,
     -        FORM='FORMATTED',ERR=900)

C----     Read the first record
          READ(4,220,END=900) IREC
 220      FORMAT(A)

C----     Check whether the first record contains the actual model
C         number
          READ(4,320,END=902) IREC
 320      FORMAT(A)
          IF (IREC(1:5).EQ.'MODEL') THEN
             READ(IREC,360,IOSTAT=IERR) IMODEL
 360         FORMAT(9X,I5)
             IF (IERR.NE.0) THEN
                IMODEL = IFILE
             ENDIF
          ELSE
             IMODEL = IFILE
          ENDIF

C----     Check whether this model is one of the ones wanted
          WANTED = INMODL(IMODEL,MODFRM,MODTO,MAXRNG,MRANGE)
          MWANT(IFILE) = WANTED
          IF (WANTED) THEN
              FOUND = .TRUE.
              IF (FFILE.EQ.0) FFILE = IFILE
              NMODEL = NMODEL + 1
              MODNUM(NMODEL) = IFILE
              ACTNUM(NMODEL) = IMODEL
              TOPMOD = MAX(IMODEL,TOPMOD)
          ENDIF

C----     Close the file
          CLOSE(4)

C---- Loop back if no hits found
 500  CONTINUE

C---- If still no hits found, then don't have anything to plot!
      IF (.NOT.FOUND) GO TO 904

C---- If have selected fewer models than there are in the ensemble, then
C     modify general plot subtitle
      IF (NMODEL.LT.NFILE) THEN
          WRITE(NUMBER,520) NMODEL
 520      FORMAT(I3)
          IPOS = 1
          IF (NUMBER(1:1).EQ.' ') IPOS = 2
          IF (NUMBER(2:2).EQ.' ') IPOS = 3
          TITLE = TITLE(1:NAMLEN + 2) // NUMBER(IPOS:) // ' of ' //
     -         TITLE(NAMLEN + 3:)
          TLEN = TLEN + 8 - IPOS
      ENDIF

C---- If there is a selection on the residue ranges, then identify this
C     also on the general plot subtitle
      IF (RSELEC) THEN
          TITLE = TITLE(1:TLEN) // '**'
          TLEN = TLEN + 2
      ENDIF

      GO TO 999

 900  CONTINUE
      PRINT*, '*** Error opening .rin file.'
      PRINT*, FILRIN(IFILE)
      GO TO 990

 902  CONTINUE
      PRINT*, '*** Data error reading .rin file.'
      PRINT*, FILRIN(IFILE)
      GO TO 990

 904  CONTINUE
      PRINT*, '*** ERROR. No models within supplied ranges found!'
      GO TO 990

990   CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION INMODL  -  Check whether supplied model-number is within one or
C                      more of the ranges defined by the user
C
C----------------------------------------------------------------------+---

      LOGICAL FUNCTION INMODL(IMODEL,MODFRM,MODTO,MAXRNG,MRANGE)

      INTEGER       MAXRNG

      INTEGER       IMODEL, IRANGE, MODFRM(MAXRNG), MODTO(MAXRNG),
     -              MRANGE
      LOGICAL       MATCHF, MATCHL, WANTED


C---- Initialise
      WANTED = .FALSE.

C---- Loop through all the defined model ranges and check model number
C     against each one
      DO 500, IRANGE = 1, MRANGE

C----     Initialise flags for this range (first and last parameters)
          MATCHF = .FALSE.
          MATCHL = .FALSE.

C----     Check for lower-bound
          IF (IMODEL.GE.MODFRM(IRANGE)) THEN
              MATCHF = .TRUE.
          ENDIF

C----     Check for upper-bound
          IF (IMODEL.LE.MODTO(IRANGE)) THEN
              MATCHL = .TRUE.
          ENDIF

C----     If model satisfies both the first and last parameters, then
C         it is wanted
          IF (MATCHF .AND. MATCHL) WANTED = .TRUE.
500   CONTINUE

C---- Return the result
      INMODL = WANTED

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION INRANG  -  Check whether supplied residue is within one or
C                      more of the ranges defined by the user
C
C----------------------------------------------------------------------+---

      LOGICAL FUNCTION INRANG(RESNUM,CHAIN,RESFRM,RESTO,MAXRNG,NRANGE)

      INTEGER       MAXRNG

      CHARACTER*1   CH, CHAIN
      CHARACTER*5   DUMMY, RESNUM
      CHARACTER*6   RESFRM(MAXRNG), RESID, RESTO(MAXRNG)
      INTEGER       ILEN, IRANGE, LENSTR, NRANGE
      LOGICAL       MATCHF, MATCHL, WANTED


C---- If all residues selected, then don't need to perform any more tests
      IF (RESFRM(1).EQ.'*ALL  ') THEN
          WANTED = .TRUE.

C---- Otherwise, test whether this residue is required
      ELSE

C----     Right-justify the entered number
          ILEN = LENSTR(RESNUM)
          DUMMY = ' '
          DUMMY(5 - ILEN + 1:) = RESNUM
          RESNUM = DUMMY

C----     If the last character is a letter of an insertion code, take
C         that to be an insertion code, so keep whole token right-justified
          CH = RESNUM(5:5)
          IF (LGE(CH,'A') .AND. LLE(CH,'Z')) THEN

C----     Otherwise, shift back a character
          ELSE
              DUMMY = RESNUM(2:)
              RESNUM = DUMMY
          ENDIF

C----     Initialise
          WANTED = .FALSE.
          RESID = CHAIN // RESNUM

C----     Loop through all the defined residue ranges and check residue
C         number and chain against each one
          DO 500, IRANGE = 1, NRANGE

C----         Initialise flags for this range (first and last parameters)
              MATCHF = .FALSE.
              MATCHL = .FALSE.

C----         Check for ALL of a given chain
              IF (RESFRM(IRANGE)(2:6).EQ.'ALL  ') THEN

C----             Check that is the correct chain
                  IF (CHAIN.EQ.RESFRM(IRANGE)(1:1)) THEN
                      MATCHF = .TRUE.
                      MATCHL = .TRUE.
                  ENDIF

C----         Otherwise, check whether range is from FIRST residue
              ELSE IF (RESFRM(IRANGE)(2:6).EQ.'FIRST') THEN

C----             Check that this residue is of a higher chain than that
C                 stated
                  IF (LGE(CHAIN,RESFRM(IRANGE)(1:1))) THEN
                      MATCHF = .TRUE.
                  ENDIF

C----         Otherwise check that the residue number is indeed higher
C             than (or equal to) the number defining the start of the
C             range
              ELSE
                  IF (LGE(RESID,RESFRM(IRANGE))) THEN
                      MATCHF = .TRUE.
                  ENDIF
              ENDIF

C----         Check whether range-end is from LAST residue
              IF (RESTO(IRANGE)(2:6).EQ.'LAST ') THEN

C----             Check that this residue is of a lower chain than that
C                 stated
                  IF (LLE(CHAIN,RESTO(IRANGE)(1:1))) THEN
                      MATCHL = .TRUE.
                  ENDIF

C----         If only a single residue defined, check whether this is the
C             one
              ELSE IF (RESTO(IRANGE)(2:6).EQ.'XXXXX') THEN
                  IF (LLE(RESID,RESFRM(IRANGE))) THEN
                      MATCHL = .TRUE.
                  ENDIF

C----         Otherwise check that the residue number is lower
C             than (or equal to) the number defining the end of the
C             range
              ELSE
                  IF (LLE(RESID,RESTO(IRANGE))) THEN
                      MATCHL = .TRUE.
                  ENDIF
              ENDIF

C----         If residue satisfies both the first and last parameters,
C             then it is wanted
              IF (MATCHF .AND. MATCHL) WANTED = .TRUE.
500       CONTINUE
      ENDIF

C---- Return the result
      INRANG = WANTED

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PINTIC  -  Show ticks every so many residues
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PINTIC(XPORIG,YPORIG,YDIFF,XPWID,MXPINS,NPIN,TSIZE,
     -    TICPOS,BRKPOS,NOTICK)

      INTEGER       IRES, IMARK, MXPINS, NPIN
      LOGICAL       BRKPOS(MXPINS), NOTICK, TICPOS(MXPINS)
      REAL          SCALEX, TSIZE, X, XGAP, XL, XPORIG, XPWID, XR, Y,
     -              YDIFF, YPORIG
 
C---- Plot the ticks
      XGAP = XPWID / MXPINS
      SCALEX = XPWID / MXPINS
      X = XPORIG
      Y = YPORIG - YDIFF
      IRES = 0
      CALL PSLWID(0.2)

C---- Plot ticks at appropriate positions
      XL = X
      XR = X
      X = X + XGAP / 2.0
      DO 100, IMARK = 1, NPIN

C----     Draw tick at this position, if required
          IF (.NOT.NOTICK .AND. TICPOS(IMARK)) THEN
              CALL PSLINE(X,Y,X,Y + TSIZE)
          ENDIF

C----     Increment x-position
          X = X + XGAP

C----     Build up line up to break
          IF (.NOT.BRKPOS(IMARK)) THEN
              XR = XR + XGAP

C----     If have come to line-break, draw line so far
          ELSE
              IF (XL.NE.XR) CALL PSLINE(XL,Y,XR,Y)
              XL = XR + XGAP
              XR = XL
          ENDIF
 100  CONTINUE

C---- Draw last bit of line
      IF (XR.GT.XL) THEN
          CALL PSLINE(XL,Y,XR,Y)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PINRNO  -  Draw residue numbers along x-axis every 5 or so
C                        residues
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE PINRNO(XLIM1,XLIM2,YLIM1,PINSNO,MXPINS,NPIN,FPIN,
     -    TICPOS,TICLAB,TICLEN,SIZLAB,TSIZE)

      INTEGER        MXPINS
      CHARACTER*5    LABEL, PINSNO(MXPINS), TICLAB(MXPINS)
      INTEGER        FPIN, IFROM, IPIN, ITO, LASTNO, LSTPOS, NBLANK,
     -               NPIN, NUMBER, PGAP, PSTEP, TICLEN(MXPINS)
      LOGICAL        LPRINT, TICPOS(MXPINS)
      REAL           SIZLAB, TSIZE, X, XGAP, XLIM1, XLIM2, Y,
     -               YLIM1

      DATA  PSTEP  / 5 /

C---- Initialise variables
      LASTNO = FPIN
      LSTPOS = -999
      NBLANK = 0
      PGAP = PSTEP
      DO 100, IPIN = 1, MXPINS
          TICLAB(IPIN) = ' '
          TICLEN(IPIN) = 0
          TICPOS(IPIN) = .FALSE.
 100  CONTINUE
      XGAP = (XLIM2 - XLIM1) / MXPINS
      X = XLIM1 + XGAP / 2.0
      Y = YLIM1 - TSIZE / 2.0

C---- Loop through all the residues in the current graph
      DO 200, IPIN = 1, NPIN

C----     Increment count of blanks
          NBLANK = IPIN - LSTPOS
          NUMBER = -999

C----     Initialise flag determining whether to print this residue number or
C         not
          LPRINT = .FALSE.

C----     If have gone given number of residues without printing, then print
          IF (NBLANK.EQ.PSTEP) THEN
              LPRINT = .TRUE.

C----     Otherwise, pick off the residue number, if can
          ELSE
              READ(PINSNO(IPIN),120,ERR=130) NUMBER
 120          FORMAT(1X,I3,1X)

C----         If sequence number is a multiple of the required gap, then
C             want to print it
              IF (MOD(NUMBER,PSTEP).EQ.0) THEN

C----             Check that not too close to last one
                  IF (NBLANK.GE.PGAP) THEN
                      LPRINT = .TRUE.

C----             If too close, then reset last position printed
                  ELSE
                      LSTPOS = IPIN
                  ENDIF

C----         Otherwise, if odd gap in sequence numbering, then print
              ELSE IF (NUMBER.NE.LASTNO + 1) THEN

C----             Check that not too close to last one
                  IF (NBLANK.GE.PGAP) THEN
                      LPRINT = .TRUE.
                  ENDIF
              ENDIF
              LASTNO = NUMBER
 130          CONTINUE
          ENDIF

C----     Check that not at chain-break
          IF (PINSNO(IPIN)(2:4).EQ.' ') THEN
              LPRINT = .FALSE.
              LASTNO = 0
          ENDIF

C----     If residue number is to be printed, print it
          IF (LPRINT) THEN
CHECK v.3.5-->
C              LABEL = PINSNO(IPIN)(2:5)
C              IFROM = 1
C              ITO = 4
C              IF (LABEL(1:1).EQ.' ') IFROM = 2
C              IF (LABEL(2:2).EQ.' ') IFROM = 3
C              IF (LABEL(4:4).EQ.' ') ITO = 3
              LABEL = PINSNO(IPIN)(1:5)
              IFROM = 1
              ITO = 5
              IF (LABEL(1:1).EQ.' ') IFROM = 2
              IF (LABEL(2:2).EQ.' ') IFROM = 3
              IF (LABEL(3:3).EQ.' ') IFROM = 4
              IF (LABEL(4:4).EQ.' ') ITO = 4
CHECK v.3.5<--
              CALL PSCTXT(X,Y,SIZLAB,LABEL(IFROM:ITO))
              CALL PSLINE(X,YLIM1,X,YLIM1 - SIZLAB / 3.0)

C----         Reset counts
              NBLANK = 0
              LSTPOS = IPIN
              TICPOS(IPIN) = .TRUE.
              TICLAB(IPIN) = LABEL(IFROM:ITO)
              TICLEN(IPIN) = ITO - IFROM + 1
          ENDIF
          X = X + XGAP
 200  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE DHELIX  -  Draw the helix segment
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE DHELIX(LENRUN,X,Y,IRES,DX,HWIDTH,HEIGHT,SHOWAC,INCOLR,
     -    MXCOLR,RGB,COLOUR)

      SAVE

      INTEGER       BLANK, COIL, HELIX, STRAND
      PARAMETER    (
     -              BLANK  =   0,
     -              COIL   =   3,
     -              HELIX  =   1,
     -              STRAND =   2
     -             )

      REAL          DSHADE, HSHADE
      PARAMETER    (
     -              DSHADE = 0.3,
     -              HSHADE = 0.7
     -             )

      INTEGER       COLOUR, IHALF, IRES, ITURN, LENRUN, MXCOLR, NHALF,
     -              UPDOWN
      LOGICAL       INCOLR, SHOWAC
      REAL          DH, DX, HWIDTH, HEIGHT, RGB(3,MXCOLR), TWIDTH, X,
     -              X1, X2, X3, X4, X5, X6, X7, X8, XM1, XM2, XTRA, Y,
     -              Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, YM

C---- Helix representation:
C
C   Half-turns                          5_____7 4
C                                     2 /\     \ 
C                                      /  \     \
C   UPDOWN = 0                     YM /_ _ \     \ 
C                                    /     /\     \
C      Y   1 ______3  6______ 8     /_____/  \_____\
C            \     \  /     /      1      3  6     8
C             \     \/ _ _ / YM     XM1   XM2
C              \     \    /                    UPDOWN = 1
C               \     \  / 
C              2 \_____\/4
C                5     7
C                   XM1   XM2
C
C  Triangle defined by XM1, XM2, YM has slightly darker shading (DSHADE)
C  than the remainder of the inside of the helix

C---- Calculate the length of each half-turn so that get approximately 3.5
C     residues per turn for given run-length
      NHALF = 2.0 * REAL(LENRUN) / 3.5 + 0.5
      IF (NHALF.LT.1) NHALF = 1
      TWIDTH = HWIDTH
      DH = (DX * REAL(LENRUN) - HWIDTH) / REAL(NHALF)
      IF (TWIDTH.GT.DH) THEN
CHECK v.3.5-->
C          DH = DX / 2.0
C          TWIDTH = DH
          IF (LENRUN.EQ.1) THEN
              TWIDTH = DX * 0.5
          ELSE
              TWIDTH = DX * 0.75
          ENDIF
          DH = (DX * REAL(LENRUN) - TWIDTH) / REAL(NHALF)
CHECK v.3.5<--
      ENDIF
      ITURN = IRES
      XTRA = HEIGHT * TWIDTH / DH
      X = X + TWIDTH / 2.0

C---- Loop through the given number of turns
      DO 500, IHALF = 1, NHALF

C----     Initialise direction of helix
          ITURN = ITURN + 1
          UPDOWN = 2 * MOD(ITURN,2) - 1

C----     Draw middle portion of turn
C----     First segment - calculate coordinates
          X1 = X - TWIDTH / 2.0
          X2 = X + (DH - TWIDTH) / 2.0
          Y1 = Y
          Y2 = Y + REAL(UPDOWN) * HEIGHT
          X3 = X + TWIDTH / 2.0
          X4 = X + DH - (DH - TWIDTH) / 2.0
          Y3 = Y
          Y4 = Y + REAL(UPDOWN) * HEIGHT

C----     Second segment - calculate coordinates
          X5 = X + (DH - TWIDTH) / 2.0
          X6 = X + DH - TWIDTH / 2.0
          Y5 = Y + REAL(UPDOWN) * HEIGHT
          Y6 = Y
          X7 = X + DH - (DH - TWIDTH) / 2.0
          X8 = X + DH + TWIDTH / 2.0
          Y7 = Y + REAL(UPDOWN) * HEIGHT
          Y8 = Y

C----     Adjust for cover-up of nearest bit of turn over furthest bit
          IF (UPDOWN.GT.0) THEN
              X4 = X + DH / 2.0
              Y4 = Y + HEIGHT - XTRA
              XM1 = X4 - TWIDTH
              XM2 = X4
              YM = Y4
          ELSE
              X5 = X + DH / 2.0
              Y5 = Y - HEIGHT + XTRA
              XM2 = X5 + TWIDTH
              XM1 = X5
              YM = Y5
          ENDIF

CTESTING
C      X1 = 6.0 * X1 - 400.0
C      X2 = 6.0 * X2 - 400.0
C      X3 = 6.0 * X3 - 400.0
C      X4 = 6.0 * X4 - 400.0
C      X5 = 6.0 * X5 - 400.0
C      X6 = 6.0 * X6 - 400.0
C      X7 = 6.0 * X7 - 400.0
C      X8 = 6.0 * X8 - 400.0
C      XM1 = 6.0 * XM1 - 400.0
C      XM2 = 6.0 * XM2 - 400.0
C      Y1 = 6.0 * Y1 - 300.0
C      Y2 = 6.0 * Y2 - 300.0
C      Y3 = 6.0 * Y3 - 300.0
C      Y4 = 6.0 * Y4 - 300.0
C      Y5 = 6.0 * Y5 - 300.0
C      Y6 = 6.0 * Y6 - 300.0
C      Y7 = 6.0 * Y7 - 300.0
C      Y8 = 6.0 * Y8 - 300.0
C      YM = 6.0 * YM - 300.0
CTESTING

C----     Blank out the region that contains this turn of the helix
          CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
          CALL PSUBOX(X1,Y1,X2,Y2,X4,Y4,X3,Y3)
          CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
          CALL PSUBOX(X5,Y5,X6,Y6,X8,Y8,X7,Y7)

C----     Shade in inside of helix
          IF (UPDOWN.GT.0) THEN
              IF (.NOT.SHOWAC) THEN
                  CALL PSHADE(HSHADE,COLOUR,RGB,MXCOLR,INCOLR)
              ELSE
                  CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
              ENDIF
              CALL PSUBOX(X1,Y1,X3,Y3,XM2,YM,XM1,YM)
              IF (.NOT.SHOWAC) THEN
                  CALL PSHADE(DSHADE,COLOUR,RGB,MXCOLR,INCOLR)
              ELSE
                  CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
              ENDIF
              CALL PSUBOX(XM1,YM,XM2,YM,X4,Y4,X2,Y2)
          ELSE
              IF (.NOT.SHOWAC) THEN
                  CALL PSHADE(DSHADE,COLOUR,RGB,MXCOLR,INCOLR)
              ELSE
                  CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
              ENDIF
              CALL PSUBOX(X5,Y5,XM1,YM,XM2,YM,X7,Y7)
              IF (.NOT.SHOWAC) THEN
                  CALL PSHADE(HSHADE,COLOUR,RGB,MXCOLR,INCOLR)
              ELSE
                  CALL PSHADE(1.0,COLOUR,RGB,MXCOLR,INCOLR)
              ENDIF
              CALL PSUBOX(XM1,YM,XM2,YM,X8,Y8,X6,Y6)
          ENDIF

C----     Draw start of helix
          IF (IHALF.EQ.1) THEN
              CALL PSLINE(X1,Y1,X3,Y3)
          ENDIF

C----     Draw helix edges
          CALL PSLINE(X1,Y1,X2,Y2)
          CALL PSLINE(X3,Y3,X4,Y4)
          CALL PSLINE(X5,Y5,X6,Y6)
          CALL PSLINE(X7,Y7,X8,Y8)
          CALL PSLINE(X2,Y2,X7,Y7)

C----     Increment x-coordinate
          X = X + DH
 500  CONTINUE

C---- Draw the helix end
      CALL PSLINE(X - TWIDTH / 2.0,Y,X + TWIDTH / 2.0,Y)
      X = X + TWIDTH / 2.0
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4<--
CHECK v.3.3.2-->
C**************************************************************************
C
C  SUBROUTINE OPEHTM  -  Open the html file, reading through to the
C                        end if the file already exists
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPEHTM(FNAME,NEWFIL,IFAIL)

      CHARACTER*80   IREC
      CHARACTER*(*)  FNAME
      LOGICAL        IFAIL, NEWFIL

C---- Initialise variables
      IFAIL = .FALSE.

C---- Open the summary file
      OPEN(UNIT=15,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
CHECK v.3.4.3-->
CCVAX     -    CARRIAGECONTROL='LIST',
CVAX     -    CARRIAGECONTROL='LIST',RECL=160,
CHECK v.3.4.3<--
     -    ERR=900)

C---- If this is an existing file, then read through all its records until
C     get to the end of file
      IF (.NOT.NEWFIL) THEN

C----     Loop through the file until reach the end
 100      CONTINUE
             READ(15,110,END=500) IREC
 110         FORMAT(A)
          GO TO 100

C---- If this is a new file, then write the header records to it
      ELSE
          WRITE(15,120)
 120      FORMAT(/,'<h2>PROCHECK statistics</h2>')
      ENDIF

C---- End of file reached
 500  CONTINUE

      GO TO 999

 900  CONTINUE
      PRINT*, '*** ERROR opening summary file'
      PRINT*, '* ', FNAME
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.
      GO TO 999

 999  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------+---
CHECK v.3.3.2<--
CHECK v.3.4.3-->
C**************************************************************************
C
C  SUBROUTINE GETDAT  -  Read in the torsion angle distributions from the
C                        prodata file.
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      SUBROUTINE GETDAT(DISTRB,TWODEE,NOBSER,NCELL,NCELL1,NAMINO,VALBEG,
C     -    VALEND,STEP,NCOUNT,NRMEAN,NRMSTD,FILNUM,IFAIL)
      SUBROUTINE GETDAT(DISTRB,TWODEE,NOBSER,MXCELL,NCELL1,NCELL2,
     -    NAMINO,VALBEG,VALEND,STEP,NCOUNT,NRMEAN,NRMSTD,FILNUM,NONORM,
     -    IFAIL)
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       NCELL, NCELL1, NAMINO
      INTEGER       MXCELL, NCELL1, NCELL2, NAMINO
CHECK v.3.5.2<--

      CHARACTER*1   CH
CTESTING
      CHARACTER*3   AMNAME(21)
CTESTING
      CHARACTER*10  CNUMB, RNUMB
      CHARACTER*52  LETTER
      CHARACTER*80  IREC
      INTEGER       DISTRB, FILNUM, IAMINO, ICELL, ILOOP, IDISTR, IPOS,
CHECK v.3.5.2-->
C     -              ISTATE, IVALUE, JCELL, LINE, LVALUE,
C     -              NCOUNT(NAMINO +1), NPOS, NVALS, NZERO
     -              ISTATE, IVALUE, JCELL, LINE, LVALUE, NCELL,
     -              NCOUNT(NAMINO +1), NPOS, NVALS, NZERO, OFFSET,
     -              TCELL, TOTCEL
CHECK v.3.5.2<--
CTESTING
C      INTEGER       IOBSER(45)
CTESTING
CHECK v.3.5.2-->
C      LOGICAL       ENDFIL, GOTLET, GOTNUM, IFAIL, INEW, INORM, TWODEE,
C     -              UPDATE
C      REAL          MEAN, NOBSER(NCELL,NCELL,NAMINO+1),
      LOGICAL       ENDFIL, GOTLET, GOTNUM, IFAIL, INEW, INORM, NONORM,
     -              TWODEE, UPDATE
      REAL          MEAN, NOBSER(MXCELL*(NAMINO+1)),
CHECK v.3.5.2<--
     -              NRMEAN(NAMINO + 1), NRMSTD(NAMINO + 1), STDEV,
     -              STEP(2), VALBEG(2), VALEND(2)

      DATA LETTER( 1:26) / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LETTER(27:52) / 'abcdefghijklmnopqrstuvwxyz'/

CTESTING
      DATA AMNAME /'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
     -             'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
     -             'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ALL' /
CTESTING

C---- Rewind the compressed distributions file, procheck.dat
      REWIND(FILNUM)

C---- Initialise variables
      IDISTR = 0
      LINE = 0

C---- Define distribution ranges
      TWODEE = .FALSE.
      IF (DISTRB.EQ.1) THEN
          VALBEG(1) = -180.0
          VALBEG(2) = -180.0
          VALEND(1) =  180.0
          VALEND(2) =  180.0
          TWODEE = .TRUE.
      ELSE IF (DISTRB.EQ.2) THEN
          VALBEG(1) =    0.0
          VALBEG(2) =    0.0
          VALEND(1) =  360.0
          VALEND(2) =  360.0
          TWODEE = .TRUE.
      ELSE IF (DISTRB.EQ.3) THEN
          VALBEG(1) =    0.0
          VALEND(1) =  360.0
      ELSE IF (DISTRB.EQ.4) THEN
          VALBEG(1) =    0.0
          VALEND(1) =  360.0
      ELSE IF (DISTRB.EQ.5) THEN
          VALBEG(1) =    0.0
          VALEND(1) =  360.0
      ELSE IF (DISTRB.EQ.6) THEN
          VALBEG(1) =    0.0
          VALEND(1) =  360.0
CHECK v.3.5.2-->
      ELSE IF (DISTRB.EQ.7) THEN
          VALBEG(1) =    0.0
          VALBEG(2) = -180.0
          VALEND(1) =  180.0
          VALEND(2) =  180.0
          TWODEE = .TRUE.
CHECK v.3.5.2<--
      ENDIF

CHECK v.3.5.2-->
C---- Define the number of cells in the distribution
      IF (TWODEE) THEN
          IF (DISTRB.EQ.7) THEN
              NCELL1 = 60
              NCELL2 = 120
          ELSE
              NCELL1 = 45
              NCELL2 = 45
          ENDIF
      ELSE
          NCELL1 = 45 * 45
          NCELL2 = 1
      ENDIF

C---- Size of this distribution
      NCELL = NCELL1 * NCELL2
      IF (NCELL.GT.MXCELL) GO TO 900
      OFFSET = NCELL * NAMINO
      TOTCEL = MXCELL * (NAMINO + 1)
CHECK v.3.5.2<--

C---- Determine cell-size in each direction
      IF (TWODEE) THEN
CHECK v.3.5.2-->
C          STEP(1) = (VALEND(1) - VALBEG(1)) / NCELL
C          STEP(2) = (VALEND(2) - VALBEG(2)) / NCELL
          STEP(1) = (VALEND(1) - VALBEG(1)) / NCELL1
          STEP(2) = (VALEND(2) - VALBEG(2)) / NCELL2
CHECK v.3.5.2<--
      ELSE
          STEP(1) = (VALEND(1) - VALBEG(1)) / NCELL1
      ENDIF

C---- Search for the start of the data for the current distribution
 100  CONTINUE
          LINE = LINE + 1
          READ(FILNUM,120,ERR=904,END=902) IREC
 120      FORMAT(A)
          IF (IREC(1:13).EQ.'Distribution ') THEN
              READ(IREC,160,ERR=904) IDISTR
 160          FORMAT(13X,I1)
          ENDIF
      IF (IDISTR.NE.DISTRB) GO TO 100

C---- Initialise variables
      ENDFIL = .FALSE.
      INEW =  .TRUE.
      INORM = .FALSE.
      ISTATE = 0
      DO 300, IAMINO = 1, NAMINO + 1
          NCOUNT(IAMINO) = 0
          NRMEAN(IAMINO) = 0.0
          NRMSTD(IAMINO) = 0.0
CHECK v.3.5.2-->
C          DO 250, JCELL = 1, NCELL
C              DO 200, ICELL = 1, NCELL
C                  NOBSER(ICELL,JCELL,IAMINO) = 0.0
C 200          CONTINUE
C 250      CONTINUE
          DO 200, ICELL = 1, TOTCEL
              NOBSER(ICELL) = 0.0
 200      CONTINUE
CHECK v.3.5.2<--
 300  CONTINUE

CHECK v.3.5.2-->
C---- Initialise cell count
      ICELL = 0
CHECK v.3.5.2<--

C---- Read through the data file
 400  CONTINUE
          LINE = LINE + 1
          READ(FILNUM,120,END=800,ERR=904) IREC
          IF (IREC(1:13).EQ.'Distribution ') GO TO 800

C----     If looking for the current amino acid's normalization factors,
C         then read them in
          IF (INORM) THEN
CHECK v.3.5.2-->
C              IF (IREC(1:4).NE.'Norm') GO TO 916
              IF (IREC(1:4).EQ.'Norm') THEN
CHECK v.3.5.2<--
                  READ(IREC,430,ERR=904) MEAN, STDEV
 430              FORMAT(15X,2F10.0)

CHECK v.3.5.2-->
C----         If no "Norm" record present, then check whether we need
C             it
              ELSE IF (NONORM) THEN
                  MEAN = 0.0
                  STDEV = 0.0
                  INEW = .TRUE.
                  INORM = .FALSE.

C----         Otherwise, we need the normalization data but it is missing
              ELSE
                  GO TO 916
              ENDIF
C              INEW = .TRUE.
C              INORM = .FALSE.
CHECK v.3.5.2<--
              NRMEAN(IAMINO) = MEAN
              NRMSTD(IAMINO) = STDEV
CHECK v.3.5.2-->
          ENDIF
CHECK v.3.5.2<--

C----     If looking for next amino acid, read its details in
CHECK v.3.5.2-->
C          ELSE IF (INEW) THEN
          IF (INEW) THEN
CHECK v.3.5.2<--
              READ(IREC,440,ERR=904) IAMINO
 440          FORMAT(5X,I2)

C----         If have a valid amino acid, then initialise variables
              IF (IAMINO.GT.0 .AND.IAMINO.LE.NAMINO) THEN
                  INEW = .FALSE.
CHECK v.3.5.2-->
C                  ICELL = NCELL
                  TCELL = 0
CHECK v.3.5.2<--
                  ISTATE = 0
CHECK v.3.5.2-->
C                  JCELL = 0
CHECK v.3.5.2<--
                  LVALUE = 0
                  NVALS = 0
                  NZERO = 0
                  UPDATE = .FALSE.
              ELSE
                  GO TO 906
              ENDIF

C----     Otherwise, process the current record
CHECK v.3.5.2-->
C          ELSE
          ELSE IF (.NOT.INORM) THEN
CHECK v.3.5.2<--

C----         Loop through the characters in this record
              DO 600, IPOS = 1, 80
                  CH = IREC(IPOS:IPOS)

C----             Check if this is a letter and, if so, extract its value
                  IF ((LGE(CH,'A') .AND. LLE(CH,'Z')) .OR.
     -                (LGE(CH,'a') .AND. LLE(CH,'z'))) THEN
                      IVALUE = INDEX(LETTER,CH)
                      GOTLET = .TRUE.
                  ELSE
                      GOTLET = .FALSE.
                  ENDIF

C----             Check for a number
                  IF (LGE(CH,'0') .AND. LLE(CH,'9')) THEN
                      GOTNUM = .TRUE.
                  ELSE
                      GOTNUM = .FALSE.
                  ENDIF

C----             Process character depending on current state

C----             State 0: Looking at first character. Check to see what
C                 it is
                  IF (ISTATE.EQ.0) THEN

C----                 If a letter, then can update straight away
                      IF (GOTLET) THEN
                          UPDATE = .TRUE.

C----                 Otherwise, if the start of a number, will need to
C                     build it up
                      ELSE

C----                     If a hash-number, prepare to build it up
                          IF (CH.EQ.'#') THEN
                              ISTATE = 1

C----                     Otherwise, if it's the start of the number of
C                         zeros, get ready to build that up
                          ELSE IF (GOTNUM) THEN
                              ISTATE = 3
                              CNUMB = ' '
                              NPOS = 0
                          ENDIF
                      ENDIF
                  ENDIF

C----             State 1: Have a hash, so skip to next state for first
C                          digit
                  IF (ISTATE.EQ.1) THEN
                      ISTATE = 2
                      CNUMB = ' '
                      NPOS = 0

C----             State 2: Bulding up the hash number
                  ELSE IF (ISTATE.EQ.2) THEN

C----                 If this is a hash, an asterisk or a letter, then
C                     have come to the end of the number
                      IF (CH.EQ.'#' .OR. CH.EQ.'*' .OR. GOTLET) THEN
                          UPDATE = .TRUE.
                          IF (NPOS.GT.0 .AND. NPOS.LE.10) THEN
                              RNUMB = ' '
                              RNUMB(10 - NPOS + 1:10) = CNUMB(1:NPOS)
                              READ(RNUMB,480,ERR=910) LVALUE
 480                          FORMAT(I10)
                              UPDATE = .TRUE.
                          ELSE
                              GO TO 910
                          ENDIF

C----                 Otherwise, it must be a continuation of the
C                     repetition number
                      ELSE IF (GOTNUM) THEN
                          NPOS = NPOS + 1
                          IF (NPOS.LE.10) CNUMB(NPOS:NPOS) = CH
                      ELSE
                          GO TO 908
                      ENDIF

C----             State 3: Bulding up the number of zeros
                  ELSE IF (ISTATE.EQ.3) THEN

C----                 If this is a hash an asterisk or a letter,
C                     then have come to the end of the number
                      IF (CH.EQ.'#' .OR. CH.EQ.'*' .OR. GOTLET) THEN
                          IF (NPOS.GT.0 .AND. NPOS.LE.10) THEN
                              RNUMB = ' '
                              RNUMB(10 - NPOS + 1:10) = CNUMB(1:NPOS)
                              READ(RNUMB,480,ERR=910) NZERO
                              UPDATE = .TRUE.
                          ELSE
                              GO TO 910
                          ENDIF

C----                 Otherwise, it must be a continuation of the
C                     number
                      ELSE IF (GOTNUM) THEN
                          NPOS = NPOS + 1
                          IF (NPOS.LE.10) CNUMB(NPOS:NPOS) = CH
                      ELSE
                          GO TO 908
                      ENDIF
                  ENDIF

C----             If have a value to update, then do so
                  IF (UPDATE) THEN
                     ISTATE = 0
                     UPDATE = .FALSE.

C----                Update number of zeros or uncompressed number
                     IF (NZERO.GT.0 .OR. LVALUE.GT.0) THEN
                         IF (NZERO.GT.0) THEN
                             LVALUE = 0
                             NVALS = NZERO
                         ELSE
                             NVALS = 1
                         ENDIF

C----                    Perform the update
                         DO 500, ILOOP = 1, NVALS

C----                        Calculate subscripts for current cell
CHECK v.3.5.2-->
C                             ICELL = ICELL + 1
C                             IF (ICELL.GT.NCELL) THEN
C                                 ICELL = 1
C                                 JCELL = JCELL + 1
C                                 IF (JCELL.GT.NCELL) GO TO 912
C                             ENDIF
C                             NOBSER(ICELL,JCELL,IAMINO) = LVALUE
C                             NOBSER(ICELL,JCELL,NAMINO + 1)
C     -                           = NOBSER(ICELL,JCELL,NAMINO + 1)
C     -                           + LVALUE
                             TCELL = TCELL + 1
                             JCELL = TCELL + OFFSET
                             IF (JCELL.GT.TOTCEL) GO TO 913
                             ICELL = NCELL * (IAMINO - 1) + TCELL
                             IF (ICELL.GT.TOTCEL) GO TO 912

C----                        Store this cell's count of observations
                             NOBSER(ICELL) = LVALUE
                             NOBSER(JCELL) = NOBSER(JCELL) + LVALUE
CHECK v.3.5.2<--
                             NCOUNT(IAMINO) = NCOUNT(IAMINO) + LVALUE
                             NCOUNT(NAMINO + 1) = NCOUNT(NAMINO + 1)
     -                           + LVALUE
 500                     CONTINUE

C----                    Reinitialise
                         LVALUE = 0
                         NZERO = 0
                     ENDIF

C----                If have a letter, then update its value
                     IF (GOTLET) THEN

C----                    Calculate subscripts for current cell
CHECK v.3.5.2-->
C                         ICELL = ICELL + 1
C                         IF (ICELL.GT.NCELL) THEN
C                             ICELL = 1
C                             JCELL = JCELL + 1
C                             IF (JCELL.GT.NCELL) GO TO 912
C                         ENDIF
C                         NOBSER(ICELL,JCELL,IAMINO) = IVALUE
C                         NOBSER(ICELL,JCELL,NAMINO + 1)
C     -                       = NOBSER(ICELL,JCELL,NAMINO + 1) + IVALUE
                         TCELL = TCELL + 1
                         JCELL = TCELL + OFFSET
                         IF (JCELL.GT.TOTCEL) GO TO 913
                         ICELL = NCELL * (IAMINO - 1) + TCELL
                         IF (ICELL.GT.TOTCEL) GO TO 912

C----                    Store this cell's count of observations
                         NOBSER(ICELL) = IVALUE
                         NOBSER(JCELL) = NOBSER(JCELL) + IVALUE
CHECK v.3.5.2<--
                         NCOUNT(IAMINO) = NCOUNT(IAMINO) + IVALUE
                         NCOUNT(NAMINO + 1) = NCOUNT(NAMINO + 1)
     -                       + IVALUE
                      ENDIF
                  ENDIF

C----             If asterisk, then have reached the end
                  IF (CH.EQ.'*') THEN
                      ISTATE = 10
                      INORM = .TRUE.
CHECK v.3.5.2-->
C                      IF (ICELL.NE.NCELL .AND. JCELL.NE.NCELL) GO TO 914
                      IF (TCELL.NE.NCELL) GO TO 914
CHECK v.3.5.2<--
                  ELSE IF (CH.EQ.'#') THEN
                      ISTATE = 2
                      CNUMB = ' '
                      NPOS = 0
                  ENDIF
 600          CONTINUE
CHECK v.3.5.2-->
C          ENDIF

C----     If last record was the normalization factor, then next is a new
C         distribution
          ELSE
              INEW = .TRUE.
              INORM = .FALSE.
          ENDIF
CHECK v.3.5.2<--

      GO TO 400

C---- End of file reached
 800  CONTINUE
CHECK v.3.5.2-->
C      IF (.NOT.INEW) GO TO 902
      IF (.NOT.INEW .AND. .NOT.NONORM) GO TO 902
CHECK v.3.5.2<--
      IF (NCOUNT(NAMINO + 1).EQ.0) GO TO 918

CTESTING
C      IF (DISTRB.EQ.1) THEN
C          DO 1300, IAMINO = 1, NAMINO + 1
C              PRINT*
C              PRINT*
C              PRINT 1001, IAMINO, AMNAME(IAMINO), NCOUNT(IAMINO)
C 1001         FORMAT(I2,'. ',A3,'   Count:',I6)
C              DO 1200, ICELL = 1, NCELL
C                  IOBSER(ICELL) = NOBSER(ICELL,IAMINO)
C 1200         CONTINUE
C              PRINT 1005, (IOBSER(ICELL), ICELL = 1, NCELL)
C 1005         FORMAT(45I5)
C 1300     CONTINUE
C      ENDIF
CTESTING

      GO TO 999


C---- Errors reading parameter file
 900  CONTINUE
      PRINT*, '*** ERROR. Program error. Distribution too large for ',
     -    'array:', NCELL
      GO TO 990

 902  CONTINUE
      PRINT*, '*** ERROR. Premature end of parameters file,',
     -        ' prodata, encountered at line', LINE
      GO TO 990

 904  CONTINUE
      PRINT*, '*** Error reading parameter file, prodata, at',
     -     ' line', LINE
      GO TO 990

 906  CONTINUE
      PRINT*, '*** ERROR. Invalid amino acid code in parameter file,',
     -    ' prodata, at line', LINE, '       :', IAMINO
      GO TO 990

 908  CONTINUE
      PRINT*, '*** ERROR. Invalid character in parameter file,',
     -    ' prodata, at line', LINE, '    Pstn:', IPOS,
     -    '   [', IREC(IPOS:IPOS), ']'
      GO TO 990

 910  CONTINUE
      PRINT*, '*** ERROR. Invalid repetition number in parameter',
     -    ' file, prodata, at line', LINE, '    Pstn:', IPOS,
     -    '   [', IREC(IPOS:IPOS), ']', CNUMB, NPOS
      GO TO 990

 912  CONTINUE
      PRINT*, '*** ERROR. Cell-size too large in parameter',
     -    ' file, prodata, at line', LINE, '    Pstn:', IPOS,
CHECK v.3.5.2-->
C     -    '   [', IREC(IPOS:IPOS), ']', ICELL, NCELL
     -    '   [', IREC(IPOS:IPOS), ']', ICELL, TOTCEL
CHECK v.3.5.2<--
      PRINT*, IREC
      GO TO 990

CHECK v.3.5.2-->
 913  CONTINUE
      PRINT*, '*** ERROR. Totals cell off end of array in ',
     -    ' file, prodata, at line', LINE, '    Pstn:', IPOS,
     -    '   [', IREC(IPOS:IPOS), ']', TCELL, OFFSET, TOTCEL
      PRINT*, IREC
      GO TO 990
CHECK v.3.5.2<--

 914  CONTINUE
      PRINT*, '*** ERROR. Array NOBSER not completed from parameter',
     -    ' file, prodata. Size:', ICELL
      GO TO 990

 916  CONTINUE
      PRINT*, '*** ERROR. Missing normalization factor in parameter',
     -    ' file, prodata, at line', LINE
      GO TO 990

 918  CONTINUE
      PRINT*, '*** ERROR. No distribution data in file prodata for ',
     -    'distribution: ', DISTRB
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE SMEAR  -  Smear out the observations over the 2D distributions
C                       to make them less clumpy
C
C----------------------------------------------------------------------+--- 

CHECK v.3.5.2-->
C      SUBROUTINE SMEAR(NOBSER,ENERGY,NCELL,NCOUNT,NAMINO)
      SUBROUTINE SMEAR(NOBSER,ENERGY,NCELL1,NCELL2,NCOUNT,NAMINO)
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       NCELL, NAMINO
      INTEGER       NCELL1, NCELL2, NAMINO
CHECK v.3.5.2<--

      INTEGER       I, IAMINO, ICELL, ICELLO, J, JCELL, JCELLO,
     -              NCOUNT(NAMINO + 1)
CHECK v.3.5.2-->
C      REAL          ENERGY(NCELL,NCELL,NAMINO + 1), FACTOR, NOBS,
C     -              NOBSER(NCELL,NCELL,NAMINO + 1)
      REAL          ENERGY(NCELL1,NCELL2,NAMINO + 1), FACTOR, NOBS,
     -              NOBSER(NCELL1,NCELL2,NAMINO + 1)
CHECK v.3.5.2<--

C---- Loop over all the residue types for this distribution
      DO 1000, IAMINO = 1, NAMINO

C----     Loop through all the grid points to initialise the energies
C         array (used here as a temporary array for storing the
C         smoothed observations)
CHECK v.3.5.2-->
C          DO 200, JCELL = 1, NCELL
C              DO 100, ICELL = 1, NCELL
          DO 200, JCELL = 1, NCELL2
              DO 100, ICELL = 1, NCELL1
CHECK v.3.5.2<--
                  ENERGY(ICELL,JCELL,IAMINO) = 0.0
 100          CONTINUE
 200      CONTINUE

C----     Process only if have data points for this residue's distribution
          IF (NCOUNT(IAMINO).GT.0) THEN

C----         Loop through all the cells for the current residue
CHECK v.3.5.2-->
C              DO 600, JCELL = 1, NCELL
C                  DO 500, ICELL = 1, NCELL
              DO 600, JCELL = 1, NCELL2
                  DO 500, ICELL = 1, NCELL1
CHECK v.3.5.2<--

C----                 Extract the number of observations
                      NOBS = NOBSER(ICELL,JCELL,IAMINO)

C----                 Loop over the 8 cells surrounding the current one
                      DO 400, J = -1, 1

C----                     Convert relative position to actual position
                          JCELLO = JCELL + J
CHECK v.3.5.2-->
C                          IF (JCELLO.LT.1) JCELLO = JCELLO + NCELL
C                          IF (JCELLO.GT.NCELL) JCELLO = JCELLO - NCELL
                          IF (JCELLO.LT.1) JCELLO = JCELLO + NCELL2
                          IF (JCELLO.GT.NCELL2) JCELLO = JCELLO - NCELL2
CHECK v.3.5.2<--

                          DO 300, I = -1, 1

C----                         Convert relative position to actual position
                              ICELLO = ICELL + I
CHECK v.3.5.2-->
C                              IF (ICELLO.LT.1) ICELLO = ICELLO + NCELL
C                              IF (ICELLO.GT.NCELL)
C     -                            ICELLO = ICELLO - NCELL
                              IF (ICELLO.LT.1) ICELLO = ICELLO + NCELL1
                              IF (ICELLO.GT.NCELL1)
     -                            ICELLO = ICELLO - NCELL1
CHECK v.3.5.2<--

C----                         Compute the smearing factor to be applied
                              FACTOR = REAL(2 - ABS(I))
     -                            * REAL(2 - ABS(J)) / 16.0

C----                         Apply the smear
                              ENERGY(ICELLO,JCELLO,IAMINO)
     -                            = ENERGY(ICELLO,JCELLO,IAMINO)
     -                            + FACTOR * NOBS
 300                      CONTINUE
 400                  CONTINUE
 500              CONTINUE
 600          CONTINUE

C----         Loop over all the cells again to transfer the smoothed
C             values back across to the original array
CHECK v.3.5.2-->
C              DO 800, JCELL = 1, NCELL
C                  DO 700, ICELL = 1, NCELL
              DO 800, JCELL = 1, NCELL2
                  DO 700, ICELL = 1, NCELL1
CHECK v.3.5.2<--
                      NOBSER(ICELL,JCELL,IAMINO)
     -                    = ENERGY(ICELL,JCELL,IAMINO)
 700              CONTINUE
 800          CONTINUE
          ENDIF
 1000 CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE CALCLO  -  Convert the raw data-counts into log-odds scores
C
C----------------------------------------------------------------------+--- 

CHECK v.3.5.2-->
C      SUBROUTINE CALCLO(NOBSER,ENERGY,NCELL,NAMINO,NCOUNT)
      SUBROUTINE CALCLO(NOBSER,ENERGY,NCELL1,NCELL2,NAMINO,NCOUNT)
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       NCELL, NAMINO
      INTEGER       NCELL1, NCELL2, NAMINO
CHECK v.3.5.2<--

      INTEGER       IAMINO, ICELL, JCELL, NCOUNT(NAMINO + 1)
CHECK v.3.5.2-->
C      REAL          ENERGY(NCELL,NCELL,NAMINO + 1), GRDMAX, GRDMIN,
C     -              NOBSER(NCELL,NCELL,NAMINO + 1), PROBAB
      REAL          ENERGY(NCELL1,NCELL2,NAMINO + 1), GRDMAX, GRDMIN,
     -              NOBSER(NCELL1,NCELL2,NAMINO + 1), PROBAB
CHECK v.3.5.2<--

C---- Loop over all the residue types for this distribution
      DO 600, IAMINO = 1, NAMINO

C----     Process only if have data points for this residue's distribution
          IF (NCOUNT(IAMINO).GT.0) THEN

C----         Initialise the minimum score for this residue-type
              GRDMAX = -999.9
              GRDMIN = 0.0

C----         Loop through all the cells for the current residue
CHECK v.3.5.2-->
C              DO 200, JCELL = 1, NCELL
C                  DO 100, ICELL = 1, NCELL
              DO 200, JCELL = 1, NCELL2
                  DO 100, ICELL = 1, NCELL1
CHECK v.3.5.2<--

C----                 If no data points at all for this cell, set
C                     energy to a high level
                      IF (NOBSER(ICELL,JCELL,IAMINO).EQ.0.0) THEN
                          ENERGY(ICELL,JCELL,IAMINO) = 999.9

C----                 Calculate log-odds score for this cell
                      ELSE
                          PROBAB = NOBSER(ICELL,JCELL,IAMINO)
     -                        / REAL(NCOUNT(IAMINO))
                          ENERGY(ICELL,JCELL,IAMINO)
     -                        = LOG(PROBAB)
                          GRDMAX
     -                        = MAX(ENERGY(ICELL,JCELL,IAMINO),GRDMAX)
                          GRDMIN
     -                        = MIN(ENERGY(ICELL,JCELL,IAMINO),GRDMIN)
                      ENDIF
 100              CONTINUE
 200          CONTINUE

C----         Loop through all the grid points setting all invalid data
C             points to the maximum energy value
CHECK v.3.5.2-->
C              DO 400, JCELL = 1, NCELL
C                  DO 300, ICELL = 1, NCELL
              DO 400, JCELL = 1, NCELL2
                  DO 300, ICELL = 1, NCELL1
CHECK v.3.5.2<--
                      IF (ENERGY(ICELL,JCELL,IAMINO).GT.999.0) THEN
                          ENERGY(ICELL,JCELL,IAMINO) = GRDMIN
                      ENDIF
 300              CONTINUE
 400          CONTINUE
          ENDIF
 600  CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION CALC1D  -  Calculate log-odds score for the properties of
C                      this residue for a 1D distribution
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      REAL FUNCTION CALC1D(IAMINO,XANG,NAMINO,NCELL,NCELL1,ENERGY,STEP,
C     -    VALBEG)
      REAL FUNCTION CALC1D(IAMINO,XANG,NAMINO,NCELL,ENERGY,STEP,VALBEG)
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       NAMINO, NCELL, NCELL1
      INTEGER       NAMINO, NCELL
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       CELL1, IAMINO, ICELL, JCELL
C
C      REAL          ENERGY(NCELL,NCELL,NAMINO + 1), STEP(2), VALBEG(2),
C     -              XANG
      INTEGER       CELL1, IAMINO

      REAL          ENERGY(NCELL,NAMINO + 1), STEP(2), VALBEG(2), XANG
CHECK v.3.5.2<--

C---- Determine in which cell the current pair of values falls and
C     calculate the log-odds score
      CELL1 = (XANG - VALBEG(1)) / STEP(1) + 1

C---- If value falls outside plot-region, then return zero
CHECK v.3.5.2-->
C      IF (CELL1.LT.1 .OR. CELL1.GT.NCELL1) THEN
      IF (CELL1.LT.1 .OR. CELL1.GT.NCELL) THEN
CHECK v.3.5.2<--
          CALC1D = 999.99

C---- Otherwise, get the corresponding energy
      ELSE

C----     Calculate cell's physical position in the array
CHECK v.3.5.2-->
C          JCELL = 1 + (CELL1 - 1) / NCELL
C          ICELL = CELL1 - (JCELL - 1) * NCELL
C
CC----     Return the energy value
C          CALC1D = ENERGY(ICELL,JCELL,IAMINO)

C----     Return the energy value
          CALC1D = ENERGY(CELL1,IAMINO)
CHECK v.3.5.2<--
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION CALC2D  -  Calculate Boltzmann energy for the properties of
C                      this residue for a 2D distribution
C
C----------------------------------------------------------------------+---

CHECK v.3.5.2-->
C      REAL FUNCTION CALC2D(IAMINO,XANG,YANG,NAMINO,NCELL,ENERGY,STEP,
C     -    VALBEG)
      REAL FUNCTION CALC2D(IAMINO,XANG,YANG,NAMINO,NCELL1,NCELL2,ENERGY,
     -    STEP,VALBEG)
CHECK v.3.5.2<--

CHECK v.3.5.2-->
C      INTEGER       NAMINO, NCELL
      INTEGER       NAMINO, NCELL1, NCELL2
CHECK v.3.5.2<--

      INTEGER       IAMINO, ICELL, JCELL
CHECK v.3.5.2-->
C      REAL          ENERGY(NCELL,NCELL,NAMINO + 1), STEP(2), VALBEG(2),
C     -              XANG, YANG
      REAL          ENERGY(NCELL1,NCELL2,NAMINO + 1), STEP(2),
     -              VALBEG(2), XANG, YANG
CHECK v.3.5.2<--

C---- Determine in which cell the current pair of values falls and
C     calculate the log-odds score
      ICELL = (XANG - VALBEG(1)) / STEP(1) + 1
      JCELL = (YANG - VALBEG(2)) / STEP(2) + 1

C---- If value falls outside plot-region, then return zero
CHECK v.3.5.2-->
C      IF (ICELL.LT.1 .OR. ICELL.GT.NCELL .OR. JCELL.LT.1 .OR.
C     -    JCELL.GT.NCELL) THEN
      IF (ICELL.LT.1 .OR. ICELL.GT.NCELL1 .OR. JCELL.LT.1 .OR.
     -    JCELL.GT.NCELL2) THEN
CHECK v.3.5.2<--
          CALC2D = 999.99

C---- Otherwise, retrieve the energy value
      ELSE

C----     Return the energy value
          CALC2D = ENERGY(ICELL,JCELL,IAMINO)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.4.3<--
