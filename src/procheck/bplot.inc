C----------------------------------------------------------------------+---
C BPLOT.INC - Include variables for BPLOT.FOR
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INTEGER       MCHANG, MCHLEN, MXBARS, MXCOLR, NACROS, NALP,
     -              NDOWN, NGAPA, NGAPL, NGAPP, NMAIN, NPLANE, NPLATM
      PARAMETER    (
     -              MCHANG =   20,
     -              MCHLEN =   11,
     -              MXBARS =  120,
CHECK v.3.4-->
C     -		    MXCOLR =   20,
     -		    MXCOLR =   50,
CHECK v.3.4<--
     -              NACROS =    3,
     -              NDOWN  =    4,
     -              NGAPA  =    6,
     -              NGAPL  =    6,
     -              NGAPP  =    3,
     -              NMAIN  =  MCHANG + MCHLEN,
     -              NPLANE =    9,
     -              NALP   =  MCHANG + MCHLEN + NPLANE,
     -              NPLATM =   10
     -             )

      REAL          ASEP, AWIDTH, FPOSA, FPOSL, LSEP, LWIDTH, PWIDTH,
     -              XDAYLT, XFQLMG, XFQRMG, XFQWID, YDAYLT, YFQBMG,
     -              YFQTMG, YFQWID
      PARAMETER    (
     -              ASEP   =    5.0,
     -              AWIDTH =    0.5,
     -              FPOSA  =    1.0,
     -              FPOSL  =    0.01,
     -              LSEP   =    0.04,
     -              LWIDTH =    0.005,
     -              PWIDTH =    0.002,
     -              XDAYLT =   40.0,
     -              XFQLMG =   50.0,
     -              XFQRMG =   20.0,
     -              XFQWID =  120.0,
     -              YDAYLT =   40.0,
     -              YFQBMG =  100.0,
     -              YFQTMG =  100.0,
     -              YFQWID =  150.0
     -             )

      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2
      PARAMETER    (
     -              BBOXX1 =  30,
     -              BBOXX2 = 550,
     -              BBOXY1 =  50,
     -              BBOXY2 = 780
     -             )

      INTEGER       NDGACA, NDGACL, NDGACP
      PARAMETER    (
     -              NDGACA = 6,
     -              NDGACL = 6,
     -              NDGACP = 6
     -             )

      REAL          DGASPA, DGASPL, DGASPP, DGLENA, DGLENL, DGLENP,
     -              DGHDOF, DGMRGX, DGMRGY, DGHDSZ, DGNOTE, DGSTDL,
     -              DGTEX, DGWIDX, DGWIDY
      PARAMETER    (
     -              DGASPA =   1.15,
     -              DGASPL =   1.5,
     -              DGASPP =   0.8,
     -              DGLENA =   4.0,
     -              DGLENL =   3.0,
     -              DGLENP =   8.0,
     -              DGHDOF =   0.3,
     -              DGMRGX =  60.0,
     -              DGMRGY =  50.0,
     -              DGHDSZ =  15.0,
     -              DGNOTE =   0.02,
     -              DGSTDL =   1.4,
     -              DGTEX  =   0.1,
     -              DGWIDX = 460.0,
     -              DGWIDY = 650.0
     -             )

      CHARACTER*1   AL(NALP), CHAIN
      CHARACTER*3   PNAME(NPLANE)
      CHARACTER*8   LENANG(NMAIN)
      CHARACTER*12  COLNAM(MXCOLR)
      CHARACTER*13  ENCODE(NMAIN)
      CHARACTER*16  ENDESC(NMAIN)
      CHARACTER*30  HEADIN
CHECK v.3.4-->
CHECK v.3.2-->
CC      CHARACTER*80  BRCODE, FILPS, FILLAN, FILPLN, PDBFIL
C      CHARACTER*80  BRCODE, FILPS, FILLAN, FILPLN, FILSUM, PDBFIL
CCHECK v.3.2<--
      CHARACTER*80  BRCODE, FILHTM, FILLAN, FILPLN, FILPS, FILSUM,
     -              PDBFIL
CHECK v.3.4<--
      COMMON/CHARS/
     -              AL, CHAIN,

     -              PNAME,

     -              LENANG,

     -		    COLNAM,

     -              ENCODE,

     -              ENDESC,

     -              HEADIN,

CHECK v.3.4-->
CCHECK v.3.2-->
CC     -              BRCODE, FILPS, FILLAN, FILPLN, PDBFIL
C     -              BRCODE, FILPS, FILLAN, FILPLN, FILSUM, PDBFIL
CCHECK v.3.2<--
     -              BRCODE, FILHTM, FILLAN, FILPLN, FILPS, FILSUM,
     -              PDBFIL
CHECK v.3.4<--


      INTEGER       BCOUNT(MXBARS+2,NALP), BLEN, COLPLT(4,4), HLEN,
     -              IEND, ILEN, INTRLX, IPLOT, ISTART, MAXFRQ(NALP),
CHECK v.3.4-->
C     -              MAXLEN, MCHNUM(NALP), MINLEN, NBARS(NALP), NFGAPY,
C     -              NPATOM(NPLANE), NRES, NTICKS(NALP), PLSTAT, PLTYPE,
C     -		    PSLEN, TODO(3)
     -              MAXLEN, MCHNUM(NALP), MINLEN, NBARS(NALP), NCOLOR,
     -              NFGAPY, NPATOM(NPLANE), NRES, NTICKS(NALP), PLSTAT,
     -		    PLTYPE, PSLEN, TODO(3)
CHECK v.3.4<--
      COMMON/INTS/
     -              BCOUNT, BLEN, COLPLT, HLEN,
     -              IEND, ILEN, INTRLX, IPLOT, ISTART, MAXFRQ,
CHECK v.3.4-->
C     -              MAXLEN, MCHNUM, MINLEN, NBARS, NFGAPY,
C     -              NPATOM, NRES, NTICKS, PLSTAT, PLTYPE,
C     -		    PSLEN, TODO
     -              MAXLEN, MCHNUM, MINLEN, NBARS, NCOLOR,
     -              NFGAPY, NPATOM, NRES, NTICKS, PLSTAT,
     -		    PLTYPE, PSLEN, TODO
CHECK v.3.4<--

CHECK v.3.2-->
C      LOGICAL       BSHADE(3), GOTOFF, IFAIL, INCOLR(4), PLOTBA, PLOTBL,
C     -              PLOTDG, PLOTPL, TOPLOT(3)
C      COMMON/LOGS/
C     -              BSHADE, GOTOFF, IFAIL, INCOLR, PLOTBA, PLOTBL,
C     -              PLOTDG, PLOTPL, TOPLOT
      LOGICAL       BSHADE(3), COMBPS, GOTOFF, IFAIL, INCOLR(4), PLABEL,
CHECK v.3.4.3-->
C     -              PLOTBA, PLOTBL, PLOTDG, PLOTPL, TOPLOT(3), WITHAN
     -              PLOTBA, PLOTBL, PLOTDG, PLOTPL, TOPLOT(3), WITHAN,
     -              WWWOUT
CHECK v.3.4.3<--
      COMMON/LOGS/
     -              BSHADE, COMBPS, GOTOFF, IFAIL, INCOLR, PLABEL,
CHECK v.3.4.3-->
C     -              PLOTBA, PLOTBL, PLOTDG, PLOTPL, TOPLOT, WITHAN
     -              PLOTBA, PLOTBL, PLOTDG, PLOTPL, TOPLOT, WITHAN,
     -              WWWOUT
CHECK v.3.4.3<--
CHECK v.3.2<--

      REAL          DGLIMA, DGLIML, DGLIMP, DGLIMR, ENGMEA(NMAIN),
     -              ENGSTD(NMAIN), GAPX, LIMHST, LIMPLN(2), LIMRES(2),
     -              MAXX(NALP), MCHMAX(NALP), MCHMIN(NALP),
     -              MCHMEA(NALP), MCHSTD(NALP), MINX(NALP), RADDEG,
     -		    RESMIN, RGB(3,MXCOLR), STDLIM
      COMMON/REALS/
     -              DGLIMA, DGLIML, DGLIMP, DGLIMR, ENGMEA,
     -              ENGSTD, GAPX, LIMHST, LIMPLN, LIMRES,
     -              MAXX, MCHMAX, MCHMIN,
     -              MCHMEA, MCHSTD, MINX, RADDEG,
     -		    RESMIN, RGB, STDLIM

CHECK v.3.2-->
      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/
CHECK v.3.2<--

C----------------------------------------------------------------------+---
C Parameters
C ----------
C
C ASEP   - Distance between numbers on x-axis for bond angles
C AWIDTH - Width of histogram bars on angle plots, in degrees
C AL(NALP) - 1-character marker indicating whether this is a bond-length
C          or a bond angle
C BCOUNT(MXBARS+2,NALP) - Count of lengths/angles in each histogram bar
C          on each plot
C BLEN   - Length of subheading on all plots
C BSHADE(3) - Flags indicating whether background shading required on
C         the plots
C CHAIN  - Chain identifier for current protein structure. If blank, then
C	   all chains are to be included
C COLNAM(MXCOLR) - Names identifying the user-defined colours
C COLPLT(3,13) - Colour nunbers for the various colours used in the 3 plot
C          types
C COMBPS - Flag indicating whether all pages of a given plot are to be
C          combined into a single PostScript file
C DGASPA, DGASPL, DGASPP - Aspect ratios of boxes containing distorted
C          geometries of angles, lengths and planes
C DGLENA, DGLENL, DGLENP - Equivalent "length" of each distorted geometry
C          box
C DGLIMA, DGLIML, DGLIMP, DGLIMR - Cut-off values for display of distorted
C          geometry
C DGHDOF - Relative offset of distorted-geometry heading in text box
C DGMRGX, DGMRGY - x- and y-margins for plot area for distorted geometry
C DGHDSZ - Size of distorted-geometry headings
C DGNOTE - Size of region for nore at the end of the distorted geometry
C          plots
C DGSTDL - "Standard" bond length used when drawing distorted bond angles
C DGTEX  - Fraction of distorted-geometry plot area occupied by heading
C DGWIDX, DGWIDY - x- and y-sizes of plot area for distorted geometry
C ENCODE(NMAIN) - Engh & Huber coding for the bond length/angle
C ENDESC(NMAIN) - Continuation of Engh & Huber coding
C ENGMEA(NMAIN) - Mean value from Engh & Huber data for this length/angle
C ENGSTD(NMAIN) - St dev from Engh & Huber data for this length/angle
C FILHTM - Name of output html file holding summary information in html
C          format
C FILLAN - Filename containing input data created by program ANGLEN. Name is
C          of the form <filename>.lan
C FILPLN - Filename containing transformed coords of all planar groups,
C          created by program ANGLEN. Name is of the form <filename>.pln
C FILPS  - Name of output PostScript files
C FILSUM - Name of output .sum summary file
C GAPX   - Gap between ticks on the x-axes of the quality plots
C GOTOFF - Flag indicating whether any plots on this page have data
C          off the ends of the graph
C HEADIN - Plot heading
C HLEN   - Length of plot heading (in characters)
C IEND   - Last character of filename holding original structure
C IFAIL  - Failure flag causing program to abort due to some error condition
C ILEN   - Length of filename holding original structure
C INCOLR(3) - Flags indicating whether colour PostScript files are required
C          for the different plot types
C INTRLX - Integer/Real for numbers on the x-axis (0/1)
C IPLOT  - Number of current PostScript output file
C ISTART - First character of filename holding original structure
C LENANG(NMAIN) - Character-string showing atoms making up the bond-length
C          or bond angle
C LIMPLN(2) - RMS distances for highlighting for ring and other groups,
C          respectively
C LIMRES(2) - No. of standard deviations for highlighting bars on plot
C LSEP   - Distance between numbers on x-axis for bond lengths
C LWIDTH - Width of histogram bars on bond length plots, in Angstroms
C MCHANG - Number of main-chain bond angle types
C MAXFRQ(NALP) - Maximum frequency count for each bond length/angle
C MCHLEN - Number of main-chain bond length types
C MCHMAX(NALP) - Maximum value for each type of main-chain bond length and
C          angle
C MCHMIN(NALP) - Minimum value for each type of main-chain bond length and
C          angle
C MCHMEA(NALP) - Mean value for each type of main-chain bond length and
C          angle
C MCHNUM(NALP) - Number of values for each  type of main-chain bond
C          length and angle
C MCHSTD(NALP) - Standard deviation for each type of main-chain bond
C          length and angle
C MINLEN, MAXLEN - Minimum and maximum sequence lengths found in the
C          .res files
C MAXX(NALP), MINX(NALP) - Maximum and minimum values for the x-axes of
C          each bond-length/angle plot
C MXBARS - Maximum number of histogram bars that can be stored
C MXCOLR - Maximum number of RGB colours that can be defined by the user
C NACROS - Number of graphs across the page for frequency distrib plot
C NALP   - Total number of bond length and angles, and planar groups
C NBARS(NALP) - Number of histogram bars on each graph
C NCOLOR - Number of user-defined colours in the parameter file
C NDGACA, NDGACL, NDGACP - Number of boxes across the page giving the
C          distorted geometry of angles, lengths and planes
C NDOWN  - Number of graphs down the page for frequency distrib plot
C NFGAPY - Number of gaps along y-axis for frequency distrib plot
C NMAIN  - Number of main-chain bond length and bond angle types
C NPLANE - Number of planar groups defined
C NPLATM - Maximum number of atoms in planar groups
C NRES   - Number of residues in the protein
C NGAPA, NGAPL, NGAPP - Number of standard deviations either side of mean
C          defining size of graph in x-direction
C NTICKS(NALP) - Number of ticks required for each graph along x-axis
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C PLABEL - Flag indicating whether plot filename to be printed on the plot
C          itself
C PLOTBA, PLOTBL, PLOTDG, PLOTPL - Flags indicating whether bond angle,
C          bond length, distorted geometry, and/or planar groups are to
C          be plotted
C PLTYPE - Number indicating which of the plot-types is being processed
C PNAME(NPLANE) - Name of Planar group
C PSLEN  - Length of PostScript filename
C PWIDTH - Width of histogram bars on planar groups plots
C RADDEG - Conversion factor between radians and degrees
C RESMIN - Minimum resolution - used for the comparison statistics
C RGB(3,MXCOLR) - User-defined RGB colour definitions
C TODO(3) - Count of each type of distorted geometry to be plotted
C TOPLOT(3) - Flag indicating whether each of the 3 histogram plot-types
C          has any data to be plotted
C WITHAN - Flag indicating whether file-handles are required for
C          PostScript output filenames
C WWWOUT - Flag indicating whether output .html file required
C XDAYLT - "Daylight" gap between plots in the x-direction
C XFQLMG, XFQRMG - Left and right margins for freq distribs plot
C XFQWID, YFQWID - Width and height of each litte graph in the freq
C          distribs plot
C YDAYLT - "Daylight" gap between plots in the y-direction
C YFQBMG, YFQTMG - Bottom and top margins for freq distribs plot
C YBAND(10)  - Band widths for shaded curves on main-chain and side-chain
C          plots
C---------------------------------------------------------------------------
