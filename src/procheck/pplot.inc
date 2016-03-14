C----------------------------------------------------------------------+--- 
C PPLOT.INC - Include variables for PPLOT.FOR
C----------------------------------------------------------------------+--- 

CVAX      IMPLICIT NONE

CHECK v.3.5.1-->
C      INTEGER       MCHANG, MCHLEN, MXCOLR, MXCOLS, MXPINS, NALP,
C     -              NCOVAL, NDIHED, NDISTR, NMAIN, NPLANE, NQUALP
      INTEGER       MCHANG, MCHLEN, MXCOLR, MXCOLS, MXPINS, MXRES, NALP,
     -              NCOVAL, NDIHED, NDISTR, NMAIN, NPLANE, NQUALP
CHECK v.3.5.1<--
      PARAMETER    (
     -              MCHANG =   20,
     -              MCHLEN =   11,
CHECK v.3.4-->
C     -              MXCOLR =  20,
     -              MXCOLR =   50,
CHECK v.3.4<--
     -              MXCOLS =   15,
     -              MXPINS =  100,
CHECK v.3.5.1-->
     -              MXRES  = 8000,
CHECK v.3.5.1<--
     -              NCOVAL =    2,
     -              NDIHED =    5,
     -              NDISTR =  NCOVAL + NDIHED + 3,
     -              NMAIN  =  MCHANG + MCHLEN,
     -              NPLANE =    9,
     -              NQUALP =   11,
     -              NALP   =  MCHANG + MCHLEN + NPLANE
     -             )

CHECK v.3.1-->
C      INTEGER       MXPLOT, PLTCH1, PLTOME, PLTCAL, PLTHB, PLTCGB,
C     -              PLTLAV, PLTLCC, PLTLPP, PLTMCB, PLTOOI, PLTSCB
C      PARAMETER    (
C     -              MXPLOT = 11,
      INTEGER       MXPLOT, PLTACC, PLTCH1, PLTOME, PLTCAL, PLTHB,
     -              PLTCGB, PLTLAV, PLTLCC, PLTLPP, PLTMCB, PLTOOI,
CHECK v.3.2-->
C     -              PLTSCB
C      PARAMETER    (
C     -              MXPLOT = 12,
     -              PLTSCB, PLTSDM, PLTSDS
      PARAMETER    (
     -              MXPLOT = 14,
CHECK v.3.2<--
CHECK v.3.1<--
     -              PLTCH1 = 1,
     -              PLTOME = 2,
     -              PLTCAL = 3,
     -              PLTHB  = 4,
     -              PLTCGB = 5,
     -              PLTMCB = 6,
     -              PLTSCB = 7,
     -              PLTLPP = 8,
     -              PLTLCC = 9,
     -              PLTLAV = 10,
CHECK v.3.1-->
C     -              PLTOOI = 11
     -              PLTOOI = 11,
CHECK v.3.2-->
C     -              PLTACC = 12
     -              PLTACC = 12,
     -              PLTSDM = 13,
     -              PLTSDS = 14
CHECK v.3.2<--
CHECK v.3.1<--
     -             )

      REAL          GLOBC2, GLOBDH, GLOBGM, GLOBGP, GLOBGT, GLOBHB,
     -              GLOBHP, GLOBOM, GLOBPH, GLOBPP, GLOBSS, GLOB3L,
     -              GLOB3R, SDEVCS, SDEVC2, SDEVDH, SDEVGM, SDEVGP,
     -              SDEVGT, SDEVHB, SDEVHP, SDEVOM, SDEVPH, SDEVPO,
     -              SDEVPP, SDEVSS, SDEV3R, SDEV3L   

      PARAMETER    (
     -              GLOBC2 = 177.4,
     -              GLOBDH =  33.9,
     -              GLOBGM =  64.1,
     -              GLOBGP = 293.3,
     -              GLOBGT = 183.6,
     -              GLOBHB =  -2.03,
     -              GLOBHP = -39.4,
     -              GLOBOM = 180.0,
     -              GLOBPH = -65.3,
     -              GLOBPP = -65.4,
     -              GLOBSS =   2.0,
     -              GLOB3L = -85.8,
     -              GLOB3R =  96.8
     -             )
      PARAMETER    (
     -              SDEVC2 =  18.5,
     -              SDEVCS =   7.8,
     -              SDEVDH =   3.5,
     -              SDEVGM =  15.7,
     -              SDEVGP =  15.0,
     -              SDEVGT =  16.8,
     -              SDEVHB =   0.75,
     -              SDEVHP =  11.3,
     -              SDEVOM =   5.8,
     -              SDEVPH =  11.9,
     -              SDEVPO =  15.7,
     -              SDEVPP =  11.2,
     -              SDEVSS =   0.1,
     -              SDEV3R =  14.8,
     -              SDEV3L =  10.7
     -             )

      REAL          VWDIST
      PARAMETER    (
     -              VWDIST = 2.6
     -             )
 
      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2
      PARAMETER    (
     -              BBOXX1 =  30,
     -              BBOXX2 = 550,
     -              BBOXY1 =  50,
     -              BBOXY2 = 780
     -             )

      REAL          HWIDX, HWIDY, SSSIZE, XGAP, XPORIG, XPWID, YHT,
     -              YPORIG, YPSEP, YPWID
      PARAMETER    (
     -              HWIDX  =   2.0,
     -              HWIDY  =   2.0,
     -              SSSIZE =   3.0,
     -              XPORIG =  95.0,
     -              XPWID  = 425.0,
     -              XGAP   = XPWID / MXPINS,
     -              YHT    =   5.0,
     -              YPORIG = 355.0,
     -              YPSEP  =  47.0,
     -              YPWID  =  77.0
     -             )

C---- Parameters for main-chain and side-chain parameters plot
      INTEGER       NGAPX
      PARAMETER    (NGAPX = 6)

      REAL          RSTART, REND, RGAP, XORIGN, XSEP, XMWID, YORIGN,
     -              YSEP, YMWID
      PARAMETER    (
     -              REND   = 4.0,
     -              RGAP   = 0.1,
     -              RSTART = 1.0,
     -              XORIGN = 114.0,
     -              XSEP   =  56.0,
     -              XMWID  = 150.0,
     -              YORIGN = 190.0,
     -              YSEP   =  40.0,
     -              YMWID  = 146.0
     -             )

      CHARACTER*1   CHAIN
      CHARACTER*12  COLNAM(MXCOLR)
      CHARACTER*30  DESC(NQUALP)
CHECK v.3.4-->
C      CHARACTER*78  BRCODE, FILAST, FILLAN, FILNB, FILOUT, FILPS,
CCHECK v.3.2-->
CC     -              FILRIN, FILRSA, FILSCO, FILSDH, PDBFIL
C     -              FILRIN, FILRSA, FILSCO, FILSDH, FILSUM, PDBFIL
CCHECK v.3.2<--
      CHARACTER*78  BRCODE, FILAST, FILHTM, FILLAN, FILNB, FILOUT,
     -              FILPS, FILRIN, FILRSA, FILSCO, FILSDH, FILSUM,
     -              PDBFIL
CHECK v.3.4<--
      COMMON/CHARS/
     -	            CHAIN,

     -              COLNAM,

     -              DESC,

CHECK v.3.4-->
C     -              BRCODE, FILAST, FILLAN, FILNB, FILOUT, FILPS,
CCHECK v.3.2-->
CC     -              FILRIN, FILRSA, FILSCO, FILSDH, PDBFIL
C     -              FILRIN, FILRSA, FILSCO, FILSDH, FILSUM, PDBFIL
CHECK v.3.2<--
     -              BRCODE, FILAST, FILHTM, FILLAN, FILNB, FILOUT,
     -              FILPS, FILRIN, FILRSA, FILSCO, FILSDH, FILSUM,
     -              PDBFIL
CHECK v.3.4<--

      INTEGER       BLEN, COLPLT(3,13), IEND, ILEN, ILINE, IPAGE, IPLOT,
CHECK v.3.2-->
C     -              ISTART, MCHNUM(NMAIN), NALLOW, NALPHC, NAROMA, NASN,
C     -              NBAD, NBADPR, NCHICH, NCHI2, NCHI3L, NCHI3R, NCIS,
C     -              NCORE, NCPOOL, NDAM, NDISUL, NENDCH, NGENER, NGLY,
C     -		     NGMINU, NGPLUS, NGTRAN, NHBOND, NHXPSI, NLINES,
C     -              NOMEGA, NONPG, NOUTSI, NPHIHX, NPLOTS, NPRINT, NPRO,
C     -              NPROPH, NRAMPL, NRESID, NRGRAF, NVALUE(NQUALP),
C     -              OOIMAX, OOIMIN, PLTNUM(3), PLTYPE
     -              ISTART, MCHNUM(NMAIN), MCLAST(3), NALLOW, NALPHC,
     -              NAROMA, NASN, NBAD, NBADPR, NCHICH, NCHI2, NCHI3L,
CHECK v.3.4-->
C     -              NCHI3R, NCIS, NCORE, NCPOOL, NDAM, NDISUL, NENDCH,
C     -              NGENER, NGLY, NGMINU, NGPLUS, NGTRAN, NHBOND,
C     -              NHXPSI, NLINES, NOMEGA, NONPG, NOUTSI, NPHIHX,
C     -              NPLOTS, NPRINT, NPRO, NPROPH, NRAMPL, NRESID,
C     -              NRGRAF, NVALUE(NQUALP), OOIMAX, OOIMIN, PLTNUM(3),
C     -              PLTYPE
     -              NCHI3R, NCIS, NCOLOR, NCORE, NCPOOL, NDAM, NDISUL,
     -		    NENDCH, NGENER, NGLY, NGMINU, NGPLUS, NGTRAN,
     -              NHBOND, NHXPSI, NLINES, NOMEGA, NONPG, NOUTSI,
     -              NPHIHX, NPLOTS, NPRINT, NPRO, NPROPH, NRAMPL,
     -              NRESID, NRGRAF, NVALUE(NQUALP), OOIMAX, OOIMIN,
     -              PLTNUM(3), PLTYPE
CHECK v.3.4<--
CHECK v.3.2<--
      COMMON/INTS/
     -              BLEN, COLPLT, IEND, ILEN, ILINE, IPAGE, IPLOT,
CHECK v.3.2-->
C     -              ISTART, MCHNUM, NALLOW, NALPHC, NAROMA, NASN,
C     -              NBAD, NBADPR, NCHICH, NCHI2, NCHI3L, NCHI3R, NCIS,
C     -              NCORE, NCPOOL, NDAM, NDISUL, NENDCH, NGENER, NGLY,
C     -		     NGMINU, NGPLUS, NGTRAN, NHBOND, NHXPSI, NLINES,
C     -              NOMEGA, NONPG, NOUTSI, NPHIHX, NPLOTS, NPRINT, NPRO,
C     -              NPROPH, NRAMPL, NRESID, NRGRAF, NVALUE,
C     -              OOIMAX, OOIMIN, PLTNUM, PLTYPE
     -              ISTART, MCHNUM, MCLAST, NALLOW, NALPHC,
     -              NAROMA, NASN, NBAD, NBADPR, NCHICH, NCHI2, NCHI3L,
CHECK v.3.4-->
C     -              NCHI3R, NCIS, NCORE, NCPOOL, NDAM, NDISUL, NENDCH,
C     -		    NGENER, NGLY, NGMINU, NGPLUS, NGTRAN, NHBOND,
C     -              NHXPSI, NLINES, NOMEGA, NONPG, NOUTSI, NPHIHX,
C     -              NPLOTS, NPRINT, NPRO, NPROPH, NRAMPL, NRESID,
C     -              NRGRAF, NVALUE, OOIMAX, OOIMIN, PLTNUM,
C     -              PLTYPE
     -              NCHI3R, NCIS, NCOLOR, NCORE, NCPOOL, NDAM, NDISUL,
     -		    NENDCH, NGENER, NGLY, NGMINU, NGPLUS, NGTRAN,
     -              NHBOND, NHXPSI, NLINES, NOMEGA, NONPG, NOUTSI,
     -              NPHIHX, NPLOTS, NPRINT, NPRO, NPROPH, NRAMPL,
     -              NRESID, NRGRAF, NVALUE, OOIMAX, OOIMIN,
     -              PLTNUM, PLTYPE
CHECK v.3.4<--
CHECK v.3.2<--

CHECK v.3.2-->
C      LOGICAL       ACCESS, ASTONL, BCLIST, BSHADE(3), HEADIN, IDITIS,
C     -              IFAIL, INCOLR(3), NONBON, NOPLOT, ONLYHI, PLOTPN,
C     -              PLOTQ1, PLOTQ2, SHOWAC
      LOGICAL       ACCESS, ASTONL, BCLIST, BSHADE(3), COMBPS, HEADIN,
     -              IDITIS, IFAIL, INCOLR(3), NONBON, NOPLOT, ONLYHI,
CHECK v.3.4.3-->
C     -              PLABEL, PLOTPN, PLOTQ1, PLOTQ2, SHOWAC, WITHAN
     -              PLABEL, PLOTPN, PLOTQ1, PLOTQ2, SHOWAC, WITHAN,
     -              WWWOUT
CHECK v.3.4.3<--
CHECK v.3.2<--
      COMMON/LOGS/
CHECK v.3.2-->
C     -              ACCESS, ASTONL, BCLIST, BSHADE, HEADIN, IDITIS,
C     -              IFAIL, INCOLR, NONBON, NOPLOT, ONLYHI, PLOTPN,
C     -              PLOTQ1, PLOTQ2, SHOWAC
     -              ACCESS, ASTONL, BCLIST, BSHADE, COMBPS, HEADIN,
     -              IDITIS, IFAIL, INCOLR, NONBON, NOPLOT, ONLYHI,
CHECK v.3.4.3-->
     -              PLABEL, PLOTPN, PLOTQ1, PLOTQ2, SHOWAC, WITHAN,
     -              WWWOUT
CHECK v.3.4.3<--
CHECK v.3.2<--

      REAL          ACCMAX, ACDEV, ACMEAN, ACSTD, ALLOWP, BADCOP,
     -              CH2DEV, CH2MEA, CH2STD, COEFF(3,NQUALP), COREPC,
     -              C3LMEA, C3LSTD, C3RMEA, C3RSTD, GENERP,
     -              GFACTR(NDISTR), GMIDEV, GMIMEA, GMISTD, GPLDEV,
     -              GPLMEA, GPLSTD, GTRDEV, GTRMEA, GTRSTD, HBDEV,
     -              HBMEAN, HBSTD, HPMEAN, HPSTD, LIMRES, MCHMAX(NMAIN),
     -              MCHMIN(NMAIN), MCHMEA(NMAIN), MCHSTD(NMAIN),
CHECK v.3.2-->
C     -              MXDEVC(MXCOLS), OMDEV, OMMEAN, OMSTD, OUTSIP,
C     -              PHMEAN, PHSTD, PINMAX(MXPLOT), PLTSTD(MXPLOT),
C     -              POOMEA, POOSTD, PPMEAN, PPSTD, RESOL,
C     -              RGB(3,MXCOLR), SFACT, STDLIM, SSMEAN, SSSTD,
C     -              SVALUE(NQUALP), YBAND(NQUALP)
     -              MXDEVC(MXCOLS), OMAXDV(2), OMDEV, OMMEAN, OMSTD,
     -              OUTSIP, PHMEAN, PHSTD, PINMAX(MXPLOT),
     -              PLTSTD(MXPLOT), POOMEA, POOSTD, PPMEAN, PPSTD,
CHECK v.3.5.1-->
C     -              RESOL, RGB(3,MXCOLR), SFACT, STDLIM, SSMEAN,
C     -              SSSTD, SVALUE(NQUALP), YBAND(NQUALP)
     -              RESGF(MXRES), RESOL, RGB(3,MXCOLR), SFACT, STDLIM,
     -              SSMEAN, SSSTD, SVALUE(NQUALP), YBAND(NQUALP)
CHECK v.3.5.1<--
CHECK v.3.2<--

      COMMON/REALS/
     -              ACCMAX, ACDEV, ACMEAN, ACSTD, ALLOWP, BADCOP,
     -              CH2DEV, CH2MEA, CH2STD, COEFF, COREPC,
     -              C3LMEA, C3LSTD, C3RMEA, C3RSTD, GENERP,
     -              GFACTR, GMIDEV, GMIMEA, GMISTD, GPLDEV,
     -              GPLMEA, GPLSTD, GTRDEV, GTRMEA, GTRSTD, HBDEV,
     -              HBMEAN, HBSTD, HPMEAN, HPSTD, LIMRES, MCHMAX,
     -              MCHMIN, MCHMEA, MCHSTD,
CHECK v.3.2-->
C     -              MXDEVC, OMDEV, OMMEAN, OMSTD, OUTSIP,
C     -              PHMEAN, PHSTD, PINMAX, PLTSTD,
C     -              POOMEA, POOSTD, PPMEAN, PPSTD, RESOL,
C     -              RGB, SFACT, STDLIM, SSMEAN, SSSTD,
C     -              SVALUE, YBAND
     -              MXDEVC, OMAXDV, OMDEV, OMMEAN, OMSTD,
     -              OUTSIP, PHMEAN, PHSTD, PINMAX,
     -              PLTSTD, POOMEA, POOSTD, PPMEAN, PPSTD,
CHECK v.3.5.1-->
C     -              RESOL, RGB, SFACT, STDLIM, SSMEAN,
C     -              SSSTD, SVALUE, YBAND
     -              RESGF, RESOL, RGB, SFACT, STDLIM,
     -              SSMEAN, SSSTD, SVALUE, YBAND
CHECK v.3.5.1<--

      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/

CHECK v.3.2<--

C----------------------------------------------------------------------+--- 
C Parameters
C ----------
C
C ACCESS - Flag indicating whether calculated accessibilities are to be
C          picked up from the .rsa file (where present), or whether
C          estimated accessibilities (based on Ooi numbers) are to be used
C          instead
C ACCMAX - Maximum accessibility value
C ACDEV  - Deviation from ideal value of alpha-carbon chirality
C ACMEAN - Mean value for alpha-carbon chirality
C ACSTD  - Standard deviation of alpha-carbon chirality
C ALLOWP - Percentage of residues in the allowed regions of the Ramachandran
C          plot
C ASTONL - Flag signifying that only asterisk-lines are to be printed (ie
C          without the data itself)
C BADCOP - Number of bad contacts per 100 residues
C BBOXX1, BBOXX2, BBOXY1, BBOXY2 - Bounding box limits for PostScript area
C          on page. Entire plot is placed weithin this region
C BCLIST - Flag indicating whether a listing of the bad contacts is required
C BLEN   - Length of subheading on all plots
C BSHADE(3) - Flag indicating whether background shading required on plots
C CGBMAX - Maximum gamma atom B-value for pins plot
C CHAIN  - Chain identifier for current protein structure. If blank, then
C	   all chains are to be included
C CH1MAX - Maximum Chi-1 deviation for pins plot
C CH2DEV - Deviation from ideal value of the Chi-2 dihedral angle
C CH2MEA - Mean value of the Chi-2 dihedral angle
C CH2STD - Standard deviation of the Chi-2 dihedral angle
C COEFF(3,NQUALP) - Coefficients for shaded curves on main-chain and
C          side-chain plots
C COLNAM(MXCOLR) - Names identifying the user-defined colours
C COLPLT(3,13) - Colour nunbers for the various colours used in the 3 plot
C          types
C COMBPS - Flag indicating whether all pages of a given plot are to be
C          combined into a single PostScript file
C COREPC - Percentage of residues in the core regions of the Ramachandran
C          plot
C C3LMEA - Mean value of the Chi-3 right-hand values
C C3LSTD - Standard deviation of the Chi-3 right-hand values
C C3RMEA - Mean value of the Chi-3 left-hand values
C C3RSTD - Standard deviation of the Chi-3 left-hand values
C DESC(NQUALP) - Descriptions of the stereochemical parameters on the
C          Quality Plots
C DISALL - Flag indicating whether the given residue is in the disallowed
C          region of the Ramachandran plot
C FILAST - Output file giving maximum asterisk values, residue by residue
C FILHTM - Name of output html file holding summary information in html
C          format
C FILLAN - Input file, <filename>.lan containing mainchain bond lengths and
C          bond angles
C FILOUT - Output print file, <filename>.out
C FILNB  - Input file, <filename>.nb containing non-bonded interactions
C FILPS  - Name of output PostScript files
C FILRIN - Filename containing input data created by program CLEAN. Name is
C          of the form <filename>.rin
C FILRSA - Input file, <filename>.rsa containing residue accessibilities
C FILSCO - Output file, <filename>.sco containing overall stereochemical 
C          parameters
C FILSDH - Input file, <filename>.sdh containing log-odds scores for
C          various properties, created by program TPLOT
C FILSUM - Name of output .sum summary file
C GENERP - Percentage of residues in the generous regions of the Ramachandran
C          plot
C GFACTR(NDISTR) - G-factors giving measures of the stereochemical quality
C          of different aspects of the protein's structure
C GLOBC2 - Global mean value for Chi-2 trans
C GLOBDH - Global mean value for alpha-carbon chirality
C GLOBGM - Global mean value for Chi-1 gauche minus
C GLOBGP - Global mean value for Chi-1 gauche plus
C GLOBGT - Global mean value for Chi-1 trans
C GLOBHB - Global mean value for hydrogen bond energy
C GLOBHP - Global mean value for helix psi
C GLOBOM - Global mean value for omega dihedral angle
C GLOBPH - Global mean value for phi helix
C GLOBPP - Global mean value for proline phi
C GLOBSS - Global mean value for disulphide bond
C GLOB3L - Global mean value for Chi-3 left-handed
C GLOB3R - Global mean value for Chi-3 right-handed
C GMIDEV - Deviation from ideal value of Chi-1 gauche-minus values
C GMIMEA - Mean of Chi-1 gauche-minus values
C GMISTD - Standard deviation of gauche-minus Chi-1 values
C GPLDEV - Deviation from ideal value of Chi-1 gauche-plus values
C GPLMEA - Mean of Chi-1 gauche-plus values
C GPLSTD - Standard deviation of gauche-plus Chi-1 values
C GTRDEV - Deviation from ideal value of Chi-1 trans values
C GTRMEA - Mean of Chi-1 trans values
C GTRSTD - Standard deviation of trans Chi-1 values
C HBDEV  - Deviation from ideal value of H-bond energy values
C HBMEAN - Mean of H-bond energy values
C HBSTD  - Standard deviation of H-bond energy values
C HEADIN - Flag indicating whether the explanatory heading is required on
C          the asterisks print-out
C HPMEAN - Mean value of the helix psi values
C HPSTD  - Standard deviation of the helix psi values
C HWIDX, HWIDY - Marker size on Quality plots
C IDITIS - Flag indicating the IDITIS table output is required
C IEND   - Last character of filename holding original structure
C IFAIL  - Failure flag causing program to abort due to some error condition
C ILEN   - Length of filename holding original structure
C ILINE  - Line number on listing
C INCOLR(3) - Flags indicating whether colour PostScript files are required
C          for the different plot types
C IPAGE  - Listing page-number
C IPLOT  - Number of current PostScript output file
C ISTART - First character of filename holding original structure
C LIMRES - Number of standard deviations for highlighting residues on the
C          Residues plot
C MCHANG - Number of main-chain bond angle types
C MCHLEN - Number of main-chain bond length types
C MCHMAX(NMAIN) - Maximum value for each type of main-chain bond length and
C          angle
C MCHMIN(NMAIN) - Minimum value for each type of main-chain bond length and
C          angle
C MCHMEA(NMAIN) - Mean value for each type of main-chain bond length and
C          angle
C MCHNUM(NMAIN) - Number of values for each  type of main-chain bond
C          length and angle
C MCHSTD(NMAIN) - Standard deviation for each type of main-chain bond
C          length and angle
C MCLAST(3) - Morris et al classification
C MXCOLR - Maximum number of RGB colours that can be defined by the user
C MXCOLS - Number of columns in the asterisks print-out
C MXDEVC(MXCOLS) - Maximum deviation value in each column
C MXPINS - Number of residues per page of the pins plot
C MXPLOT - Maximum number of options for the graphs that can be plotted
C          on the Residues plot
C MXRES  - Maximum number of residues
C NALLOW - Number of residues in the allowed regions of the Ramachandran plot
C NALP   - Total number of angle/length/plane types held in the .lan file
C NALPHC - Number of alpha-carbon chirality values
C NAROMA - Number of aromatic + Aps residues in the protein
C NASN   - Number of asparagine residues in the protein
C NBAD   - Number of bad contacts
C NBADPR - Number of bad protein-protein contacts (ie not including hetatoms)
C NCHICH - Number of points plotted on the Chi-1 vs Chi-2 plot
C NCHI2  - Number of Chi-2 dihedral angles
C NCHI3L - Number of Chi-3 left-hand dihedral angles
C NCHI3R - Number of Chi-3 right-hand dihedral angles
C NCIS   - Number of cis-peptides in the structure
C NCOLOR - Number of user-defined colours in the parameter file
C NCORE  - Number of residues in the core regions of the Ramachandran plot
C NCOVAL - Number of covalent log-odds scores in .sdh file
C NCPOOL - Number of Chi-1 pooled values
C NDAM   - Number of D-amino acids in the structure
C NDIHED - Number of dihedral angle scores in the .sdh file
C NDISTR - Total number of log-odds scores in the .sdh file
C NDISUL - Number of disulphide bond values
C NENDCH - Number of end residues having meaningless phi-psi values
C NGAPX  - Number of ticks on x-axis of properties plots
C NGENER - Number of residues in the generous regions of the Ramachandran
C          plot
C NGLY   - Number of glycine residues in the protein
C NGMINU - Number of Chi-1 gauche-minus values
C NGPLUS - Number of Chi-1 gauche-plus values
C NHBOND - Number of H-bond energy values
C NHXPSI - Number of helix psi values
C NGTRAN - Number of Chi-1 trans values
C NMAIN  - Number of main-chain bond length and bond angle types
C NOMEGA - Number of omega dihedral angles
C NONBON - Flag signifying whether the non-bonded interactions file,
C          <filename>.nb, exists
C NONPG  - Number of non-Gly and non-Pro residues in the protein
C NOPLOT - Flag indicating whether or not any Simpleplot plots are required
C NOUTSI - Number of residues in the disallowed regions of the Ramachandran
C          plot
C NPHIHX - Number of phi helix values
C NPLANE - Number of types of planar groups held in the .lan file
C NPRINT - Number of the printout currently being produced - identifies
C          which column headings are required at the top of the page
C NPRO   - Number of proline residues
C NPROPH - Number of proline phi values
C NQUALP - Number of plots shown on the Quality Plots pages
C NRAMPL - Number of points plotted on the Ramachandran plot
C NRESID - Number of residues in protein
C NRGRAF - Number of graphs (0-3) to be plotted on the Residues plot
C NVALUE(NQUALP) - Number of values of the stereochemical parameters
C          shown on the Quality plots
C OMAXDV(2) - Maximum deviation score for the residue properties and
C          bond lengths/bond angles
C OMDEV  - Deviation from ideal value of the omega dihedral angle
C OMMEAN - Mean value of the omega dihedral angle
C OMSTD  - Standard deviation of the omega dihedral angle
C ONLYHI - Flag indicating whether only highlighted residues are to be
C          printed
C OOMAX, OOIMIN - Maximum and minimum Ooi numbers
C OUTSIP - Percentage of residues in the disallowed regions of the
C          Ramachandran plot
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C PHMEAN - Mean value of the phi helix values
C PHSTD  - Standard deviation of the phi helix values
C PINMAX(MXPLOT) - Maximum value for each of the Residue plots
C PLABEL - Flag indicating whether plot filename to be printed on the plot
C          itself
C PLOTPN - Flag signifying whether pins plot is required
C PLOTQ1 - Flag signifying whether main-chain quality plot is required
C PLOTQ2 - Flag signifying whether side-chain quality plot is required
C PLTNUM(3) - User-defined selection of which polts are to appear on
C          the Residues plot
C PLTSTD(MXPLOT) Standard deviation for each plot
C PLTYPE - Number indicating which of the plot-types is being processed
C POOMEA - Mean of Chi-1 pooled values
C POOSTD - Standard deviation of pooled Chi-1 values
C PPMEAN - Mean value of the proline phi values
C PPSTD  - Standard deviation of the proline phi values
C REND, RSTART - Resolution end and start for properties plots
C RESGF(MXRES) - Residue G-factor
C RESOL  - Protein's resolution
C RGAP   - Step-size for shading region on properties plots
C RGB(3,MXCOLR) - User-defined RGB colour definitions
C SDEVC2 - Global standard deviation for Chi-2 trans
C SDEVCS - Global standard deviation for cis-peptide omega dihedral angle
C SDEVDH - Global standard deviation for alpha-carbon chirality
C SDEVGM - Global standard deviation for Chi-1 gauche minus
C SDEVGP - Global standard deviation for Chi-1 gauche plus
C SDEVGT - Global standard deviation for Chi-1 trans
C SDEVHB - Global standard deviation for hydrogen bond energy
C SDEVHP - Global standard deviation for helix psi
C SDEVOM - Global standard deviation for omega dihedral angle
C SDEVPH - Global standard deviation for phi helix
C SDEVPO - Global standard deviation for pooled Chi-1 torsion angle
C SDEVPP - Global standard deviation for proline phi
C SDEVSS - Global standard deviation for disulphide bond
C SDEV3L - Global standard deviation for Chi-3 left-handed
C SDEV3R - Global standard deviation for Chi-3 right-handed
C SFACT  - Stereochemical R-factor, as calculated by program TPLOT and
C          held at the end of the .sdh file
C SHOWAC - Flag indicating whether accessibility shading is required
C SSMEAN - Mean value of the disulphide bond values
C SSSIZE - Half-height of wiring diagram figures
C SSSTD  - Standard deviation of the disulphide bond values
C STDLIM - Lower limit for standard deviation figures to be shown on
C          asterisks listing
C SVALUE(NQUALP) - Values of the stereochemical parameters shown on
C          the Quality plots
C VWDIST - van der Waals distance defining bad contacts
C WITHAN - Flag indicating whether file-handles are required for
C          PostScript output filenames
C WWWOUT - Flag indicating whether output .html file required
C XGAP   - Width of each residue on the Residues plot
C YBAND(NQUALP)  - Band widths for shaded curves on main-chain and
C          side-chain plots
C YHT    - Height of squares on the PROCHECKer board
C
C---------------------------------------------------------------------------
