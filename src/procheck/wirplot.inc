C----------------------------------------------------------------------+---
C WIRPLOT.INC - Include variables for WIRPLOT.FOR
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INTEGER       FILEN, MAXPDB, MINHEL, MINSTR, MXCOBJ, MXCOLR,
     -              MXCYS, MXDATA, MXPINS, MXPINW, MXPIN2, MXRES,
     -              MXSITE, NPPAGE
      PARAMETER    (
     -              FILEN  =   512,
     -              MAXPDB =    27,
     -              MINHEL =     3,
     -              MINSTR =     3,
     -              MXCOBJ =     8,
     -              MXCOLR =    20,
     -              MXCYS  =   200,
     -              MXDATA =  1000,
     -              MXPINS =   100,
     -              MXPINW =    75,
     -              MXPIN2 =    60,
     -              MXRES  = 10000,
     -              MXSITE =   500,
     -              NPPAGE =    11
     -             )

      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2
      PARAMETER    (
     -              BBOXX1 =    -5,
     -              BBOXX2 =   690,
     -              BBOXY1 =    -5,
     -              BBOXY2 =   880
     -             )

      REAL          HWIDX, HWIDY, SSSIZ, SSSIZW, XPORG, XPORW, XPWIDE,
     -              XPWIDW, XPWID2
      PARAMETER    (
     -              HWIDX  =   2.0,
     -              HWIDY  =   2.0,
     -              SSSIZ  =   3.0,
     -              SSSIZW =   4.0,
     -              XPORG  =  75.0,
     -              XPORW  =  10.0,
     -              XPWIDE = 465.0,
     -              XPWIDW = 580.0,
     -              XPWID2 = (XPWIDW * MXPIN2) / MXPINW
     -             )

      REAL          DASCOY, MAXDVY, NOTEX1, NOTEX2, NOTEY1, NOTEY2,
     -              SECSTY, SEQRMY, STSIZ, STSIZW, YPCOMP, YPDIFF, YPTOP
      PARAMETER    (
     -              DASCOY = 189.0,
     -              MAXDVY = 144.0,
     -              NOTEX1 = 112.0,
     -              NOTEX2 = 157.0,
     -              NOTEY1 =  22.0,
     -              NOTEY2 =  30.0,
     -              SECSTY =  20.0,
     -              SEQRMY =  17.0,
     -              STSIZ  =   7.5,
     -              STSIZW =  10.0,
     -              YPCOMP =  35.0,
     -              YPDIFF =  60.0,
     -              YPTOP  = 730.0
     -             )

      REAL          SIZLB, SIZLBW, TSIZE
      PARAMETER    (
     -              SIZLB  =   8.0,
     -              SIZLBW =  10.0,
     -              TSIZE  =  15.0
     -             )
      CHARACTER*1   CHAIN, CYSUPD(MXCYS), RESBRK(MXRES), RESCIS(MXRES),
     -              RESDAM(MXRES), RESCHN(MXRES), RESSST(MXRES)
      CHARACTER*3   CYSNUM(MXCYS), DATINF(MXDATA), DTYPE(MXDATA),
     -              RESNAM(MXRES), SITNAM(MXSITE)
      CHARACTER*5   RANGE(2,MXDATA), RESSEQ(MXRES)
      CHARACTER*12  COLNAM(MXCOLR)
      CHARACTER*78  BRCODE, FILPS, FILRIN, FILRSA, PDBFIL
      COMMON/CHARS/
     -              CHAIN, CYSUPD, RESBRK, RESCIS,
     -              RESDAM, RESCHN, RESSST,

     -              CYSNUM, DATINF, DTYPE,
     -              RESNAM, SITNAM,

     -              RANGE, RESSEQ,

     -              COLNAM,

     -              BRCODE, FILPS, FILRIN, FILRSA, PDBFIL

      INTEGER       BLEN, COLPLT(4), FLEN, IEND, ILEN, IPLOT, ISTART,
     -              NCIS, NCOLOR, NCYS, NDAM, NDATA, NDATAR, NDISUL,
     -              NDOM, NENDCH, NHELIX, NPAGE, NPINS, NRESID, NSITE,
     -              NSTRND, NUSED, OOIMAX, OOIMIN, RESDOM(MXRES),
     -              RESDNA(MXRES), RESITE(2,MXRES), RESLIG(MXRES),
     -              RESMET(MXRES)
      COMMON/INTS/
     -              BLEN, COLPLT, FLEN, IEND, ILEN, IPLOT, ISTART,
     -              NCIS, NCOLOR, NCYS, NDAM, NDATA, NDATAR, NDISUL,
     -              NDOM, NENDCH, NHELIX, NPAGE, NPINS, NRESID, NSITE,
     -              NSTRND, NUSED, OOIMAX, OOIMIN, RESDOM,
     -              RESDNA, RESITE, RESLIG,
     -              RESMET

      LOGICAL       ACCESS, CGIRUN, HAVCON, HAVDOM, HAVINT, HAVSIT,
     -              HAVSST, IFAIL, INCOLR, NOLABL, NOPRNT, NORINF,
     -              RESINT(3,MXRES), SHOWAC, SHOWDS, SHOWKY, SHOWRS,
     -              SUSED(MXSITE), USEPRM, WFORM
      COMMON/LOGS/
     -              ACCESS, CGIRUN, HAVCON, HAVDOM, HAVINT, HAVSIT,
     -              HAVSST, IFAIL, INCOLR, NOLABL, NOPRNT, NORINF,
     -              RESINT, SHOWAC, SHOWDS, SHOWKY, SHOWRS,
     -              SUSED, USEPRM, WFORM

      REAL          ACCMAX, CYSPOS(2,MXCYS), CYSRAD, RESACC(MXRES),
     -              RGB(3,MXCOLR), SSSIZE, XPORIG, XPWID

      COMMON/REALS/
     -              ACCMAX, CYSPOS, CYSRAD, RESACC,
     -              RGB, SSSIZE, XPORIG, XPWID

C----------------------------------------------------------------------+---
C Parameters
C ----------
C
C ACCMAX - Maximum accessibility value
C CGIRUN - Flag indicating that running in a CGI script
C CHAIN  - Chain identifier for current protein structure. If blank, then
C          all chains are to be included
C COLPLT(4) - Colours for: 1. Background colour, 2. secondary structure,
C          3. minimum accessibility, 4. maximum accessibility
C CYSNUM(MXCYS) - Identifying number for each disulphide bridge pair
C CYSPOS(2,MXCYS) - Postscript coordinates of each cysteine in a disulphide
C          bond
C CYSRAD - Radius of circles for cysteines in disulphide bridges
C CYSUPD(MXCYS) - Up-down marker indicating whether cysteine is above or
C          below the coil line
C DATINF(MXDATA) - Class/ID associated with each PROMOTIF data item
C DTYPE(MXDATA) - PROMOTIF data item type (ie beta-turn, gamma-turn, etc)
C FILEN  - Filename lengths
C FILPS  - Name of output PostScript files
C FILRIN - Filename containing input data created by program CLEAN. Name is
C          of the form <filename>.rin
C FILRSA - Input file, <filename>.rsa containing residue accessibilities
C FILSEQ - Sequence file, sequence.dat, created by pdbsplit.f (used for
C          CA-only files)
C FLEN   - Length of PostScript filename
C HAVCON - Flag indicating whether contacts data present
C HAVDOM - Flag indicating whether sequence contains any domains
C HAVINT - Flag indicating whether SITE records present in original PDB file
C HAVSIT - Flag indicating whether SITE records present in original PDB file
C HAVSST - Flag denoting whether structure contains any secondary structure
C          elements
C IEND   - Last character of filename holding original structure
C IFAIL  - Failure flag causing program to abort due to some error condition
C ILEN   - Length of filename holding original structure
C INCOLR - Flag denoting whether the plot is to be in colour
C IPLOT  - Number of current PostScript output file
C ISTART - First character of filename holding original structure
C MAXPDB - Maximum number of alternative locations for the PDB files
C MINHEL - Minimum helix length
C MINSTR - Minimum strand length
C MXCOBJ - Maximum number of different colours for the different objects
C MXCOLR - Maximum number of RGB colours that can be defined by the user
C MXCYS  - Maximum number of cysteine for plotting disulphide bonds
C MXDATA - Maximum number of data items that can be read in from the
C          PROMOTIF analysis files
C MXPINS - Number of residues per page of the pins plot
C MXPINW - Number of residues per page of the pins plot when done for
C          WWW pages
C MXPIN2 - Number of residues per page of the pins plot when done for
C          the second pass of the WWW pages
C MXRES  - Maximum number of residues for which secondary structure
C          information can be stored
C MXSITE - Maximum number of SITE details that can be stored
C NCIS   - Number of cis-peptides in the structure
C NCOLOR - Number of the available colours that are actually made use of
C NCYS   - Number of cysteines on the current PostScript page
C NDAM   - Number of D-amino acids in the structure
C NDATA  - Number of PROMOTIF data items read in
C NDATAR - Number of secondary structure elements read in from .rin file
C NDOM   - Number of domains read in from the domains file
C NGAPX  - Number of ticks on x-axis of properties plots
C NHELIX - Helix number
C NOLABL - Flag indicating that PROMOTIF helix and strand labelling not
C          to be applied
C NOPRNT - Flag indicating that no PRINT statements to be issued at start
C NORINF - Flag indicating that no .rin file found. Program will assume
C          that the file is a CA-only file.
C NPAGE  - PostScript page number
C NPINS  - Number os residues across the page
C NPPAGE - Number of diagrams fitted down the page
C NRESID - Number of residues in protein
C NSITE  - Number of different sites picked up from the SITE records
C NSTRND - Strand number
C NUSED  - Number of active sites actually used for this chain
C OOMAX, OOIMIN - Maximum and minimum Ooi numbers
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C PLTYPE - Number indicating which of the plot-types is being processed
C RANGE(2,MXDATA) - Residue range corresponding to each PROMOTIF data item
C RESACC(MXRES) - Estimated accessibility for each residue
C RESCHN(MXRES) - Chain for each residue
C RESCIS(MXRES) - Cis-residue marker for each residue
C RESDAM(MXRES) - D-amino acid marker for each residue
C RESDNA(MXRES) - Number of H-bonds and contacts to DNA
C RESDOM(MXRES) - Domain to which residue belongs
C RESINT(3,MXRES) - Intron flags for current residue position (phases 1-3)
C RESITE(2,MXRES) - Site identifiers for each residue, as read in from the
C                 SITE records
C RESLIG(MXRES) - Number of H-bonds and contacts to ligands
C RESMET(MXRES) - Number of H-bonds and contacts to metals
C RESNAM(MXRES) - 3-letter residue name for each residue
C RESSEQ(MXRES) - Residue sequence number for each residue
C RESSST(MXRES) - Secondary structure assignments for each residue
C RGB(3,MXCOLR) - User-defined RGB colour definitions
C SHOWAC - Flag indicating whether accessibility shading is required
C SHOWDS - Flag indicating whether disulphide bonds to be plotted
C SHOWKY - Flag indicating whether key to be plotted
C SHOWRS - Flag indicating whether residue numbers are required in plot
C SITNAM(MXSITE) - 3-letter identifier for each site described by the SITE
C          records in the PDB file
C SSSIZ, SSSIZE, SSSIZW - Sizes of secondary structure plot graphics
C SUSED(MXSITE) - Flags indicating which active sites are included in the
C          plot 
C USEPRM - Flag indicating whether PROMOTIF data to be used for the
C          secondary structure plot
C WFORM  - Flag indicating whether WWW format required
C XGAP   - Width of each residue on the Residues plot
C XPORG, XPORIG, XPORW - PostScript x-coord of left-hand margin of plot
C XPWID, XPWIDE, XPWIDW, XPWID2 - Different PostScript widths of plot
C
C---------------------------------------------------------------------------
