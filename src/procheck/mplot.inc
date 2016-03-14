C----------------------------------------------------------------------+--- 
C MPLOT.INC - Include variables for MPLOT.FOR
C----------------------------------------------------------------------+--- 

CVAX      IMPLICIT NONE

CHECK v.3.4-->
C      INTEGER       MAXARR, MINOBS, MXCOLR, MXCONS, MXFILE, MXLEND,
C     -              MXMARK, MXPINS, MXRCON, MXRES, MXTORS
      INTEGER       MAXARR, MINOBS, MXCOLR, MXFILE, MXLEND,
     -              MXMARK, MXPINS, MXRCON, MXRES, MXTORS
CHECK v.3.4<--
      PARAMETER    (
     -              MAXARR  = 300000,
     -              MINOBS  =   10,
CHECK v.3.4-->
C     -              MXCOLR  =   20,
C     -              MXCONS  = 10000,
     -              MXCOLR  =   50,
CHECK v.3.4<--
     -              MXFILE  =   60,
     -              MXLEND  =    8,
     -              MXMARK  =    9,
     -              MXPINS  =   60,
     -              MXRCON  =  500,
CHECK v.3.5.1-->
C     -              MXRES   = 1000,
     -              MXRES   = 5000,
CHECK v.3.5.1<--
     -              MXTORS  =    4
     -             )

CHECK v.3.5.2-->
C      INTEGER       NACROS, NAMINO, NCELL, NCELL1, NCTYPE, NDISTR,
C     -              NDOWN, NDOWN2, NGFACT, NHIST, NPLOT, NPROPT,
C     -              NQUALP
      INTEGER       MXCELL, NACROS, NAMINO, NCTYPE, NDISTR, NDOWN,
     -              NDOWN2, NGFACT, NHIST, NPLOT, NPROPT, NQUALP
CHECK v.3.5.2<--
      PARAMETER    (
CHECK v.3.5.2-->
     -              MXCELL  =  7200,
CHECK v.3.5.2<--
     -              NACROS  =    4,
     -              NAMINO  =   20,
CHECK v.3.5.2-->
C     -              NCELL   =   45,
C     -              NCELL1  =   NCELL * NCELL,
CHECK v.3.5.2<--
CHECK v.3.4-->
C     -              NCTYPE  =    2,
     -              NCTYPE  =    3,
CHECK v.3.4<--
     -              NDISTR  =    6,
     -              NDOWN   =   10,
     -              NDOWN2  =    5,
     -              NGFACT  =    3,
     -              NHIST   =  180,
CHECK v.3.4-->
     -              NPLOT   = NDISTR + 4,
CHECK v.3.4<--
CHECK v.3.4.3-->
C     -              NPROPT  =    6,
     -              NPROPT  =    8,
CHECK v.3.4.3<--
CHECK v.3.3<--
C     -              NQUALP  =    4
     -              NQUALP  =    5
CHECK v.3.3-->
     -             )

CHECK v.3.4-->
      INTEGER       MAXRNG
      PARAMETER    (MAXRNG = 200)
CHECK v.3.4<--

CHECK v.3.4.3-->
C      INTEGER       LOCHI1, LOCHI2, LOCOME, LOCPHI, LOCPSI, LOCZET
      INTEGER       LOCHI1, LOCHI2, LOCHI3, LOCHI4, LOCOME, LOCPHI,
     -              LOCPSI, LOCZET
CHECK v.3.4.3<--
      PARAMETER    (
     -              LOCHI1  =    3,
     -              LOCHI2  =    4,
CHECK v.3.4.3-->
     -              LOCHI3  =    7,
     -              LOCHI4  =    8,
CHECK v.3.4.3<--
     -              LOCOME  =    5,
     -              LOCPHI  =    1,
     -              LOCPSI  =    2,
     -              LOCZET  =    6
     -             )

      INTEGER       MXPLOT, PLTCH1, PLTOME, PLTCAL
      PARAMETER    (
     -              MXPLOT =  3,
     -              PLTCH1 =  1,
     -              PLTOME =  2,
     -              PLTCAL =  3
     -             )

      REAL          DANGLE
      PARAMETER    (
     -              DANGLE = 360.0 / NHIST
     -             )

      REAL          DIARAD, MSIZE, PHIGHT, XFQLMG, XFQWID, YDAYLT,
     -              YDAYL2, YFQTMG, YFQTM2, YFQWID, YFQWD2
      PARAMETER    (
     -              DIARAD =   3.0,
     -              PHIGHT =   3.0,
     -              MSIZE  =   4.0,
     -              XFQLMG =  50.0,
     -              XFQWID = 100.0,
     -              YDAYLT =   4.0,
     -              YDAYL2 =  18.0,
     -              YFQTMG = 100.0,
     -              YFQTM2 =  85.0,
     -              YFQWID =  55.0,
     -              YFQWD2 = 118.0
     -             )

      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2
      PARAMETER    (
     -              BBOXX1 =  30,
     -              BBOXX2 = 550,
     -              BBOXY1 =  50,
     -              BBOXY2 = 780
     -             )

      REAL          HWIDX, HWIDY, SSSIZE, XGAP, XPORIG, XPWID, YCHQHT,
     -              YPORIG, YPSEP, YPWID
      PARAMETER    (
     -              HWIDX  = (BBOXX2 - BBOXX1) / 200.0,
     -              HWIDY  = (BBOXY2 - BBOXY1) / 200.0,
     -              SSSIZE =   4.0,
     -              XPORIG =  95.0,
     -              XPWID  = 425.0,
     -              XGAP   = XPWID / MXPINS,
     -              YCHQHT =   5.0,
     -              YPORIG = 355.0,
     -              YPSEP  =  47.0,
     -              YPWID  =  77.0
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

CHECK v.3.4-->
CHECK v.3.4.3-->
C      REAL          MINSHD, VALLOW, VALUPP
      REAL          MAXSHD, MINSHD, VALLOW, VALUPP
CHECK v.3.4.3<--
      PARAMETER    (
     -              MINSHD =  0.995,
CHECK v.3.4.3-->
     -              MAXSHD =  0.4,
CHECK v.3.4.3<--
     -              VALLOW =  0.00,
     -              VALUPP = -3.00
     -             )
CHECK v.3.4<--

CHECK v.3.4-->
C      CHARACTER*1   CHAIN, CONSST(MXRES), RESTYP(MXCONS),
      CHARACTER*1   CHAIN, CHNBRK(MXRES), CONSST(MXRES),
CHECK v.3.4<--
     -              STOSST(MXRES,MXFILE)
      CHARACTER*3   AMINO(NAMINO+1), VALRES(MXRES)
CHECK v.3.4-->
      CHARACTER*6   RESFRM(MAXRNG), RESTO(MAXRNG)
CHECK v.3.4<--
      CHARACTER*9   VALSNO(MXRES)
      CHARACTER*10  MKNAME(MXFILE)
      CHARACTER*12  COLNAM(MXCOLR)
      CHARACTER*30  DESC(NQUALP)
CHECK v.3.4-->
C      CHARACTER*40  TITLE
      CHARACTER*78  TITLE
CHECK v.3.4<--
      CHARACTER*80  BRCALL, BRCODE(MXFILE), FILPS, FILCNS, FILEDT,
CHECK v.3.4-->
C     -              FILRIN(MXFILE), FILRMS, FILTOR, FILVIO(NCTYPE),
C     -              PDBFIL
     -              FILRIN(MXFILE), FILRMS, FILRNG, FILTOR, PDBFIL
CHECK v.3.4<--
      CHARACTER*(MXLEND) FILID(MXFILE)
      COMMON/CHARS/
CHECK v.3.4-->
C     -              CHAIN, CONSST, RESTYP,
     -              CHAIN, CHNBRK, CONSST,
CHECK v.3.4<--
     -              STOSST,

     -              AMINO, VALRES,

CHECK v.3.4-->
     -              RESFRM, RESTO,
CHECK v.3.4<--

     -              VALSNO,

     -              MKNAME,

     -              COLNAM,

     -              DESC,

     -              TITLE,

     -              BRCALL, BRCODE, FILPS, FILCNS, FILEDT,
CHECK v.3.4-->
C     -              FILRIN, FILRMS, FILTOR, FILVIO,
C     -              PDBFIL,
     -              FILRIN, FILRMS, FILRNG, FILTOR, PDBFIL,
CHECK v.3.4<--

     -              FILID

CHECK v.3.4-->
C      INTEGER       ACONST, ASTSTR, AVESTR, COLPLT(13,NPLOT), DATLOC(3),
C     -              DVIOL(MXTORS,MXRES), ENASTR, FSTCON(MXRES), FSTRES,
C     -              IEND, ILEN, INTRLX, IPLOT, ISTART, LENID(MXFILE),
C     -              MAXFRQ, MAXLEN, MINLEN, MKTYPE(MXFILE), NALLOW,
C     -              NALPHC, NAROMA, NASN, NBAD, NCHICH, NCHI2, NCHI3L,
C     -              NCHI3R, NCIS, NCONST(NCTYPE), NCORE,
C     -              NCOUNT(NAMINO+1), NCPOOL, NDAM, NDISUL, NENDCH,
C     -              NFGAPY, NFOUND, NGENER, NGLY, NGMINU, NGPLUS,
C     -              NGTRAN, NHBOND, NHXPSI, NLABEL, NMODEL, NOMEGA,
C     -              NONPG, NOOI(MXRES), NOUTSI, NPHIHX, NPRO, NPROPH,
C     -              NUMSCO(MXRES + 1,MXFILE + 1), NRAMPL, NRES, NRGRAF,
C     -              NTICKX, NVALUE(NQUALP,MXFILE), NXTRES(MXRES),
C     -              OOIMAX, OOIMIN, PLTNUM(3), PLSTAT, PSLEN,
C     -              SSTCNT(3,MXRES), TCONST, TLEN, VIORES(2,MXCONS)
      INTEGER       ACONST, ACTNUM(MXFILE), ASTSTR, AVESTR,
     -              COLPLT(13,NPLOT), DATLOC(3), DVIOL(MXTORS,MXRES),
     -              ENASTR, FFILE, FSTCON(MXRES), FSTRES, IEND, ILEN,
     -              INTRLX, IPLOT, ISTART, LENID(MXFILE), MAXFRQ,
     -              MAXLEN, MINLEN, MKTYPE(MXFILE), MODFRM(MAXRNG),
     -              MODNUM(MXFILE), MODTO(MAXRNG), MRANGE, NALLOW,
     -              NALPHC, NAMLEN, NAROMA, NASN, NBAD, NCHICH, NCHI2,
     -              NCHI3L, NCHI3R, NCIS, NCOLOR, NCONST(NCTYPE), NCORE,
     -              NCOUNT(NAMINO+1), NCPOOL, NDAM, NDISUL, NENDCH,
     -              NFGAPY, NFOUND, NGENER, NGLY, NGMINU, NGPLUS,
     -              NGTRAN, NHBOND, NHXPSI, NLABEL, NMODEL, NOMEGA,
     -              NONPG, NOOI(MXRES), NOUTSI, NPHIHX, NPRO, NPROPH,
CHECK v.3.4.3-->
C     -              NUMSCO(MXRES + 1,MXFILE + 1), NRAMPL, NRANGE, NRES,
     -              NUMSCO(NDISTR + 1,MXRES + 1,MXFILE + 1),
     -              NRAMPL, NRANGE, NRES,
CHECK v.3.4.3<--
     -              NRGRAF, NTICKX, NVALUE(NQUALP,MXFILE), NWANTD,
     -              NXTRES(MXRES), OOIMAX, OOIMIN, PLTNUM(3), PLSTAT,
     -              PSLEN, SSTCNT(3,MXRES), TLEN, TOPMOD
CHECK v.3.4<--

CHECK v.3.4-->
C      COMMON/INTS/
C     -              ACONST, ASTSTR, AVESTR, COLPLT, DATLOC,
C     -              DVIOL, ENASTR, FSTCON, FSTRES,
C     -              IEND, ILEN, INTRLX, IPLOT, ISTART, LENID,
C     -              MAXFRQ, MAXLEN, MINLEN, MKTYPE, NALLOW,
C     -              NALPHC, NAROMA, NASN, NBAD, NCHICH, NCHI2, NCHI3L,
C     -              NCHI3R, NCIS, NCONST, NCORE,
C     -              NCOUNT, NCPOOL, NDAM, NDISUL, NENDCH,
C     -              NFGAPY, NFOUND, NGENER, NGLY, NGMINU, NGPLUS,
C     -              NGTRAN, NHBOND, NHXPSI, NLABEL, NMODEL, NOMEGA,
C     -              NONPG, NOOI, NOUTSI, NPHIHX, NPRO, NPROPH,
C     -              NUMSCO, NRAMPL, NRES, NRGRAF,
C     -              NTICKX, NVALUE, NXTRES,
C     -              OOIMAX, OOIMIN, PLTNUM, PLSTAT, PSLEN,
C     -              SSTCNT, TCONST, TLEN, VIORES
      COMMON/INTS/
     -              ACONST, ACTNUM, ASTSTR, AVESTR,
     -              COLPLT, DATLOC, DVIOL,
     -              ENASTR, FFILE, FSTCON, FSTRES, IEND, ILEN,
     -              INTRLX, IPLOT, ISTART, LENID, MAXFRQ,
     -              MAXLEN, MINLEN, MKTYPE, MODFRM,
     -              MODNUM, MODTO, MRANGE, NALLOW,
     -              NALPHC, NAMLEN, NAROMA, NASN, NBAD, NCHICH, NCHI2,
     -              NCHI3L, NCHI3R, NCIS, NCOLOR, NCONST, NCORE,
     -              NCOUNT, NCPOOL, NDAM, NDISUL, NENDCH,
     -              NFGAPY, NFOUND, NGENER, NGLY, NGMINU, NGPLUS,
     -              NGTRAN, NHBOND, NHXPSI, NLABEL, NMODEL, NOMEGA,
     -              NONPG, NOOI, NOUTSI, NPHIHX, NPRO, NPROPH,
     -              NUMSCO, NRAMPL, NRANGE, NRES,
     -              NRGRAF, NTICKX, NVALUE, NWANTD,
     -              NXTRES, OOIMAX, OOIMIN, PLTNUM, PLSTAT,
     -              PSLEN, SSTCNT, TLEN, TOPMOD
CHECK v.3.4<--

CHECK v.3.4-->
C      LOGICAL       BSHADE(NPLOT), CALALL, COMBPS, DONPAG,
C     -              DOUBL(MXRES), ENSEMB, HAVANG, HAVCON, HAVEG, HAVRMS,
C     -              HISTOG, IFAIL, INCOLR(NPLOT), NMR, PLABEL,
C     -              PLOTDS(NDISTR), PLOTMD, PLOTPN, PLOTQY, PLOTSS,
C     -              PLOTVI, RESMOD(MXRES,MXFILE), SHOWAC,
C     -              SHOWMA, SHOWSA, TWODEE, WITHAN
      LOGICAL       BSHADE(NPLOT), CALALL, COMBPS, DONPAG, DOUBL(MXRES),
     -              ENSEMB, HAVANG, HAVCON, HAVEG, HAVRAN, HAVRMS,
     -              HISTOG, IFAIL, INCOLR(NPLOT), MWANT(MXFILE), NMR,
     -              OMITAS, OMITCV, OMITGF,PLABEL, PLOTDS(NDISTR),
     -              PLOTPN, PLOTQY, PLOTSS, RESMOD(MXRES,MXFILE),
CHECK v.3.4.3-->
C     -              RSELEC, SHOWAC, SHOWMA, SHOWSA, TWODEE, WITHAN
     -              RSELEC, SHOWAC, SHOWMA, SHOWSA, SPLAYN, TWODEE,
     -              WITHAN
CHECK v.3.4.3<--
CHECK v.3.4<--

CHECK v.3.4-->
C      COMMON/LOGS/
C     -              BSHADE, CALALL, COMBPS, DONPAG,
C     -              DOUBL, ENSEMB, HAVANG, HAVCON, HAVEG, HAVRMS,
C     -              HISTOG, IFAIL, INCOLR, NMR, PLABEL,
C     -              PLOTDS, PLOTMD, PLOTPN, PLOTQY, PLOTSS,
C     -              PLOTVI, RESMOD, SHOWAC,
C     -              SHOWMA, SHOWSA, TWODEE, WITHAN
      COMMON/LOGS/
     -              BSHADE, CALALL, COMBPS, DONPAG, DOUBL,
     -              ENSEMB, HAVANG, HAVCON, HAVEG, HAVRAN, HAVRMS,
     -              HISTOG, IFAIL, INCOLR, MWANT, NMR,
     -              OMITAS, OMITCV, OMITGF,PLABEL, PLOTDS,
     -              PLOTPN, PLOTQY, PLOTSS, RESMOD,
CHECK v.3.4.3-->
C     -              RSELEC, SHOWAC, SHOWMA, SHOWSA, TWODEE, WITHAN
     -              RSELEC, SHOWAC, SHOWMA, SHOWSA, SPLAYN, TWODEE,
     -              WITHAN
CHECK v.3.4.3<--
CHECK v.3.4<--

      REAL          ACDEV,
     -              ACMEAN, ACSTD,
CHECK v.3.4.3-->
C     -              ARRAY(NCELL,NCELL,NAMINO+1), AVRACC(MXRES), BADCOP,
C     -              CH2DEV, CH2MEA, CH2STD, CIRCVR(NDISTR,MXRES),
C     -              COEFF(3,NQUALP),
CCHECK v.3.2.2-->
CC     -              COREPC, C3LMEA, C3LSTD, C3RMEA,
C     -              COREPC, CVCUT, C3LMEA, C3LSTD, C3RMEA,
CCHECK v.3.2.2<--
C     -              C3RSTD, DIHCON(2,MXTORS,MXRES),
C     -              DVIOLM(MXTORS,MXRES,MXFILE), GAPX,
C     -              GFACT(NGFACT + 1, MXRES), GFACTR(NGFACT + 1),
C     -              GMIDEV, GMIMEA, GMISTD, GPLDEV, GPLMEA, GPLSTD,
C     -              GTRDEV, GTRMEA, GTRSTD, HBDEV, HBMEAN, HBSTD,
C     -              HPMEAN, HPSTD, LIMRES, MAXX, MEAVAL(NQUALP), MINX,
C     -              NOETOP, NRMEAN(NAMINO+1), NRMSTD(NAMINO+1), OMDEV,
C     -              OMMEAN, OMSTD, PHMEAN, PHSTD, PINDMX, PINMAX(MXPLOT)
     -              AVRACC(MXRES), BADCOP, CH2DEV, CH2MEA, CH2STD,
     -              CIRCVR(NDISTR,MXRES), COEFF(3,NQUALP), COREPC,
     -              CVCUT, C3LMEA, C3LSTD, C3RMEA, C3RSTD,
     -              DIHCON(2,MXTORS,MXRES), DVIOLM(MXTORS,MXRES,MXFILE),
CHECK v.3.5.2-->
C     -              ENERGY(NCELL,NCELL,NAMINO+1), GAPX,
     -              GAPX,
CHECK v.3.5.2<--
CHECK v.3.4.3-->
C     -              GFACT(NGFACT + 1, MXRES), GFACTR(NGFACT + 1),
     -              GFACT(NDISTR + 1, MXRES), GFACTR(NDISTR + 1),
CHECK v.3.4.3<--
     -              GMIDEV, GMIMEA, GMISTD, GPLDEV, GPLMEA, GPLSTD,
     -              GTRDEV, GTRMEA, GTRSTD, HBDEV, HBMEAN, HBSTD,
     -              HPMEAN, HPSTD, LIMRES, MAXX, MEAVAL(NQUALP), MINX,
CHECK v.3.5.2-->
C     -              NOBSER(NCELL,NCELL,NAMINO+1),  NOETOP,
     -              NOETOP,
CHECK v.3.5.2<--
     -              NRMEAN(NAMINO+1), NRMSTD(NAMINO+1), OMDEV, OMMEAN,
     -              OMSTD, PHMEAN, PHSTD, PINDMX, PINMAX(MXPLOT)
CHECK v.3.4.3<--
      REAL
     -              PLTSTD(MXPLOT), POOMEA, POOSTD, PPMEAN, PPSTD,
     -              RESMIN, RESOL(MXFILE), RGB(3,MXCOLR),
     -              RMSALL(3,MXRES), RMSDEV(3,MXRES,MXFILE),
CHECK v.3.4.3-->
C     -              SAVOOI(MXRES,MXFILE), SCOMOD(MXRES + 1,MXFILE + 1),
     -              SAVOOI(MXRES,MXFILE),
     -              SCOMOD(NDISTR + 1,MXRES + 1,MXFILE + 1),
CHECK v.3.4.3<--
     -              SSMEAN, SSSTD, STATVL(NQUALP), STEP(2),
     -              SVALUE(NQUALP,MXFILE), VALBEG(2), VALEND(2),
     -              VALOOI(MXRES), VALUE(NPROPT,MXRES,MXFILE),
     -              XVALUE(MXFILE), YBAND(NQUALP)
      COMMON/REALS1/
     -              ACDEV,
     -              ACMEAN, ACSTD,
CHECK v.3.4.3-->
C     -              ARRAY, AVRACC, BADCOP,
C     -              CH2DEV, CH2MEA, CH2STD, CIRCVR,
C     -              COEFF,
CCHECK v.3.2.2-->
CC     -              COREPC, C3LMEA, C3LSTD, C3RMEA,
C     -              COREPC, CVCUT, C3LMEA, C3LSTD, C3RMEA,
CCHECK v.3.2.2<--
C     -              C3RSTD, DIHCON,
C     -              DVIOLM, GAPX,
C     -              GFACT, GFACTR,
C     -              GMIDEV, GMIMEA, GMISTD, GPLDEV, GPLMEA, GPLSTD,
C     -              GTRDEV, GTRMEA, GTRSTD, HBDEV, HBMEAN, HBSTD,
C     -              HPMEAN, HPSTD, LIMRES, MAXX, MEAVAL, MINX,
C     -              NOETOP, NRMEAN, NRMSTD, OMDEV,
C     -              OMMEAN, OMSTD, PHMEAN, PHSTD, PINDMX, PINMAX
     -              AVRACC, BADCOP, CH2DEV, CH2MEA, CH2STD,
     -              CIRCVR, COEFF, COREPC,
     -              CVCUT, C3LMEA, C3LSTD, C3RMEA, C3RSTD,
     -              DIHCON, DVIOLM,
CHECK v.3.5.2-->
C     -              ENERGY, GAPX,
     -              GAPX,
CHECK v.3.5.2<--
     -              GFACT, GFACTR,
     -              GMIDEV, GMIMEA, GMISTD, GPLDEV, GPLMEA, GPLSTD,
     -              GTRDEV, GTRMEA, GTRSTD, HBDEV, HBMEAN, HBSTD,
     -              HPMEAN, HPSTD, LIMRES, MAXX, MEAVAL, MINX,
CHECK v.3.5.2-->
C     -              NOBSER,  NOETOP,
     -              NOETOP,
CHECK v.3.5.2<--
     -              NRMEAN, NRMSTD, OMDEV, OMMEAN,
     -              OMSTD, PHMEAN, PHSTD, PINDMX, PINMAX
CHECK v.3.4.3<--
      COMMON/REALS2/
     -              PLTSTD, POOMEA, POOSTD, PPMEAN, PPSTD,
     -              RESMIN, RESOL, RGB,
     -              RMSALL, RMSDEV,
     -              SAVOOI, SCOMOD,
     -              SSMEAN, SSSTD, STATVL, STEP,
     -              SVALUE, VALBEG, VALEND,
     -              VALOOI, VALUE,
     -              XVALUE, YBAND

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

      SAVE /CHARS/, /INTS/, /LOGS/, /REALS1/, /REALS2/

C----------------------------------------------------------------------+--- 
C Parameters
C ----------
C
C ACDEV  - Deviation from ideal value of alpha-carbon chirality
C ACMEAN - Mean value for alpha-carbon chirality
C ACONST - Total number of dihedral angle constraints
C ACTNUM(MXFILE) - Actual model number corresponding to every entry in
C          MODNUM(MXFILE)
C ACSTD  - Standard deviation of alpha-carbon chirality
C AMINO(NAMINO) - Three-letter residue code
C ASTSTR - File-number of structure corresponding to the asterisked
C          structure
C AVRACC(MXRES) - Average accessibility for each residue (averaged over
C          the ensemble)
C AVESTR - File-number of structure corresponding to the average structure
C BADCOP - Number of bad contacts per 100 residues
C BBOXX1, BBOXX2, BBOXY1, BBOXY2 - Bounding box limits for PostScript area
C          on page. Entire plot is placed weithin this region
C BRCALL - Brookhaven code for ensemble
C BRCODE(MXFILE) - Individual Brookhaven codes for each member of the
C          ensemble
C BSHADE(NPLOT) - Flag indicating whether background shading required on plots
C CALALL - Flag indicating whether the G-factors and circular
C          variances, etc have to be calculated for all the distributions
C CHAIN  - User-defined chain - only this chain in the structure is to be
C          analysed (SPACE means all chains)
C CHNBRK(MXRES) - Markers indicating breaks in protein chain
C CH2DEV - Deviation from ideal value of the Chi-2 dihedral angle
C CH2MEA - Mean value of the Chi-2 dihedral angle
C CH2STD - Standard deviation of the Chi-2 dihedral angle
C CIRCVR(NDISTR,MXRES) - Circular variance for each frequency distribution
C COEFF(3,NQUALP) - Coefficients for shaded curves on main-chain and
C          side-chain plots
C COLPLT(5,NPLOT) - Colours for the different parts of each plot
C COLNAM(MXCOLR) - Names identifying the user-defined colours
C COMBPS - Flag indicating whether all pages of a given plot are to be
C          combined into a single PostScript file
C CONSST(MXRES) - Consensus secondary structure for each residue
C COREPC - Percentage of residues in the core regions of the Ramachandran
C          plot
C CVCUT  - Cut-off circular variance value for determing which frequency
C          plots are to be excluded
C C3LMEA - Mean value of the Chi-3 right-hand values
C C3LSTD - Standard deviation of the Chi-3 right-hand values
C C3RMEA - Mean value of the Chi-3 left-hand values
C C3RSTD - Standard deviation of the Chi-3 left-hand values
C DANGLE - Angular width of histogram bars in torsion-angle plots
C DATLOC(3) - Location of the data (in the VALUE array) corresponding to
C          each of the 3 graphs on the residue properties plot
C DESC(NQUALP) - Descriptions of the stereochemical parameters plotted on
C          the Quality Plots
C DIARAD - Radius of circular variance dials
C DIHCON(2,MXTORS,MXRES) - Dihedral angle constraints
C DONPAG - Flag indicating whether any duplicate residues have been printed
C          on the current page of Chi-1 values
C DOUBL(MXRES) - Flag indicating whether a given residue position is taken
C          by different residue types in the structures read in
C DVIOL(MXTORS,MXRES) - Number of constraint violations for the dihedral
C          angle constraints
C DVIOLM(MXTORS,MXRES,MXFILE) - Dihedral-angle violation for each residue
C          in each model
C ENASTR - File-number of structure corresponding to the energy-minimized
C          average structure
C ENERGY(MXCELL,NAMINO+1) - Log-odds scores for the different dihedral
C          angle distributions, obtained from the counts stored in NOBSER
C ENSEMB - Flag indicating whether Chi-1 plots are to be produced as
C          ensembles or as progressions of the values
C FFILE  - First file in a NMR ensemble that contains a wanted model
C FILPS  - Name of output PostScript files
C FILCNS - Filename for constraints and constraint-violations file
C FILEDT - Filename for output .edt file
C FILID(MXFILE) - Identifying name of each file in the ensemble
C FILRIN - Filename containing input data created by program SECSTR. Name is
C          of the form <filename>.rin
C FILRMS - Filename for RMS deviations file, created by RMSDEV.
C FILRNG - Name of input ranges file
C FILTOR - Filename for torsion angle constraints file, created by AQUA.
C FSTCON(MXRES) - Pointer to first constraint for the given residue
C FSTRES - Pointer to the first residue of the sequence
C GAPX   - Gap between ticks on the x-axes of the quality plots
C GFACT(NDISTR + 1,MXRES) - Averaged G-factors for each property for each
C          residue
C GFACTR(NDISTR + 1) - Overall G-factors for ensemble as a whole
C GMIDEV - Deviation from ideal value of Chi-1 gauche-minus values
C GMIMEA - Mean of Chi-1 gauche-minus values
C GMISTD - Standard deviation of gauche-minus Chi-1 values
C GPLDEV - Deviation from ideal value of Chi-1 gauche-plus values
C GPLMEA - Mean of Chi-1 gauche-plus values
C GPLSTD - Standard deviation of gauche-plus Chi-1 values
C GTRDEV - Deviation from ideal value of Chi-1 trans values
C GTRMEA - Mean of Chi-1 trans values
C GTRSTD - Standard deviation of trans Chi-1 values
C HAVANG - Flag indicating whether have torsion angle constraints
C HAVCON - Flag indicating whether have a constraints file
C HAVRAN - Flag indicating whether have a ranges file
C HAVRMS - Flag indicating whether have a RMS deviations file
C HAVEG  - Flag indicating whether have a G-factor distribution for the
C          current property being plotted
C HBDEV  - Deviation from ideal value of H-bond energy values
C HBMEAN - Mean of H-bond energy values
C HBSTD  - Standard deviation of H-bond energy values
C HISTOG - Flag indicating whether progressions required on distribution
C          plots, rather than histograms
C HWIDX, HWIDY - Half-widths of markers for points in Ramachandran plot
C IEND   - Last character of filename holding original structure
C IFAIL  - Failure flag causing program to abort due to some error condition
C ILEN   - Length of filename holding original structure
C INCOLR(NPLOT) - Flag indicating whether each plot is to be in colour or
C          black-and-white
C INTRLX - Integer/Real for numbers on the x-axis (0/1)
C IPLOT  - Number of current PostScript output file
C ISTART - First character of filename holding original structure
C LENID(MXFILE) - Length of each identifying name, FILID, for each file
C          in the ensemble
C LIMRES - Number of standard deviations for highlighting residues on the
C          Residues plot
C LOCHI1, LOCHI2, LOCPHI, LOCPSI - Array locations for stored properties
C MAXARR - Maximum size of constraints array
C MAXFRQ - Maximum frequency count for residues in the torsion-angles plot
C MAXRNG - Maximum number of model- or residue-ranges entered in the
C          ranges file
C MAXSHD - Maximum shading corresponding to a unfavourable G-factor
C MAXSTD - Maximum number of standard deviation values for which standard
C          deviations have to be calculated for the Quality Plots
C MINLEN, MAXLEN - Minimum and maximum sequence lengths found in the
C          .rin files
C MINSHD - Minimum shading corresponding to a favourable G-factor
C MAXX, MINX - Maximum and minimum values for the x-axes of the quality plots
C MEAVAL(NQUALP) - Mean value of each stereochemical parameter across all
C          the files read in
C MINOBS - Minimum number of observations in an array for calculation of
C          probabilities
C MKNAME(MXFILE) - Label to be attached to each point on the quality plots
C MKTYPE(MXFILE) - Marker-type for each point plotted
C MODFRM(MAXRNG), MODTO(MAXRNG) - Model-number ranges read in from the
C          ranges file
C MODNUM(MXFILE) - Model numbers (as defined by number of .rin file) of
C          those members of the ensemble that have been selected
C MRANGE - Number of model-number ranges in the ranges file
C MSIZE  - Marker size on 2D plots
C MWANT(MXFILE) - Flags indicating which models are required
C MXCELL - Maximum number of cells stored for each distribution
C MXCOLR - Maximum number of RGB colours that can be defined by the user
C MXRCON - Maximum number of distance restraints per residue
C MXDEVC(NQUALPS) - Maximum deviation value in each column
C MXFILE - Maximum number of .rin files allowed in the ensemble
C MXLEND - Maximum length of ensemble file-id that can be stored
C MXMARK - Maximum number of marker types that can appear on the quality
C          plots
C MXPINS - Number of residues per page of the pins plot
C MXPLOT - Maximum number of options for the graphs that can be plotted
C          on the Residues plot
C MXRES  - Maximum number of residues allowed for the protein structure
C          being plotted
C MXTORS - Maximum number of torsion-angle restraints per residue
C NACROS - Number of graphs across the page for frequency distrib plot
C NALLOW - Number of residues in the allowed regions of the Ramachandran plot
C NALPHC - Number of alpha-carbon chirality values
C NAMINO - Number of amino acid types
C NAMLEN - Length of the PDB code plus chain ID (if relevant)
C NAROMA - Number of aromatic + Aps residues in the protein
C NASN   - Number of asparagine residues in the protein
C NBAD   - Number of bad contacts
C NCELL1, NCELL2 - Dimensions of array holding distribution
C NCOUNT(NAMINO+1) - Count of residues of each residue type in the dataset
C NCHI2  - Number of Chi-2 dihedral angles
C NCHI3L - Number of Chi-3 left-hand dihedral angles
C NCHI3R - Number of Chi-3 right-hand dihedral angles
C NCIS   - Number of cis-peptides in the structure
C NCOLOR - Number of user-defined colours in the parameter file
C NCONRS(MXRES) - Number of constraints per residue
C NCONST(NCTYPE) - Total number of distance constraints of each type picked
C          up from the .rvio file
C NCORE  - Number of residues in the core regions of the Ramachandran plot
C NCPOOL - Number of Chi-1 pooled values
C NDAM   - Number of D-amino acids in the structure
C NDISTR - Number of different frequency distributions to be plotted
C NDISUL - Number of disulphide bond values
C NDOWN  - Number of graphs down the page for frequency distrib plot
C NDOWN2 - Number of graphs down the page for 2D plots
C NENDCH - Number of end residues having meaningless phi-psi values
C NFGAPY - Number of gaps along y-axis for frequency distrib plot
C NFILE  - Number of input .rin files
C NFOUND - Number of x-axis values found in the mplot.in file
C NGAPX  - Number of ticks on x-axis of properties plots
C NGENER - Number of residues in the generous regions of the Ramachandran
C          plot
C NGFACT - Number of calculated residue-by-residue G-factors
C NGLY   - Number of glycine residues in the protein
C NGMINU - Number of Chi-1 gauche-minus values
C NGPLUS - Number of Chi-1 gauche-plus values
C NGTRAN - Number of Chi-1 trans values
C NHBOND - Number of H-bond energy values
C NHIST  - Number of histogram bars in the frequency distrib plot
C NHXPSI - Number of helix psi values
C NLABEL - Number of point-labels found in the mplot.in file
C NMODEL - Number of model structures in the ensemble
C NMR    - Flag indicating that the ensemble is a set of NMR models coming
C          from a single file
C NOBSER(MXCELL,NAMINO+1) - Residue-by-residue counts of numbers of
C          points falling in each of the defined cells in the plot. In
C          routine CALCLO, these counts are converted into log-odds
C          scores, stored in array ENERGY
C NOETOP - Maximum y-value for upper-bound NOE violation graphs
C NOMEGA - Number of omega dihedral angles
C NONPG  - Number of non-Gly and non-Pro residues in the protein
C NOOI(MXRES) - Number of Ooi values stored at each residue position
C NOUTSI - Number of residues in the disallowed regions of the Ramachandran
C          plot
C NPHIHX - Number of phi helix values
C NPLOT  - Total number of plots available (= the frequency and 2D distribs
C          plus the two quality plots
C NPRO   - Number of proline residues
C NPROPH - Number of proline phi values
C NQUALP - Number of plots shown on the Quality Plots pages
C NRAMPL - Number of points plotted on the Ramachandran plot
C NRANGE - Number of residue-number ranges in the ranges file
C NRES   - Number of residues in the protein
C NRGRAF - Number of graphs (0-3) to be plotted on the Residues plot
C NTICKX - Number of ticks on the x-axes of the quality plots
C NPROPT - Number of values to be stored from the .rin files
C NUMSCO(NDISTR + 1,MXRES + 1,MXFILE + 1) - Number of values contributing to
C          the mean overall G-factor for each residue for each model in the
C          ensemble
C NVALUE(NQUALP) - Number of values of the stereochemical parameters
C          shown on the Quality plots for each file in the ensemble
C NWANTD - Number of models wanted
C NXTRES(MXRES) - Pointer to next residue in sequence (necessary to cope
C          with situation where have slight differences in the sequences
C          of the .rin files
C OMDEV  - Deviation from ideal value of the omega dihedral angle
C OMITAS, OMITCV, OMITGF - Flags indicating whether accessibility symbols,
C          circular variance dials and G-factors are to be omitted from
C          the dihedral angle distribution plots
C OMMEAN - Mean value of the omega dihedral angle
C OMSTD  - Standard deviation of the omega dihedral angle
C OOMAX, OOIMIN - Maximum and minimum Ooi numbers
C PDBFIL - Name of original file holding the coordinates of the structure
C          being processed
C PHIGHT - Size of point-marker on main graphs of residues plot
C PHMEAN - Mean value of the phi helix values
C PHSTD  - Standard deviation of the phi helix values
C PINDMX - Maximum RMS deviation read in
C PLABEL - Flag indicating whether plot filename to be printed on the plot
C          itself
C PLTNUM(3) - User-defined selection of which polts are to appear on
C          the Residues plot
C PLOTDS(NDISTR) - Flags indicating which of the frequency distributions
C          are required
C PLOTPN - Flag signifying whether pins plot is required
C PLOTQY - Flag indicating whether main-chain parameters plot to be produced
C PLOTSS - Flag signifying whether model-by-model secondary-structure
C          wiring diagram is required
C PLSTAT - Number indicating what the plot statistics for main-chain and
C          side-chain plots are to be based on: 0=No statistics, 1=Ave. of
C          structures, 2="A" structure, 3="E" structure, 4="*" structure
C PLTSTD(MXPLOT) - Standard deviation for each plot
C POOMEA - Mean of Chi-1 pooled values
C POOSTD - Standard deviation of pooled Chi-1 values
C PPMEAN - Mean value of the proline phi values
C PPSTD  - Standard deviation of the proline phi values
C PSLEN  - Length of PostScript filename
C REND, RSTART - Resolution end and start for properties plots
C RESFRM(MAXRNG), RESTO(MAXRNG) - Residue-number ranges read in from the
C          ranges file
C RESMOD(MXRES,MXFILE) - Flags indicating which residues are present in
C          each model (necessary for insertions and deletions)
C RESMIN - Minimum resolution - used for the comparison statistics
C RESOL(MXFILE) - Resolution given in each .sco file
C RGAP   - Step-size for shading region on properties plots
C RGB(3,MXCOLR) - User-defined RGB colour definitions
C RMS(3,MXRES,MXFILE) - Residue-by-residue deviations for each model,
C          giving the main-chain, side-chain and all-atom deviations
C RMSALL(3,MXRES) - Overall residue-by-reside RMS deviations
C RSELEC - Flag indicating whether residue-range selection on
C SAVOOI(MXRES,MXFILE) - Individual Ooi numbers
C SCOMOD(NDISTR + 1,MXRES + 1,MXFILE + 1) - Mean overall G-factor for
C          each residue for each model in the ensemble
C SHOWAC - Flag indicating whether accessibility shading is required
C SHOWMA - Flag indicating whether accessibility shading is required on the
C          model-by-model G-factors plot
C SHOWSA - Flag indicating whether accessibility shading is required on the
C          model-by-model secondary structures plot
C SPLAYN - Flag indicating whether splay-numbering of the models on the
C          2D distribution plots is required
C SSMEAN - Mean value of the disulphide bond values
C SSSIZE - Half-height of wiring diagram figures
C SSTCNT(3,MXRES) - Count of the Kabsch & Sander assignments for each
C          of the 3 different secondary structure types assigned for
C          each instance of this residue in the ensemble
C STATVL(NQUALP) - Values used for plot statistics on main-chain and
C          side-chain plots
C STEP(2) - Size of each cell in each direction
C STOSST(MXRES,MXFILE) - Stored secondary structures for each model
C SVALUE(NQUALP,MXFILE) - Values of the 10 stereochemical parameters shown on
C          the plots, and others, for each .rin file read in
C TITLE  - Title for plot
C TLEN   - Length of title
C TOPMOD - Highest model-number in the ensemble
C TWODEE - Flag indicating whether this is a 2D plot rather than a freq
C          distrib
C VALBEG(2), VALEND(2) - Start and end limits for the value pair selected
C VALLOW, VALUPP - Lower and upper G-factors defining the colour range
C          (values outside this range get the same colour as at the
C          range boundary
C VALOOI(MXRES) - Average Ooi numbers for each residue
C VALUE(NPROPT,MXRES,MXFILE) - Stored values for each residue in each
C          .rin file
C VALRES(MXRES) - Residue name for each residue in the sequence
C VALSNO(MXRES) - Sequence number for each residue in the sequence
C WITHAN - Flag indicating whether file-handles are required for
C          PostScript output filenames
C XFQLMG - Left margin for freq distribs plot
C XFQWID, YFQWID - Width and height of each litte graph in the freq
C          distribs plot
C XGAP   - Width of each residue on the Residues plot
C XMWID, YMWID - Width/height of graphs on properties plots
C XORIGN, YORIGN - Origin for properties plots
C XVALUE(MXFILE) - X-axis value assigned to the points from the given file
C YFQWD2 - Height of each 2D plot graph
C YDAYLT - "Daylight" gap between plots in the y-direction
C YDAYL2 - "Daylight" gap between 2D plots in the y-direction
C YFQBMG, YFQTMG - Bottom and top margins for freq distribs plot
C YBAND(NQUALP)  - Band widths for shaded curves on main-chain and side-chain
C          plots
C YCHQHT - Height of squares on the PROCHECKer board
C---------------------------------------------------------------------------
