C-----------------------------------------------------------------------------
C
C CLEAN.F - Amended version of brkcln.f for use with structure-checking suite
C           of programs, PROCHECK
C
C           Amendments identified by CHECK v.m.n--> and CHECK v.m.n<--
C           comment lines, where m.n is the version number.
C
C v.1.0   - Original amendments
C                                                Roman Laskowski (Apr 1992)
C
C v.2.0   - Rewrite to allow entry of any filename for the coordinates file
C           rather than just the 4-letter Brookhaven code
C                                            Roman Laskowski (Jul/Aug 1992)
C
C v.2.0.1 - Minor amendment to highlight bad chirality of isoleucine
C           residues
C                                             Gail Hutchinson (14 Aug 1992)
C
C v.2.1   - Write out zero-occupancy atoms to .new file marked with ATZERO
C           and HEZERO, rather than ATOM and HETATM.
C         - Write out alternate occupancies to .new file marked with ATALT
C           and HEALT, rather than ATOM and HETATM.
C                                           Roman Laskowski (4-9 Dec 1992)
C
C v.2.1.1 - Amendment to treat case where occupancy field is blank (rather
C           than explicitly zero) as being equivalent to occupancy of 1.0
C                                             Roman Laskowski (1 Feb 1993)
C
C v.2.1.3 - Addition of asterisks to error prints so that can be plucked
C           out of the log files and displayed by the script file. Removal
C           of .clnlog file and transfer of all its output messages to the
C           standard output
C           Display of counts of number of zero occupancy and alternate
C           position atoms written out
C                                            Roman Laskowski (26 Mar 1993)
C
C v.2.1.4 - Minor change to error message where chirality of ILE residue is
C           wrong. Crude change to limit number of 'unknown amino acid' error
C           messages displayed (particularly irritating for lots of waters)
C                                            Roman Laskowski (12 May 1993)
C
C v.2.3   - Tiny change so that VAX .COM routines don't show error if no
C           asterisks found in log file
C                                            Roman Laskowski (14 Nov 1993)
C
C v.3.0.1 - Changed order of DATA statements in one routine so that they
C           follow the variable declarations (was causing problems when
C           compiling with the f2c compiler).
C           Commented out unreferenced labels (producing compiler warnings
C           when compiling on Convex).
C                                         Roman Laskowski (29/30 Mar 1994)
C
C v.3.1   - Amendments for PROCHECK-NMR to deal with each MODEL in an
C           ensemble separately.
C                                         Roman Laskowski (18/19 Apr 1994)
C
C v.3.2   - Residue sequence number replaced by residue ID in display of
C           which residues have had side-chain atoms swapped.
C                                            Roman Laskowski (30 Apr 1994)
C         - Amendment to make NMR pseudo-atoms (labelled with a Q in
C           pstn 2 of the atom name) be written straight out to the .new
C           file (eg like hydrogens are).
C         - Amendment to deal with slightly non-standard atom names
C           output by X-PLOR.
C                                            Roman Laskowski (13 May 1994)
C         - Minor amendments to various statements to make them
C           acceptable to f2c, and to deal with various uninitialised
C           variables. SAVE statements for variables moved prior to
C           declaration of those variables in the routine. (Amendments
C           supplied by Dave Love at Daresbury).
C                                 David Love/Roman Laskowski (13 Oct 1994)
C         - Amendment for case where one one model in the PDB file and
C           hence no MODEL and ENDMDL records
C                                            Roman Laskowski (25 Oct 1994)
C         - Reinstatement of facility to process a list of PDB files in
C           a single run. (For use with PROCOMP when comparing a set of
C           similar structures).
C                                           Roman Laskowski (4-7 Nov 1994)
C         - Amendment so that zero-occupancy atoms are not rejected for
C           NMR structures.
C                                            Roman Laskowski (28 Nov 1994)
C
C v.3.3.1 - Increased number of chain-breaks allowed, and checked for
C           overflow. (Can get a large number of apparent chain-breaks
C           in structures where one chain is normal, and the other has
C           C-alpha coordinates only)
C                                             Roman Laskowski (5 Jun 1995)
C
C v.3.3.3 - Bug-fix to stop program writing first atom out twice.
C                                             Roman Laskowski (1 Feb 1996)
C
C v.3.4   - Amendment for NMR structures to show which atoms are
C           incorrectly labelled, but not to actually perform the swap.
C                                             Roman Laskowski (1 Feb 1996)
C
C v.3.4.3 - Addition of warning message when an unrecognised atom-type is
C           encountered.
C           Correction of non-standard FORMAT statement.
C                                             Roman Laskowski (5 Jun 1996)
C
C v.3.4.4 - Increase in MAXRES to 10000 residues.
C                                            Roman Laskowski (15 Oct 1996)
C
C v.3.5   - Addition of Residue Definitions file and removal of hard-coded
C           definitions.
C                                            Roman Laskowski (31 Oct 1996)
C           Check for H-atom where H is in the 1st character position.
C                                            Roman Laskowski ( 6 Nov 1996)
C           Replacement of initialisation DATA statements for the large
C           by initilising loops. The DATA statements were resulting in 
C           the size of the executable file being enormous in, in the case
C           of Bernard Murray's Linux g77 compiler, making the program
C           impossible to compile.
C           Reinstatement of header-record reading routines which were
C           being bypassed since GETMOD routine added.
C                                            Roman Laskowski (20 Nov 1996)
C           Removal of condition where program gives up reading the PDB
C           file if it hasn't hit a protein residue by line 1000 (hence
C           fails on structures such as 1ber).
C                                            Roman Laskowski ( 8 May 1997)
C
C v.3.5.3 - Increase of maxalt parameter from 600 to 4000 to deal with
C           structures with a very large number of alternate atom
C           positions (eg 1egw).
C                                            Roman Laskowski (14 Nov 2000)
C
C           New flag to omit alternate atom location markers as some
C           programs (eg JESS) are thrown by them.
C                                            Roman Laskowski (26 Sep 2002)
C
C           Bug fix on overrun of aacod3 array for non-standard amino
C           acids.
C                                            Roman Laskowski (30 Sep 2002)
C v.3.6.2 - Increase size of altmax as program failing on some PDB entries.
C                                            Roman Laskowski (15 Mar 2004)
C
C----------------------------------------------------------------------+--- 
C
C Compilation
C
C f77 -Wimplicit -fbounds-check -o clean clean.f
C
C----------------------------------------------------------------------+--- 
C
C Program variables
C -----------------
C
C aacnt  - Count of residues read in
C aain(aapln) - 3-letter residue id-s read in from current SEQRES record
C aapln  - Number of residue names in each SEQRES record
C aastd(maxres) - Marker indicating whether this is a standard amino acid
C          or not
C abetsz - Maximum number of residue types that can be defined in resdefs
C          file
C brcode - PDB code, picked up from the HEADER record
C brkrec - Input line read in from PDB file
C class  - Protein class name from the HEADER record
C clslen - Length of protein class name in class
C compnd - Protein compound name from COMPND record
C cmplen - Length of protein compound name in compnd
C contch - Continuation character in pstn 10 of PDB record
C extdat(maxlin)*66 - Extra lines of data associated with residues (ie
C          additional to the identitfied atom records)
C extres(maxlin) - Residue number for to which corresponding extdat
C          record applies.
C inseqr*(3 * maxtyp) - Stored residue types, as read in from the
C          SEQRES records
C maxlin - Maximum total number of lines of extra records that can be
C          stored. Each line contains extra data associated with residues
C          (eg H-atom records, TER records, unidentified atoms, etc).
C maxtyp - Maximum number of different residue types that can be stored
C          from those given in the SEQRES records
C mxhtyp - Maximum number of HETATM types that can be defined in resdefs
C          file
C NMR    - Flag indicating whether this is an NMR structure, with MODEL
C          and ENDMDL records
C readok - Flag indicating whether file-reads successful
C resnam - Residue name read in from file
C restot(sideaa) - Total number of each type of residue in the file
C setloc - Alternate location identifier
C sideaa - Maximum number of all residue types that can be defined
C          (= abetsz + mxhtyp)
C soulen - Length of protein source name in source
C source - Protein source name from SOURCE record
C sqrsct - Number of resides from the SEQRES records stored
C stdaa  - Number of standard amino acids
C
C----------------------------------------------------------------------+--- 



      PROGRAM CLEAN

CHECK v.1.0-->
C      program brkcln
CHECK v.1.0<--
c     copyright by David Keith Smith, 1989
c*******************************************************************************
c     filename:      brkcln.f
c     last modified: 26th april 1991
c*******************************************************************************
c     global variables
c
      include  'brkcln.par'
c
CHECK v.3.2-->
C      character aacod1*(abetsz), aacod3*(3*abetsz), 
C     +          atname*((4*mnchna)-4), seqbcd(maxres)*6, 
C     +          ssord*21, brksym(maxres)*1, headln*132, 
C     +          altcod(maxres)*1, aastd(maxres)*1, atpres(3,maxres)*1,
C     +          thmprs(maxres)*1, hetcod*63, outsty*1, dudcod*105, 
C     +          modcod*(nmod*3), stdmod*(nmod+1), nstdmd*(nmod+1), 
C     +          modpl*(nmod), scname(sideaa)*44, altnam(3)*44, cbet*4, 
      character aacod3*(3*abetsz), 
     +          atname*((4*mnchna)-4), seqbcd(maxres)*6, 
     +          brksym(maxres)*1,
     +          altcod(maxres)*1, aastd(maxres)*1, atpres(3,maxres)*1,
CHECK v.3.5-->
C     +          thmprs(maxres)*1, hetcod*63, outsty*1, dudcod*105, 
C     +          modcod*(nmod*3), stdmod*(nmod+1), nstdmd*(nmod+1), 
C     +          modpl*(nmod), scname(sideaa)*44, altnam(3)*44,
CCHECK v.3.2<--
C     +          allspc*44, lab(maxres)*1, extdat(maxlin)*66, 
C     +          tphnam*((4*mnchna)-4), olenam*((4*mnchna)-4), 
C     +          pphnam*((4*mnchna)-4), chstop(maxres)*1,
C     +          chstrt(maxres)*1, mcalt(mnchna,maxres)*1, 
C     +          scalt(sdchna,maxres)*1, altmcc(maxres)*1, 
C     +          altscc(maxres)*1, brcode*4,
C     +          pyrnam*((4*mnchna)-4), altdat(maxalt)*66, 
C     +          batnam*((4*mnchna)-4)
     +          thmprs(maxres)*1, hetcod*(3*MXHTYP), outsty*1,
     -          dudcod*105, MCNAME(sideaa)*44,
     +          modcod*(nmod*3), stdmod*(nmod+1), nstdmd*(nmod+1), 
     +          modpl*(nmod), scname(sideaa)*44, altnam(3)*44,
     +          allspc*44, lab(maxres)*1, extdat(maxlin)*66, 
     +          chstop(maxres)*1, chstrt(maxres)*1,
     +          mcalt(mnchna,maxres)*1,  scalt(sdchna,maxres)*1,
     +          altmcc(maxres)*1, altscc(maxres)*1, brcode*4,
     +          altdat(maxalt)*66
CHECK v.3.5<--



c
c     character autsum(maxres)*1, dsdef(maxrse)*1, 
c    +          ssbond(nstruc,maxres)*1, summpl(maxres)*1,
c    +          sssumm(maxres)*1
c   
      integer   seqcod(maxres), seqlen, chnsz(maxchn), chnct, hedlen, 
     +          chnend(maxchn), nscats(sideaa), mcser(mnchna,maxres), 
     +          scser(sdchna,maxres), naltmc, naltsc, strlen, 
     +          altsci(3,altmax), altmci(3,altmax), extres(maxlin), 
     +          nline, nrec, totswp(sideaa), ntot, restot(sideaa), 
     +          hist(17), hismin(17), chnlt, altlin, altres(maxalt)

CHECK v.3.5-->
      INTEGER   altpos(sideaa), extprc(sideaa), j, k, natsdc(sideaa),
     -          sswap(2,sideaa)
CHECK v.3.5<--
c
c     integer   hbond(maxbnd,maxres), bridpt(2,maxres), dspart(maxres),
c    +          oois(2,maxres), j
c
      real      mcaaat(ncoord,mnchna,maxres), 
     +          scaaat(ncoord,sdchna,maxres), mcangs(nmnang,maxres), 
     +          scangs(nsdang,maxres), mcocc(mnchna,maxres), 
     +          mcthrm(mnchna,maxres), scocc(sdchna,maxres), 
     +          scthrm(sdchna,maxres), altmcr(2,altmax), 
     +          altscr(2,altmax), altmca(3,altmax), altsca(3,altmax), 
     +          sx2, sx, ave, sd
c
c     real      hbonde(maxbnd,maxres), cadsts(11,maxres), 
c    +          dsdsts(3,maxres)
c
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres), caonly, 
     +          readok, astruc, altern
c
      parameter (
     +          allspc = '                                            ')
c
c     local variables
c
      integer   finish, opnera, opnerb, opnerc, opnerd, 
     +          opnere, opnerf, opnerg, opnerh, fiend, ok, endnfi
c
      character indirn*1, logfil*80
c
      parameter (finish = 2, opnera = 1, opnerb = 3, fiend = 4, ok = 0, 
     +           indirn = '@', opnerc = 5, endnfi = 6, opnerd=7,
     +           opnere=9, opnerf=11, opnerg=13, opnerh=15)
c
      logical   finuse
c
c     integer   filety
c
CHECK v.2.1-->
C      integer   fistat, namlen, iocode, i, corepl, corend,
C     +          logpl, dirpl, nodepl, outlen
      integer   fistat, namlen, iocode, i,
     +          outlen
CHECK v.2.1<--
c
      character finam*(fnamln), brknam*(fnamln), outnam*(fnamln),
CHECK v.2.1-->
C     +          altnm*(fnamln), outdir*(fnamln)
     +          altnm*(fnamln)
CHECK v.2.1<--
CHECK v.2.1-->
C      integer   ifail, outln, brklen, lastsl
CHECK v.2.1<--
      character wkres*3
      character pronam*11
      character chnnam*9
CHECK v.2.1-->
C      character brkid*4
CHECK v.2.1<--

CHECK v.2.0-->
      CHARACTER*78  PDBFIL
      INTEGER       IEND, ILEN, ISTART
CHECK v.3.5-->
C      LOGICAL       IERROR
      LOGICAL       IERROR, IFAIL
      INTEGER       AAPLN, INSLEN, SQRSCT
      PARAMETER    (AAPLN = 13)
      CHARACTER*3   AAIN(AAPLN)
CHECK v.3.5<--
CHECK v.2.0<--

CHECK v.2.1.3-->
      INTEGER       NATALT, NATZER, NHEALT, NHEZER
CHECK v.2.1.3<--

CHECK v.3.1-->
      INTEGER       IMODEL, NMODEL
CHECK v.3.2-->
C      LOGICAL       ENDFIL, FIRST, NMR
CHECK v.3.5.3-->
C      LOGICAL       ENDFIL, ENSEMB, FIRST, NMR
      LOGICAL       ENDFIL, ENSEMB, FIRST, NMR, NOALT
CHECK v.3.5.3<--
CHECK v.3.2<--
CHECK v.3.1<--
      CHARACTER inseqr*(3*MAXTYP)

CHECK v.3.2-->
      integer mcsern, scsern, sctn
      parameter (mcsern=(mnchna*maxres), scsern=(sdchna*maxres),
     +     sctn=(ncoord*sdchna*maxres))
CHECK v.3.2<--
c
c     input file control loop
c
      data hismin(1) / -60 /
      data hismin(2) / -50 /
      data hismin(3) / -40 /
      data hismin(4) / -30 /
      data hismin(5) / -20 /
      data hismin(6) / -10 /
      data hismin(7) /   0 /
      data hismin(8) /  10 /
      data hismin(9) /  20 /
      data hismin(10) / 30 /
      data hismin(11) / 40 /
      data hismin(12) / 50 /
      data hismin(13) / 60 /
      data hismin(14) / 70 /
      data hismin(15) / 80 /
      data hismin(16) / 90 /
      data hismin(17) / 100 /
      data sx2 / 0.0 /
      data sx / 0.0 /
      data ntot / 0 /
CHECK v.3.2-->
C      data aacod1(1:29) / 'ABCDEFGHIXKLMNXPQRSTXVWXYZX--' /
CHECK v.3.2<--
CHECK v.3.5-->
C      data aacod3(1:39)  / 'ALAASXCYSASPGLUPHEGLYHISILEXXXLYSLEUMET' /
C      data aacod3(40:78) / 'ASNXXXPROGLNARGSERTHRXXXVALTRPXXXTYRGLX' /
C      data aacod3(79:87) / 'UNKPCAINI' /
C      data hetcod(1:39)  / 'AIBPHLSECALMMPRFRDPYRLYMGLMPPHPGLOLETPH' /
C      data hetcod(40:63) / 'ABANLEB2VB2IB1FBNOB2AB2F' /
CHECK v.3.5<--
      data modcod(1:39)  / 'ACEANIBZOCBZCLTFORNH2PHORHATFATOS NHMYR' /
      data modcod(40:42) / 'BOC' /
CDEBUG
C      data dudcod(1:39) / '  A  C  G  T  U1MA5MCOMC1MG2MGM2G7MGOMG' /
C      data dudcod(40:78) / ' YG  I +UH2U5MUPSU +C +AGALAGLGCUASGMAN' /
C      data dudcod(79:105) / 'GLCCEGG4SAGSNAGGLSNGSNAMEXC' /
      data dudcod(1:39) / '                                       ' /
      data dudcod(40:78) / '                                       ' /
      data dudcod(79:105) / '                           ' /
CDEBUG
      data stdmod / 'ABCDEFGHIJKLZ' /
      data nstdmd / 'abcdefghijklz' /
CHECK v.3.5-->
C      data modpl  / '+0++++0++++0' /
      data modpl  / '+0++++0++++0++' /
CHECK v.3.5<--
CHECK v.3.2-->
C      data ssord / 'HEheBbGgIiJKLMNOPTtS ' /
CHECK v.3.2<--
      data atname / ' N   CA  C   O      ' /
CHECK v.3.5-->
C      data tphnam / ' N   CA  C   O1  O2 ' /
C      data pphnam / ' N   CA  P   OP1 OP2' /
C      data olenam / ' ON  CA  C   O      ' /
C      data pyrnam / ' ON  CA  C   O      ' /
C      data batnam / ' N   CA  B   O1  O2 ' /
CHECK v.3.5<--
CHECK v.3.4.4-->
C      data atomin / 34200*.false. /
C      data mcaaat / 102600*0 /
CHECK v.3.5-->
C      data atomin / 60000*.false. /
C      data mcaaat / 180000*0 /
CHECK v.3.5<--
CHECK v.3.4.4<--
CHECK v.3.5-->
C      data nscats / 1, 4, 2, 4, 5, 7, 0, 10, 4, -1, 5, 3*4, -1, 3, 5, 
C     +              7, 2, 3, -1, 3, 10, -1, 8, 5, -1, 4, 12, 3, 7, 4, 
C     +              2*2, 7, 8, 6, 2, 10, 4, 5, 9, 2, 4, 3, 4, 7, 4, 
C     +              1, 7 /
C      data scname(1) / allspc /
C      data scname(2) / ' CG  AD1 AD2                                ' /
C      data scname(3) / ' SG                                         ' /
C      data scname(4) / ' CG  OD1 OD2                                ' /
C      data scname(5) / ' CG  CD  OE1 OE2                            ' /
C      data scname(6) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(7) / allspc /
C      data scname(8) / ' CG  ND1 CD2 CE1 NE2 AD1 AD2 AE1 AE2        ' /
C      data scname(9) / ' CG1 CG2 CD1                                ' /
C      data scname(10) / allspc /
C      data scname(11) / ' CG  CD  CE  NZ                             ' /
C      data scname(12) / ' CG  CD1 CD2                                ' /
C      data scname(13) / ' CG  SD  CE                                 ' /
C      data scname(14) / ' CG  OD1 ND2                                ' /
C      data scname(15) / allspc /
C      data scname(16) / ' CG  CD                                     ' /
C      data scname(17) / ' CG  CD  OE1 NE2                            ' /
C      data scname(18) / ' CG  CD  NE  CZ  NH1 NH2                    ' /
C      data scname(19) / ' OG                                         ' /
C      data scname(20) / ' OG1 CG2                                    ' /
C      data scname(21) / allspc /
C      data scname(22) / ' CG1 CG2                                    ' /
C      data scname(23) / ' CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2        ' /
C      data scname(24) / allspc /
C      data scname(25) / ' CG  CD1 CD2 CE1 CE2 CZ  OH                 ' /
C      data scname(26) / ' CG  CD  AE1 AE2                            ' /
C      data scname(27) / allspc /
C      data scname(28) / ' CG  CD  OE                                 ' /
C      data scname(29) / ' CG  CD  CE  NZ  CI1 CI2 CI3 CI4 NI2 CI5 CI6' /
C      data scname(30) / ' CB1 CB2                                    ' /
C      data scname(31) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(32) / 'SEG  OD1 OD2                                ' /
C      data scname(33) / ' CM                                         ' /
C      data scname(34) / ' SG                                         ' /
C      data scname(35) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(36) / '                                            ' /
C      data scname(37) / ' CG  CD  CE  NZ  CM                         ' /
C      data scname(38) / ' CM                                         ' /
C      data scname(39) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(40) / ' P   O1P O2P                                ' /
C      data scname(41) / ' CG  CD1 CD2                                ' /
C      data scname(42) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(43) / ' CG                                         ' /
C      data scname(44) / ' CG  CD  CE                                 ' /
C      data scname(45) / ' CG1 CG2                                    ' /
C      data scname(46) / ' CG1 CG2 CD1                                ' /
C      data scname(47) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
C      data scname(48) / ' CG  CD  CE                                 ' /
C      data scname(49) / allspc /
C      data scname(50) / ' CG  CD1 CD2 CE1 CE2 CZ                     ' /
CHECK v.3.5<--
      data altnam(1) / ' CG  AD1 AD2                                ' /
      data altnam(2) / ' CG  AD1 AD2                                ' /
      data altnam(3) / ' CG  CD  AE1 AE2                            ' /
CHECK v.3.2-->
C      data cbet / ' CB ' /
CHECK v.3.2<--
      data fistat / 0 /
      data finuse / .false. /
      data pronam / 'ldamino.dat' /
      data chnnam / 'chain.dat' /
CHECK v.3.2-->
CHECK v.3.5-->
C      data scaaat /sctn*0.0/,  mcser /mcsern*0/, scser/scsern*0/,
C     +     mcalt /mcsern*' '/, mcthrm /mcsern*0.0/, hist /17*0/,
C     +     scthrm /scsern*0.0/, extdat /maxlin*' '/,
C     +     scalt /scsern*' '/, nrec/0/
      data hist /17*0/, nrec /0/
CHECK v.3.5<--
CHECK v.3.2<--
c
CHECK v.3.0.1-->
C  209 format (a4, 1x, i2, 1x, i2, 1x, f5.2, 1x, f5.2)
C  210 format (a4, 1x, i4)
C  211 format (a4, 1x, i2, 1x, i2, 1x, f5.2, 1x, f5.2, 1x, a2)
CHECK v.3.0.1<--

CHECK v.1.0-->
C 1000 format (' Brookhaven clean-up program BRKCLN - copyright by',
C     + ' David Keith Smith, 1989')
C 1010 format (/,' Enter the name of the next file to be processed.',/,
C     + '  Just a file name for a Brookhaven format file,',/,
C     + '  ', A1, 'filename for a file holding names of Brookhaven',
C     + ' format files,',/,
C     + '  <RETURN> to exit.')
CHECK v.1.0<--

CHECK v.3.0.1-->
C 1020 format (A)
C 1030 format (' Working on Brookhaven file ',A,' . . . ')
CHECK v.3.0.1<--

CHECK v.3.2-->
 1020 format (A)
 1030 format (' Working on Brookhaven file ',A,' . . . ')
CHECK v.3.2<--

CHECK v.2.1.3-->
C 1040 format (' Unable to open output file "',A, '"') 
C 1050 format (' Unable to open the list file "', A, '"')
C 1060 format (' Unable to open the Brookhaven file "', A, '"')
C 1070 format (' Unable to open output file "', A, '"')
 1040 format (' **** Unable to open output file "',A, '"') 
CHECK v.3.2-->
C 1050 format (' **** Unable to open the list file "', A, '"')
CHECK v.3.2<--
 1060 format (' **** Unable to open the Brookhaven file "', A, '"')
 1070 format (' **** Unable to open output file "', A, '"')
CHECK v.2.1.3<--
CHECK v.3.4-->
 1080 format (' Total number of residues changed and total number',
     + ' of residues:')
 1082 format (' Total number of residues that need to be changed and t',
     + 'otal number of residues:')
CHECK v.3.4<--
 1090 format (5X, A3, 2i6)
CHECK v.2.1.3-->
C 1100 format (' C-alpha only file')
C 1110 format (' Error while reading file')
C 3001 format (' Error reading input file on unit ',i3,
C     +        ' after ',i5,' lines')
C 7000 format ('Error opening file "',A,'" for input.')
C 7001 format ('Error opening file "',A,'" for output.')
 1100 format (' **** C-alpha only file')
 1110 format (' **** Error while reading file')
CHECK v.3.0.1-->
C 3001 format (' **** Error reading input file on unit ',i3,
C     +        ' after ',i5,' lines')
C 7000 format ('**** Error opening file "',A,'" for input.')
CHECK v.3.0.1<--
 7001 format ('**** Error opening file "',A,'" for output.')
CHECK v.2.1.3<--
c
      write (*, 1000)

CHECK v.2.1.3-->
C---- Initialise counts of zero- and alternate-occupancy atoms written out
      NATZER = 0
      NHEZER = 0
      NATALT = 0
      NHEALT = 0
CHECK v.2.1.3<--
CHECK v.3.1-->
      NMODEL = 0
CHECK v.3.1<--

CHECK v.3.5-->
C---- Initialise all the variables previously initialised in DATA
C     statements
      DO 30, k = 1, maxres
          DO 20, j = 1, mnchna
              atomin(j,k) = .FALSE.
              mcalt(j,k) = ' '
              mcser(j,k) = 0
              mcthrm(j,k) = 0.0
              DO 10, i = 1, ncoord
                  mcaaat(i,j,k) = 0.0
 10           CONTINUE
 20       CONTINUE
 30   CONTINUE
      DO 60, k = 1, maxres
          DO 50, j = 1, sdchna
              scalt(j,k) = ' '
              scser(j,k) = 0
              scthrm(j,k) = 0.0
              DO 40, i = 1, ncoord
                  scaaat(i,j,k) = 0.0
 40           CONTINUE
 50       CONTINUE
 60   CONTINUE
      DO 80, i = 1, maxlin
          extdat(i) = ' '
 80   CONTINUE
      inslen = 0
      inseqr = ' '
      sqrsct = 0
CHECK v.3.5<--


CHECK v.1.0-->
C---- Print initial heading
 1000 FORMAT (//,
     -    ' Clean-up program CLEAN - Written by David Keith Smith,',
     -    ' 1989.',/,
     -    '                          Amended by Roman Laskowski,',
     -    ' 1992.',//)

CHECK v.3.5-->
C---- Read in the Residue Definitions from the resdefs.dat file
      CALL GETDEF(AACOD3,ABETSZ,HETCOD,MXHTYP,NSCATS,SCNAME,MCNAME,
     -    SIDEAA,NATSDC,EXTPRC,ALTPOS,SSWAP,IFAIL)
      IF (IFAIL) GO TO 990
CHECK v.3.5<--

CHECK v.2.0-->
C---- Accept 4-letter Brookhaven code
C      WRITE (*, 1010)
C 1010 FORMAT (' Enter 4-letter Brookhaven code of structure:')
C      READ (*,1020) BRCODE

C---- Accept full path to input .pdb file
C      WRITE (*, 1015)
C 1015 FORMAT (//,
C     -    ' Enter full path to .pdb file holding coordinates.',/,
C     -    ' Note: File name must be of the form [path]p<code>.pdb')
C      READ (*,1020) FINAM
C      NAMLEN = INDEX(FINAM,SPACE) - 1
C      IF (NAMLEN.LE.0) THEN
C          BRKNAM = 'p' // BRCODE // '.pdb'
C          NAMLEN = 0
C      ELSE
C          BRKNAM = FINAM(1:NAMLEN) // 'p' // BRCODE // '.pdb'
C      ENDIF
C      NAMLEN = NAMLEN + 9
C      WRITE (*,1018) BRKNAM(1:NAMLEN)
C 1018 FORMAT(//,
C     -    ' Processing file: ',A,/)

C---- Form the names of the output files
C      OUTNAM = 'p' // BRCODE // '.new'
C      ALTNM  = 'p' // BRCODE // '.alt'
C      LOGFIL = 'p' // BRCODE // '.clnlog'
C      OUTLEN = 9
CHECK v.2.0<--
CHECK v.1.0<--
 
CHECK v.2.0-->
C---- Accept name of original .pdb file holding the structure
      PRINT*, 'Enter filename containing coordinates of structure'
CHECK v.3.2-->
      PRINT*, '  (for file containing ensemble of NMR structures en',
     -    'ter @filename;'
      PRINT*, '   for set of separate PDB files to be processed, en',
     -    'ter %filelist,'
      PRINT*, '   where filelist contains a list of PDB files to be',
     -    ' cleaned up)'
CHECK v.3.2<--
CHECK v.3.5-->
C      READ(*,10) PDBFIL
C10    FORMAT(A)
      READ(*,170) PDBFIL
 170  FORMAT(A)
CHECK v.3.5<--

CHECK v.3.5.3-->
C---- Initialise flag
      NOALT = .FALSE.
CHECK v.3.5.3<--

CHECK v.3.1-->
C---- If the very first character of the file-name is an @, then the
C     PDB file is an NMR file containing several models of the structure
      IF (PDBFIL(1:1).EQ.'@') THEN
          NMR = .TRUE.
          PDBFIL = PDBFIL(2:)
          FIRST = .TRUE.
CHECK v.3.2-->
          ENSEMB = .FALSE.
      ELSE IF (PDBFIL(1:1).EQ.'%') THEN
          ENSEMB = .TRUE.
          FIRST = .TRUE.
          NMR = .FALSE.
CHECK v.3.2<--
CHECK v.3.5.3-->
      ELSE IF (PDBFIL(1:1).EQ.'!') THEN
          PDBFIL = PDBFIL(2:)
          ENSEMB = .FALSE.
          FIRST = .TRUE.
          NMR = .FALSE.
          NOALT = .TRUE.
CHECK v.3.5.3<--
      ELSE
CHECK v.3.2-->
          ENSEMB = .FALSE.
          FIRST = .TRUE.
CHECK v.3.2<--
          NMR = .FALSE.
      ENDIF
CHECK v.3.1<--

CHECK v.3.2-->
C---- If processing a set of PDB files listed in the %filelist specified
C     by the user, open the file-list
      fistat = ok
      IF (ENSEMB) THEN
          finam = PDBFIL
          finuse = .true.
          namlen = INDEX(finam,space)
          do 180 i = 1, namlen - 1
              finam(i:i) = finam(i+1:i+1)
  180     continue
          finam(namlen:namlen) = space
          open (namfi, file=finam, status='old',
CVAX     -        READONLY,
     -        iostat=iocode) 
          if (iocode .ne. 0) fistat = opnera
          if (fistat.eq.opnera) GO TO 902
      ENDIF
CHECK v.3.2<--

CHECK v.3.2-->
C---- If processing a list of structures, get the name of the next one
C     from the file list
 200  CONTINUE
      IF (ENSEMB .AND. fistat .eq. ok) then
          read (namfi, 1020, iostat=iocode) brknam
          if (iocode .lt. 0) then
             close (namfi) 
             finuse = .false.
             fistat = endnfi
          else
             namlen = index(brknam,space) - 1
             write (*, 1030) brknam(1:namlen)
          end if
          PDBFIL = brknam
          NMR = .FALSE.
      ENDIF
CHECK v.3.2<--

C---- Peel off directory path and extension
      CALL GETNAM(PDBFIL,ISTART,IEND,IERROR)
      IF (IERROR) GO TO 990

C---- Form names of other files that will be required in default directory
      ILEN = IEND - ISTART + 1
      BRCODE = ' '
      brknam = PDBFIL
      NAMLEN = 78
      OUTNAM = PDBFIL(ISTART:IEND) // '.new'
CHECK v.3.2-->
      outlen = ILEN + 4
CHECK v.3.2<--
      ALTNM = PDBFIL(ISTART:IEND) // '.alt'
      LOGFIL = PDBFIL(ISTART:IEND) // '.clnlog'
CHECK v.2.0<--

c
c     write (*, *) 'Do you require Standard or Full output (<S> or F)'
c     read (*, '(A1)') outsty
c     if ((outsty .eq. 'F') .or. (outsty .eq. 'f')) then
c        outsty = 'F'
c     else
c        outsty = 'S'
c      end if
c
CHECK v.3.2-->
C      fistat = ok
CHECK v.3.2<--

CHECK v.1.0-->
C      call getstr ('Enter output directory', 
C     +             '/home/bsm/naylor/data/', outdir, outln, ifail)
C      if (ifail .ne. 0)  then
C         write(6,*)'Error in GETSTR: Ifail = ', ifail
C         stop
C      end if
CHECK v.1.0<--

      outsty = 'F'

CHECK v.1.0-->
C      open (protfi, file=pronam, status='unknown', iostat=iocode) 
C      if (iocode .ne. 0) fistat = opnerd
C      open (chnfi,  file=chnnam, status='unknown', iostat=iocode) 
C      if (iocode .ne. 0) fistat = opnere
CHECK v.1.0<--

c
c     Get the name of the next Brookhaven file or file of names from the user,
c     finish is set when they've had enough.
c

CHECK v.1.0-->
C  100 continue
C      write (*, 1010) indirn
C      read  (*, 1020) finam
C      namlen = index(finam,space) - 1
C      if (namlen .le. 0) then
C         fistat = finish
C      else
C         if (finam(1:1) .eq. indirn) then
Cc           open the list file if necessary
C            finuse = .true.
C            do 180 i = 1, namlen - 1
C            finam(i:i) = finam(i+1:i+1)
C  180       continue
C            finam(namlen:namlen) = space
C            open (namfi, file=finam, status='old', iostat=iocode) 
C            if (iocode .ne. 0) fistat = opnera
C         end if
C      end if
CHECK v.1.0<--

c
c     Get the next Brookhaven file name, either from the list file or from 
c     what the user just typed in, and load it into BRKNAM.
c

CHECK v.1.0-->
C  200 continue
C      if (fistat .eq. ok) then
C         if (finuse) then
C            read (namfi, 1020, iostat=iocode) brknam
C            if (iocode .lt. 0) then
C               close (namfi) 
C               finuse = .false.
C               fistat = endnfi
C            else
C               namlen = index(brknam,space) - 1
C               write (*, 1030) brknam(1:namlen)
C            end if
C         else
C            brknam = finam
C         end if
C      end if
CHECK v.1.0<--

c
c     Open the Brookhaven file if everything ok
c
      if (fistat .eq. ok) then
CHECK v.1.0-->
C         open (brkfi, file=brknam, status='old', iostat=iocode) 
CHECK v.1.0<--
         open (brkfi, file=brknam, status='old', iostat=iocode
CVAX     -        , READONLY
     -        ) 
         if (iocode .ne. 0) fistat = opnerb
      end if
c
c     Parse the file name to strip off nodenames, logical names, directory 
c     names and file extensions. Then open the corresponding ".NEW" and ".ALT"
c     output files. NB: This code seems still to be VAX specific for files 
c     NOT in the user's current working directory.
c
      if (fistat .eq. ok) then

CHECK v.1.0-->
C         nodepl = index(brknam,'::')
C         if (nodepl .ne. 0) then
C            logpl = (nodepl+1) + index(brknam(nodepl+2:),':')
C         else
C            logpl = index(brknam,':')
C         end if
C         dirpl = index(brknam,']')
C         corepl = max(logpl,dirpl) + 1
C         corepl = lastsl(brknam) + 1
C         corend = ( (corepl-1) + index(brknam(corepl:),'.') ) - 1
C         if (corend .eq. (corepl-2)) corend = namlen
C         if (corend .lt. corepl) then
C            write (6,*) brknam(1:namlen)
C            fistat = opnerc
C         end if
C         brklen = corend - corepl +1
C         outnam = outdir(1:outln)//brknam(corepl:corend)// '.new'
C         altnm  = outdir(1:outln)//brknam(corepl:corend)// '.alt'
C         logfil = outdir(1:outln)//brknam(corepl:corend)// '.clnlog'
C         outlen = outln + brklen + 4
CHECK v.1.0<--

CHECK v.2.1.3-->
C         open (logfi, file=logfil, status='unknown', iostat=iocode) 
C         if (iocode .ne. 0) fistat = opnerf
CHECK v.2.1.3<--
CHECK v.1.0-->
C         open (outfi, file=outnam, status='unknown', iostat=iocode) 

C----    Open the .new output file
         open (outfi, file=outnam, status='unknown', iostat=iocode
CVAX     -       , CARRIAGECONTROL='LIST'
     -   ) 
CHECK v.1.0<--
         if (iocode .ne. 0) fistat = opnerg
      end if

CHECK v.3.1-->
C---- If this is an NMR structure, then will loop back to here until
C     all models have been individually processed
 100  CONTINUE

CHECK v.3.2-->
C      IF (NMR) THEN
      if (fistat.eq.ok) then
CHECK v.3.2<--

C----     If this is not the first model encountered, write out the
C         ENDMDL record for the previous one
CHECK v.3.2-->
C          IF (.NOT.FIRST) THEN
          IF (NMR .AND. .NOT.FIRST) THEN
CHECK v.3.2<--
              WRITE(OUTFI,110) ENDMDL
 110          FORMAT(A)
          ENDIF
          FIRST = .FALSE.

C----     Loop until first ATOM record encountered, or a MODEL record
C         is found (indicating an NMR structure) or end-of-file encountered
CHECK v.3.2-->
C          CALL GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR)
CHECK v.3.5-->
C          CALL GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR)
          CALL GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR,
     -        BRCODE,AAPLN,AAIN,INSEQR,MAXTYP,SQRSCT,INSLEN,AUTHFI,
     -        ASTRUC)
CHECK v.3.5<--
CHECK v.3.2<--

C----     Check whether end of file or data error encountered
          IF (ENDFIL) GO TO 500
          IF (IERROR) GO TO 990
CHECK v.3.2-->
          IF (NMR) THEN
CHECK v.3.2<--
              NMODEL = NMODEL + 1
              PRINT*
              PRINT*, '   Processing NMR model', IMODEL
          ENDIF
CHECK v.3.1<--
CHECK v.3.2-->
      endif
CHECK v.3.2<--

c
c     If we have survived this far, we have a file to process otherwise
c     provide error messages.
c
      if (fistat .eq. ok) then
c        process file
c        open (errfi, status='SCRATCH') 
         readok = .true.
         caonly = .true.
         astruc = .false.
         altern = .false.
         nline = 0
         ave    = -9.99
         sd     = -9.99
         do 345 i = 1,sideaa
         totswp(i) = 0
         restot(i) = 0
  345    continue

C----    Read in the ATOM and HETATM records from the PDB file
         call rdbrk (mcaaat, atname, atomin, aacod3, hetcod, seqcod,
     +               seqbcd, scaaat, scatin, altcod, aastd,  atpres, 
CHECK v.3.2-->
C     +               thmprs, caonly, astruc, headln, hedlen, readok,
     +               thmprs, caonly, astruc, hedlen, readok,
CHECK v.3.2<--
     +               dudcod, modcod, stdmod, nstdmd, modpl,  seqlen, 
     +               nscats, mcser,  scser,  mcocc,  mcthrm, scocc, 
     +               scthrm, scname, altnam, naltmc, naltsc, altmci, 
     +               altsci, altmcr, altscr, altmca, altsca, lab, 
CHECK v.3.5-->
C     +               extdat, extres, nline,  tphnam, olenam, pphnam, 
C     +               nrec,   mcalt,  scalt,  altmcc, altscc, pyrnam, 
CCHECK v.3.1-->
CC     +               restot, altdat, altlin, altres, batnam, brcode)
C     +               restot, altdat, altlin, altres, batnam, brcode,
C     -               NMR)
CCHECK v.3.1<--
     +               extdat, extres, nline,  nrec,   mcalt,  scalt,
     +               altmcc, altscc, restot, altdat, altlin, altres,
     +               brcode, altpos, MCNAME, NMR,    FIRST,  inseqr,
     -               sqrsct)
CHECK v.3.5<--
         if (naltmc .gt. 0 .or. naltsc .gt. 0 .or. altlin .gt. 0) 
     +          altern = .true.
c        open file for alternates if there are any
         if (altern) open (altfi, file=altnm, status='unknown', 
     +                     iostat=iocode) 
         if (iocode .ne. 0) fistat = opnerh
         if (readok .and. (.not. caonly)) then
c           look for breaks across the peptide bonds
            call chnbrk (mcaaat, atomin, seqbcd, brksym, chnsz,  chnct, 
     +                   chnend, caonly, seqlen, chstop, chstrt, chnlt)
c           find main chain torsion angles
CHECK v.3.2-->
C            call mkangl (mcaaat, mcangs, atomin, chnsz,  chnct,  caonly, 
C     +                   seqlen, scaaat, scatin, seqcod, aacod3, hetcod,
C     +                   sx2,    sx,     ntot,   hist,   hismin, ave, 
C     +                   sd)

C----       Calculate all the main-chain torsion angles
            call mkangl (mcaaat, mcangs, atomin, chnsz,  chnct,  caonly, 
     +                   scaaat, scatin, seqcod, aacod3, hetcod,
     +                   sx2,    sx,     ntot,   hist,   hismin, ave, 
     +                   sd)
CHECK v.3.2<--

CHECK v.1.0-->
C            write (protfi, 209) brcode, chnct, chnlt, ave, sd
C            do 220 i = 1, chnct
C            write (chnfi, 210) brcode, chnsz(i)
C  220       continue
CHECK v.1.0<--

c           find sidechain chi angles
            if (outsty .eq. 'F') call mksang (mcaaat, atomin, scaaat, 
     +                                        scangs, scatin, seqcod, 
     +                                        aacod3, hetcod, seqlen, 
     +                                        nscats, chstop, chstrt, 
     +                                        seqbcd, mcser,  scser, 
     +                                        mcocc,  mcthrm, scocc, 
CHECK v.3.2-->
C     +                                        scthrm, scname, altnam, 
     +                                        scthrm, scname,
CHECK v.3.2<--
CHECK v.3.5-->
C     +                                        naltmc, naltsc, altmci, 
C     +                                        altsci, altmcr, altscr, 
C     +                                        altmca, altsca, lab, 
C     +                                        extdat, extres, nline, 
C     +                                        tphnam, olenam, pphnam, 
C     +                                        atname, mcangs, nrec, 
C     +                                        mcalt,  scalt,  altmcc,  
C     +                                        altscc, pyrnam, totswp, 
C     +                                        altdat, altlin, altres, 
CCHECK v.2.1.3-->
CC     +                                        batnam)
C     +                                        batnam, NATALT, NATZER,
CCHECK v.3.2-->
CC     +                                        NHEALT, NHEZER)
C     +                                        NHEALT, NHEZER,NMR)
CCHECK v.3.2<--
CCHECK v.2.1.3<--
     -                                        natsdc, extprc, sswap,
     +                                        naltmc, naltsc, altmci, 
     +                                        altsci, altmcr, altscr, 
     +                                        altmca, altsca, lab, 
     +                                        extdat, extres, nline, 
     +                                        atname, mcangs, nrec, 
     +                                        mcalt,  scalt,  altmcc,  
     +                                        altscc, totswp, 
     +                                        altdat, altlin, altres, 
     +                                        MCNAME, NATALT, NATZER,
CHECK v.3.5.3-->
C     +                                        NHEALT, NHEZER, NMR)
     +                                        NHEALT, NHEZER, NMR, 
     +                                        NOALT)
CHECK v.3.5.3<--
CHECK v.3.5<--
CHECK v.3.1-->
C            close (outfi) 
            IF (.NOT.NMR) close (outfi) 
CHECK v.3.1<--
            if (altern) close (altfi) 
         else

C----       CA-only file, write out all the atom records
            if (readok .and. caonly .and. outsty .eq. 'F') then
               if (fistat .ne. ok) then
                  write (*, 1040) outnam(1:outlen)
c                 close (errfi) 
               else

C----             Compute any chain breaks
                  call chnbrk (mcaaat, atomin, seqbcd, brksym, chnsz, 
     +                         chnct,  chnend, caonly, seqlen, chstop, 
     +                         chstrt, chnlt)

C----             Write atom records to output .new file
                  call calout (mcaaat, atomin, seqcod, aacod3, hetcod, 
     +                         seqlen, chstop, chstrt, seqbcd, mcser, 
     +                         mcocc,  mcthrm, extdat, nline,  extres, 
CHECK v.3.5-->
C     +                         nrec)
     +                         nrec,   NMR)
CHECK v.3.5<--

CHECK v.1.0-->
C                  write (protfi, 211) brcode, chnct, chnlt, ave, 
C     +                                sd, 'CA'
C                  do 230 i = 1, chnct
C                  write (chnfi, 210) brcode, chnsz(i)
C  230             continue
CHECK v.1.0<--

CHECK v.3.1-->
C                  close (outfi) 
                  IF (.NOT.NMR) close (outfi) 
CHECK v.3.1<--
                  if (altern) close (altfi) 
               end if
c              close (errfi) 
            else
               close (outfi, status='delete') 
               if (altern) close (altfi, status='delete') 
c              close (errfi) 
            end if
            if (astruc) close (authfi) 
CHECK v.3.5-->
C            if (caonly .and. readok) write (*, 1100) 
            if (caonly .and. readok .and. seqlen.gt.0)
     -          write (*, 1100) 
CHECK v.3.5<--
            if (.not. readok) then
               write (*, 1110) 
            end if
         end if
CHECK v.2.1.3-->
C         write (logfi, 1080) 
CHECK v.3.4-->
C         write (*, 1080) 
         IF (.NOT.NMR) THEN
             write (*, 1080)
         ELSE
             write (*, 1082)
         ENDIF
CHECK v.3.4<--
CHECK v.2.1.3<--
         do 300 i = 1, sideaa
         if (i .le. abetsz) then
            wkres = aacod3((i*3)-2:i*3)
         else
            wkres = hetcod(((i-abetsz)*3)-2:(i-abetsz)*3)
         end if

CHECK v.1.0-->
C         write (logfi, 1090) wkres, totswp(i), restot(i)
         IF (RESTOT(I).GT.0) THEN
CHECK v.2.1.3-->
C             write (logfi, 1090) wkres, totswp(i), restot(i)
             write (*, 1090) wkres, totswp(i), restot(i)
CHECK v.2.1.3<--
         ENDIF
CHECK v.1.0<--

  300    continue
CHECK v.2.1.3-->
C         close (logfi)
CHECK v.2.1.3<--
      else if (fistat .eq. opnera) then
CHECK v.3.2-->
C         write (*, 1050) finam(1:namlen)
CHECK v.3.2<--
         finuse = .false.
      else if (fistat .eq. opnerb) then
         write (*, 1060) brknam(1:namlen)
      else if (fistat .eq. opnerc) then
         write (*, 1070) outnam(1:outlen)
      else if (fistat .eq. opnerd) then
         write (*, 7001) pronam(1:strlen(pronam))
      else if (fistat .eq. opnere) then
         write (*, 7001) chnnam(1:strlen(chnnam))
      else if (fistat .eq. opnerf) then
         write (*, 7001) logfil(1:strlen(logfil))
      else if (fistat .eq. opnerg) then
         write (*, 7001) outnam(1:strlen(outnam))
      else if (fistat .eq. opnerh) then
         write (*, 7001) altnm(1:strlen(altnm))
      end if
c
c     if they want more we'll cooperate
c

CHECK v.1.0-->
C      if (fistat .ne. finish) then
C         fistat = ok
C         if (finuse) then
C            goto 200
C         else
C            goto 100
C         end if
C      end if
C      close (protfi) 
C      close (chnfi) 
CHECK v.1.0<--

c
c     compute total standard deviation and average dihedral angle
c     call tsdev(sx2,sx,ntot)
c     write (*,*)' distribution of dihedral angles'
c     do 400 i=1,17
c     write (*,*) hismin(i),'-',hismin(i)+10, hist(i)
c 400 continue
c

CHECK v.2.1.3-->
C---- Print numbers of alternate and zero occupancies written out
      IF (NATZER.GT.0) THEN
          PRINT*, ' * Number of zero-occupancy atoms (labelled',
     -        ' ATZERO)   ', NATZER
      ENDIF
      IF (NHEZER.GT.0) THEN
          PRINT*, ' * Number of zero-occupancy hetatoms (labelled',
     -        ' HEZERO)', NHEZER
      ENDIF
      IF (NATALT.GT.0) THEN
          PRINT*, ' * Number of alternate atoms (labelled',
     -        ' ATALT)      ', NATALT
      ENDIF
      IF (NHEALT.GT.0) THEN
          PRINT*, ' * Number of alternate hetatoms (labelled',
     -        ' HEALT)   ', NHEALT
      ENDIF
CHECK v.2.1.3<--

CHECK v.3.1-->
C---- If this is an NMR structure, then loop back for next model
      IF (NMR) GO TO 100

C---- If this is an NMR structure, show how many models have been processed
 500  CONTINUE
      IF (NMR) THEN
          PRINT*
          PRINT*, '* NMR ensemble comprises', NMODEL, ' model ',
     -        'structures'
      ENDIF

CHECK v.3.2-->
      if (fistat .ne. finish) then
         fistat = ok
         if (finuse) then
            goto 200
         end if
      end if
      close (protfi) 
      close (chnfi) 
CHECK v.3.2<--

      GO TO 990

CHECK v.3.2-->
C---- Error opening mplot.in file
C 900  CONTINUE
C      PRINT*, '*** ERROR opening mplot.in file. Program aborted'
CHECK v.3.2<--
CHECK v.3.1<--

CHECK v.3.2-->
 902  CONTINUE
      PRINT*, '*** ERROR opening file:', finam(1:namlen)
CHECK v.3.2<--

CHECK v.2.0-->
990   CONTINUE
CHECK v.2.0<--
CHECK v.2.3-->
CHECK v.3.5-->
C      PRINT*, '* Program completed'
      IF (IFAIL) THEN
          PRINT*, '* Program terminated with error'
      ELSE
          PRINT*, '* Program completed'
      ENDIF
CHECK v.3.5<--
CHECK v.2.3<--
      end
c
CHECK v.3.5-->
C**************************************************************************
C
C  SUBROUTINE GETDEF  -  Read in the Residue Definitions from the
C                        resdefs.dat file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETDEF(AACOD3,ABETSZ,HETCOD,MXHTYP,NSCATS,SCNAME,
     -    MCNAME,SIDEAA,NATSDC,EXTPRC,ALTPOS,SSWAP,IFAIL)

      INTEGER       ABETSZ, MXHTYP, SIDEAA

      CHARACTER*(*) AACOD3, HETCOD

      CHARACTER*3   RESNAM
      CHARACTER*6   IDENT
      CHARACTER*44  SCNAME(SIDEAA), MCNAME(SIDEAA)
      CHARACTER*80  IREC
      INTEGER       ALTPOS(SIDEAA), EXTPRC(SIDEAA), HRES, IPOS, IRES,
     -              ISIDE, ISWAP, LINE, NATSDC(SIDEAA), NSCATS(SIDEAA),
     -              SSWAP(2,SIDEAA)
      LOGICAL       IFAIL

C---- Initialise variables
      DO 50, ISIDE = 1, SIDEAA
          ALTPOS(ISIDE) = 0
          EXTPRC(ISIDE) = 0
          MCNAME(ISIDE) = ' '
          NATSDC(ISIDE) = 0
          NSCATS(ISIDE) = -1
          SCNAME(ISIDE) = ' '
          SSWAP(1,ISIDE) = 0
          SSWAP(2,ISIDE) = 0
 50   CONTINUE
      AACOD3 = ' '
      HETCOD = ' '
      IFAIL = .FALSE.
      HRES = 0
      IRES = 0
      ISIDE = 0
      LINE = 0

C---- Open the data file
      OPEN(UNIT=1, FILE='resdefs.dat', STATUS='OLD', FORM='FORMATTED',
     -    ACCESS='SEQUENTIAL',
CVAX     -    READONLY,
     -    ERR=900)

C---- Read through the Residue Definitions file, storing the data
 100  CONTINUE

C----     Read in the next record
          READ(1,120,ERR=902,END=500) IREC
 120      FORMAT(A)
          LINE = LINE + 1

C----     Get the record identifier
          IDENT = IREC(1:6)

C----     Check for residue-start
          IF (IDENT.EQ.'RESDUE') THEN

C----         Extract the residue name
              RESNAM = IREC(8:)
              IRES = IRES + 1
              IF (IRES.GT.ABETSZ) GO TO 904
              ISIDE = ISIDE + 1
              IF (ISIDE.GT.SIDEAA) GO TO 908

C----         Store the name
              IPOS = 3 * (IRES - 1) + 1
              AACOD3(IPOS:IPOS + 2) = RESNAM

C----         Get the atom parameters
              READ(IREC,140,ERR=910) NATSDC(ISIDE), EXTPRC(ISIDE),
     -            ALTPOS(ISIDE)
 140          FORMAT(10X,3I2)

C----     Check for HETATM residue type
          ELSE IF (IDENT.EQ.'HETNAM') THEN

C----         Extract the residue name
              RESNAM = IREC(8:)
              HRES = HRES + 1
              IF (HRES.GT.MXHTYP) GO TO 906
              ISIDE = ISIDE + 1
              IF (ISIDE.GT.SIDEAA) GO TO 908

C----         Store the name
              IPOS = 3 * (HRES - 1) + 1
              HETCOD(IPOS:IPOS + 2) = RESNAM

C----         Get the atom parameters
              READ(IREC,140,ERR=910) NATSDC(ISIDE), EXTPRC(ISIDE)

C----     If this record defines the non-standard mainchain atoms, then store
          ELSE IF (IDENT.EQ.'MAINAT') THEN

C----         Store the mainchain atom names
              MCNAME(ISIDE) = IREC(8:)

C----     If this record defines the sidechain atoms, then store
          ELSE IF (IDENT.EQ.'SIDEAT') THEN

C----         Initialise the count of sidechain atoms
              NSCATS(ISIDE) = 0
              IREC = IREC(8:)

C----         Remove the CB, if present
              IF (IREC(1:4).EQ.' CB ') THEN
                  NSCATS(ISIDE) = NSCATS(ISIDE) + 1
                  IREC = IREC(5:)
              ENDIF

C----         Store the sidechain atom names
              SCNAME(ISIDE) = IREC

C----         Count the number of atom names
              DO 200, IPOS = 1, 44 - 3, 4
                  IF (IREC(IPOS:IPOS+3).NE.'    ') THEN
                      NSCATS(ISIDE) = NSCATS(ISIDE) + 1
                  ENDIF
 200          CONTINUE

C----         Check for special XXX residues
              IF (RESNAM.EQ.'XXX' .OR. RESNAM.EQ.'UNK') THEN
                  NSCATS(ISIDE) = -1
              ENDIF

C----     If this is a record indicating which atoms may need to be
C         swapped, then read in the swap positions
          ELSE IF (IDENT.EQ.'ATSWAP') THEN

C----         Get the atom-swap markers
              IREC = IREC(8:)
              ISWAP = 0

C----         Loop through the positions to pick up the swap atoms
              DO 300, IPOS = 1, 44 - 3, 4
                  IF (IREC(IPOS:IPOS+3).EQ.'^^^^') THEN
                      ISWAP = ISWAP + 1
                      IF (ISWAP.LE.2) THEN
                          SSWAP(ISWAP,ISIDE) = (IPOS - 1) / 4 + 1
                      ENDIF
                  ENDIF
 300          CONTINUE

          ENDIF

C---- Loop back for next record
      GO TO 100

C---- End of file reached
 500  CONTINUE
      CLOSE(1)

      GO TO 999

 900  CONTINUE
      PRINT*, '*** ERROR. Unable to open parameter file, resdefs.dat'
      GO TO 990

 902  CONTINUE
      PRINT*, '*** ERROR. Read error on file, resdefs.dat, at line',
     -    LINE + 1
      GO TO 990

 904  CONTINUE
      PRINT*, '*** ERROR. Maximum number of definable residues in file',
     -     ' resdefs.dat exceeded at residue [', RESNAM, ']'
      PRINT*, '***        Number allowed, ABETSZ, =', ABETSZ
      GO TO 990

 906  CONTINUE
      PRINT*, '*** ERROR. Maximum number of definable HETATM residues ',
     -     'in file resdefs.dat exceeded at residue [', RESNAM, ']'
      PRINT*, '***        Number allowed, MXHTYP, =', MXHTYP
      GO TO 990

 908  CONTINUE
      PRINT*, '*** ERROR. Maximum number of definable residue types in',
     -     'file resdefs.dat exceeded at residue [', RESNAM, ']'
      PRINT*, '***        Number allowed, SIDEAA, =', SIDEAA
      GO TO 990

 910  CONTINUE
      PRINT*, '*** ERROR reading atom numbers in file, resdefs.dat, at',
     -    ' line', LINE + 1
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CHECK v.3.5<--
c*******************************************************************************
c
      subroutine rdbrk (mcaaat, atname, atomin, aacod3, hetcod, seqcod, 
     +                  seqbcd, scaaat, scatin, altcod, aastd,  atpres, 
CHECK v.3.2-->
C     +                  thmprs, caonly, astruc, headln, hedlen, readok, 
     +                  thmprs, caonly, astruc, hedlen, readok, 
CHECK v.3.2<--
     +                  dudcod, modcod, stdmod, nstdmd, modpl,  seqlen, 
     +                  nscats, mcser,  scser,  mcocc,  mcthrm, scocc, 
     +                  scthrm, scname, altnam, naltmc, naltsc, altmci, 
     +                  altsci, altmcr, altscr, altmca, altsca, lab, 
CHECK v.3.5-->
C     +                  extdat, extres, nline,  tphnam, olenam, pphnam, 
C     +                  nrec,   mcalt,  scalt,  altmcc, altscc, pyrnam, 
CCHECK v.3.1-->
CC     +                  restot, altdat, altlin, altres, batnam, brcode)
C     +                  restot, altdat, altlin, altres, batnam, brcode,
C     -                  NMR)
CCHECK v.3.1<--
     +                  extdat, extres, nline,  nrec,   mcalt,  scalt,
     +                  altmcc, altscc, restot, altdat, altlin, altres,
     +                  brcode, altpos, MCNAME, NMR,    FIRST,  inseqr,
     -                  sqrsct)
CHECK v.3.5<--
c
      include 'brkcln.par'
c
      integer   seqlen, seqcod(maxres), hedlen, mcser(mnchna,maxres), 
     +          scser(sdchna,maxres), naltmc, naltsc, 
     +          altsci(3,altmax), altmci(3,altmax), extres(maxlin), 
     +          nline, nrec, altres(maxalt), altlin
c 
      character aacod3*(abetsz*3), atname*((mnchna*4)-4), 
CHECK v.3.2-->
C     +          seqbcd(maxres)*6, headln*132, altcod(maxres)*1, 
     +          seqbcd(maxres)*6, altcod(maxres)*1, 
CHECK v.3.2<--
     +          aastd(maxres)*1, atpres(3,maxres)*1, thmprs(maxres)*1,
CHECK v.3.5-->
C     +          hetcod*63, modcod*(nmod*3), stdmod*(nmod+1), 
C     +          nstdmd*(nmod+1), modpl*(nmod), dudcod*105, 
C     +          lab(maxres)*1, extdat(maxlin)*66, 
C     +          tphnam*((mnchna*4)-4), olenam*((mnchna*4)-4), 
C     +          pphnam*((mnchna*4)-4), mcalt(mnchna,maxres)*1, 
C     +          scalt(sdchna,maxres)*1, altmcc(maxres)*1, 
C     +          altscc(maxres)*1, 
C     +          pyrnam*((mnchna*4)-4), altdat(maxalt)*66, 
C     +          batnam*((4*mnchna)-4)
     +          hetcod*(3*MXHTYP), modcod*(nmod*3), stdmod*(nmod+1), 
     +          nstdmd*(nmod+1), modpl*(nmod), dudcod*105, 
     +          lab(maxres)*1, extdat(maxlin)*66, 
     +          mcalt(mnchna,maxres)*1, MCNAME(SIDEAA)*44,
     +          scalt(sdchna,maxres)*1, altmcc(maxres)*1, 
     +          altscc(maxres)*1, altdat(maxalt)*66
CHECK v.3.5<--
c
      real      mcaaat(ncoord,mnchna,maxres), 
     +          scaaat(ncoord,sdchna,maxres), mcocc(mnchna,maxres), 
     +          mcthrm(mnchna,maxres), 
     +          scocc(sdchna,maxres), scthrm(sdchna,maxres), 
     +          altmcr(2,altmax), altscr(2,altmax), altmca(3,altmax),
     +          altsca(3,altmax)
c
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres), caonly, 
     +          astruc, readok
c
c     local variables
c
CHECK v.3.5-->
C      integer   aaconv, aapln, num
C      parameter (aapln = 13)
      integer   aaconv, num
CHECK v.3.5<--
c
      real      acoord(ncoord), oldocc, occup, therm, oldthm, extocc, 
     +          altocc
c
      logical   altset, xney
c
c     integer   nrecs, aathln, rderr, altype, nxaa, nxrec, skipaa, 
c    +          hedend
c
CHECK v.3.5-->
C      integer   atcnt, attype, i, j, aacnt, wkpl, clslen, cmplen, 
C     +          soulen, sctype, altpos(sideaa), respl, scatct, thmcnt, 
C     +          aacode, nscats(sideaa), sqrsct, restot(sideaa)
      integer   atcnt, attype, i, j, aacnt, wkpl,
     +          sctype, altpos(sideaa), respl, scatct, thmcnt, 
     +          aacode, nscats(sideaa), sqrsct, restot(sideaa)
CHECK v.3.5<--
c
CHECK v.3.5-->
C      character aain(maxres)*3, atmnam*4, resnam*3, altloc*1, 
C     +          brkrec*80, oldchn*1, seqnum*6, keywd*(keylen),
C     +          oldatm*4, setloc*1, class*40, compnd*50, source*50, 
C     +          brcode*4, cbet*4, scname(sideaa)*44, altnam(3)*44, 
C     +          contch*1, incode*4, wkstr*50, inseqr*150, setatm*4, 
      character atmnam*4, resnam*3, altloc*1, 
     +          brkrec*80, oldchn*1, seqnum*6, keywd*(keylen),
     +          oldatm*4, setloc*1,
     +          brcode*4, cbet*4, scname(sideaa)*44, altnam(3)*44, 
     +          setatm*4, 
CHECK v.3.5<--
     +          wkres*3, extatm*4, extnum*6, string*66, altatm*4, 
     +          altnum*6
CHECK v.3.2-->
      character rname*3
CHECK v.3.2<--
CHECK v.3.5-->
      CHARACTER*(*) inseqr
CHECK v.3.5<--
c
c     character chid*1, chnid(maxchn)*1, innum*6
c
CHECK v.3.5-->
C      integer   strlen, resfnd, kount
      integer   NREQRD, resfnd
CHECK v.3.5<--

CHECK v.3.1-->
      LOGICAL   NMR
CHECK v.3.1<--
CHECK v.3.5-->
      LOGICAL   FIRST
CHECK v.3.5<--
c
c     Routine to read the brookhaven data file and extract the sequence
c     and the atom names and coordinates. Checking is also done for
c     residue compatibility between sequence and atom. C-alpha only
c     files are flagged. The sequence residues are expected together,
c     then the secondary structure records (sheets together), and then
c     the atom records which are expected in residue order.
c
CHECK v.3.5-->
C      save altpos
C      data altpos / 7*0, 1, 5*0, 2, 2*0, 3, 33*0 /
CHECK v.3.5<--
      data cbet / ' CB ' /
c
CHECK v.3.5-->
CHECK v.3.4.3-->
CC  220 format (18x, 13(x,a3))
C  220 format (18x, 13(1x,a3))
CCHECK v.3.4.3<--
C  520 format (9x, a1, a40, 12x, a4)
C  620 format (9x, a1, a50)
C  720 format (9x, a1, a50)
CHECK v.3.5<--
 1080 format (21x, a1)
 1120 format (12x, a4, a1, a3, 1x, a6, 27x, 2f6.2)
 1121 format (12x, a4, 5x,         a6, 27x, f6.2)
 1320 format (6x, i5, 1x, a4, a1, 4x, a6, 3x, 3f8.3, 2f6.2)
CHECK v.2.1.3-->
C 1374 format('WARNING: This side chain atom (',a4,') not found for',
C     + ' this residue ',a3, 1x, i4)
CHECK v.3.0.1-->
C 1374 format('* WARNING: This side chain atom (',a4,') not found for',
C     + ' this residue ',a3, 1x, i4)
CHECK v.3.0.1<--
CHECK v.2.1.3<--
 6000 format (A)
CHECK v.2.1.3-->
C 6010 format (' Too many amino acids - input terminated.')
 6010 format (' **** Too many amino acids - input terminated.')
CHECK v.2.1.3<--
c

C---- Initialise variables
CHECK v.3.5-->
C      kount  = 0
CHECK v.3.5<--
      altlin = 0
      naltmc = 0
      naltsc = 0
      seqlen = 0
      hedlen = 0
CHECK v.3.5-->
C      inseqr = space
C      sqrsct = 0
      NREQRD = 0
CHECK v.3.5<--
      brcode = '    '
      astruc = .false.
      do 50 i = 1, mnchna
      do 49 j = 1, maxres
      mcocc(i,j) = 0.0
   49 continue
   50 continue
      do 52 i = 1, sdchna
      do 51 j = 1, maxres
      scocc(i,j) = 0.0
   51 continue
   52 continue
c

C---- Loop while reading through the PDB file
  100 continue

C----     Read in the next record from the PDB file
          read (brkfi, 6000, end=3000) brkrec
CHECK v.3.5-->
C          kount = kount + 1
CHECK v.3.5<--
          keywd = brkrec(1:keylen)

C----     Go to the appropriate routine, depending on the record type
CHECK v.3.5-->
C          if (keywd .ne. atmkey .and. 
          if (keywd .ne. atmkey .and. keywd .ne. modkey .and.
CHECK v.3.5<--
     +        keywd .ne. htakey) write (outfi, 6000) brkrec
CHECK v.3.5-->
C          if (keywd .eq. seqkey) goto 200
C          if (keywd .eq. hedkey) goto 500
C          if (keywd .eq. cmpkey) goto 600
C          if (keywd .eq. soukey) goto 700
CHECK v.3.5<--
          if (keywd .eq. atmkey) goto 900
          if (keywd .eq. htakey) goto 900
CHECK v.3.1-->
          IF (KEYWD .EQ. ENDMDL) GO TO 3000
CHECK v.3.1<--
          if (keywd .eq. modkey) then
CHECK v.3.5-->
C             readok = .false.
CHECK v.3.5<--
             seqlen = 0
CHECK v.3.1-->
CHECK v.3.5-->
C             PRINT*, '*** ERROR: MODEL record encountered unexpecte',
C     -           'dly in PDB file'
             PRINT*, '*** Warning. MODEL contains no atoms'
             WRITE(OUTFI,120)
 120         FORMAT('ENDMDL')
             write (outfi, 6000) brkrec
             FIRST = .TRUE.
CHECK v.3.5<--
CHECK v.3.1<--
             goto 3000
          end if
CHECK v.3.5-->
C          if (keywd .eq. hlxkey .or. keywd .eq. shtkey .or. 
C     +        keywd .eq. trnkey .or. keywd .eq. dsfkey) then
C             if (.not. astruc) then
C                open (authfi, status='SCRATCH') 
C                astruc = .true.
C             end if
C             write (authfi, 6000) brkrec
C          end if
CHECK v.3.5<--

C---- Loop back for next record
      goto 100
c
CHECK v.3.5-->
C  200 continue
C      read (brkrec, 220) (aain(i), i=1,aapln)
C      do 250 i = 1, aapln
C      if (aain(i) .ne. '   ') then
C         if (resfnd(inseqr,aain(i)) .eq. 0) then
C            sqrsct = sqrsct + 1
C            inseqr(((sqrsct-1)*3)+1:sqrsct*3) = aain(i)
C         end if
C      end if
C  250 continue
C      goto 100
Cc
Cc     set up header records
Cc
C  500 continue
C      read (brkrec, 520) contch, wkstr, incode
C      if (contch .eq. space) then
C         brcode = incode
C         class = wkstr
C         clslen = strlen(class)
C      end if
C      goto 100
Cc
C  600 continue
C      read (brkrec, 620) contch, wkstr
C      if (contch .eq. space) then
C         compnd = wkstr
C         cmplen = strlen(compnd)
C      end if
C      goto 100
Cc
C  700 continue
C      read (brkrec, 720) contch, wkstr
C      if (contch .eq. space) then
C         source = wkstr
C         soulen = strlen(source)
C      end if
C      goto 100
CHECK v.3.5<--

c
c     look for next atom records to process
c
  900 continue

 1000 continue
      if (brkrec(1:keylen) .ne. atmkey .and. 
     +    brkrec(1:keylen) .ne. htakey) then
         read (brkfi, '(A)', end=3000) brkrec
CHECK v.3.1-->
         IF (BRKREC(1:KEYLEN) .EQ. ENDMDL) GO TO 3000
CHECK v.3.1<--
         goto 1000
      end if

C---- First atom of a new residue

c
c     Arrival here implies we have just read in the record for the first atom 
c     of a particular amino acid. For each aa get the seqnum and resnam and 
c     check against the read in sequence. If ok initialise to no atoms present
c     and get the coords.
c

C---- Initialise variables
      aacnt = 0
      aastd(1) = space

C---- Pick up this atom's details
      read (brkrec, 1080) oldchn
 1100 continue
      keywd  = brkrec(1:keylen)
      aacode = 0
      altset = .false.
      setatm = '    '
      read (brkrec, 1120) oldatm, setloc, resnam, seqnum, oldocc, oldthm
CHECK v.2.1.1-->
      IF (brkrec(57:60).EQ.'    ') oldocc = 1.0
CHECK v.2.1.1<--

C---- Check whether this residue is present in the SEQRES details
      if (sqrsct .ne. 0 .and. resfnd(inseqr,resnam) .eq. 0) then
         if (brkrec(1:keylen) .eq. atmkey .or. 
     +       brkrec(1:keylen) .eq. htakey) then
            if (setloc .ne. space .and. nline .ne. 0) then
               string = extdat(nline)
               read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
               IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
               if (extatm .eq. oldatm .and. extnum .eq. seqnum) then
                  if (extocc .ge. oldocc) then
                     altlin = altlin + 1
                     altdat(altlin) = brkrec(1:66)
                  else
                     altlin = altlin + 1
                     altdat(altlin) = string
                     extdat(nline) = brkrec(1:66)
                  end if
                  altres(altlin) = aacnt
               else
CHECK v.3.5-->
C                  nline = nline + 1
C                  extdat(nline) = brkrec(1:66)
C                  extres(nline) = aacnt
                   CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -                 brkrec(1:66),NREQRD)
CHECK v.3.5<--
               end if
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         else
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         end if
         goto 1400
      end if
c
      if (resfnd(dudcod,resnam) .ne. 0) then
         if (brkrec(1:keylen) .eq. atmkey .or. 
     +       brkrec(1:keylen) .eq. htakey) then
            if (setloc .ne. space .and. nline .ne. 0) then
               string = extdat(nline)
               read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
               IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
               if (extatm .eq. oldatm .and. extnum .eq. seqnum) then
                  if (extocc .ge. oldocc) then
                     altlin = altlin + 1
                     altdat(altlin) = brkrec(1:66)
                  else
                     altlin = altlin + 1
                     altdat(altlin) = string
                     extdat(nline) = brkrec(1:66)
                  end if
                  altres(altlin) = aacnt
               else
CHECK v.3.5-->
C                  nline = nline + 1
C                  extdat(nline) = brkrec(1:66)
C                  extres(nline) = aacnt
                   CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -                 brkrec(1:66),NREQRD)
CHECK v.3.5<--
               end if
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         else
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         end if
         goto 1400
      end if
c
      respl = resfnd(modcod,resnam)
      if (respl .ne. 0) then

C----    Modified residue code
         if (brkrec(1:keylen) .eq. atmkey .or. 
     +       brkrec(1:keylen) .eq. htakey) then
            if (setloc .ne. space .and. nline .ne. 0) then
               string = extdat(nline)
               read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
               IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
               if (extatm .eq. oldatm .and. extnum .eq. seqnum) then
                  if (extocc .ge. oldocc) then
                     altlin = altlin + 1
                     altdat(altlin) = brkrec(1:66)
                  else
                     altlin = altlin + 1
                     altdat(altlin) = string
                     extdat(nline) = brkrec(1:66)
                  end if
                  altres(altlin) = aacnt
               else
CHECK v.3.5-->
C                  nline = nline + 1
C                  extdat(nline) = brkrec(1:66)
C                  extres(nline) = aacnt
                   CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -                 brkrec(1:66),NREQRD)
CHECK v.3.5<--
               end if
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         else
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         end if
         respl = ((respl - 1) / 3) + 1
         if (modpl(respl:respl) .eq. '+') then
            aastd(aacnt+1) = stdmod(respl:respl)
         else
            if (aastd(aacnt) .eq. 'S') then
               aastd(aacnt) = stdmod(respl:respl)
            else
               if (aastd(aacnt) .eq. 'O') then
                  aastd(aacnt) = nstdmd(respl:respl)
               else
                  if (llt(aastd(aacnt),'a')) then
                     aastd(aacnt) = stdmod(nmod+1:nmod+1)
                  else
                     aastd(aacnt) = nstdmd(nmod+1:nmod+1)
                  end if
               end if
            end if
         end if
         goto 1400
      end if
c
      aacode = aaconv(resnam,aacod3,hetcod)
      if (aacode .gt. 0) restot(aacode) = restot(aacode) + 1
      if (aacode .le. 0) then
         if (brkrec(1:keylen) .eq. atmkey .or. 
     +       brkrec(1:keylen) .eq. htakey) then
            if (setloc .ne. space .and. nline .ne. 0) then
               string = extdat(nline)
               read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
               IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
               if (extatm .eq. oldatm .and. extnum .eq. seqnum) then
                  if (extocc .ge. oldocc) then
                     altlin = altlin + 1
                     altdat(altlin) = brkrec(1:66)
                  else
                     altlin = altlin + 1
                     altdat(altlin) = string
                     extdat(nline) = brkrec(1:66)
                  end if
                  altres(altlin) = aacnt
               else
CHECK v.3.5-->
C                  nline = nline + 1
C                  extdat(nline) = brkrec(1:66)
C                  extres(nline) = aacnt
                   CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -                 brkrec(1:66),NREQRD)
CHECK v.3.5<--
               end if
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         else
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         end if
         goto 1400
      end if
c
      aacnt = aacnt + 1
      if (aastd(aacnt) .eq. space) then
         aastd(aacnt) = 'S'
         if (aacode .gt. stdaa) aastd(aacnt) = 'O'
      else
         if (aacode .gt. stdaa) then
            wkpl = index(stdmod,aastd(aacnt))
            aastd(aacnt) = nstdmd(wkpl:wkpl)
         end if
      end if
c
      seqcod(aacnt) = aacode
      seqbcd(aacnt) = seqnum
      atcnt  = 0
      scatct = 0
      thmcnt = 0
      thmprs(aacnt) = space
      do 1200 i = 1, mnchna
      atomin(i,aacnt) = .false.
 1200 continue
      do 1220 i = 1, sdchna
      scatin(i,aacnt) = .false.
 1220 continue
      do 1240 i = 1, 3
      atpres(i,aacnt) = space
 1240 continue
 1300 continue
c
c     process cordinates find atom type and if main chain set in coords
c
      read (brkrec, 1320) num, atmnam, altloc, altnum, 
     +                   (acoord(i), i=1,3), occup, therm
CHECK v.2.1.1-->
      IF (brkrec(57:60).EQ.'    ') occup = 1.0
CHECK v.2.1.1<--
      if (atmnam      .eq. ' OXT' .or. atmnam      .eq. ' NXT' .or. 
CHECK v.3.2-->
C     +    atmnam(2:2) .eq. 'H'    .or. atmnam(2:2) .eq. 'D') then
     +    atmnam(2:2) .eq. 'H'    .or. atmnam(2:2) .eq. 'D' .or.
CHECK v.3.5-->
     +    atmnam(1:1) .eq. 'H'    .or.
CHECK v.3.5<--
     +    atmnam(2:2) .eq. 'Q') then
CHECK v.3.2<--
         if (altloc .ne. space .and. nline .ne. 0) then
            string = extdat(nline)
            read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
            IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
            if (extatm .eq. atmnam .and. extnum .eq. altnum) then
               if (extocc .ge. occup) then
                  altlin = altlin + 1
                  altdat(altlin) = brkrec(1:66)
               else
                  altlin = altlin + 1
                  altdat(altlin) = string
                  extdat(nline) = brkrec(1:66)
               end if
               altres(altlin) = aacnt
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         else
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         end if
         goto 1400
      end if
c
      attype = index(atname,atmnam)
      if (attype .eq. 0) then
CHECK v.3.5-->
C         if (aacode .gt. 44 .and. aacode .le. 50) then
C            attype = index(batnam,atmnam)
C         else if (aacode .eq. 42) then
C            attype = index(tphnam,atmnam)
C         else if (aacode .eq. 41) then
C            attype = index(olenam,atmnam)
C         else if (aacode .eq. 36) then
C            attype = index(pyrnam,atmnam)
C         else
C            if (aacode .eq. 39) then
C               attype = index(pphnam,atmnam)
C            end if
C         end if

C----     Check for unusual mainchain atoms as identified in the MAINAT
C         record of the resdefs.dat file
          IF (MCNAME(AACODE).NE.' ') THEN
              attype = index(MCNAME(AACODE),atmnam)
          ENDIF
CHECK v.3.5<--
      end if
c
      if (altloc .ne. space) then
         if (setatm .eq. '    ') then
            setatm = atmnam
            setloc = altloc
            oldocc = occup
            oldthm = therm
         end if
         if (.not. altset) then
            if (setatm .eq. atmnam) then
               if (occup .gt. oldocc) then
                  setloc = altloc
                  oldocc = occup
                  if (xney(oldthm,0.0)) thmcnt = thmcnt - 1
                  oldthm = therm
                  if (attype .ne. 0) then
                     atcnt = atcnt - 1
                  else
                     scatct = scatct - 1
                  end if
               end if
            else
               altset = .true.
            end if
         end if
      end if
c
      if (xney(therm,0.0) .and. (altloc .eq. setloc .or. 
     +                           altloc .eq. space)) thmcnt = thmcnt + 1
      if (attype .ne. 0 .and. (altloc .eq. setloc .or. 
     +                         altloc .eq. space)) then 
         attype = (attype / 4) + 1
         if (brkrec(1:keylen) .eq. atmkey) then
            lab(aacnt) = 'A'
         else
            lab(aacnt) = 'H'
         end if
c        next section applies to cases where first atom has higher occupancy
         if (xney(mcocc(attype,aacnt),0.0)) then
            if (mcocc(attype,aacnt) .lt. occup) then
               naltmc = naltmc + 1
               altmci(1,naltmc) = aacnt
               altmci(2,naltmc) = attype
               altmci(3,naltmc) = mcser(attype,aacnt)
               altmcr(1,naltmc) = mcocc(attype,aacnt)
               altmcr(2,naltmc) = mcthrm(attype,aacnt)
               altmcc(naltmc)   = mcalt(attype,aacnt)
               do 1349 i = 1, ncoord
               altmca(i,naltmc) = mcaaat(i,attype,aacnt)
 1349          continue
            end if
         end if
         mcser(attype,aacnt)  = num
         mcocc(attype,aacnt)  = occup
         mcthrm(attype,aacnt) = therm
         mcalt(attype,aacnt)  = altloc
         do 1350 i = 1, ncoord
         mcaaat(i,attype,aacnt) = acoord(i)
 1350    continue
         atomin(attype,aacnt) = .true.
         if (attype .ne. calph) caonly = .false.
         atcnt = atcnt + 1
      else
         if (attype .ne. 0) then
c           alternate main chain locations
            attype = (attype / 4) + 1
            naltmc = naltmc + 1
            altmci(1,naltmc) = aacnt
            altmci(2,naltmc) = attype
            altmci(3,naltmc) = num
            altmcr(1,naltmc) = occup
            altmcr(2,naltmc) = therm
            altmcc(naltmc) = altloc
            do 1351 i = 1, ncoord
            altmca(i,naltmc) = acoord(i)
 1351       continue
         end if
      end if
c
      if (attype .eq. 0 .and. (altloc .eq. setloc .or. 
     +                         altloc .eq. space)) then
         scatct = scatct + 1
         if (aacode .le. sideaa) then
            if (aacode .ne. 27) then

CHECK v.3.2-->
C----          Check for special cases of non-standard atom-labelling
C              (eg if file created by X-PLOR)
CHECK v.3.5.3-->
C               rname = aacod3(3*(aacode-1)+1:3*aacode)
               IF (aacode.LE.abetsz) THEN
                   rname = aacod3(3*(aacode-1)+1:3*aacode)
               ELSE
                   rname
     -                 = hetcod(3*(aacode-abetsz-1)+1:3*(aacode-abetsz))
               ENDIF
CHECK v.3.5.3<--
               if (rname.eq.'ILE' .and. atmnam.eq.' CD ')
     -             atmnam = ' CD1'
               if (rname.eq.'GLN' .and. atmnam.eq.' OE ')
     -             atmnam = ' OE1'
               if (rname.eq.'GLU' .and. atmnam.eq.' OE ')
     -             atmnam = ' OE1'
CHECK v.3.2<--
               if (atmnam .eq. cbet) then
                  sctype = 1
               else
                  sctype = index(scname(aacode),atmnam)
                  if (sctype .eq. 0 .and. altpos(aacode) .ne. 0)
     +                sctype = index(altnam(altpos(aacode)),atmnam)
                  if (sctype .ne. 0) sctype = (sctype / 4) + 2
               end if
c              next section applies to cases where first atom has higher 
c              occupancy
               if (sctype .ne. 0) then
                  if (xney(scocc(sctype,aacnt),0.0)) then
                     if (scocc(sctype,aacnt) .lt. occup) then
                        naltsc = naltsc + 1
                        altsci(1,naltsc) = aacnt
                        altsci(2,naltsc) = sctype
                        altsci(3,naltsc) = scser(sctype,aacnt)
                        altscr(1,naltsc) = scocc(sctype,aacnt)
                        altscr(2,naltsc) = scthrm(sctype,aacnt)
                        altscc(naltsc)   = scalt(sctype,aacnt)
                        do 1353 i = 1, ncoord
                        altsca(i,naltsc) = scaaat(i,sctype,aacnt)
1353                    continue
                     end if
                  end if
                  scatin(sctype,aacnt) = .true.
                  scocc(sctype,aacnt)  = occup
                  scthrm(sctype,aacnt) = therm
                  scser(sctype,aacnt)  = num
                  scalt(sctype,aacnt)  = altloc
                  do 1375 i = 1, ncoord
                  scaaat(i,sctype,aacnt) = acoord(i)
 1375             continue
               else
CHECK v.3.4.3-->
                  PRINT*, '** Unrecognized atom type [', atmnam, ']',
     -                ' in residue ', resnam, ' ', seqbcd(aacnt)
CHECK v.3.4.3<--
                  if (seqcod(aacnt) .le. abetsz) then
                     wkres = aacod3((seqcod(aacnt)*3)-2:
     +                               seqcod(aacnt)*3)
                  else
                     wkres = hetcod(((seqcod(aacnt)-abetsz)*3)-2:
     +                               (seqcod(aacnt)-abetsz)*3)
                  end if
               end if
            else
CHECK v.3.5-->
C               nline = nline + 1
C               extdat(nline) = brkrec(1:66)
C               extres(nline) = aacnt
                CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -              brkrec(1:66),NREQRD)
CHECK v.3.5<--
            end if
         end if
      else
         if (atmnam .eq. cbet) then
            sctype = 1
         else
CHECK v.3.2-->
C----        Check for special cases of non-standard atom-labelling
C            (eg if file created by X-PLOR)
CHECK v.3.5.3-->
C             rname = aacod3(3*(aacode-1)+1:3*aacode)
             IF (aacode.LE.abetsz) THEN
                 rname = aacod3(3*(aacode-1)+1:3*aacode)
             ELSE
                 rname = hetcod(3*(aacode-abetsz-1)+1:3*(aacode-abetsz))
             ENDIF
CHECK v.3.5.3<--
             if (rname.eq.'ILE' .and. atmnam.eq.' CD ')
     -           atmnam = ' CD1'
             if (rname.eq.'GLN' .and. atmnam.eq.' OE ')
     -           atmnam = ' OE1'
             if (rname.eq.'GLU' .and. atmnam.eq.' OE ')
     -           atmnam = ' OE1'
CHECK v.3.2<--
            sctype = index(scname(aacode),atmnam)
            if (sctype .eq. 0 .and. altpos(aacode) .ne. 0)
     +          sctype = index(altnam(altpos(aacode)),atmnam)
            if (sctype .ne. 0) sctype = (sctype / 4) + 2
         end if
c        alternate side chain locations
         if (sctype .ne. 0) then
            naltsc = naltsc + 1
            altsci(1,naltsc) = aacnt
            altsci(2,naltsc) = sctype
            altsci(3,naltsc) = num
            altscr(1,naltsc) = occup
            altscr(2,naltsc) = therm
            altscc(naltsc)   = altloc
            do 1352 i = 1, ncoord
            altsca(i,naltsc) = acoord(i)
 1352       continue
         end if
      end if
c
c     search for the next atom record and treat as appropriate
c     process coords only if same aa 
c     if new aa start aa process again after setting the count fields
c
 1400 continue
      read (brkfi, '(a)', end=2000) brkrec
CHECK v.3.1-->
      IF (BRKREC(1:KEYLEN) .EQ. ENDMDL) GO TO 2000
CHECK v.3.1<--
      if (brkrec(1:keylen) .ne. endkey) then
         if (brkrec(1:keylen) .eq. terkey) then
CHECK v.3.5-->
C            nline = nline + 1
C            extdat(nline) = brkrec(1:66)
C            extres(nline) = aacnt
             CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -           brkrec(1:66),NREQRD)
CHECK v.3.5<--
         else if (brkrec(1:keylen) .eq. maskey) then
            read (brkrec, '(50X,I5)') nrec
         else
            if (brkrec(1:keylen) .eq. atmkey .or. 
     +          brkrec(1:keylen) .eq. htakey) then
               read (brkrec, '(12x,a4,a1,4x,a6,27x,f6.2)') altatm, 
     +                                     altloc, altnum, altocc
CHECK v.2.1.1-->
               IF (brkrec(57:60).EQ.'    ') altocc = 1.0
CHECK v.2.1.1<--
               if (altnum .eq. seqnum .and. 
     +              keywd .eq. brkrec(1:keylen)) then
                  if (aacode .le. 0) then
CHECK v.3.5-->
C                     if (nline .gt. 1000 .and. aacnt .eq. 0) then
                     if (nline .gt. maxlin .and. aacnt .eq. 0) then
CHECK v.3.5<--
                        seqlen = 0
                        goto 3000
                     end if
                     if (altloc .ne. space .and. nline .ne. 0) then
                        string = extdat(nline)
                        read (string, 1121) extatm, extnum, extocc
CHECK v.2.1.1-->
                        IF (string(57:60).EQ.'    ') extocc = 1.0
CHECK v.2.1.1<--
                        if (extatm .eq. altatm .and. 
     +                      extnum .eq. altnum) then
                           if (extocc .ge. altocc) then
                              altlin = altlin + 1
                              altdat(altlin) = brkrec(1:66)
                           else
                              altlin = altlin + 1
                              altdat(altlin) = string
                              extdat(nline) = brkrec(1:66)
                           end if
                           altres(altlin) = aacnt
                        else
CHECK v.3.5-->
C                           nline = nline + 1
C                           extdat(nline) = brkrec(1:66)
C                           extres(nline) = aacnt
                            CALL STOEXT(extdat,maxlin,nline,extres,
     -                          aacnt,brkrec(1:66),NREQRD)
CHECK v.3.5<--
                        end if
                     else
CHECK v.3.5-->
C                        nline = nline + 1
C                        extdat(nline) = brkrec(1:66)
C                        extres(nline) = aacnt
                         CALL STOEXT(extdat,maxlin,nline,extres,aacnt,
     -                       brkrec(1:66),NREQRD)
CHECK v.3.5<--
                     end if
                     goto 1400
                  end if
                  goto 1300
               else
                  if (aacode .gt. 0) call cntset (atpres, thmprs, 
     +                                            altcod, nscats, 
     +                                            aacode, aacnt, 
     +                                            atcnt,  scatct, 
     +                                            thmcnt, setloc)
                  if (aacnt .ge. maxres) then
CHECK v.1.0-->
C                     write (*, 6010) 
CHECK v.1.0<--
CHECK v.2.1.3-->
C                     write (logfi, 6010) 
                     write (*, 6010) 
CHECK v.2.1.3<--
                     goto 2100
                  else
                     if (aacode .gt. 0) aastd(aacnt+1) = space
                     goto 1100
                  end if
               end if
            end if
         end if
         goto 1400
      end if
c
c     end of input file reached
c
 2000 continue
      if (aacode .gt. 0) call cntset (atpres, thmprs, altcod, nscats, 
     +                                aacode, aacnt,  atcnt,  scatct, 
     +                                thmcnt, setloc)
c
 2100 continue
      seqlen = aacnt
c
 3000 continue
CHECK v.3.5-->
C      if (seqlen .eq. 0) readok = .false.
      if (seqlen .eq. 0 .AND. .NOT.NMR) readok = .false.

C---- If the number of extra residue records has been exceeded, then
C     show how many are required
      IF (NREQRD.GT.MAXLIN) THEN
          PRINT*, '*** Total number of extra lines in file', NREQRD
          PRINT*, '*** Maximum allowed                    ', MAXLIN
          PRINT*, '*** Number of lines lost               ',
     -       NREQRD - MAXLIN
      ENDIF
CHECK v.3.5<--
CHECK v.3.1-->
C      close (brkfi) 
      IF (.NOT.NMR) THEN
          close (brkfi) 
      ENDIF
CHECK v.3.1<--
c 901 write(*, 3001) brkfi, kount
c3001 format (' Error reading input file on unit 'i2,
c    +        ' after 'i5,' lines')
      end
c
c*******************************************************************************
c
      subroutine cntset (atpres, thmprs, altcod, nscats, aacode, aacnt, 
     +                   atcnt,  scatct, thmcnt, setloc)
c     perform end of amino acid flag setting
c
      include 'brkcln.par'
c
      character atpres(3,maxres)*1, thmprs(maxres)*1, altcod(maxres)*1,
     +          setloc*1
c
      integer   nscats(sideaa), aacode, aacnt, scatct, thmcnt, atcnt
      integer   totats
c
      character atmchk*1
c
      altcod(aacnt)   = setloc
      atpres(1,aacnt) = space
      atpres(2,aacnt) = space
      atpres(3,aacnt) = space
      thmprs(aacnt)   = space
      totats          = atcnt + scatct
      if (aacode .le. sideaa) then
         if (nscats(aacode) .ge. 0) then
            atpres(1,aacnt) = atmchk((mnchna+nscats(aacode))-1,totats)
            atpres(2,aacnt) = atmchk(mnchna-1,atcnt)
            atpres(3,aacnt) = atmchk(nscats(aacode),scatct)
            thmprs(aacnt)   = atmchk((mnchna+nscats(aacode))-1,thmcnt)
         end if
      end if
      end
c
c*******************************************************************************
c
      character*1 function atmchk (total, count)
c     set 'present' fields to all none or some as appropriate
c
      include 'brkcln.par'
c
      integer  total, count
c
      if (total .eq. 0) then
         atmchk = space
      else
         if (count .eq. 0) then
            atmchk = 'N'
         else
            if (count .ge. total) then
               atmchk = 'A'
            else
               atmchk = 'S'
            end if
         end if
      end if
      end
c
c*******************************************************************************
c
      integer function strlen (string)
c     find out the real length of the string
c
      include 'brkcln.par'
c
      character string*(*)
c
      integer nxch
c
      do 200 nxch = len(string), 1, -1
      if (string(nxch:nxch) .ne. space) goto 300
  200 continue
      nxch = 0
  300 continue
      strlen = nxch
      end
c
c*******************************************************************************
c
      integer function aaconv (aacode, aacod3, hetcod)
c     find the numerical position of an amino acid code
c
      include 'brkcln.par'
c
CHECK v.3.2-->
      save lastcd
CHECK v.3.2<--
CHECK v.3.5-->
C      character aacode*3, aacod3*(3*abetsz), hetcod*63
      character aacode*3, aacod3*(3*abetsz), hetcod*(3*MXHTYP)
CHECK v.3.5<--
CHECK v.3.0.1--> (Order of statements changed)
      integer resfnd, aaval
CHECK v.3.0.1<--
CHECK v.2.1.4-->
      character lastcd*3
      data lastcd / '   ' /

CHECK v.3.2-->
C      save lastcd
CHECK v.3.2<--
CHECK v.2.1.4<--
c
c
  530 format (' "',a3, '":', a)
c
      aaval = resfnd(aacod3,aacode)
      if (aaval .eq. 0) then
         aaval = resfnd(hetcod,aacode)
         if (aaval .eq. 0) then
            aaval = unkcod
CHECK v.2.1.3-->
C            write (logfi, 530) aacode,' unknown amino acid code'
CHECK v.2.1.4-->
C            write (*, 530) aacode,' unknown amino acid code'
            if (aacode.ne.lastcd) then
                write (*, 530) aacode,' unknown amino acid code'
                lastcd = aacode
            endif
CHECK v.2.1.4<--
CHECK v.2.1.3<--
CHECK v.1.0-->
C            write (*,     530) aacode,' unknown amino acid code'
CHECK v.1.0<--
         else
            aaval = (abetsz + (aaval / 3)) + 1
         end if
      else
         aaval = (aaval / 3) + 1
      end if
      aaconv = aaval
      end
c
c*******************************************************************************
c
      integer function resfnd (soustr, questr)
c
c     implicit none
c     implicit complex(a-z)
c
      character soustr*(*), questr*3
c
      integer respl, nxpl
c
      respl = index(soustr,questr)
  100 continue
      if (mod(respl,3) .ne. 1 .and. respl .ne. 0) then
         nxpl = index(soustr(respl+1:),questr)
         if (nxpl .ne. 0) then
            respl = respl + nxpl
         else
            respl = nxpl
         end if
         goto 100
      end if
      resfnd = respl
      end
c
c*******************************************************************************
c
      subroutine chnbrk (mcaaat, atomin, seqbcd, brksym, chnsz,  chnct, 
     +                   chnend, caonly, seqlen, chstop, chstrt, chnlt)
c     search through looking for distant peptide bonds. if any over
c     the limit are found or if a c or n atom doesn't exist then a
c     chain break is deemed to occur
c
      include 'brkcln.par'
c
      real      mcaaat(ncoord,mnchna,maxres)
c
      logical   atomin(mnchna,maxres), caonly
c
      integer   chnsz(maxchn), chnct, seqlen, chnend(maxchn), chnlt
c
      character seqbcd(maxres)*6, brksym(maxres)*1, chstop(maxres)*1, 
     +          chstrt(maxres)*1
c
      real      atdist
c
      integer   nxres, ninchn, i
c
  500 format (' Chain letter change but no chain break between ', i4,
     + 2x, a6, ' and ', i4, 2x, a6)
c
      do 100 i = 1, maxchn
      chnsz(i) = 0
  100 continue
      do 200 nxres = 1, seqlen
      brksym(nxres) = space
      chstop(nxres) = space
      chstrt(nxres) = space
  200 continue
c
      ninchn = 1
      chnct  = 1
      chnlt  = 1
      do 1000 nxres = 2, seqlen
      if (caonly) then
         if (atdist(mcaaat(1,calph,nxres-1),mcaaat(1,calph,nxres)) 
     +      .gt. cadbnd) then
            call chnset (brksym, chnsz,  chnct,  chnend, ninchn, seqbcd, 
     +                   nxres,  seqlen, chstop, chstrt, chnlt)
         else
            ninchn = ninchn + 1
            if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) then
CHECK v.2.1.3-->
C               write (logfi, 500) nxres-1, seqbcd(nxres-1), 
C     +                            nxres,   seqbcd(nxres)
               write (*, 500) nxres-1, seqbcd(nxres-1), 
     +                            nxres,   seqbcd(nxres)
CHECK v.2.1.3<--
CHECK v.1.0-->
C               write (*,     500) nxres-1, seqbcd(nxres-1), 
C     +                            nxres,   seqbcd(nxres)
CHECK v.1.0<--
            end if
         end if
      else
         if (atomin(carb,nxres-1) .and. atomin(nmain,nxres)) then
            if (atdist(mcaaat(1,carb,nxres-1),
     +                 mcaaat(1,nmain,nxres)) .gt. pepbnd) then
               call chnset (brksym, chnsz, chnct,  chnend, ninchn, 
     +                      seqbcd, nxres, seqlen, chstop, chstrt, 
     +                      chnlt)
            else
               ninchn = ninchn + 1
               if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) then
CHECK v.2.1.3-->
C                  write (logfi, 500) nxres-1, seqbcd(nxres-1), 
C     +                               nxres,   seqbcd(nxres)
                  write (*, 500) nxres-1, seqbcd(nxres-1), 
     +                               nxres,   seqbcd(nxres)
CHECK v.2.1.3<--
CHECK v.1.0-->
C                  write (*,     500) nxres-1, seqbcd(nxres-1), 
C     +                               nxres,   seqbcd(nxres)
CHECK v.1.0<--
                  chstop(nxres-1) = '>'
                  chstrt(nxres)   = '<'
               end if
            end if
         else
            call chnset (brksym, chnsz,  chnct,  chnend, ninchn, seqbcd,
     +                   nxres,  seqlen, chstop, chstrt, chnlt)
         end if
      end if
 1000 continue
      call chnset (brksym,   chnsz,  chnct,  chnend, ninchn, seqbcd, 
     +             seqlen+1, seqlen, chstop, chstrt, chnlt)
      end
c
c*******************************************************************************
c
      subroutine chnset (brksym, chnsz,  chnct,  chnend, ninchn, seqbcd, 
     +                   nxres,  seqlen, chstop, chstrt, chnlt)
c
      include 'brkcln.par'
c
      integer   chnsz(maxchn), chnct, ninchn, seqlen, nxres, 
     +          chnend(maxchn), chnlt
c
      character seqbcd(maxres)*6, brksym(maxres)*1, chstop(maxres)*1, 
     +          chstrt(maxres)*1
c
  100 format (' Chain break between ',i4, ' (', a6, ') and ', i4, ' (',
     + a6, ')')
c
      chnsz(chnct) = ninchn
      if (chnct .eq. 1) then
         chnend(chnct) = ninchn
      else
         chnend(chnct) = chnend(chnct-1) + ninchn
      end if
      ninchn = 1
      if (nxres .le. seqlen) then
         brksym(nxres-1) = brkch
         brksym(nxres)   = brkch
CHECK v.2.1.3-->
C         write (logfi, 100) nxres-1, seqbcd(nxres-1), 
C     +                      nxres,   seqbcd(nxres)
         write (*, 100) nxres-1, seqbcd(nxres-1), 
     +                      nxres,   seqbcd(nxres)
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     100) nxres-1, seqbcd(nxres-1), 
C     +                      nxres,   seqbcd(nxres)
CHECK v.1.0<--
         if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) chnlt = 
     +                                                     chnlt + 1
         chnct = chnct + 1
CHECK v.3.3.1-->
         if (chnct.gt.maxchn) then
             chnct = maxchn
             print*, '*** Maximum number of chain-breaks exceeded',
     -           maxchn
         endif
CHECK v.3.3.1<--
         chstop(nxres-1) = '>'
         chstrt(nxres)   = '<'
      end if
      end
c
c*******************************************************************************
c
      real function atdist (coorda, coordb)
c     routine to determine the distance in space between the two
c     atoms represented by the coordinate arrays.
c
      include 'brkcln.par'
c
      real    coorda(ncoord), coordb(ncoord), dist
c
      integer i
c
      dist = 0.0
      do 200 i = 1, ncoord
      dist = dist + ((coorda(i) - coordb(i)) ** 2)
  200 continue
      atdist = sqrt(dist)
      end
c
c*******************************************************************************
c
      subroutine mksang (mcaaat, atomin, scaaat, scangs, scatin, seqcod, 
     +                   aacod3, hetcod, seqlen, nscats, chstop, chstrt, 
     +                   seqbcd, mcser,  scser,  mcocc,  mcthrm, scocc, 
CHECK v.3.2-->
C     +                   scthrm, scname, altnam, naltmc, naltsc, altmci, 
CHECK v.3.5-->
C     +                   scthrm, scname, naltmc, naltsc, altmci, 
CCHECK v.3.2<--
C     +                   altsci, altmcr, altscr, altmca, altsca, lab, 
C     +                   extdat, extres, nline,  tphnam, olenam, pphnam, 
C     +                   atname, mcangs, nrec,   mcalt,  scalt,  altmcc, 
C     +                   altscc, pyrnam, totswp, altdat, altlin, altres, 
CCHECK v.2.1.3-->
CC     +                   batnam)
CCHECK v.3.2-->
CC     +                   batnam, NATALT, NATZER, NHEALT, NHEZER)
C     +                   batnam, NATALT, NATZER, NHEALT, NHEZER, NMR)
CCHECK v.3.2<--
CCHECK v.2.1.3<--
     +                   scthrm, scname, natsdc, extprc, sswap, naltmc,
     -                   naltsc, altmci,
     +                   altsci, altmcr, altscr, altmca, altsca, lab, 
     +                   extdat, extres, nline,
     +                   atname, mcangs, nrec,   mcalt,  scalt,  altmcc, 
     +                   altscc, totswp, altdat, altlin, altres, 
CHECK v.3.5.3-->
C     +                   MCNAME, NATALT, NATZER, NHEALT, NHEZER, NMR)
     +                   MCNAME, NATALT, NATZER, NHEALT, NHEZER,
     +                   NMR, NOALT)
CHECK v.3.5.3<--
CHECK v.3.5<--
c     calculate chi angles from atom data. check for that alternates
c     are closest to c alpha
c
      include 'brkcln.par'
c
      real      mcaaat(ncoord,mnchna,maxres), 
     +          scaaat(ncoord,sdchna,maxres), scangs(nsdang,maxres), 
     +          mcocc(mnchna,maxres), mcthrm(mnchna,maxres), 
     +          scocc(sdchna,maxres), scthrm(sdchna,maxres), 
     +          altmcr(2,altmax), altscr(2,altmax), 
     +          altmca(3,altmax), altsca(3,altmax), 
     +          mcangs(nmnang,maxres)
c
CHECK v.3.5-->
C      character aacod3*(abetsz*3), hetcod*63, seqbcd(maxres)*6, 
      character aacod3*(abetsz*3), hetcod*(3*MXHTYP), seqbcd(maxres)*6, 
CHECK v.3.5<--
CHECK v.3.2-->
C     +          scname(sideaa)*44, altnam(3)*44, lab(maxres)*1, 
     +          scname(sideaa)*44, lab(maxres)*1, 
CHECK v.3.2<--
CHECK v.3.5-->
C     +          extdat(maxlin)*66, tphnam*((4*mnchna)-4), 
C     +          olenam*((4*mnchna)-4), pphnam*((4*mnchna)-4), 
C     +          atname*((4*mnchna)-4), chstop(maxres)*1, 
C     +          chstrt(maxres)*1, mcalt(mnchna,maxres)*1, 
C     +          scalt(sdchna,maxres)*1, altmcc(maxres)*1, 
C     +          altscc(maxres)*1, pyrnam*((4*mnchna)-4), 
C     +          altdat(maxalt)*66, batnam*((4*mnchna)-4)
     +          extdat(maxlin)*66,
     +          atname*((4*mnchna)-4), chstop(maxres)*1, 
     +          chstrt(maxres)*1, mcalt(mnchna,maxres)*1, 
     +          scalt(sdchna,maxres)*1, altmcc(maxres)*1, 
     +          altscc(maxres)*1, MCNAME(SIDEAA),
     +          altdat(maxalt)*66
CHECK v.3.5<--
c
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres)
c
      integer   seqcod(maxres), seqlen, mcser(mnchna,maxres), 
     +          scser(sdchna,maxres), naltmc, naltsc, 
     +          altmci(3,altmax), altsci(3,altmax), 
     +          extres(maxlin), nline, nrec, altres(maxalt), 
     +          altlin, natsdc(sideaa), nxang, i, j, nxchi, nxres, nchi, 
     +          aacod, extprc(sideaa), swatm, extatm, chipl, nswap, 
     +          wkswap, nold(sdchna,maxres), nscats(sideaa), k, 
     +          change(maxres), n, ii, istart, nathet, 
CHECK v.3.5-->
C     +          totswp(sideaa), astart
     +          totswp(sideaa), astart, sswap(2,sideaa)
CHECK v.3.5<--
c
c     integer   altref
c
CHECK v.3.2-->
C      logical   athere(maxchi+3), alt, xney
CHECK v.3.5.3-->
C      logical   athere(maxchi+3), NMR, xney
      logical   athere(maxchi+3), NMR, NOALT, xney
CHECK v.3.5.3<--
CHECK v.3.2<--
c
      real      relats(ncoord,maxchi+3), distmn, distal, ang1, ang2, 
     +          distg1, distg2
c
CHECK v.3.5.3-->
      CHARACTER ALTCH*1
CHECK v.3.5.3<--
      character outerr*80, wkres*3, scdata*44, ld*1
c     character altps1, altps2
      character cbet*4, scatom*4, scold*4, mcatom*4
      character label*6
CHECK v.2.1-->
      character altlab*6
CHECK v.2.1<--
c
      real      atdist, dihed
CHECK v.2.0.1-->
      real      dihang
CHECK v.2.0.1<--
c

CHECK v.2.1.3-->
      INTEGER       NATALT, NATZER, NHEALT, NHEZER
CHECK v.2.1.3<--

CHECK v.3.5-->
C      save extprc, natsdc
Cc
C      data natsdc / 2*0, 1, 2, 3, 2, 0, 2*2, 0, 4, 2, 3, 2, 0, 1, 3, 5, 
C     +              2*1, 0, 1, 2, 0, 2, 2*0, 3, 4, 0, 2, 3*1, 2, 0, 4, 
C     +              0, 2, 0, 2*2, 1, 3, 1, 2*2, 3, 0, 2 /
Cc
C      data extprc / 3*0, 3*4, 2*0, 5, 2*0, 2, 5*0, 4, 0, 3, 0, 2, 2*0, 
C     +              4, 5*0, 4, 3*0, 4, 3*0, 4, 0, 2, 4, 2*0, 2, 5, 4, 
C     +              2*0, 4 /
CHECK v.3.5<--
c
      data cbet / ' CB ' /
c
 6000 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2)
CHECK v.3.0.1-->
C 6001 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, a2)
C 6002 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, i1)
CHECK v.3.0.1<--
 6003 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, a1,
     + a1, a4, a1)
CHECK v.3.0.1-->
C 6004 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, a1,
C     + a1)
CHECK v.3.0.1<--
 6005 format (a66)
 6006 format ('MASTER',i6)
 6007 format ('END')
 6008 format (a6, i5, 1x, a4, a1, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, a1,
     + a1, 4x, a1)
 6010 format (' distg1 is', 2f10.2)
 6020 format (' distg2 is', 2f10.2)
CHECK v.2.0.1-->
C 6030 format (' ILE residue still wrong!!')
CHECK v.2.1.3-->
C 6030 format (' ILE residue still wrong!!',1x,i4)
CHECK v.2.1.4-->
 6030 format (' **** ILE residue has wrong chirality at CB',1x,i4)
CHECK v.2.1.4<--
CHECK v.2.1.3<--
CHECK v.2.0.1<--
CHECK v.2.1.3-->
C 6040 format (' side chain atoms swapped for residues:')
 6040 format (' * Side chain atoms swapped for residues:')
CHECK v.3.4-->
 6042 format (' * Side chain atoms should be swapped for the following',
     -    ' residues:',/,
     -    ' * As this is an NMR ensemble the swap has not been made')
CHECK v.3.4<--
CHECK v.2.1.3<--
CHECK v.3.2-->
C 6050 format (a3,1x,i4,2x)
 6050 format (a3,1x,a6,1x)
CHECK v.3.2<--
 6060 format (1x,'*',1X, a)
CHECK v.3.2-->
C 6070 format (1x,'#',1X, a)
CHECK v.3.2<--
 6100 format (' WARNING: Number of atoms/hetatms in new file differs',
     + ' from old file - ',/,
     + 10X, ' new number = ',I5,/,
     + 10X, ' old number = ',I5)
 6120 format (' Number of alternate mainchain atoms is ', i5)
 6130 format (' Number of alternate sidechain atoms is ', i5)
c
      astart = 1
      istart = 1
      nswap  = 0
      outerr = space
      nathet = 0
      do 10 i = 1, sdchna
      do 9 j = 1, maxres
      nold(i,j) = 0
    9 continue
c     write out any extra data before residue 1
   10 continue
      if (nline .gt. 0) then
         do 50 ii = 1, nline
         if (extres(ii) .eq. 0) then
            write (outfi, 6005) extdat(ii)
            if (extdat(ii)(1:keylen) .eq. atmkey .or. 
     +          extdat(ii)(1:keylen) .eq. htakey) nathet = nathet + 1
         else
            if (extres(ii) .gt. 0) then
               istart = ii
               goto 51
            end if
         end if
   50    continue
      end if
   51 continue
      if (altlin .gt. 0) then
         do 60 ii = 1, altlin
         if (altres(ii) .eq. 0) then
            write (altfi, 6005) altdat(ii)
            if (altdat(ii)(1:keylen) .eq. atmkey .or. 
     +          altdat(ii)(1:keylen) .eq. htakey) nathet = nathet + 1
         else
            if (altres(ii) .gt. 0) then
               astart = ii
               goto 61
            end if
         end if
   60    continue
      end if
   61 continue
      do 5000 nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         wkres = aacod3((seqcod(nxres)*3)-2:seqcod(nxres)*3)
      else
         wkres = hetcod(((seqcod(nxres)-abetsz)*3)-2:
     +                   (seqcod(nxres)-abetsz)*3)
      end if
      aacod = seqcod(nxres)
      nchi = 0
      nchi = natsdc(aacod)
      do 100 nxang = nchi+1, nsdang
      scangs(nxang,nxres) = null
  100 continue
      if (nchi .gt. 0) then
         athere(1) = atomin(nmain,nxres)
         athere(2) = atomin(calph,nxres)
         do 200 i = 1, ncoord
         relats(i,1) = mcaaat(i,nmain,nxres)
         relats(i,2) = mcaaat(i,calph,nxres)
  200    continue
         extatm = 0
         if (extprc(aacod) .gt. 1) extatm = 1
         do 400 i = 1, nchi+1+extatm
         athere(2+i) = scatin(i,nxres)
         do 300 j = 1, ncoord
         relats(j,2+i) = scaaat(j,i,nxres)
  300    continue
  400    continue
c        ILE CASE
         swatm = 0
         if (extprc(aacod) .eq. 5) then
            if (scatin(nchi,nxres) .and. scatin(nchi+2,nxres)) then
               distg1 = atdist(scaaat(1,nchi,nxres),
     +                         scaaat(1,nchi+2,nxres))
               if (distg1 .gt. ccdst) then
CHECK v.2.1.3-->
C                  write (logfi, 6010) distg1, ccdst
                  write (*, 6010) distg1, ccdst
CHECK v.2.1.3<--
CHECK v.1.0-->
C                  write (*,     6010) distg1, ccdst
CHECK v.1.0<--
                  if (scatin(nchi+1,nxres)) then
                     distg2 = atdist(scaaat(1,nchi+1,nxres),
     +                               scaaat(1,nchi+2,nxres))
                     if (distg2 .lt. ccdst) then
CHECK v.2.1.3-->
C                        write (logfi, 6020) distg2, ccdst
                        write (logfi, 6020) distg2, ccdst
CHECK v.2.1.3<--
CHECK v.1.0-->
C                        write (*,     6020) distg2, ccdst
CHECK v.1.0<--
                        swatm = 1
CHECK v.3.2-->
C                        call swap (aacod, nxres, scaaat, scocc, scthrm, 
C     +                             Scatin, nold, totswp)
                        call swap (aacod, nxres, scaaat, scocc, scthrm, 
CHECK v.3.4-->
C     +                             nold, totswp)
CHECK v.3.5-->
C     +                             nold, totswp,NMR)
     +                             nold, totswp,sswap,NMR)
CHECK v.3.5<--
CHECK v.3.4<--
CHECK v.3.2<--
                        change(nxres) = 1
                     end if
                  end if
               end if
            end if
CHECK v.2.0.1-->
            dihang=dihed(1,scaaat(1,1,nxres),scaaat(1,2,nxres),
     -          scaaat(1,3,nxres),mcaaat(1,calph,nxres))
CHECK v.2.1.3-->
C            if (dihang.lt.0) write(logfi,6030) nxres
CHECK v.2.1.3<--
            if (dihang.lt.0) write(*,6030) nxres
CHECK v.2.0.1<--
         end if
         if (extprc(aacod) .eq. 1) then
            if (atomin(calph,nxres) .and. scatin(nchi+1,nxres) .and. 
     +                                    scatin(nchi+2,nxres)) then
               distmn = atdist(mcaaat(1,calph,nxres),
     +                         scaaat(1,nchi+1,nxres))
               distal = atdist(mcaaat(1,calph,nxres),
     +                         scaaat(1,nchi+2,nxres))
               if (distal .lt. distmn) then
                  swatm = 1
CHECK v.3.2-->
C                  call swap (aacod, nxres, scaaat, scocc, scthrm, 
C     +                       scatin, nold, totswp)
                  call swap (aacod, nxres, scaaat, scocc, scthrm, 
CHECK v.3.4-->
C     +                       nold, totswp)
CHECK v.3.5-->
C     +                       nold, totswp,NMR)
     +                       nold, totswp,sswap,NMR)
CHECK v.3.5<--
CHECK v.3.4<--
CHECK v.3.2<--
                  change(nxres) = 1
                  athere(nchi+3) = scatin(nchi+2,nxres)
                  do 500 i = 1, ncoord
                  relats(i,nchi+3) = scaaat(i,nchi+2,nxres)
  500             continue
               end if
            end if
         end if
         do 1000 nxchi = 1, nchi
         if (athere(nxchi)   .and. athere(nxchi+1) .and. 
     +       athere(nxchi+2) .and. athere(nxchi+3)) then
            scangs(nxchi,nxres) = dihed(nmnang+nxchi,
     +                                  relats(1,nxchi),
     +                                  relats(1,nxchi+1),
     +                                  relats(1,nxchi+2),
     +                                  relats(1,nxchi+3))
         else
            scangs(nxchi,nxres) = null
         end if
 1000    continue
         if (extprc(aacod) .gt. 1 .and. xney(scangs(nchi,nxres),null)) 
     +       then
            if (athere(nchi)   .and. athere(nchi+1) .and. 
     +          athere(nchi+2) .and. athere(nchi+3+extatm)) then
               scangs(nchi+1,nxres) = dihed(nmnang+nchi+extatm,
     +                                      relats(1,nchi),
     +                                      relats(1,nchi+1),
     +                                      relats(1,nchi+2),
     +                                      relats(1,nchi+3+extatm))
            else
               scangs(nchi+1,nxres) = null
            end if
            if (xney(scangs(nchi+1,nxres),null)) then
               if (extprc(aacod) .eq. 2 .or. 
     +             extprc(aacod) .eq. 3 .or.
     +             extprc(aacod) .eq. 5) then
                  ang1 = scangs(nchi,nxres)
                  ang2 = scangs(nchi+1,nxres)
                  if (ang1 .lt. 0.0) ang1 = 360.0 + ang1
                  if (ang2 .lt. 0.0) ang2 = 360.0 + ang2
                  chipl = nchi
                  if (abs(ang1-ang2) .gt. 180.0) then
                     if (ang2 .gt. ang1 .and. extprc(aacod) .eq. 2)
     +                  chipl = nchi + 1
                     if (ang1 .gt. ang2 .and. extprc(aacod) .eq. 3)
     +                  chipl = nchi + 1
                  else
                     if (ang2 .lt. ang1 .and. extprc(aacod) .eq. 2)
     +                  chipl = nchi + 1
                     if (ang1 .lt. ang2 .and. extprc(aacod) .eq. 3)
     +                  chipl = nchi + 1
                  end if
                  if (chipl .ne. nchi) then
CHECK v.2.1.3-->
C                     if (extprc(aacod) .eq. 5) write (logfi, 6030) 
                     if (extprc(aacod) .eq. 5) write (*, 6030) 
CHECK v.2.1.3<--
CHECK v.1.0-->
C                     if (extprc(aacod) .eq. 5) write (*,     6030) 
CHECK v.1.0<--
                     swatm = 1
CHECK v.3.2-->
C                     call swap (aacod, nxres, scaaat, scocc, scthrm, 
C     +                          scatin, nold, totswp)
                     call swap (aacod, nxres, scaaat, scocc, scthrm, 
CHECK v.3.4-->
C     +                          nold, totswp)
CHECK v.3.5-->
C     +                          nold, totswp,NMR)
     +                          nold, totswp,sswap,NMR)
CHECK v.3.5<--
CHECK v.3.4<--
CHECK v.3.2<--
                     change(nxres) = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               if (extprc(aacod) .eq. 4) then
                  if (abs(scangs(nchi+1,nxres)) .lt. 
     +                abs(scangs(nchi,nxres))) then
                     swatm = 1
CHECK v.3.2-->
C                     call swap (aacod, nxres, scaaat, scocc, scthrm, 
C     +                          scatin, nold, totswp)
                     call swap (aacod, nxres, scaaat, scocc, scthrm, 
CHECK v.3.4-->
C     +                          nold, totswp)
CHECK v.3.5-->
C     +                          nold, totswp,NMR)
     +                          nold, totswp,sswap,NMR)
CHECK v.3.5<--
CHECK v.3.4<--
CHECK v.3.2<--
                     change(nxres) = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               scangs(nchi+1,nxres) = null
            end if
         end if
         if (swatm .eq. 1) then
CHECK v.2.1.3-->
C            if (nswap .eq. 0) write (logfi, 6040) 
CHECK v.3.4-->
C            if (nswap .eq. 0) write (*, 6040) 
            if (nswap .eq. 0) THEN
                IF (.NOT.NMR) THEN
                    write (*, 6040)
                ELSE
                    write (*, 6042)
                ENDIF
            endif
CHECK v.3.4<--
CHECK v.2.1.3<--
CHECK v.1.0-->
C            if (nswap .eq. 0) write (*,     6040) 
CHECK v.1.0<--
CHECK v.3.2-->
C            wkswap = mod(nswap,8) + 1
C            write (outerr(((wkswap-1)*10)+1:wkswap*10), 6050) wkres, 
C     +                                                        nxres
C            if (wkswap .eq. 8) then
            wkswap = mod(nswap,7) + 1
            write (outerr(((wkswap-1)*11)+1:wkswap*11), 6050) wkres, 
     +          seqbcd(nxres)
            if (wkswap .eq. 7) then
CHECK v.3.2<--
CHECK v.2.1.3-->
C               write (logfi, 6060) outerr
               write (*, 6060) outerr
CHECK v.2.1.3<--
CHECK v.1.0-->
C               write (*,     6060) outerr
CHECK v.1.0<--
               outerr = space
            end if
            nswap = nswap + 1
         end if
      end if
      if (lab(nxres) .eq. 'A') then
         label = 'ATOM  '
      else
         label = 'HETATM'
      end if
c     write data out to new brookhaven file
      do 3000 j = 1, 5
CHECK v.3.5-->
C      if (aacod .gt. 44 .and. aacod .le. 50) then
C         mcatom = batnam((4*j)-3:4*j)
C      else if (aacod .eq. 42) then
C         mcatom = tphnam((4*j)-3:4*j)
C      else if (aacod .eq. 41) then
C         mcatom = olenam((4*j)-3:4*j)
C      else if (aacod .eq. 39) then
C         mcatom = pphnam((4*j)-3:4*j)
C      else if (aacod .eq. 36) then
C         mcatom = pyrnam((4*j)-3:4*j)
C      else
C         mcatom = atname((4*j)-3:4*j)
C      end if
       IF (MCNAME(AACOD).NE.' ') THEN
           MCATOM = MCNAME(AACOD)((4*j)-3:4*j)
       ELSE
           MCATOM = ATNAME((4*j)-3:4*j)
       ENDIF
CHECK v.3.5<--
      if ((mcangs(1,nxres) .lt. 23.0 .or. 
     +     mcangs(1,nxres) .gt. 45.0) .and.
     +      xney(mcangs(1,nxres),null)) then
         ld = '*'
      else
         ld = ' '
      end if
      if (atomin(j,nxres)) then
CHECK v.2.1-->
         altlab = label
CHECK v.3.2-->
C         if (mcocc(j,nxres).eq.0.0) then
         if (.NOT.NMR .AND. mcocc(j,nxres).eq.0.0) then
CHECK v.3.2<--
             altlab = label(1:2) // 'ZERO'
CHECK v.2.1.3-->
             IF (label(1:2).EQ.'AT') THEN
                 NATZER = NATZER + 1
             ELSE
                 NHEZER = NHEZER + 1
             ENDIF
CHECK v.2.1.3<--
         endif
C         write (outfi, 6008) label, mcser(j,nxres), mcatom, 
CHECK v.3.5.3-->
         ALTCH = mcalt(j,nxres)
         IF (NOALT) ALTCH = ' '
CHECK v.3.5.3<--
         write (outfi, 6008) altlab, mcser(j,nxres), mcatom, 
CHECK v.2.1<--
CHECK v.3.5.3-->
C     +                       mcalt(j,nxres), wkres, seqbcd(nxres), 
     +                       ALTCH, wkres, seqbcd(nxres), 
CHECK v.3.5.3<--
     +                      (mcaaat(i,j,nxres), i=1,ncoord), 
     +                       mcocc(j,nxres), mcthrm(j,nxres), 
     +                       chstop(nxres), chstrt(nxres), ld
         nathet = nathet + 1
      end if
CHECK v.3.5.3-->
C      if (naltmc .gt. 0 .and. mcalt(j,nxres) .ne. ' ') then
      if (naltmc .gt. 0 .and. mcalt(j,nxres) .ne. ' ' .AND.
     -    .NOT.NOALT) then
CHECK v.3.5.3<--
         do 2900 i = 1, naltmc
         if (altmci(2,i) .eq. j .and. altmci(1,i) .eq. nxres) then
            write (altfi, 6000) label, altmci(3,i), mcatom, altmcc(i), 
     +                          wkres, seqbcd(nxres), 
     +                         (altmca(k,i), k=1,3), altmcr(1,i), 
     +                          altmcr(2,i)

CHECK v.2.1-->
               altlab = label(1:2) // 'ALT'
               write (outfi, 6000) altlab, altsci(3,j), scatom, 
     +                             altscc(j), wkres, seqbcd(nxres), 
     +                            (altsca(k,j), k=1,3), altscr(1,j), 
     +                             altscr(2,j)
CHECK v.2.1<--

CHECK v.2.1.3-->
             IF (label(1:2).EQ.'AT') THEN
                 NATALT = NATALT + 1
             ELSE
                 NHEALT = NHEALT + 1
             ENDIF
CHECK v.2.1.3<--
            nathet = nathet + 1
         end if
 2900    continue
      end if
c     side chain atoms
 3000 continue
      scdata = scname(aacod)
      do 4000 i = 1, nscats(aacod)
      if (scatin(i,nxres)) then
         if (i .eq. 1) then
            scatom = cbet
         else
            scatom = scdata((4*(i-2))+1:4*(i-1))
c          determine the original atom name
         end if
         if (nold(i,nxres) .ne. 0) then
            n = nold(i,nxres)
            scold = scdata((4*(n-2))+1:4*(n-1))
         else
            scold = '    '
         end if
CHECK v.2.1-->
         altlab = label
CHECK v.3.2-->
C         if (scocc(i,nxres).eq.0.0) then
         if (.NOT.NMR .AND. scocc(i,nxres).eq.0.0) then
CHECK v.3.2<--
             altlab = label(1:2) // 'ZERO'
CHECK v.2.1.3-->
             IF (label(1:2).EQ.'AT') THEN
                 NATZER = NATZER + 1
             ELSE
                 NHEZER = NHEZER + 1
             ENDIF
CHECK v.2.1.3<--
         endif
C         write (outfi, 6003) label, scser(i,nxres), scatom, 
CHECK v.3.5.3-->
         ALTCH = scalt(i,nxres)
         IF (NOALT) ALTCH = ' '
CHECK v.3.5.3<--
         write (outfi, 6003) altlab, scser(i,nxres), scatom,
CHECK v.2.1<--
CHECK v.3.5.3-->
C     +                       scalt(i,nxres), wkres, seqbcd(nxres), 
     +                       ALTCH, wkres, seqbcd(nxres), 
CHECK v.3.5.3<--
     +                      (scaaat(j,i,nxres), j=1,ncoord), 
     +                       scocc(i,nxres), scthrm(i,nxres), 
     +                       chstop(nxres), chstrt(nxres), scold, ld
         nathet = nathet + 1
CHECK v.3.5.3-->
C         if (naltsc .gt. 0 .and. scalt(i,nxres) .ne. ' ') then
         if (naltsc .gt. 0 .and. scalt(i,nxres) .ne. ' ' .AND.
     -       .NOT.NOALT) then
CHECK v.3.5.3<--
            do 3960 j = 1, naltsc
            if (altsci(1,j) .eq. nxres .and. 
     +          altsci(2,j) .eq. i) then
               write (altfi, 6000) label, altsci(3,j), scatom, 
     +                             altscc(j), wkres, seqbcd(nxres), 
     +                            (altsca(k,j), k=1,3), altscr(1,j), 
     +                             altscr(2,j)

CHECK v.2.1-->
               altlab = label(1:2) // 'ALT'
               write (outfi, 6000) altlab, altsci(3,j), scatom, 
     +                             altscc(j), wkres, seqbcd(nxres), 
     +                            (altsca(k,j), k=1,3), altscr(1,j), 
     +                             altscr(2,j)
CHECK v.2.1<--

CHECK v.2.1.3-->
               IF (label(1:2).EQ.'AT') THEN
                   NATALT = NATALT + 1
               ELSE
                   NHEALT = NHEALT + 1
               ENDIF
CHECK v.2.1.3<--
               nathet = nathet + 1
            end if
 3960       continue
         end if
      end if

c     now write out any extra data associated with this residue
 4000 continue
      if (nline .gt. 0) then
         do 4500 ii = istart, nline
         if (extres(ii) .eq. nxres) then
            write (outfi, 6005) extdat(ii)
            if (extdat(ii)(1:keylen) .eq. atmkey .or. 
     +          extdat(ii)(1:keylen) .eq. htakey)
     +         nathet = nathet + 1
         else
            if (extres(ii) .gt. nxres) then
               istart = ii
               goto 4510
            end if
         end if
 4500    continue
      end if
 4510 continue
      if (altlin .gt. 0) then
         do 4600 ii = astart, altlin
         if (altres(ii) .eq. nxres) then
            write (altfi, 6005) altdat(ii)
            if (altdat(ii)(1:keylen) .eq. atmkey .or. 
     +          altdat(ii)(1:keylen) .eq. htakey)
     +         nathet = nathet + 1
         else
            if (altres(ii) .gt. nxres) then
               astart = ii
               goto 4610
            end if
         end if
 4600    continue
      end if
 4610 continue
 5000 continue
CHECK v.3.5<--
      IF (.NOT.NMR) THEN
CHECK v.3.5-->
          write (outfi, 6006) nathet
CHECK v.3.5-->
      ENDIF
CHECK v.3.5<--
      if (nathet .ne. nrec) then
CHECK v.2.1.3-->
C         write (logfi, 6100) nathet, nrec
         write (*, 6100) nathet, nrec
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     6100) nathet, nrec
CHECK v.1.0<--
      end if
CHECK v.3.2-->
C      if (alt) write (altfi, 6006) 
CHECK v.3.2<--
CHECK v.3.5<--
      IF (.NOT.NMR) THEN
CHECK v.3.5-->
          write (outfi, 6007) 
CHECK v.3.5-->
      ENDIF
CHECK v.3.5<--
CHECK v.3.2-->
C      if (mod(nswap,8) .ne. 0) then
      if (mod(nswap,7) .ne. 0) then
CHECK v.3.2<--
CHECK v.2.1.3-->
C         write (logfi, 6070) outerr
CHECK v.3.2-->
C         write (*, 6070) outerr
         write (*, 6060) outerr
CHECK v.3.2<--
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     6070) outerr
CHECK v.1.0<--
      end if
      if (naltmc .gt. 0) then
CHECK v.2.1.3-->
C         write (logfi, 6120) naltmc
         write (*, 6120) naltmc
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     6120) naltmc
CHECK v.1.0<--
      end if
      if (naltsc .gt. 0) then
CHECK v.2.1.3-->
C         write (logfi, 6130) naltsc
         write (*, 6130) naltsc
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     6130) naltsc
CHECK v.1.0<--
      end if
c     call msgprt
      end
c
c*******************************************************************************
c
CHECK v.3.2-->
C      subroutine swap (aacod, nxres, scaaat, scocc, scthrm, scatin, 
C     +                 nold, totswp)
      subroutine swap (aacod, nxres, scaaat, scocc, scthrm,
CHECK v.3.4-->
C     +                 nold, totswp)
CHECK v.3.5-->
C     +                 nold, totswp,NMR)
     +                 nold, totswp,sswap,NMR)
CHECK v.3.5<--
CHECK v.3.4<--
CHECK v.3.2<--
c     subroutine to swap round side chain coordinates for misnamed 
c     residues and record the original atom corresponding to the 
c     coords
c
      include 'brkcln.par'
c
      real    scaaat(ncoord,sdchna,maxres), scocc(sdchna,maxres), 
     +        scthrm(sdchna,maxres), dumaat(ncoord), dumocc, dumthm
c
CHECK v.3.2-->
C      logical scatin(sdchna,maxres)
CHECK v.3.2<--
CHECK v.3.4-->
      LOGICAL       NMR
CHECK v.3.4<--
c
      integer sswap(2,sideaa), aacod, nold(sdchna,maxres), i, j, 
     +        nside, nxres, totswp(sideaa)
c
CHECK v.3.5-->
C      data sswap(1, 1) / 0 /
C      data sswap(1, 2) / 0 /
C      data sswap(1, 3) / 0 /
C      data sswap(1, 4) / 3 /
C      data sswap(1, 5) / 4 /
C      data sswap(1, 6) / 3 /
C      data sswap(1, 7) / 0 /
C      data sswap(1, 8) / 0 /
C      data sswap(1, 9) / 2 /
C      data sswap(1, 10) / 0 /
C      data sswap(1, 11) / 0 /
C      data sswap(1, 12) / 3 /
C      data sswap(1, 13) / 0 /
C      data sswap(1, 14) / 0 /
C      data sswap(1, 15) / 0 /
C      data sswap(1, 16) / 0 /
C      data sswap(1, 17) / 0 /
C      data sswap(1, 18) / 6 /
C      data sswap(1, 19) / 0 /
C      data sswap(1, 20) / 2 /
C      data sswap(1, 21) / 0 /
C      data sswap(1, 22) / 2 /
C      data sswap(1, 23) / 0 /
C      data sswap(1, 24) / 0 /
C      data sswap(1, 25) / 3 /
C      data sswap(1, 26) / 0 /
C      data sswap(1, 27) / 0 /
C      data sswap(1, 28) / 0 /
C      data sswap(1, 29) / 0 /
C      data sswap(1, 30) / 0 /
C      data sswap(1, 31) / 3 /
C      data sswap(1, 32) / 0 /
C      data sswap(1, 33) / 0 /
C      data sswap(1, 34) / 0 /
C      data sswap(1, 35) / 3 /
C      data sswap(1, 36) / 0 /
C      data sswap(1, 37) / 0 /
C      data sswap(1, 38) / 0 /
C      data sswap(1, 39) / 3 /
C      data sswap(1, 40) / 0 /
C      data sswap(1, 41) / 3 /
C      data sswap(1, 42) / 3 /
C      data sswap(1, 43) / 0 /
C      data sswap(1, 44) / 0 /
C      data sswap(1, 45) / 2 /
C      data sswap(1, 46) / 2 /
C      data sswap(1, 47) / 3 /
C      data sswap(1, 48) / 0 /
C      data sswap(1, 49) / 0 /
C      data sswap(1, 50) / 3 /
C      data sswap(2, 1) / 0 /
C      data sswap(2, 2) / 0 /
C      data sswap(2, 3) / 0 /
C      data sswap(2, 4) / 0 /
C      data sswap(2, 5) / 0 /
C      data sswap(2, 6) / 5 /
C      data sswap(2, 7) / 0 /
C      data sswap(2, 8) / 0 /
C      data sswap(2, 9) / 0 /
C      data sswap(2, 10) / 0 /
C      data sswap(2, 11) / 0 /
C      data sswap(2, 12) / 0 /
C      data sswap(2, 13) / 0 /
C      data sswap(2, 14) / 0 /
C      data sswap(2, 15) / 0 /
C      data sswap(2, 16) / 0 /
C      data sswap(2, 17) / 0 /
C      data sswap(2, 18) / 0 /
C      data sswap(2, 19) / 0 /
C      data sswap(2, 20) / 0 /
C      data sswap(2, 21) / 0 /
C      data sswap(2, 22) / 0 /
C      data sswap(2, 23) / 0 /
C      data sswap(2, 24) / 0 /
C      data sswap(2, 25) / 5 /
C      data sswap(2, 26) / 0 /
C      data sswap(2, 27) / 0 /
C      data sswap(2, 28) / 0 /
C      data sswap(2, 29) / 0 /
C      data sswap(2, 30) / 0 /
C      data sswap(2, 31) / 5 /
C      data sswap(2, 32) / 0 /
C      data sswap(2, 33) / 0 /
C      data sswap(2, 34) / 0 /
C      data sswap(2, 35) / 5 /
C      data sswap(2, 36) / 0 /
C      data sswap(2, 37) / 0 /
C      data sswap(2, 38) / 0 /
C      data sswap(2, 39) / 5 /
C      data sswap(2, 40) / 0 /
C      data sswap(2, 41) / 0 /
C      data sswap(2, 42) / 5 /
C      data sswap(2, 43) / 0 /
C      data sswap(2, 44) / 0 /
C      data sswap(2, 45) / 0 /
C      data sswap(2, 46) / 0 /
C      data sswap(2, 47) / 5 /
C      data sswap(2, 48) / 0 /
C      data sswap(2, 49) / 0 /
C      data sswap(2, 50) / 5 /
CHECK v.3.5<--
c
          do 100 i = 1, 2
              if (sswap(i,aacod) .ne. 0) then
                 nside = sswap(i,aacod)
CHECK v.3.4-->
      IF (.NOT.NMR) THEN
CHECK v.3.4<--
                 do 10 j = 1, ncoord
                     dumaat(j) = scaaat(j,nside,nxres)
                     scaaat(j,nside,nxres) = scaaat(j,nside+1,nxres)
                     scaaat(j,nside+1,nxres) = dumaat(j)
   10            continue
                 dumocc = scocc(nside,nxres)
                 scocc(nside,nxres) = scocc(nside+1,nxres)
                 scocc(nside+1,nxres) = dumocc
                 dumthm = scthrm(nside,nxres)
                 scthrm(nside,nxres) = scthrm(nside+1,nxres)
                 scthrm(nside+1,nxres) = dumthm
CHECK v.3.4-->
      ENDIF
CHECK v.3.4<--
                 nold(nside,nxres) = nside + 1
                 nold(nside+1,nxres) = nside
              end if
  100     continue
      totswp(aacod) = totswp(aacod) + 1
      end
c
c*******************************************************************************
c
      subroutine calout (mcaaat, atomin, seqcod, aacod3, hetcod, seqlen, 
     +                   chstop, chstrt, seqbcd, mcser,  mcocc,  mcthrm, 
CHECK v.3.5-->
C     +                   extdat, nline,  extres, nrec)
     +                   extdat, nline,  extres, nrec,   NMR)
CHECK v.3.5<--
c     generates output file for c-alpha only files
c
      include 'brkcln.par'
c
      real      mcaaat(ncoord,mnchna,maxres), mcocc(mnchna,maxres), 
     +          mcthrm(mnchna,maxres)
c
CHECK v.3.5-->
C      character aacod3*(abetsz*3), hetcod*63, chstop(maxres)*1, 
      character aacod3*(abetsz*3), hetcod*(3*MXHTYP), chstop(maxres)*1, 
CHECK v.3.5<--
     +          seqbcd(maxres)*6, wkres*3, chstrt(maxres)*1, 
     +          extdat(maxlin)*66
c
      logical   atomin(mnchna,maxres)
CHECK v.3.5-->
      LOGICAL   NMR
CHECK v.3.5<--
c
      integer   seqcod(maxres), seqlen, mcser(mnchna,maxres), nxres, i, 
     +          aacod, nathet, ii, nline, extres(maxlin), istart, nrec
c
 6000 format (a6, i5, 1x, a4, 1x, a3, 1x, a6, 3x, 3f8.3, 2f6.2, 1x, a1,
     + a1)
 6005 format (a66)
 6006 format ('MASTER',i6)
 6007 format ('END')
 6100 format (' WARNING: Number of atoms/hetatms in new file differs',
     + ' from old file - ',/,
     + 10X, ' new number = ',I5,/,
     + 10X, ' old number = ',I5)
c
      nathet = 0
      if (nline .gt. 0) then
         do 50 ii = 1, nline
         if (extres(ii) .eq. 0) then
            write (outfi, 6005) extdat(ii)
            if (extdat(ii)(1:keylen) .eq. atmkey .or. 
     +          extdat(ii)(1:keylen) .eq. htakey) nathet = nathet + 1
         else
            if (extres(ii) .gt. 0) then
               istart = ii
               goto 51
            end if
         end if
   50    continue
      end if
   51 continue
      do 5000 nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         wkres = aacod3((seqcod(nxres)*3)-2:seqcod(nxres)*3)
      else
         wkres = hetcod(((seqcod(nxres)-abetsz)*3)-2:
     +                   (seqcod(nxres)-abetsz)*3)
      end if
      aacod = seqcod(nxres)
      if (atomin(2,nxres)) then
         write (outfi, 6000) 'ATOM  ', mcser(calph,nxres), ' CA ', 
     +                        wkres, seqbcd(nxres), 
     +                       (mcaaat(i,calph,nxres), i=1,ncoord), 
     +                        mcocc(calph,nxres), mcthrm(calph,nxres), 
     +                        chstop(nxres), chstrt(nxres)
         nathet = nathet + 1
      end if
      if (nline .gt. 0) then
         do 4500 ii = istart, nline
         if (extres(ii) .eq. nxres) then
            write (outfi, 6005) extdat(ii)
            if (extdat(ii)(1:keylen) .eq. atmkey .or. 
     +          extdat(ii)(1:keylen) .eq. htakey) nathet = nathet + 1
         else
            if (extres(ii) .gt. nxres) then
               istart = ii
               goto 4510
            end if
         end if
 4500    continue
      end if
 4510 continue
 5000 continue
CHECK v.3.5-->
      IF (.NOT.NMR) THEN
CHECK v.3.5<--
          write (outfi, 6006) nathet
          if (nathet .ne. nrec) then
CHECK v.2.1.3-->
C         write (logfi, 6100) nathet, nrec
             write (*, 6100) nathet, nrec
CHECK v.2.1.3<--
CHECK v.1.0-->
C         write (*,     6100) nathet, nrec
CHECK v.1.0<--
          end if
          write (outfi, 6007) 
CHECK v.3.5-->
      ENDIF
CHECK v.3.5<--
c     call msgprt
      end
c
c*******************************************************************************
c
      real function dihed (angnum, atoma, atomb, atomc, atomd)
c     routine to calculate the dihedral angle between the 4 atoms
c     given
c
      include 'brkcln.par'
c
      real    atoma(ncoord), atomb(ncoord), atomc(ncoord), atomd(ncoord)
c
      integer angnum
c
      real    dihatm(ncoord,dihdat), codist(ncoord,dihdat-1), 
     +        atmdst(dihdat-1), dotprd(dihdat-1,dihdat-1), ang, 
     +        cosang, detant, atdist
c
      integer i, j, k
c
      logical xeqy
c
      do 100 i = 1, ncoord
      dihatm(i,1) = atoma(i)
      dihatm(i,2) = atomb(i)
      dihatm(i,3) = atomc(i)
      dihatm(i,4) = atomd(i)
CHECK v.3.2-->
      do 99 j=5,dihdat
        dihatm (i,j) = 0.0
 99   continue
CHECK v.3.2<--
  100 continue
      do 200 i = 1, dihdat-1
      atmdst(i) = atdist(dihatm(1,i),dihatm(1,i+1))
      do 150 j = 1, ncoord
      codist(j,i) = dihatm(j,i+1) - dihatm(j,i)
  150 continue
  200 continue
      do 300 i = 1, dihdat-1
      do 300 j = 1, dihdat-1
      dotprd(i,j) = 0.0
  300 continue
      do 400 i = 1, dihdat-2
      do 400 j = i+1, dihdat-1
      do 400 k = 1, ncoord
      dotprd(i,j) = dotprd(i,j) + (codist(k,i) * codist(k,j))
  400 continue
      do 500 i = 1, dihdat-2
      do 500 j = i+1, dihdat-1
      if (xeqy(atmdst(i),0.0) .or. xeqy(atmdst(j),0.0)) then
         dotprd(i,j) = 1.0
      else
         dotprd(i,j) = ( dotprd(i,j) / atmdst(i) ) / atmdst(j)
      end if
  500 continue
      if (angnum .gt. ndihed .and. angnum .le. nmnang) then
         cosang = dotprd(1,3)
      else
         if (xeqy(dotprd(1,2),1.0) .or. xeqy(dotprd(2,3),1.0)) then
            cosang = 1.0
         else
            cosang = ((dotprd(1,2) * dotprd(2,3)) - dotprd(1,3)) / 
     +               (sqrt(1.0 - (dotprd(1,2) ** 2)) * 
     +               sqrt(1.0 - (dotprd(2,3) ** 2)))
         end if
      end if
      if (abs(cosang) .gt. 1.0) cosang = anint(cosang / abs(cosang))
      ang = acos(cosang) * radian
      if (angnum .le. ndihed .or. angnum .gt. nmnang) then
         detant = ((((((codist(1,1) * codist(2,2)) * codist(3,3)) - 
     +                ((codist(1,1) * codist(2,3)) * codist(3,2))) + 
     +                ((codist(1,2) * codist(2,3)) * codist(3,1))) - 
     +                ((codist(1,2) * codist(2,1)) * codist(3,3))) + 
     +                ((codist(1,3) * codist(2,1)) * codist(3,2))) -
     +                ((codist(1,3) * codist(2,2)) * codist(3,1))
         if (detant .lt. 0.0) ang = - ang
      end if
      dihed = ang
      end
c
c*******************************************************************************
c
c     subroutine msgprt()
c     write the error messages to the output file
c
c     include 'brkcln.par'
c
c     character wkrec*80
c
c     integer   errflg
c
c 300 format (a)
c
c     rewind (errfi) 
c     errflg = 0
c 200 continue
c     read (errfi, 300, end=1000) wkrec
c     write (*, *) wkrec
c     goto 200
c1000 continue
c     close (errfi) 
c     end
c
c*******************************************************************************
c
CHECK v.3.2-->
C      subroutine mkangl (mcaaat, mcangs, atomin, chnsz,  chnct,  caonly, 
C     +                   seqlen, scaaat, scatin, seqcod, aacod3, hetcod, 
C     +                   sx2,    sx,     ntot,   hist,   hismin, ave, 
C     +                   sd)
      subroutine mkangl (mcaaat, mcangs, atomin, chnsz,  chnct,  caonly, 
     +                   scaaat, scatin, seqcod, aacod3, hetcod, 
     +                   sx2,    sx,     ntot,   hist,   hismin, ave, 
     +                   sd)
CHECK v.3.2<--
c     Procedure to calculate the main chain dihedral angles
c     if all the atoms are present.
c
      include 'brkcln.par'
c
CHECK v.3.5-->
C      character aacod3*(3*abetsz), hetcod*63, wkres*3
      character aacod3*(3*abetsz), hetcod*(3*MXHTYP), wkres*3
CHECK v.3.5<--
c
      real      mcaaat(ncoord,mnchna,maxres), mcangs(nmnang,maxres), 
     +          scaaat(ncoord,sdchna,maxres), sx2, sx
c
      logical   atomin(mnchna,maxres), caonly, scatin(sdchna,maxres)
c
CHECK v.3.2-->
C      integer   seqlen, chnsz(maxchn), chnct, seqcod(maxres), ntot, 
C     +          hist(17)
      integer   chnsz(maxchn), chnct, seqcod(maxres), ntot, 
     +          hist(17)
CHECK v.3.2<--
c
c     local variables
c
      real      dihed
c
      integer   angrng(3), angatm(4), angoff(4), angchn, i, nxres, 
     +          chnst, nxchn, hismin(17)
c
c     integer   j, stang, finang
c
      real      ave, sd
c
      save angchn, angoff, angatm, angrng
c
      data angrng(1) / 1 /
      data angrng(2) / 0 /
      data angrng(3) / impld /
      data angatm(1) / calph /
      data angatm(2) / nmain /
      data angatm(3) / carb /
      data angatm(4) / cbeta /
      data angoff(1) / 0 /
      data angoff(2) / 0 /
      data angoff(3) / 0 /
      data angoff(4) / 0 /
      data angchn / 1 /
c
 2000 format (' Average value of CA-N-C-CB angle is ',F6.2,/,
     +        ' Standard deviation is               ',F6.2)
CHECK v.3.0.1-->
C 2001 format (' Unusual CA-N-C-CB torsion angle ', f5.1, ' residue ', 
C     + 1x, i4, 1x, a3)
CHECK v.3.0.1<--
c
      if (.not. caonly) then
         chnst = 0
         do 3000 nxchn = 1,chnct
         if (chnsz(nxchn) .lt. angchn) then
            do 1000 nxres = chnst+1, chnst+chnsz(nxchn)
            mcangs(angrng(3),1) = null
 1000       continue
         else
            do 2100 nxres = chnst+angrng(1), 
     +                      chnst+chnsz(nxchn)+angrng(2)
            if (atomin(angatm(1),nxres+angoff(1)) .and. 
     +          atomin(angatm(2),nxres+angoff(2)) .and. 
     +          atomin(angatm(3),nxres+angoff(3)) .and. 
     +          scatin(angatm(4),nxres+angoff(4))) then
               mcangs(angrng(3),nxres) = dihed(1,
     +                              mcaaat(1,angatm(1),nxres+angoff(1)),
     +                              mcaaat(1,angatm(2),nxres+angoff(2)),
     +                              mcaaat(1,angatm(3),nxres+angoff(3)),
     +                              scaaat(1,angatm(4),nxres+angoff(4)))
               do 2005 i = 1, 17
               if ( mcangs(1,nxres) .ge. real(hismin(i)) .and. 
     +              mcangs(1,nxres) .lt. real(hismin(i) + 10) ) then
                  hist(i) = hist(i) + 1
                  goto 2010
               end if
 2005          continue
 2010          continue
               if (mcangs(1,nxres) .lt. 23.0 .or. 
     +             mcangs(1,nxres). gt. 45.0) then
                  if (seqcod(nxres) .le. abetsz) then
                     wkres = aacod3((seqcod(nxres)*3)-2:
     +                               seqcod(nxres)*3)
                  else
                     wkres = hetcod(((seqcod(nxres)-abetsz)*3)-2:
     +                               (seqcod(nxres)-abetsz)*3)
                  end if

CHECK v.1.0-->
C                  write (logfi, 2001) mcangs(1,nxres), nxres, wkres
CHECK v.1.0<--

               end if
            else
               mcangs(angrng(3),nxres) = null
            end if
 2100       continue
         end if
         do 2400 nxres = chnst+1, chnst+angrng(1)-1
         mcangs(angrng(3),nxres) = null
 2400    continue
         do 2600 nxres = chnst+chnsz(nxchn)+angrng(2)+1, 
     +                   chnst+chnsz(nxchn)
         mcangs(angrng(3),nxres) = null
 2600    continue
         chnst = chnst + chnsz(nxchn)
 3000    continue
         call avevar (mcangs, chnst, ave, sd)
CHECK v.2.1.3-->
C         write (logfi, 2000) ave, sd
         write (*, 2000) ave, sd
CHECK v.2.1.3<--
         call totave (mcangs, chnst, sx2, sx, ntot)
      end if
      end
c
c*******************************************************************************
c
      subroutine avevar (mcangs, chnst, ave, sd)
c
      include 'brkcln.par'
c
      real    mcangs(nmnang,maxres), ave, sd, var, s
c
      integer j, n, chnst
c
      logical xney
c
      ave = 0.0
      sd  = 0.0
      var = 0.0
      n   = 0
      do 11 j = 1, chnst
      if (xney(mcangs(1,j),null)) then
         ave = ave + mcangs(1,j)
         n = n + 1
      end if
   11 continue
      if (n .gt. 0) then
         ave = ave / real(n)
      else
         ave = 0.0
      end if
      do 12 j = 1, chnst
      if (xney(mcangs(1,j),null)) then
         s = mcangs(1,j) - ave
         var = var + (s*s)
      end if
   12 continue
      if (n .gt. 1) then
         var = var / real(n-1)
         sd = sqrt(var)
      else
         var = 0.0
         sd = 0.0
      end if
      end
c
c*******************************************************************************
c
      subroutine totave (mcangs, chnst, sx2, sx, ntot)
c
      include 'brkcln.par'
c
      real     mcangs(nmnang,maxres), sx2, sx
c
c     real     ave
c
      integer  j, ntot, chnst
c
      logical xney
c
      sx  = 0.0
      sx2 = 0.0
      do 11 j = 1, chnst
      if (xney(mcangs(1,j),null)) then
         sx2 = sx2 + (mcangs(1,j) * mcangs(1,j))
         sx = sx + mcangs(1,j)
         ntot = ntot + 1
      end if
   11 continue
      end
c
c*******************************************************************************
c
c     subroutine tsdev (sx2, sx, ntot)
c     computes standard deviation and mean from total x and total x**
c     values
c     implicit none
c     implicit complex(a-z)
c
c     real     sx2, sx, mean, var, totsd, x
c
c     integer  ntot
c
c     x = (sx ** 2) / real(ntot)
c     var = (sx2 - x) / (ntot - 1)
c     totsd = sqrt(var)
c     mean = sx / ntot
c     write (*, *) ' overall mean value of ca-n-c-cb angle is', 
c    +      mean
c     write (*, *) 
c    +   ' overall standard deviation of this angle is', totsd
c     write (*, *) ' number of angle recorded', ntot
c     end
c
c*******************************************************************************
c
      logical  function xeqy (x, y)
c
c     implicit none
c     implicit complex (a-z)
c
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xeqy = (abs(x-y) .lt. acurcy) 
      return
      end
c
c*******************************************************************************
c
      logical  function xney (x, y)
c
c     implicit none
c     implicit complex (a-z)
c
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are not "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xney = (abs(x-y) .gt. acurcy) 
      return
      end
c
c*******************************************************************************
c 
c 
      subroutine getstr (prompt, defans, string, ilen, ifail)
c
c     implicit complex(a-z)
c     implicit none
c===============================================================================
c     General purpose routine to prompt the user to type in a character
c     string. The prompt that will be written out to the user should be passed
c     to the routine in variable PROMPT, along with a default response in
c     DEFANS. The variable STRING will return the user's input to the calling
c     program, and ILEN will hold its length. The flag IFAIL will have one of
c     the following values:
c        IFAIL = 0 : everything OK;
c        IFAIL = 1 : variable PROMPT is empty;
c        IFAIL = 2 : variable DEFANS is empty;
c        IFAIL = 3 : user's response is too long to fit in STRING;
c        IFAIL = 4 : error occurred writing out prompt to terminal;
c        IFAIL = 5 : error occurred reading in reply from keyboard;
c        IFAIL = 6 : EOF detected reading in reply from keyboard.
c===============================================================================
      character*(*) prompt, defans, string
      integer       ilen, ifail
c
      character     line*255
      integer       lenp, lens, strlen, lend
c
 1000 format (A, ' <', A, '> - ')
 1010 format (A)
c
      ilen  = 0
      ifail = 0
c
      lenp  = strlen(prompt)
      if (lenp .le. 0)  then
         ifail = 1
         return
      end if
c
      lend  = strlen(defans)
      if (lend .le. 0)  then
         ifail = 2
         return
      end if
c
      lens  = len(string)
c
CHECK v.3.0.1-->
C   10 write (6, 1000, err = 904) prompt(1:lenp), defans(1:lend)
      write (6, 1000, err = 904) prompt(1:lenp), defans(1:lend)
CHECK v.3.0.1<--
      read  (5, 1010, err = 905, end = 906) line
      ilen = strlen(line)
      if (ilen .eq. 0) then
         line = defans(1:lend)
         ilen = lend
      end if
      if (ilen .gt. lens)  then
         ifail = 3
         return
      end if
      string = line(1:ilen)
      ifail = 0
      return
c
  904 ifail = 4
      return
  905 ifail = 5
      return
  906 ifail = 6
      return
c
      end 
c
c*******************************************************************************
c
      integer function lastsl (string)
c
c     implicit complex(a-z)
c     implicit none
c
      character *(*) string
      integer   i, lens
c
      lens = len(string)
      do 10 i = lens, 1, -1
      if (string(i:i).eq.'/')  then
         lastsl = i
         return
      end if
   10 continue
      lastsl = 0
      return
      end
c
c*****************************************************************************

CHECK v.2.0-->
C*****************************************************************************
C
C  SUBROUTINE GETNAM  -  Peel off the directory path and extension from the
C                        full name of the .pdb file
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE GETNAM(PDBFIL,ISTART,IEND,IERROR)

      CHARACTER*1   PCHAR
      CHARACTER*78  PDBFIL
      INTEGER       IEND, IPOS, ISTART, ISTATE
      LOGICAL       FINISH, GOTDOT, IERROR

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
              IF (PCHAR.EQ.'/' .OR. PCHAR.EQ.'\\' .OR.
     -            PCHAR.EQ.']') THEN
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
              ELSE IF (PCHAR.EQ.'/' .OR. PCHAR.EQ.'\\'
     -            .OR. PCHAR.EQ.']') THEN
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

C---------------------------------------------------------------------------
CHECK v.2.0<--

CHECK v.3.1-->
C*****************************************************************************
C
C  SUBROUTINE GETMOD  -  Read through the PDB file until hit the next
C                        MODEL record
C
C----------------------------------------------------------------------+--- 

CHECK v.3.2-->
C      SUBROUTINE GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR)
CHECK v.3.5-->
C      SUBROUTINE GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR)
      SUBROUTINE GETMOD(BRKFI,OUTFI,IMODEL,MXMODL,ENDFIL,IERROR,NMR,
     -    BRCODE,AAPLN,AAIN,INSEQR,MAXTYP,SQRSCT,INSLEN,AUTHFI,ASTRUC)
CHECK v.3.5<--
CHECK v.3.2<--

      SAVE

      CHARACTER*80  IREC
      INTEGER       BRKFI, IMODEL, LINE, MXMODL, OUTFI
CHECK v.3.2-->
C      LOGICAL       ENDFIL, FOUND, IERROR
      LOGICAL       ENDFIL, FOUND, IERROR, NMR
CHECK v.3.2<--
CHECK v.3.5-->
      INTEGER       AAPLN

      CHARACTER*1   CONTCH
      CHARACTER*3   AAIN(AAPLN)
      CHARACTER*4   BRCODE, INCODE
      CHARACTER*6   KEYWD
      CHARACTER*50  CLASS, COMPND, SOURCE, WKSTR
      CHARACTER*(*) INSEQR
      INTEGER       AUTHFI, CLSLEN, CMPLEN, I, INSLEN, MAXTYP, RESFND,
     -              SOULEN, SQRSCT, STRLEN
      LOGICAL       ASTRUC
CHECK v.3.5<--

C---- Initialise variables
      ENDFIL = .FALSE.
      FOUND = .FALSE.
      IERROR = .FALSE.
      LINE = 0

C---- Loop through the file until hit a MODEL record or the EOF
 100  CONTINUE

C----     Read in the next record
          READ(BRKFI,120,ERR=900,END=500) IREC
 120      FORMAT(A)

C----     Write it out to the .new file
CHECK v.3.3.3-->
CHECK v.3.5-->
C          IF (IREC(1:6).NE.'ATOM  ') THEN
          KEYWD = IREC(1:6)
          IF (KEYWD.NE.'ATOM  ' .AND. KEYWD.NE.'HETATM') THEN
CHECK v.3.5<--
CHECK v.3.3.3<--
              WRITE(OUTFI,120) IREC
CHECK v.3.3.3-->
          ENDIF
CHECK v.3.3.3<--

C----     Test if this is a MODEL record
          IF (IREC(1:5).EQ.'MODEL') THEN
              LINE = LINE + 1
              READ(IREC,140,ERR=902) IMODEL
 140          FORMAT(11X,I3)
              FOUND = .TRUE.
              IF (IMODEL.GT.MXMODL) GO TO 904
CHECK v.3.2-->
              NMR = .TRUE.

CHECK v.3.5-->
C----     HELIX, SHEET, TURN or SSBOND record: Store the authors'
C----     secondary structure assignments, then write out to
C         temporary file
          ELSE IF (KEYWD.EQ.'HELIX ' .OR. KEYWD.EQ.'SHEET ' .OR. 
     -        KEYWD.EQ.'TURN  ' .OR. KEYWD.EQ.'SSBOND') THEN
              IF (.NOT. ASTRUC) THEN
                  OPEN (AUTHFI, STATUS='SCRATCH') 
                  ASTRUC = .TRUE.
              ENDIF
              WRITE(AUTHFI,120) IREC

C----     SEQRES record: store the details
          ELSE IF (KEYWD.EQ.'SEQRES') THEN
              READ(IREC,220) (AAIN(I), I = 1, AAPLN)
 220          FORMAT(18X,13(1X,A3))
              DO 250 I = 1, AAPLN
                  IF (AAIN(I).NE.'   ') THEN
                      IF (RESFND(INSEQR,AAIN(I)).EQ.0) THEN
                          SQRSCT = SQRSCT + 1
                          IF (SQRSCT.GT.MAXTYP) THEN
                              PRINT*, '*** Warning. Maximum number of ',
     -                            'different types of residue in SEQRE',
     -                            'S records, MAXTYP, exceeded', MAXTYP
                              SQRSCT = MAXTYP
                          ENDIF
                          INSEQR(((SQRSCT-1)*3)+1:SQRSCT*3) = AAIN(I)
                      ENDIF
                   ENDIF
 250          CONTINUE

C----     HEADER record: store just the first record
          ELSE IF (KEYWD.EQ.'HEADER') THEN
              READ(IREC,520) CONTCH, WKSTR, INCODE
 520          FORMAT(9X,A1,A40,12X,A4)
              IF (CONTCH.EQ.' ') THEN
                  BRCODE = INCODE
                  CLASS = WKSTR
                  CLSLEN = STRLEN(CLASS)
              ENDIF

C----     COMPND record: store just the first record
          ELSE IF (KEYWD.EQ.'COMPND') THEN
              READ(IREC,620) CONTCH, WKSTR
 620          FORMAT(9X,A1,A50)
              IF (CONTCH.EQ.' ') THEN
                  COMPND = WKSTR
                  CMPLEN = STRLEN(COMPND)
              END IF
c
C----     SOURCE record: store just the first record
          ELSE IF (KEYWD.EQ.'SOURCE') THEN
              READ(IREC,720) CONTCH, WKSTR
 720          FORMAT(9X,A1,A50)
              IF (CONTCH.EQ.' ') THEN
                  SOURCE = WKSTR
                  SOULEN = STRLEN(SOURCE)
              ENDIF
CHECK v.3.5<--

C----     If this is an ATOM or HETATM record, then can start processing
C         the structure
CHECK v.3.5-->
C          ELSE IF (IREC(1:6).EQ.'ATOM  ') THEN
          ELSE IF (KEYWD.EQ.'ATOM  ' .OR. KEYWD.EQ.'HETATM') THEN
CHECK v.3.5<--
              IMODEL = 1
              FOUND = .TRUE.
              BACKSPACE(BRKFI)
CHECK v.3.2<--
          ENDIF

C---- Loop back for next record
      IF (.NOT.FOUND) GO TO 100
      GO TO 999

C---- End of file reached
 500  CONTINUE
      ENDFIL = .TRUE.

      GO TO 999

C---- Error in file
900   CONTINUE
      PRINT*, ' *** ERROR reading NMR file while looking for MODEL',
     -    ' record'
      IF (IMODEL.EQ.0) THEN
          PRINT*, ' *** Error at line', LINE
      ELSE
          PRINT*, ' *** Error at line', LINE, '  after ENDMDL record',
     -        ' for model', IMODEL
      ENDIF
      GO TO 990

 902  CONTINUE
      PRINT*, ' *** Data error in MODEL number'
      IF (IMODEL.EQ.0) THEN
          PRINT*, ' *** Error at line', LINE
      ELSE
          PRINT*, ' *** Error at line', LINE, '  after ENDMDL record',
     -        ' for model', IMODEL
      ENDIF
      GO TO 990

 904  CONTINUE
      PRINT*, ' *** Maximum number of models exceeded'
      PRINT*, ' *** Model found:', IMODEL, '     Maxmimum allowed:',
     -    MXMODL
      GO TO 990


 990  CONTINUE
      IERROR = .TRUE.

999   CONTINUE
      RETURN
      END

C---------------------------------------------------------------------------
CHECK v.3.1<--
CHECK v.3.5-->
C*****************************************************************************
C
C  SUBROUTINE STOEXT  -  Store the current line read in from the PDB file
C                        in the list of "extra" residue details
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE STOEXT(EXTDAT,MAXLIN,NLINE,EXTRES,AACNT,INREC,NREQRD)

      SAVE

      INTEGER       MAXLIN

      CHARACTER*66  EXTDAT(MAXLIN), INREC
      INTEGER       AACNT, EXTRES(MAXLIN), NLINE, NREQRD
      LOGICAL       WARNED

      DATA WARNED / .FALSE. /

C---- Increment line count
      NLINE = NLINE + 1
      NREQRD = NREQRD + 1

C---- Check that maximum not exceeded
      IF (NLINE.GT.MAXLIN) THEN

C----     If haven't yet issued warning, then do so now
          IF (.NOT.WARNED) THEN
              PRINT*, '*** Warning. Maximum number of extra residue de',
     -            'tails exceeded. Some data lost!'
              WARNED = .TRUE.
          ENDIF
          NLINE = MAXLIN

C---- If can store this line, then do so
      ELSE
          EXTDAT(NLINE) = INREC
          EXTRES(NLINE) = AACNT
      ENDIF

      RETURN
      END

C---------------------------------------------------------------------------
CHECK v.3.5<--
