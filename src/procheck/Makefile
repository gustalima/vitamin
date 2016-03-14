# A simple Makefile for Procheck which should be a little more portable
# than the old one
#
# Andrew C.R. Martin, University of Reading
#
F77   = f77  # Your FORTRAN compiler
CC    = cc   # Your C compiler
CLIBS = -lm  # Your C maths library
COPTS = -g   # Compile/link options for C
FOPTS = -g   # Compile/link options for FORTRAN

##########################################################################
# Shouldn't need to alter anything below here
##########################################################################

OFILES = anglen.o clean.o pplot.o rmsdev.o tplot.o vplot.o \
         bplot.o gfac2pdb.o ps.o secstr.o viol2pdb.o wirplot.o
PCEXE = anglen tplot pplot bplot clean secstr nb wirplot
NMREXE = clean secstr rmsdev tplot mplot vplot gfac2pdb viol2pdb

# Main rules
# ----------
all : procheck nmr

procheck : $(PCEXE)

nmr : $(NMREXE)

cleanup : 
	\rm -f $(OFILES)

distrib : 
	\rm -f $(OFILES) $(PCEXE) $(NMREXE)

# Individual executables
# ----------------------
anglen : anglen.o
	$(F77) $(FOPTS) -o $@ anglen.o
clean : clean.o
	$(F77) $(FOPTS) -o $@ clean.o
rmsdev : rmsdev.o
	$(F77) $(FOPTS) -o $@ rmsdev.o
secstr : secstr.o
	$(F77) $(FOPTS) -o $@ secstr.o
gfac2pdb : gfac2pdb.o ps.o
	$(F77) $(FOPTS) -o $@ gfac2pdb.o ps.o
pplot : pplot.o ps.o
	$(F77) $(FOPTS) -o $@ pplot.o ps.o
bplot : bplot.o ps.o
	$(F77) $(FOPTS) -o $@ bplot.o ps.o
tplot : tplot.o ps.o
	$(F77) $(FOPTS) -o $@ tplot.o ps.o
mplot : mplot.o ps.o
	$(F77) $(FOPTS) -o $@ mplot.o ps.o
vplot : vplot.o ps.o
	$(F77) $(FOPTS) -o $@ vplot.o ps.o
viol2pdb : viol2pdb.o ps.o
	$(F77) $(FOPTS) -o $@ viol2pdb.o ps.o
wirplot : wirplot.o ps.o
	$(F77) $(FOPTS) -o $@ wirplot.o ps.o
nb : nb.c
	$(CC) $(COPTS) -o nb nb.c $(CLIBS)

# Individual rules for FORTRAN files with .inc files
# --------------------------------------------------
anglen.o : anglen.f anglen.inc
	$(F77) $(FOPTS) -c anglen.f
bplot.o : bplot.f bplot.inc
	$(F77) $(FOPTS) -c bplot.f
mplot.o : mplot.f mplot.inc
	$(F77) $(FOPTS) -c mplot.f
pplot.o : pplot.f pplot.inc
	$(F77) $(FOPTS) -c pplot.f
tplot.o : tplot.f tplot.inc
	$(F77) $(FOPTS) -c tplot.f
vplot.o : vplot.f vplot.inc
	$(F77) $(FOPTS) -c vplot.f
gfac2pdb.o : gfac2pdb.f gfac2pdb.inc
	$(F77) $(FOPTS) -c gfac2pdb.f
rmsdev.o : rmsdev.f rmsdev.inc
	$(F77) $(FOPTS) -c rmsdev.f
viol2pdb.o : viol2pdb.f viol2pdb.inc
	$(F77) $(FOPTS) -c viol2pdb.f
wirplot.o : wirplot.f wirplot.inc
	$(F77) $(FOPTS) -c wirplot.f

# Generic rules
# -------------

.f.o :
	$(F77) $(FOPTS) -c $<


