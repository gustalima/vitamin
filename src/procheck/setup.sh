# Insert the following environments and aliases into your .cshrc file

# Note: Edit the strings @@_procheck_program_directory_@@ and
#       @@_aqua_program_directory_@@ below to point to the
#       directories containing the procheck and aqua executables
#       and associated data and parameter files.


# PROCHECK environments and aliases
# ---------------------------------
prodir=@@_procheck_program_directory_@@
export prodir
alias procheck= $prodir'/procheck.scr'
alias procheck_comp=$prodir'/procheck_comp.scr'
alias procheck_nmr=$prodir'/procheck_nmr.scr'
alias proplot=$prodir'/proplot.scr'
alias proplot_comp=$prodir'/proplot_comp.scr'
alias proplot_nmr=$prodir'/proplot_nmr.scr'
alias aquapro=$prodir'/aquapro.scr'
alias gfac2pdb=$prodir'/gfac2pdb.scr'
alias viol2pdb=$prodir'/viol2pdb.scr'
alias wirplot=$prodir'/wirplot.scr'

# AQUA environment and aliases (for use with PROCHECK-NMR)
# --------------------------------------------------------
# Aliases are initialised by typing 'aqua'
#

if [ -z "$aquaroot" ]; then
   aquaroot=@@_aqua_program_directory_@@
   export aquaroot
fi

if [ "`alias aqua`" = "" ]; then
   alias aqua='source $aquaroot/aqsetupi'  # This needs a SH equivalent !!!
fi
