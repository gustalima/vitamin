from rpy import *        
import sys      

pdb = sys.argv[1]


rama_GENERAL = "General"
rama_GLYCINE = "Glycine"
rama_PROLINE = "Proline"
rama_PRE_PRO = "Pre-Pro"
ramachandran_types = [rama_GENERAL,rama_GLYCINE,rama_PROLINE,rama_PRE_PRO]

# I have used the same colours as RAMPAGE
# http://raven.bioc.cam.ac.uk/rampage.php
rama_settings = {"General" : ([0, 0.0005, 0.02, 1],
                              ['#FFFFFF','#B3E8FF','#7FD9FF'],
                              "pct/rama/rama500-general.data"),
                              # or rama500-general-nosec.data
                 "Glycine" : ([0, 0.002,  0.02, 1],
                              ['#FFFFFF','#FFE8C5','#FFCC7F'],
                              "pct/rama/rama500-gly-sym.data"),
                              # or rama500-gly-sym-nosec.data
                 "Proline" : ([0, 0.002,  0.02, 1],
                              ['#FFFFFF','#D0FFC5','#7FFF8C'],
                              "pct/rama/rama500-pro.data"),
                 "Pre-Pro" : ([0, 0.002,  0.02, 1],
                              ['#FFFFFF','#B3E8FF','#7FD9FF'],
                              "pct/rama/rama500-prepro.data")}
                 #P.S. Also rama500-ala-nosec.data

import numpy as Numeric
def load_data_file(filename) :
    STEP=2
    HALF_STEP=1
    STEP = HALF_STEP*2
    lower_bounds = range(-180, 180, STEP)
    mid_points = range(-180+HALF_STEP, 180+HALF_STEP, STEP)
    upper_bounds = range(-180+STEP, 180+STEP, STEP)

    data = Numeric.array([[0.0 for x in mid_points] for y in mid_points],
                         Numeric.dtype('<f'))

    """
    # Table name/description: "Top500 General case (not Gly, Pro, or pre-Pro) B<30"
    # Number of dimensions: 2
    # For each dimension, 1 to 2: lower_bound  upper_bound  number_of_bins  wrapping
    #   x1: -180.0 180.0 180 true
    #   x2: -180.0 180.0 180 true
    # List of table coordinates and values. (Value is last number on each line.)
    -179.0 -179.0 0.0918642445114388
    -179.0 -177.0 0.07105717866463215
    ...
    """
    input_file = open(filename,"r")
    for line in input_file :
        #Strip the newline character(s) from the end of the line
        if line[-1]=="\n" : line = line[:-1]
        if line[-1]=="\r" : line = line[:-1]
        if line[0]=="#" :
            #comment
            pass
        else :
            #data
            parts = line.split()
            assert len(parts)==3
            
            x1 = float(parts[0]) #phi
            x2 = float(parts[1]) #psi
            value = float(parts[2])
            
            assert x1 == float(int(x1))
            assert x2 == float(int(x2))
            i1 = mid_points.index(int(x1))
            i2 = mid_points.index(int(x2))
            
            data[i1,i2]=value
    input_file.close()
    return (data, lower_bounds, mid_points, upper_bounds)

#filename = "stat/rama/rama500-general.data"
#data, lower_bounds, mid_points, upper_bounds = load_data_file(filename)
#print sum(sum(data))


r.library("MASS")


print "Creating R function",
r("""
ramachandran.plot <- function(x.scatter, y.scatter,
    x.grid = seq(0, 1, len = nrow(z)), y.grid = seq(0, 1, len = ncol(z)), z.grid,
    xlim = range(x.grid, finite = TRUE), ylim = range(y.grid, finite = TRUE),
    zlim = range(z.grid, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) -
        1), plot.title="", plot.axes, key.title, key.axes, asp = NA,
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes,
    ...)
{
    if (missing(z.grid)) {
        stop("no 'z.grid' matrix specified")
    }
    else if (is.list(x.grid)) {
        y.grid <- x.grid$y
        x.grid <- x.grid$x
    }
    if (any(diff(x.grid) <= 0) || any(diff(y.grid) <= 0))
        stop("increasing 'x.grid' and 'y.grid' values expected")

    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)

    if (!is.matrix(z.grid) || nrow(z.grid) <= 1 || ncol(z.grid) <= 1)
        stop("no proper 'z.grid' matrix specified")
    if (!is.double(z.grid))
        storage.mode(z.grid) <- "double"
    .filled.contour(as.double(x.grid), as.double(y.grid), z.grid, as.double(levels),
        col = col)

    if (!(missing(x.scatter)) && !(missing(y.scatter))) {
        plot.xy(xy.coords(x.scatter,y.scatter,NULL,NULL,NULL,NULL),
                xlim=xlim, ylim=ylim, xlab="Phi", ylab="Psi", asp=asp,
                type="p", pch=20, cex=0.1)
    }
        
    if (missing(plot.axes)) {
        if (axes) {
            title(main=plot.title, xlab=expression(phi), ylab=expression(psi))
            axis(1, at=c(-180,-90,0,90,180))
            axis(2, at=c(-180,-90,0,90,180))
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}
""")
print "Done"


import math
def degrees(rad_angle) :
    """Converts and angle in radians to degrees, mapped to the range [-180,180]"""
    angle = rad_angle * 180 / math.pi
    #Note this assume the radians angle is positive as that's what MMTK does
    while angle > 180 :
        angle = angle - 360
    return angle

def next_residue(residue) :
    """Expects an MMTK residue, returns the next residue in the chain, or None"""
    #Proteins go N terminal --> C terminal
    #The next reside is bonded to the C of this atom...
    for a in residue.peptide.C.bondedTo():
        if a.parent.parent != residue:
            return a.parent.parent
    return None


def residue_amino(residue) :
    """Expects an MMTK residue, returns the three letter amino acid code in upper case"""
    if residue :
        return residue.name[0:3].upper()
    else :
        return None

def residue_ramachandran_type(residue) :
    """Expects an MMTK residue, returns ramachandran 'type' (General, Glycine, Proline or Pre-Pro)"""
    if residue_amino(residue)=="GLY" :
        return rama_GLYCINE
    elif residue_amino(residue)=="PRO" :
        return rama_PROLINE
    elif residue_amino(next_residue(residue))=="PRO" :
        #exlcudes those that are Pro or Gly
        return rama_PRE_PRO
    else :
        return rama_GENERAL

scatter_phi = dict()
scatter_psi = dict()
for ramachandran_type in ramachandran_types :
    scatter_phi[ramachandran_type]=[]
    scatter_psi[ramachandran_type]=[]
    

pdb
pdb_filename = "%s.pdb" % pdb

print "Loading PDB file: " + pdb_filename
import MMTK.PDB
import MMTK.Proteins
#protein = MMTK.Proteins.Protein("1HMP.pdb", model="no_hydrogens")
# Load the PDB file, ignore the hydrogrens, and then build a model of the peptides:
configuration = MMTK.PDB.PDBConfiguration(pdb_filename)
configuration.deleteHydrogens()
protein = MMTK.Proteins.Protein(configuration.createPeptideChains(model = "no_hydrogens"))
for chain in protein :
    print chain.name
    for residue in chain :
        phi, psi = residue.phiPsi()
        #print residue.name, phi, psi
        if phi and psi :
            ramachandran_type = residue_ramachandran_type(residue)
            assert ramachandran_type in ramachandran_types
            scatter_phi[ramachandran_type].append(degrees(phi))
            scatter_psi[ramachandran_type].append(degrees(psi))
        assert len(scatter_phi) == len(scatter_psi)
    
print "Done"

#pdf_filename = "../%s.pdf" % pdb
#r.pdf(pdf_filename)

#png_filename = "../%s.png" % pdb
#r.png(png_filename)

#To get four plots on one page, you could use :
#
#r.split_screen([2,2]) #split into two by two screen
#
#Or:
#
#r.layout(Numeric.array([[1,2],[3,4]]), respect=True)
#
#But I went for simply:

#r.par(mfrow=[2,2])
types_rama = ['General', 'Glycine', 'Proline', 'Pre-Pro']





for (i,ramachandran_type) in enumerate(ramachandran_types) :
    png_filename = "%s_%s.png"% (pdb,types_rama[i])
    r.png(png_filename, width=800, height=800, res=250, pointsize = 4)
    #pdf_filename = "../%s_%s.pdf" % (pdb, ramachandran_type)
    (rama_levels, rama_colors, rama_filename) = rama_settings[ramachandran_type]

    print "Loading data file: " + rama_filename,
    data, lower_bounds, mid_points, upper_bounds = load_data_file(rama_filename)
    print "Done"

    #print "Creating PDF output file: " + pdf_filename,
    #r.pdf(pdf_filename)
    #r.plot(scatter_phi, scatter_psi)

    print "Generating quadrant %i, %s" % (i+1, ramachandran_type)
    #r.screen(i+1)

    #Use small margins to make the plots nice and big,
    #and specify a SQUARE plot area (to go with aspect ratio, asp=1)
    r.par(mar = [4, 4, 2, 1], pty="s")

    #This function will do a Ramachandran plot in the next quadrant
    #which we setup using par(mfrow-...)
    r.ramachandran_plot(x_scatter=scatter_phi[ramachandran_type],
                        y_scatter=scatter_psi[ramachandran_type], 
                        x_grid=mid_points, y_grid=mid_points, z_grid=data,
                        xlim=[-180,180], ylim=[-180,180], asp=1.0,
                        plot_title=ramachandran_type, drawlabels=True,
                        levels=rama_levels, col=rama_colors)
    print ramachandran_type + " Done"



    
    
r.dev_off()
print "Done"


