# -*- coding: utf-8 -*-
## Got this from the princeton instruments home page on 9 Mar 2010.
import os, os.path
import sys
import piUtils


def usage():
  print "must specify a file containing SPE files to convert to FITS"
  print "One SPE file per line"
  print " e.g. "
  print "python "+sys.argv[0]+" spefiles.txt"
  print ""
  print "Or, if you specify 'all' instead of a filelist, then "
  print "all SPE files in this directory will be converted"
  print ""
  print "Returning..."
  

doUsage = False
if len(sys.argv) < 2:
  doUsage = True
elif sys.argv[1].lower() == '-h' or sys.argv[1].lower() == 'help' or sys.argv[1].lower() == '--help':
  doUsage = True

if doUsage:
  usage()
  sys.exit()

# read in the user-specified filelist
filelist = sys.argv[1]

if filelist == 'all':
  # User wants to convert all SPE files in this directory
  # Make a list of all the spe files in the directory to convert
  dirContents = os.listdir(".")
  print dirContents
  spelist = []
  for fname in dirContents:
    junk, ext = os.path.splitext(fname)
    ext = ext.upper()
    if ext == ".SPE":
      spelist.append(fname)

else:
  # verify filelist is a valid file
  filelistValid = os.access(filelist, os.R_OK)
  if not filelistValid:
    print ''
    print "invalid filelist: ["+filelist+"]"
    print ''
    usage()
    sys.exit()
  
  fh = open(filelist, "r")
  spelist = fh.readlines()
  for ii in range(len(spelist)):
    spelist[ii] = spelist[ii].strip()

print "List of SPE files to be converted to FITS:"
print "spelist = "
print spelist

# loop over that list
print 'Looping over SPE files...'
for spe in spelist:
  print " spe = ", spe
  piUtils.speToFits(spe, fitsfile=None, clobber=True)
  #piUtils.speToFits(spe, fitsfile="fitsfilename.fits", clobber=True)
