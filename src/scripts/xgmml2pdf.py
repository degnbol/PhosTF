#!/usr/bin/env python3
# Make .pdf export from .xgmml using cytoscape.
# REQUIRES:
# - py4cytoscape installed, e.g. with pip
# - cyREST installed in cytoscape app manager.
# - cytoscape running to be controlled. 
#   It will be interactive so replacing files will trigger a prompt etc.
# USE: ./xgmml2pdf.py INFILE.xgmml OUTFILE.pdf
import sys
import py4cytoscape as p4c
infile, outfile = sys.argv[1:3]
p4c.import_network_from_file(infile)
p4c.export_image(outfile, type="pdf")
