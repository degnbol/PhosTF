#!/usr/bin/env python3
# Make .png export(s) from .xgmml file(s) using cytoscape.
# REQUIRES:
# - py4cytoscape installed, e.g. with pip
# - cyREST installed in cytoscape app manager.
# - cytoscape running to be controlled
# USE: ./xgmml2png.py INFILE.xgmml [INFILE2.xgmml ...]
# RETURN: writes to INFILE.png, ...
import sys, os
import py4cytoscape as p4c

for infile in sys.argv[1:]:
    outfile = os.path.splitext(infile)[0] + '.png'

    for attempt in range(50):
        try:
            p4c.cytoscape_ping()
        except:
            time.sleep(5)
            continue
        else:
            break
        print("Cytoscape not ready")
        exit(1)

    p4c.import_network_from_file(infile)

    # fit zoom to content,
    # then let content fill 90% of view so we have room for edges pertruding from nodes.
    p4c.fit_content()
    p4c.set_network_zoom_bypass(p4c.get_network_zoom() * 0.9)

    p4c.export_image(outfile, overwrite_file=True, type="png", units="cm", width=4, height=4, resolution=150)
    print(outfile)

