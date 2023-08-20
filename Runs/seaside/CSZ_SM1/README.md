
# Seaside, OR modeling using CSZ SM1 event

Set environment variable COPES to top level of this repository.

Requires topo and dtopo files that can be obtained via:

    cd $COPES/topo
    python fetch_etopo_csz.py
    bash wget_seaside.sh

    cd $COPES/dtopo
    bash wget_csz.sh

Then

    make data

will produce .data files and also .kml files to view on Google Earth to see
the extents of topo and dtopo files used, AMR flagregions, gauge locations, etc.

To run the code:

    make .output

Then to make gauge plots:

    python plot_gauges.py

To make an animation using data in the fgout files,

    plot_fgout03T_animate.py

