
# Seaside, OR modeling using CSZ SM1 event

Set environment variable COPES to top level of this repository.

Running the code requires topo and dtopo files that can be obtained via:

    cd $COPES/topo
    python fetch_etopo_csz.py
    bash wget_seaside.sh

    cd $COPES/dtopo
    bash wget_csz.sh

Note that some topo or dtopo files use longitude E and need to be shifted by
subtracting 360 from the xlower value in the header.  For this problem the
topo files are fine but the dtopo file
    $COPES/dtopo/dtopofiles/CSZ_SM1.tt3
must be modified by changing the fourth line from
    2.33000000000000e+02   xlower
to
    -1.27000000000000e+02   xlower
Otherwise the seafloor deformation will be outside of the computational
domain and no tsunami will be generated.

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
