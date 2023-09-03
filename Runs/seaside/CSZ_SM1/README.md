
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

    cd $COPES/Runs/seaside/CSZ_SM1
    make data

will produce .data files and also .kml files to view on Google Earth to see
the extents of topo and dtopo files used, AMR flagregions, gauge locations, etc.

To run the code:

    make .output

Then to make gauge plots:

    cd $COPES/Runs/seaside
    python plot_gauges.py

To make an animation using data in the fgout files,

    python plot_fgout_animate.py
    
To make fgmax plots and kmz file:

    python process_fgmax.py
    
These scripts in `$COPES/Runs/seaside` will create plots and an animation in

    $COPES/Runs/seaside/CSZ_SM1/_plots
    
The scripts can be modified by changing `event` to be
used for additional runs using different dtopo files.

To adapt to a different location, copy the script, change `location`,
and you may also need to make some other modifications to better
suit the new location. (e.g. figsize based on the aspect ratio of the fgmax
plot.

You will also have to redefine what points to consider `onshore`,
since for Seaside we set it to consider any point with x > -123.93 onshore
in order to capture the harbor and rivers also in the plot of "depth".

What we really plot for "depth" is what we call zeta, equal to:

    zeta = h "onshore"
         = h + B0 "offshore"
         
Normally we define

    onshore = where B0 > 0
    
This is all based on the pre-seismic topography B0 (elevation relative to MHW).

We always define the water elevation (relative to MHW) at any time before or
during the tsunami as

    eta = h + B
    
as is standard in tsunami models, and is based on the elevation B at that time
(which possibly includes subsidence or uplift after the earthquake happens),
but using h+B0 for the plots together with h where the land was originally dry
has the advantage that zeta is continuous at the
original shoreline (where B0=0).

zeta has the interpretation of the apparent water elevation in a harbor
or river to an observer standing on shore, since the observer would drop
with co-seismic subsidence along with the land and the initial water
surface before the tsunami arrives. So, for example, if the tsunami
raises the water level in a harbor by 1 meter relative to the docks at
some time (i.e. the depth h has increased by 1 meter at a point where
initially h0+B0=0), then zeta would be 1 meter, as makes sense to an
onshore observer, whereas if the harbor and land all dropped by
B - B0 = -3 meters of subsidence, then eta = h+B = -2 m in the harbor,
which is less useful to know.
