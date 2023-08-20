"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""
import sys
import os
import numpy as np
import datetime

from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgmax_tools
from clawpack.geoclaw import fgout_tools
#from clawpack.geoclaw.data import ForceDry
  

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")


try:
    COPES = os.environ['COPES']
except:
    raise Exception("*** Set COPES enviornment variable to repository top")

# Set these directories where input data is found: 

root_dir = COPES
topo_dir = root_dir + '/topo/topofiles'
dtopo_dir = root_dir + '/dtopo/dtopofiles'

# load list of virtual gauges to use, with columns 
#      gaugeno, x, y
# x,y should be in decimal form, preferably cell centered on finest grid
gauges_file = root_dir + '/gauges/VGListSeaside.csv' 
gauges = np.loadtxt(gauges_file, delimiter=',')
    

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:

    # shift so that cell centers on finest grid align with DEM:
    sec16 = 1./(6*3600.)  # one sixth arcsecond, if finest grid is 1/3"

    # Using bigger domain offset by sec16

    clawdata.lower[0] = -137.0 - sec16       # west longitude
    clawdata.upper[0] = -123.0 - sec16       # east longitude

    clawdata.lower[1] = 38.0 - sec16         # south latitude
    clawdata.upper[1] = 54.0 - sec16         # north latitude

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 14*15
    clawdata.num_cells[1] = 16*15

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------

    t0 = 0            # Start at time 0
    clawdata.t0 = t0  # Start time in seconds

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False    # True to restart from prior results
    clawdata.restart_file = ''  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 4
        clawdata.tfinal = 1*3600. 
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        tfinal = 1*3600.
        dtout = 15*60.
        output_times = [0.,60.] + list(np.arange(dtout, tfinal+1, dtout))


    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True


    clawdata.output_format = 'binary'        # 'ascii' or 'binary'

    clawdata.output_q_components = 'all'     # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False     # output aux arrays each frame

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    # clawdata.source_split = 'godunov'
    clawdata.source_split = 1

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2 

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # negative checkpoint_style means alternate between aaaaa and bbbbb files
    # so that at most 2 checkpoint files exist at any time, useful when
    # doing frequent checkpoints of large problems.

    clawdata.checkpt_style = -2

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = 3600.*np.arange(1,6,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 7

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 4', 1', 30", 6", 2", 1" 1/3":
    amrdata.refinement_ratios_x = [4,2,5,3,2,3]
    amrdata.refinement_ratios_y = [4,2,5,3,2,3]
    amrdata.refinement_ratios_t = [4,2,5,3,2,3]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True 

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0

    geo_data.manning_coefficient = .025

    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.friction_depth = 200

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.02

    # ---------------
    # TOPO:
    # ---------------
    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]

    topofiles = topo_data.topofiles

    # 1-minute topo:
    topo_file = os.path.join(topo_dir, 'etopo1_-138_-122_37_55_1min.asc')
    topofiles.append([3, topo_file])

    # 3 sec topo: 
    topo_file = os.path.join(topo_dir, 'crm_vol8_3sec_cropped_for_seaside.asc')
    topofiles.append([3, topo_file])

    #2 sec is level 5, 1 sec is level 6, and 1/3sec is level 7
    topo_file = os.path.join(topo_dir, 'astoria_2s_mhw.asc')
    topofiles.append([3, topo_file])

    # 1/9 arc sec topo, will be used to support 1/3" and 1" calculations:
    # These are the newest tiles along the Oregon Coast. Randy converted
    # them from NAVD88 to MHW using factor 2.225m (MHW is 2.225m higher
    # than NAVD88, so more water above the NAVD88 datum by 2.225m).
    # Also made 1/3 arc sec topo from this to cover a larger region.
    # including the 1 arc sec region

    if 1: #1/3 arc sec
        topo_file = os.path.join(topo_dir, 'SeasideN_13s_mhw.asc')
        topofiles.append([3, topo_file])
        topo_file = os.path.join(topo_dir, 'SeasideS_13s_mhw.asc')
        topofiles.append([3, topo_file])

    if 0: #1/9 arc sec  (taking out for a test)
        topo_file = os.path.join(topo_dir, 'SeasideN_19s_mhw.asc')
        topofiles.append([3, topo_file])
        topo_file = os.path.join(topo_dir, 'SeasideS_19s_mhw.asc')
        topofiles.append([3, topo_file])

    # ---------------
    # DTOPO:
    # ---------------
    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [dtopotype, fname]
    dtopofile = dtopo_dir + '/CSZ_SM1.tt3'
    dtopo_data.dtopofiles = [[3, dtopofile]]
    dtopo_data.dt_max_dtopo = 0.2

    # ---------------
    # qinit:
    # ---------------
    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    rundata.qinit_data.variable_eta_init = True  # for subsidence


    # ---------------
    # Force dry:
    # ---------------
    #None for this project


    # ---------------
    # REGIONS:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    flagregions = rundata.flagregiondata.flagregions  # initialized to []

    # Computational domain Variable Region - 4 minutes to 1 minute:
    # Level 1 is 4 minutes, level 2 is 1 minute
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [clawdata.lower[0]-0.2,
                                 clawdata.upper[0]+0.2,
                                 clawdata.lower[1]-0.2,
                                 clawdata.upper[1]+0.2]
    flagregions.append(flagregion)

    # Region1min - fixed 1 minute:
    # Level 2 is 1 minute
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_1min'
    flagregion.minlevel = 2
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    #flagregion.spatial_region = [-127.1,-123.6,40.0,49.2] 
    flagregion.spatial_region = [-127.1,-123.6,40.0,50.0]  #for bigd
    flagregions.append(flagregion)

    # Region30sec - fixed at 30 sec :
    # Level 3 is 30 sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_30sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    #flagregion.spatial_region = [-126.6,-123.6,45.88,47.31] 
    #flagregion.spatial_region = [-126.6,-123.6,45.20,47.31] 
    #flagregion.spatial_region = [-126.6,-123.6,44.66,47.31] 
    flagregion.spatial_region = [-127.1,-123.6,44.0,47.65] 
    flagregions.append(flagregion)

    # Region6a - fixed at 6 sec :
    # Level 4 is 6 sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_6sec_A'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.58,-123.82,45.21,46.77] 
    flagregions.append(flagregion)

    # Region 2" - fixed at 2" sec, Making bigger than initially:
    # Level 5 is 2" sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Seaside_2sec'
    flagregion.minlevel = 5
    flagregion.maxlevel = 5
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    #flagregion.spatial_region = [-124.1,-123.9025,45.95,46.065] 
    flagregion.spatial_region = [-124.25,-123.9025,45.9025,46.32] 
    flagregions.append(flagregion)

    # Region 1" - fixed at 1" sec (almost old 2" becoming new 1" region):
    # Level 6 is 1" sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Seaside_1sec'
    flagregion.minlevel = 6
    flagregion.maxlevel = 6
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    #flagregion.spatial_region = [-123.9975,-123.9025,45.9625,46.0475] 
    flagregion.spatial_region = [-124.03,-123.9025,45.945,46.08] 
    #flagregion.spatial_region = [-124.0,-123.9025,45.975,46.04] 
    flagregions.append(flagregion)

    # Region 1/3" - fixed at 1/3" sec (old 1" becoming new 1/3" region, west edge slightly changed) :
    # Level 7 is 1/3" sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Seaside_onethird'
    flagregion.minlevel = 7
    flagregion.maxlevel = 7
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-123.9975,-123.9025,45.97,46.0275] 
    flagregions.append(flagregion)

    # ---------------
    # GAUGES:
    # ---------------

    rundata.gaugedata.gauges = []
    for k in range(gauges.shape[0]):
        gaugeno = int(gauges[k,0])
        xg = gauges[k,1]
        yg = gauges[k,2]
        rundata.gaugedata.gauges.append([gaugeno,xg,yg,0,1e9])
    

    # -----------------------------
    # FGMAX GRIDS:
    # NEW STYLE STARTING IN v5.7.0
    # ------------------------------
    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 5

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    if 1:
        fgmax_extent = [-123.94,-123.9025,45.972,46.02]
        # Points on a uniform 2d grid:
        dx_fine = 1./(3*3600.)  # grid resolution at finest level
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = fgmax_extent[0] #+ dx_fine/2.
        fg.x2 = fgmax_extent[1] #- dx_fine/2.
        fg.y1 = fgmax_extent[2] #+ dx_fine/2.
        fg.y2 = fgmax_extent[3] #- dx_fine/2.
        fg.dx = dx_fine
        fg.tstart_max =  t0 + 120   # when to start monitoring max values
        fg.tend_max = 1.e10         # when to stop monitoring max values
        #fg.dt_check = 20.          # target time (sec) increment between updating
        fg.dt_check = 0.            # target time (sec) increment between updating
                                    # max values
                                    # setting to 0. forces updating every time step
                                    # hopefully will now match the max of the gauge at the point.

                                    # which levels to monitor max on
        fg.min_level_check = amrdata.amr_levels_max 
        fg.arrival_tol = 1.e-2      # tolerance for flagging arrival

        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data

    # ======================================
    # fgout grids for frequent output on uniform grids
    # written to fgout_grids.data
    # NEW STYLE STARTING IN v5.9.0

    fgout_grids = rundata.fgout_data.fgout_grids

    if 1:
        # 1/3" grid around Seaside
        fgout_extent = [-123.96,-123.9025,45.972,46.0275]
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 3
        fgout.point_style = 2  # uniform rectangular x-y grid
        fgout.output_format = 'binary32'
        fgout.x1 = fgout_extent[0]
        fgout.x2 = fgout_extent[1]
        fgout.y1 = fgout_extent[2]
        fgout.y2 = fgout_extent[3]
        fgout.nx = int(round((fgout.x2 - fgout.x1)*3*3600.))
        fgout.ny = int(round((fgout.y2 - fgout.y1)*3*3600.))
        fgout.tstart = 10*60.
        fgout.tend = 60*60.
        fgout_dt = 10.  # seconds between fgout frames
        fgout.nout = int(round(((fgout.tend - fgout.tstart)/fgout_dt))) + 1
        fgout_grids.append(fgout)


    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools
    rundata = setrun(*sys.argv[1:])
    rundata.write()

#   To create kml files of inputs: 
    kmltools.make_input_data_kmls(rundata)
