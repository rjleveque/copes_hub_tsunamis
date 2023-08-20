"""
Plot fgout frames and transects
"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot, gridtools
from clawpack.visclaw import animation_tools, colormaps
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta

sys.path.insert(0,'.')
event = 'CSZ_SM1'


if 1:
    from clawpack.geoclaw import fgout_tools
    graphics_dir = os.path.abspath('../graphics')
else:
    # local versions for self-contained directory:
    import fgout_tools
    graphics_dir = './'
    
outdir = os.path.abspath('_output')
fgno = 3

nout = 301
fgframes = range(1,nout+1)
fgframes = [1,100,200,nout]  # test

if 'mmfs1' in outdir:
    # on hyak:
    outdir = outdir.replace('mmfs1/home','gscratch/tsunami')

print('Looking for output in ',outdir)

output_format = 'binary32'

#GE_image = imread(graphics_dir + '/fgout01GE.png')
#GE_extent = [-123.96,-123.9025,45.975,46.0275]
GE_image = imread(graphics_dir + '/seaside_fgout0003GE.png')
GE_extent = [-123.96,-123.9025,45.972,46.0275]

# Instantiate object for reading fgout frames:
fgout_grid1 = fgout_tools.FGoutGrid(fgno, outdir, output_format)


fgframes1 = array(fgframes)
fgframes1 = [int(i) for i in fgframes1]  # convert from numpy.int64 to int

print('fgframes1 = ',fgframes1)

fgframe1 = fgframes1[0] # start with first frame

fgout1 = fgout_grid1.read_frame(fgframe1)


# ----------
# Plotting:

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

# Note that if you change what is plotted for each frame
# (e.g. the transect locations) then you also need to change the
# update function to update with the proper data each frame.

#cmap_water = geoplot.tsunami_colormap  # opaque
alpha = 0.8  # transparency
cmap_water = colormaps.make_colormap({-0.6:[0,0,1,alpha],
                                       0.0:[0,1,1,alpha],
                                       0.6:[1,0,0,alpha]})

fig = figure(1, figsize=(14,9))
clf()


plot_extent = fgout1.extent_edges

ax = axes([.05,.05,.5,.85])
if 1:
    ax.imshow(GE_image, extent=GE_extent)
    B = nan*fgout1.B
else:
    B = fgout1.B
    
B_plot1 = ax.imshow(flipud(B.T), extent=fgout1.extent_edges,
       #cmap=geoplot.land1_colormap)
       cmap=geoplot.googleearth_transparent)

B_plot1.set_clim(0,1500)


eta_water = where(fgout1.h > 0.01, fgout1.eta, nan)
h_water = where(fgout1.h > 0.01, fgout1.h, nan)
eta_water = where(fgout1.B > 0, h_water, eta_water) # for zeta on GE
eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
       cmap=cmap_water)
       


climits = (-20,20)

eta_plot1.set_clim(climits)
axis(plot_extent)
title_text = ax.set_title('%s\nSurface/Depth at time %s  (frame %i)' \
            % (event, timedelta(seconds=fgout1.t), fgframe1))
            
ax.set_aspect(1/cos(46*pi/180))
ticklabel_format(useOffset=False)


if 1:
    cb = colorbar(eta_plot1, extend='both', shrink=0.5,
                  #orientation='vertical',anchor=(0,0))
                  orientation='horizontal', anchor=(0.4,1))
    cb.set_label('meters')
    
# Add transects to planview plot:
yt1 = 46.0015; Ttitle1 = '(12th Ave)'
yt2 = 45.9931; Ttitle2 = '(Broadway)'
yt3 = 45.9894; Ttitle3 = '(Avenue G)'

x1trans, x2trans = -123.95,  -123.91

plot([x1trans,x2trans], [yt1,yt1],'k-',linewidth=0.8)
text(x1trans-0.005,yt1+0.0005,'Transect 1 %s' % Ttitle1, fontsize=8)
plot([x1trans,x2trans], [yt2,yt2],'k-',linewidth=0.8)
text(x1trans-0.005,yt2+0.0005,'Transect 2 %s' % Ttitle2, fontsize=8)
plot([x1trans,x2trans], [yt3,yt3],'k-',linewidth=0.8)
text(x1trans-0.005,yt3+0.0005,'Transect 3 %s' % Ttitle3, fontsize=8)

# =========
# transects:
# =========


def extract_transect(fgout_soln,xtrans,ytrans):

    eta1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                   fgout_soln.eta, xtrans, ytrans)
    B1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                 fgout_soln.B, xtrans, ytrans)
    return B1d, eta1d

def annotate_transect(axtrans):
    dxkm = 1
    dxlong = dxkm/(111.* cos(pi*46/180))
    axtrans.plot([-123.936,-123.936+dxlong], [-12,-12],'k')
    axtrans.text(-123.936+dxlong/2, -13, '%i km' % dxkm, \
                 ha='center',va='top')

    axtrans.grid(True)
    axtrans.ticklabel_format(useOffset=False)
    axtrans.set_xlabel('longitude')
    axtrans.set_ylabel('meters')
    axtrans.set_xlim(x1trans,x2trans)
    xt = axtrans.get_xticks()
    #axtrans.set_xticks(xt,rotation=20)
    axtrans.set_xticks(arange(x1trans,x2trans+1e-6,.01))
    
ylimtr = (-20,20)  # ylimits for transect plots
xtrans = linspace(x1trans, x2trans, 1000)  # x points on transects

# Transect 1 (top)

y1trans, y2trans = 2*[yt1]; Ttitle = Ttitle1
ytrans = linspace(y1trans, y2trans, 1000)

axtrans = axes([.55,.7,.35,.2])
axtrans.set_title('Transect 1 at y = %.5f %s' % (y1trans,Ttitle1))

axtrans.set_ylim(ylimtr)

Btrans1, etatrans1 = extract_transect(fgout1,xtrans,ytrans)
#import pdb; pdb.set_trace()

Btrans, etatrans = Btrans1, etatrans1

# filled regions:
Bfill_plot = axtrans.fill_between(xtrans, Btrans-1e4, Btrans, 
                                  color=[.5,1,.5,1])
etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, 
                                  color=[.5,.5,1,1])

# surface and topo plots:
etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')

annotate_transect(axtrans)

# Transect 2 (middle)

y1trans, y2trans = 2*[yt2]; Ttitle = Ttitle2

ytrans = linspace(y1trans, y2trans, 1000)

axtrans2 = axes([.55,.4,.35,.2],sharex=axtrans)
axtrans2.set_title('Transect at y = %.5f %s' % (y1trans,Ttitle2))

#axtrans.set_xlim(x1trans,x2trans)
#axtrans.sharex(ax)

axtrans2.set_ylim(ylimtr)

Btrans2, etatrans2 = extract_transect(fgout1,xtrans,ytrans)

# filled regions:
Bfill_plot2 = axtrans2.fill_between(xtrans, Btrans2-1e4, Btrans2, 
                                  color=[.5,1,.5,1])
etafill_plot2 = axtrans2.fill_between(xtrans, Btrans2, etatrans2, 
                                  color=[.5,.5,1,1])

# surface and topo plots:
etatrans_plot2, = axtrans2.plot(xtrans, etatrans2, 'b')
Btrans_plot2, = axtrans2.plot(xtrans, Btrans2, 'g')

annotate_transect(axtrans2)


# Transect 3 (bottom)

y1trans, y2trans = 2*[yt3]; Ttitle = Ttitle3

ytrans = linspace(y1trans, y2trans, 1000)

axtrans3 = axes([.55,.1,.35,.2],sharex=axtrans)
axtrans3.set_title('Transect at y = %.5f %s' % (y1trans,Ttitle3))

#axtrans.set_xlim(x1trans,x2trans)
#axtrans.sharex(ax)

axtrans3.set_ylim(ylimtr)

Btrans3, etatrans3 = extract_transect(fgout1,xtrans,ytrans)

# filled regions:
Bfill_plot3 = axtrans3.fill_between(xtrans, Btrans3-1e4, Btrans3, 
                                  color=[.5,1,.5,1])
etafill_plot3 = axtrans3.fill_between(xtrans, Btrans3, etatrans3, 
                                  color=[.5,.5,1,1])

# surface and topo plots:
etatrans_plot3, = axtrans3.plot(xtrans, etatrans3, 'b')
Btrans_plot3, = axtrans3.plot(xtrans, Btrans3, 'g')

annotate_transect(axtrans3)


# The artists that will be updated for subsequent frames:

update_artists = (B_plot1, eta_plot1, \
                Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,\
                Bfill_plot2, etafill_plot2, Btrans_plot2, etatrans_plot2,\
                Bfill_plot3, etafill_plot3, Btrans_plot3, etatrans_plot3,\
                title_text)    
                    
figdummy,axdummy = subplots()

def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
        
    fgout1 = fgout_grid1.read_frame(fgframeno)

    print('Updating plot at time %s' % timedelta(seconds=fgout1.t))
    
    # unpack update_artists (must agree with definition above):
    B_plot1, eta_plot1, \
            Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot, \
            Bfill_plot2, etafill_plot2, Btrans_plot2, etatrans_plot2,\
            Bfill_plot3, etafill_plot3, Btrans_plot3, etatrans_plot3,\
            title_text = update_artists
        
    # reset title to current time:
    title_text.set_text('%s\nSurface/Depth at time %s  (frame %i)' \
            % (event,timedelta(seconds=fgout1.t), fgframeno))

    # reset eta and B in plan-view plots to current state:

    eta_water = where(fgout1.h > 0.01, fgout1.eta, nan)
    h_water = where(fgout1.h > 0.01, fgout1.h, nan)
    eta_water = where(fgout1.B > 0, h_water, eta_water) # for zeta on GE
    eta_plot1.set_array(flipud(eta_water.T))
    #B_plot1.set_array(flipud(fgout1.B.T))

    # update transect:
    
    y1trans, y2trans = 2*[yt1]

    xtrans = linspace(x1trans, x2trans, 1000)
    ytrans = linspace(y1trans, y2trans, 1000)
    
    Btrans1, etatrans1 = extract_transect(fgout1,xtrans,ytrans)
    Btrans, etatrans = Btrans1, etatrans1
    Btrans_plot.set_data(xtrans,Btrans)
    etatrans_plot.set_data(xtrans,etatrans)

    #update the PolyCollections for fill_between plots:             
    dummy = axdummy.fill_between(xtrans, Btrans-1e4, Btrans, 
                                      color=[.5,1,.5,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    Bfill_plot.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans, etatrans, 
                                      color=[.5,.5,1,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    etafill_plot.set_paths([dp.vertices])

    # update transect 2:
    y1trans, y2trans = 2*[yt2]

    ytrans = linspace(y1trans, y2trans, 1000)
    
    Btrans2, etatrans2 = extract_transect(fgout1,xtrans,ytrans)
    Btrans_plot2.set_data(xtrans,Btrans2)
    etatrans_plot2.set_data(xtrans,etatrans2)        


    #update the PolyCollections for fill_between plots:             
    dummy = axdummy.fill_between(xtrans, Btrans2-1e4, Btrans2, 
                                      color=[.5,1,.5,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    Bfill_plot2.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans2, etatrans2, 
                                      color=[.5,.5,1,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    etafill_plot2.set_paths([dp.vertices])
                      
                          
    # update transect 3:

    y1trans, y2trans = 2*[yt3]

    ytrans = linspace(y1trans, y2trans, 1000)
    
    Btrans3, etatrans3 = extract_transect(fgout1,xtrans,ytrans)
    Btrans_plot3.set_data(xtrans,Btrans3)
    etatrans_plot3.set_data(xtrans,etatrans3)        


    #update the PolyCollections for fill_between plots:             
    dummy = axdummy.fill_between(xtrans, Btrans3-1e4, Btrans3, 
                                      color=[.5,1,.5,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    Bfill_plot3.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans3, etatrans3, 
                                      color=[.5,.5,1,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    etafill_plot3.set_paths([dp.vertices])
                      
    update_artists = (B_plot1, eta_plot1, \
                    Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,\
                    Bfill_plot2, etafill_plot2, Btrans_plot2, etatrans_plot2,\
                    Bfill_plot3, etafill_plot3, Btrans_plot3, etatrans_plot3,\
                    title_text)
                    
    return update_artists

def plot_fgframe(fgframeno, save_png=False):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
    
    if save_png:
        fname = 'fgout_frame%s.png' % str(fgframeno).zfill(4)
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':
    
    anim = make_anim()
    
    rundir = os.getcwd()
    if 'mmfs1' in rundir:
        # on hyak:
        scrdir = rundir.replace('mmfs1/home','gscratch/tsunami')
        outdir = os.path.join(scrdir, '_output')
    plotdir = outdir.replace('_output','_plots')

    if 0:
        fgout_plotdir = plotdir + '/animations'
        os.system('mkdir -p %s' % fgout_plotdir)
    else:
        fgout_plotdir = '.'

    # Output files:
    name = '%s_animation' % event

    fname_mp4 = os.path.join(fgout_plotdir, name + '.mp4')
    #fname_html = os.path.join(fgout_plotdir, name + '.html')
    fname_html = None
    
    
    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        fname_html = name + '.html'
        animation_tools.make_html(anim, file_name=fname_html, title=name)
