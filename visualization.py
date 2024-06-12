import numpy as np, matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from kepler import kep2ellipse

plot_dir = './results/plots/'
vid_dir = './results/videos/'

colormaps = {'solar':  ['yellow', 'silver', 'gold', 'deepskyblue', 'tomato', 'peru', 'wheat', 'paleturquoise', 'cornflowerblue', 'pink', 'orange', 'lime'],
             'ksp':    ['yellow', 'tan', 'mediumslateblue', 'lightskyblue', 'salmon', 'darkgoldenrod', 'limegreen', 'silver'],
             'jool':   ['limegreen', 'steelblue', 'deepskyblue', 'peachpuff', 'burlywood', 'yellowgreen'],
             'default':['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']}

# Animate the orbits
def animate(results, set_range=True, plot_range=6e12, COM=False, cm='solar', kepler=False,
            lw=0.7, wspace=0.05, tmargin=0.1, bmargin=0.07, lmargin=0.2, rmargin=0.07, show=False):
    name     = results['id']
    pos      = results['x']
    sim_time = results['t']
    Nsnap    = len(sim_time) # number of snapshots, which is 1 (initial) + the number of steps
    epoch    = results['epoch']

    plt.style.use('dark_background')

    fig1 = plt.figure(figsize=(8,5))
    ax1 = p3.Axes3D(fig1, auto_add_to_figure=False)
    fig1.add_axes(ax1)

    fig2, (ax2, ax3) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    fig2.subplots_adjust(wspace=wspace, left=lmargin, right=1-rmargin, top=1-tmargin, bottom=bmargin)

    if COM: # plot in COM frame
        X_com = results['x_com']
        pos = pos - X_com[:,np.newaxis,:]

    positions = np.moveaxis(pos, 0, -1) # list of paths; the indices are: [object,coord(x,y,z),step]

    paths3d = []; points3d = []
    paths2d_z = []; points2d_z = []
    paths2d_x = []; points2d_x = []

    for i, p in enumerate(positions):
        if i < len(colormaps[cm]):
            color = colormaps[cm][i]
        else:
            color = colormaps['default'][i % len(colormaps['default'])]

        if kepler and COM: # plot Kepler orbit based on initial orbital parameters
            obj = name[i]
            if obj[0] != '*':
                ecc, a, inc, Ω, ω, _, _ = results['kepler'][obj]
                kepler_orbit = kep2ellipse(ecc, a, inc, Ω, ω) # origin at central body
                # Translate orbit since origin of plot is at COM
                kepler_orbit_x = kepler_orbit[0] - X_com[0, 0]
                kepler_orbit_y = kepler_orbit[1] - X_com[0, 1]
                kepler_orbit_z = kepler_orbit[2] - X_com[0, 2]

                ax1.plot(kepler_orbit_x, kepler_orbit_y, kepler_orbit_z, linestyle='--', linewidth=lw/2, color=color)
                ax2.plot(kepler_orbit_x, kepler_orbit_y, linestyle='--', linewidth=lw/2, color=color)
                ax3.plot(kepler_orbit_y, kepler_orbit_z, linestyle='--', linewidth=lw/2, color=color)

        paths3d.append(ax1.plot(p[0, 0:1], p[1, 0:1], p[2, 0:1], linewidth=lw, color=color)[0])
        points3d.append(ax1.plot([], [], [], marker='.', color=color, label=name[i])[0])
        
        paths2d_z.append(ax2.plot(p[0, 0:1], p[1, 0:1], linewidth=lw, color=color)[0]) # polar view
        points2d_z.append(ax2.plot([], [], marker='.', color=color, label=name[i])[0])
        paths2d_x.append(ax3.plot(p[1, 0:1], p[2, 0:1], linewidth=lw, color=color)[0]) # ecliptic view
        points2d_x.append(ax3.plot([], [], marker='.', color=color)[0])
    
    if set_range:
        ax1.set_xlim3d([-plot_range, plot_range])
        ax1.set_ylim3d([-plot_range, plot_range])
        ax1.set_zlim3d([-plot_range, plot_range])

        ax2.set_xlim([-plot_range, plot_range]); ax2.set_ylim([-plot_range, plot_range])
        ax2.set_aspect('equal')
        ax3.set_xlim([-plot_range, plot_range]); ax3.set_ylim([-plot_range, plot_range])
        ax3.set_aspect('equal')
        ax3.yaxis.tick_right(); ax3.yaxis.set_label_position("right")

    ax1.set_xlabel('X [m]'); ax1.set_ylabel('Y [m]'); ax1.set_zlabel('Z [m]')
    ax2.set_xlabel('X [m]'); ax2.set_ylabel('Y [m]')
    ax3.set_xlabel('Y [m]'); ax3.set_ylabel('Z [m]')
    
    ax1.legend(loc='center left', frameon=False, ncol=1, bbox_to_anchor=(-0.25, 0.5))
    ax2.legend(loc='center left', frameon=False, ncol=1, bbox_to_anchor=(-2*lmargin/(1-lmargin), 0.5))

    title3d = ax1.text2D(0.5, 0.95, '', horizontalalignment='center',
                         fontsize=14, transform=ax1.transAxes)

    ax2.set_title('Polar view'); ax3.set_title('Ecliptic view')
    title2d = ax2.text(1, 1+tmargin/(1-tmargin), '', horizontalalignment='center',
                       fontsize=14, transform=ax2.transAxes)

    anim3d = animation.FuncAnimation(fig1, update_paths, Nsnap,
              fargs=(positions, paths3d, points3d, title3d, sim_time, epoch), interval=1, blit=False)
    anim2d = animation.FuncAnimation(fig2, update_paths_2d, Nsnap,
              fargs=(positions, paths2d_z, paths2d_x, points2d_z, points2d_x, title2d, sim_time, epoch), interval=1, blit=False)

    if show:
        plt.show()

    return anim3d, anim2d

# Update paths for animation
# FuncAnimation takes individual iterations of range(Nsnap) and feeds them into 'step'
def update_paths(step, positions, paths, points, title, sim_time, epoch):
    trail_tail = 0 # set to 0 to plot the entire trail
    for path, point, position in zip(paths, points, positions): # for each object
        path.set_data(position[0:2, trail_tail:step])
        path.set_3d_properties(position[2, trail_tail:step])

        point.set_data(position[0, step:step+1], position[1, step:step+1])
        point.set_3d_properties(position[2, step:step+1])
        
    t = sim_time[step]
    if epoch == 0:
        title.set_text('t = %.1f ks'%(t / 1000))
    else:
        JD = epoch + sim_time[step] / 86400
        title.set_text('t = %.1f ks (JD = %.1f)'%(t / 1000, JD))

    return paths, points, title

def update_paths_2d(step, positions, paths2d_z, paths2d_x, points2d_z, points2d_x, title, sim_time, epoch):
    trail_tail = 0 # set to 0 to plot the entire trail
    for path_z, path_x, point_z, point_x, position in zip(paths2d_z, paths2d_x, points2d_z, points2d_x, positions): # for each object
        path_z.set_data(position[0:2, trail_tail:step])
        path_x.set_data(position[1:3, trail_tail:step])
        
        point_z.set_data(position[0, step:step+1], position[1, step:step+1])
        point_x.set_data(position[1, step:step+1], position[2, step:step+1])

    t = sim_time[step]
    if epoch == 0:
        title.set_text('t = %.1f ks'%(t / 1000))
    else:
        JD = epoch + sim_time[step] / 86400
        title.set_text('t = %.1f ks (JD = %.1f)'%(t / 1000, JD))

    return paths2d_z, paths2d_x, points2d_z, points2d_x, title

# plot final state of system
def plot(results, set_range=True, plot_range=6e12, COM=False, cm='solar', kepler=False,
         lw=0.7, wspace=0.05, tmargin=0.1, bmargin=0.07, lmargin=0.2, rmargin=0.07):
    name     = results['id']
    pos      = results['x']
    sim_time = results['t']
    epoch    = results['epoch']

    plt.style.use('dark_background')

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(projection='3d')
    fig1.add_axes(ax1)

    fig2, (ax2, ax3) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    fig2.subplots_adjust(wspace=wspace, left=lmargin, right=1-rmargin, top=1-tmargin, bottom=bmargin)

    X_com = results['x_com']
    if COM: # plot in COM frame
        pos = pos - X_com[:,np.newaxis,:]
        #ax1.plot([0], [0], [0], color='r', marker='+', markersize=8, mew=2)
        #ax2.plot([0], [0], color='r', marker='+', markersize=8, mew=2)
        ax3.plot([0], [0], color='r', marker='x', markersize=8, mew=2)
    else:
        #ax1.plot(X_com[-1:,0], X_com[-1:,1], X_com[-1:,2], color='r', marker='+', markersize=8, mew=2)
        #ax2.plot(X_com[-1:,0], X_com[-1:,1], color='r', marker='+', markersize=8, mew=2)
        ax3.plot(X_com[-1:,1], X_com[-1:,2], color='r', marker='x', markersize=8, mew=2)

    for i in range(len(name)):
        if i < len(colormaps[cm]):
            color = colormaps[cm][i]
        else:
            color = colormaps['default'][i % len(colormaps['default'])]

        if kepler and COM: # plot Kepler orbit based on initial orbital parameters
            obj = name[i]
            if obj[0] != '*':
                ecc, a, inc, Ω, ω, _, _ = results['kepler'][obj]
                kepler_orbit = kep2ellipse(ecc, a, inc, Ω, ω=ω) # origin at central body
                # Translate orbit since origin of plot is at COM
                kepler_orbit_x = kepler_orbit[0] - X_com[0, 0]
                kepler_orbit_y = kepler_orbit[1] - X_com[0, 1]
                kepler_orbit_z = kepler_orbit[2] - X_com[0, 2]

                ax1.plot(kepler_orbit_x, kepler_orbit_y, kepler_orbit_z, linestyle='--', linewidth=lw/2, color=color)
                ax2.plot(kepler_orbit_x, kepler_orbit_y, linestyle='--', linewidth=lw/2, color=color)
                ax3.plot(kepler_orbit_y, kepler_orbit_z, linestyle='--', linewidth=lw/2, color=color)

        ax1.plot(pos[:,i,0], pos[:,i,1], pos[:,i,2], linewidth=lw, color=color)
        ax1.plot(pos[-1,i,0], pos[-1,i,1], pos[-1,i,2], marker='.', color=color, label=name[i])
        
        ax2.plot(pos[:,i,0], pos[:,i,1], linewidth=lw, color=color)
        ax2.plot(pos[-1,i,0], pos[-1,i,1], marker='.', color=color, label=name[i])
        ax3.plot(pos[:,i,1], pos[:,i,2], linewidth=lw, color=color)
        ax3.plot(pos[-1,i,1], pos[-1,i,2], marker='.', color=color)

    if set_range:
        ax1.set_xlim3d([-plot_range, plot_range])
        ax1.set_ylim3d([-plot_range, plot_range])
        ax1.set_zlim3d([-plot_range, plot_range])

        ax2.set_xlim([-plot_range, plot_range]); ax2.set_ylim([-plot_range, plot_range])
        ax2.set_aspect('equal')
        ax3.set_xlim([-plot_range, plot_range]); ax3.set_ylim([-plot_range, plot_range])
        ax3.set_aspect('equal')
        ax3.yaxis.tick_right(); ax3.yaxis.set_label_position("right")

    ax1.set_xlabel('X [m]'); ax1.set_ylabel('Y [m]'); ax1.set_zlabel('Z [m]')
    ax2.set_xlabel('X [m]'); ax2.set_ylabel('Y [m]')
    ax3.set_xlabel('Y [m]'); ax3.set_ylabel('Z [m]')
    
    ax1.legend(loc='center left', frameon=False, ncol=1, bbox_to_anchor=(-0.4, 0.5))
    ax2.legend(loc='center left', frameon=False, ncol=1, bbox_to_anchor=(-2*lmargin/(1-lmargin), 0.5))

    if epoch == 0:
        title = 't = %.1f ks'%(sim_time[-1] / 1000)
    else:
        title = 't = %.1f ks (JD = %.1f)'%(sim_time[-1] / 1000, epoch + sim_time[-1] / 86400)

    title3d = ax1.text2D(0.5, 0.95, title, horizontalalignment='center',
                         fontsize=14, transform=ax1.transAxes)

    ax2.set_title('Polar view'); ax3.set_title('Ecliptic view')
    title2d = ax2.text(1, 1+tmargin/(1-tmargin), title, horizontalalignment='center',
                       fontsize=14, transform=ax2.transAxes)

    return fig1, fig2

def make_animation(input_path, output_dir, set_range=True, plot_range=6e12, COM=False, cm='solar', kepler=False,
                   lw=0.7, wspace=0.05, tmargin=0.1, bmargin=0.07, lmargin=0.2, rmargin=0.07):
    results = np.load(input_path, allow_pickle=True).item()
    output_name = input_path.split('/')[-1][:-4]

    anim3d, anim2d = animate(results, set_range=set_range, plot_range=plot_range, COM=COM, cm=cm, kepler=kepler, 
                     lw=lw, wspace=wspace, tmargin=tmargin, bmargin=bmargin, lmargin=lmargin, rmargin=rmargin, show=False)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=12, metadata=dict(artist='Me'), bitrate=1800)
    anim3d.save(output_dir+output_name+'_3d.mp4', writer=writer)#, fps=15, extra_args=['-vcodec', 'libx264'])
    anim2d.save(output_dir+output_name+'_2d.mp4', writer=writer)#, fps=15, extra_args=['-vcodec', 'libx264'])

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        assert len(sys.argv) >= 3

        # Example1: python3 visualization.py results/SolarSystem_10steps_3652.5d_n_2.0.npy 8e11 solar 1
        # Example2: python3 visualization.py results/KSPsystem_10steps_3652.5d_n_2.0.npy 1e11 ksp 1
        input_path = sys.argv[1]
        plot_range = float(sys.argv[2])
        if len(sys.argv) >= 4:
            colormap = sys.argv[3]
        else:
            colormap = 'default'

        if len(sys.argv) >= 5:
            kepler = bool(int(sys.argv[4]))
        else:
            kepler = False

        make_animation(input_path, vid_dir, set_range=True, plot_range=plot_range, COM=True, cm=colormap, kepler=kepler)
