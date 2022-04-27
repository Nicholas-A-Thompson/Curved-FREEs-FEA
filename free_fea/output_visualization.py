# output plots for curved FREE FEA

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
import seaborn as sns
from free_fea.geometry_setup import *

# 3D plot of centerlines and fibers
# plot_frame = 0 / 1: do not / do display coordinate system at FREE end
def plot_3D_FREE(ax, centerlines, alpha, beta, frames=[0, -1]):
    num_fibers = np.shape(alpha)[2]

    # plot display settings
    linewidth = 0.75
    # plot colors
    color_centerlines = ['black', 'gray']
    color_fiber_a = ['blue', 'lightcoral']
    color_fiber_b = ['green', 'orange']

    local_extrema_x = np.array([])
    local_extrema_y = np.array([])
    local_extrema_z = np.array([])

    # plot centerline and fibers for each specified frame
    for f in frames:
        # update list of local extrema for plot limits
        local_extrema_x = np.append(local_extrema_x, [alpha[f, :, :, 0].min(),
                                                      alpha[f, :, :, 0].min(),
                                                      beta[f, :, :, 0].max(),
                                                      beta[f, :, :, 0].max()])
        local_extrema_y = np.append(local_extrema_y, [alpha[f, :, :, 1].min(),
                                                      alpha[f, :, :, 1].min(),
                                                      beta[f, :, :, 1].max(),
                                                      beta[f, :, :, 1].max()])
        local_extrema_z = np.append(local_extrema_z, [alpha[f, :, :, 2].min(),
                                                      alpha[f, :, :, 2].min(),
                                                      beta[f, :, :, 2].max(),
                                                      beta[f, :, :, 2].max()])

        # plot centerlines
        x = centerlines[f, :, 0]
        y = centerlines[f, :, 1]
        z = centerlines[f, :, 2]
        ax.plot3D(x, y, z, color=color_centerlines[f], linewidth=linewidth)  # plot fiber

        # plot fibers
        for j in range(num_fibers):
            x = alpha[f, :, j, 0]
            y = alpha[f, :, j, 1]
            z = alpha[f, :, j, 2]
            ax.plot3D(x, y, z, color=color_fiber_a[f], linewidth=linewidth)  # plot fiber

            x = beta[f, :, j, 0]
            y = beta[f, :, j, 1]
            z = beta[f, :, j, 2]
            ax.plot3D(x, y, z, color=color_fiber_b[f], linewidth=linewidth)  # plot fiber

    # axis labels
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_zlabel('z (mm)')

    # set plot axes to be equal
    max_range = np.array([local_extrema_x.max() - local_extrema_x.min(),
                          local_extrema_y.max() - local_extrema_y.min(),
                          local_extrema_z.max() - local_extrema_z.min()]).max() / 2.0
    mid_x = (local_extrema_x.min()+ local_extrema_x.max()) / 2.0
    mid_y = (local_extrema_y.min()+ local_extrema_y.max()) / 2.0
    mid_z = (local_extrema_z.min()+ local_extrema_z.max()) / 2.0
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    ax.set_box_aspect((1, 1, 1))


# 3D plot of centerlines and fibers
# plot_frame = 0 / 1: do not / do display coordinate system at FREE end
def plot_2D_FREE(ax, centerlines, alpha, beta, frames=[0, -1]):
    num_fibers = np.shape(alpha)[2]

    # plot display settings
    linewidth = 0.75
    # plot colors
    color_centerlines = ['black', 'gray']
    color_fiber_a = ['blue', 'lightcoral']
    color_fiber_b = ['green', 'orange']

    local_extrema_x = np.array([])
    local_extrema_z = np.array([])

    # plot centerline and fibers for each specified frame
    for f in frames:
        # update list of local extrema for plot limits
        local_extrema_x = np.append(local_extrema_x, [alpha[f, :, :, 0].min(),
                                                      alpha[f, :, :, 0].min(),
                                                      beta[f, :, :, 0].max(),
                                                      beta[f, :, :, 0].max()])
        local_extrema_z = np.append(local_extrema_z, [alpha[f, :, :, 2].min(),
                                                      alpha[f, :, :, 2].min(),
                                                      beta[f, :, :, 2].max(),
                                                      beta[f, :, :, 2].max()])

        # plot centerlines
        x = centerlines[f, :, 0]
        z = centerlines[f, :, 2]
        ax.plot(x, z, color=color_centerlines[f], linewidth=linewidth)  # plot fiber

        # plot fibers
        for j in range(num_fibers):
            x = alpha[f, :, j, 0]
            z = alpha[f, :, j, 2]
            ax.plot(x, z, color=color_fiber_a[f], linewidth=linewidth)  # plot fiber

            x = beta[f, :, j, 0]
            z = beta[f, :, j, 2]
            ax.plot(x, z, color=color_fiber_b[f], linewidth=linewidth)  # plot fiber

    # axis labels
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('z (mm)')

    # set plot axes to be equal
    ax.set_aspect('equal', adjustable="datalim")


# plot one set of mean fiber angles
# pressure[f] = pressure at frame f
# fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0), midline (t = 1) and outside (t = 2) of the curvature
def plot_fiber_angles_vs_pressure(ax, pressure, fiber_angles, group=0):
    # plot display settings
    linewidth = 1.0

    # plot colors
    if group == 1:
        colors = ['lightblue', 'blue', 'midnightblue']
    elif group == 2:
        colors = ['lightgreen', 'green', 'darkgreen']
    else:
        colors = ['gold', 'orange', 'orangered']

    for t in range(np.shape(fiber_angles)[2]):
        ax.plot(pressure, fiber_angles[:, :, t], color=colors[t], linewidth=linewidth)  # plot inside angles

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('fiber angle (' + u'\N{DEGREE SIGN}' + ')')

    ax.set_xlim(pressure.min(), pressure.max())


# plot FREE radius vs fiber angle(s)
# radii[f][i] = radius at s = s[i] in frame f
# fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0), midline (t = 1) and outside (t = 2) of the curvature
# circ_position: circumferential position, 0 = inside, 1 = midline, 2 = outside
def plot_radius_vs_fiber_angle(ax, radii, fiber_angles, bulge_index, circ_position=1):
    # plot display settings
    linewidth = 1.0

    # list of indices between the endcap bulge borders
    bulge_free_zone = list(range(bulge_index, np.shape(radii)[1] - bulge_index))

    # initialize mean arrays
    mean_radii = np.zeros(np.shape(radii)[0])
    mean_fiber_angles = np.zeros(np.shape(radii)[0])

    for f in range(np.shape(radii)[0]):
        mean_radii[f] = np.mean(radii[f, bulge_free_zone])
        mean_fiber_angles[f] = np.mean(fiber_angles[f, bulge_free_zone, circ_position])

    ax.plot(mean_fiber_angles, mean_radii, color='blue', linewidth=linewidth)  # plot inside angles

    # axis labels
    ax.set_xlabel('mean fiber angle (' + u'\N{DEGREE SIGN}' + ')')
    ax.set_ylabel('mean FREE radius (mm)')

    ax.set_xlim(mean_fiber_angles.min(), mean_fiber_angles.max())


# plot fiber angles vs normalized position
# ax = matplotlib axes
# norm_positions[f][i] = position s along path with 0 <= s <= 1
# fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0), midline (t = 1) and outside (t = 2) of the curvature
# frames and colors must be same length
def plot_angle_vs_norm_position(ax, norm_positions, fiber_angles, group=0, frames=[0, -1]):
    # plot display settings
    linewidth = 1.0

    for f in frames:    # plot colors
        for t in range(np.shape(fiber_angles)[2]):
            label = ''
            if group == 1 and f == 0:
                colors = ['lightblue', 'darkblue', 'lightblue']
                if t == 1: label = r'$\alpha$' + ' undeformed'
            elif group == 1 and f == -1:
                colors = ['lightcoral', 'darkred', 'lightcoral']
                if t == 1: label = r'$\alpha$' + ' deformed'
            elif group == 2 and f == 0:
                colors = ['lightgreen', 'darkgreen', 'lightgreen']
                if t == 1: label = r'$\beta$' + ' undeformed'
            elif group == 2 and f == -1:
                colors = ['orange', 'orangered', 'orange']
                if t == 1: label = r'$\beta$' + ' deformed'
            else:
                colors = ['gray', 'black', 'gray']

            ax.plot(norm_positions[f, :], fiber_angles[f, :, t], color=colors[t], linewidth=linewidth, label=label)  # plot fiber angles

    # axis labels
    ax.set_xlabel('normalized position')
    ax.set_ylabel('fiber angle (' + u'\N{DEGREE SIGN}' + ')')

    # ax.legend()
    ax.set_xlim(0, 1)


# plot radius of curvature vs pressure
# pressure[f] = pressure at frame f
# mean_RoC[f] = mean centerline radius of curvature at in frame f
def plot_RoC_vs_pressure(ax, pressure, mean_RoC, label, linestyle, color='blue'):
    # plot display settings
    linewidth = 1.0

    # ax.plot(pressure, centerline_RoCs[:, :], color=color, linewidth=linewidth, label=label, linestyle=linestyle)  # plot radius of curvature
    ax.plot(pressure, mean_RoC, color=color, linewidth=linewidth, label=label, linestyle=linestyle)  # plot radius of curvature

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('radius of curvature (m)')

    ax.set_xlim(pressure.min(), pressure.max())


# plot lambda_1 and lamda_2 vs pressure
# pressure[f] = pressure at frame f
def plot_stretch_ratios_vs_pressure(ax, pressure, lambda_1, lambda_2, colors=sns.color_palette()):
    # plot display settings
    linewidth = 1.0

    ax.plot(pressure, lambda_1, color=colors[1], linewidth=linewidth, label='axial, ' + r'$\lambda_1$')  # plot axial stretch ratio
    ax.plot(pressure, lambda_2, color=colors[2], linewidth=linewidth, label='radial, ' + r'$\lambda_2$')  # plot radial stretch ratio

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('stretch ratio (-)')

    ax.set_xlim(pressure.min(), pressure.max())


# plot delta vs pressure
# pressure[f] = pressure at frame f
# delta[f] = axial twist at the endpoint at frame f
def plot_twist_vs_pressure(ax, pressure, delta, colors=sns.color_palette()):
    # plot display settings
    linewidth = 1.0

    ax.plot(pressure, delta, color=colors[3], linewidth=linewidth, label='twist, ' + r'$\delta$')  # plot twist

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('twist (' + r'$\degree$' + ')')

    ax.set_xlim(pressure.min(), pressure.max())


# plot radius of curvature vs normalized position
# ax = matplotlib axes
# norm_positions[f][i] = position s along path with 0 <= s <= 1
# centerline_RoC[f][i] = centerline radius of curvature at s = s[i] in frame f
# frames and colors must be same length
def plot_RoC_vs_norm_position(ax, norm_positions, centerline_RoCs, frames = [0, -1], colors=['blue', 'red']):
    # plot display settings
    linewidth = 1.0

    i = 0
    for f in frames:
        color = colors[i]
        if f == 0:
            label = 'undeformed'
        else:
            label = 'deformed'
        ax.plot(norm_positions[f, :], centerline_RoCs[f, :], color=color, linewidth=linewidth, label=label)  # plot radius of curvature
        i += 1

    # format plot
    ax.set_xlabel('normalized position')
    ax.set_ylabel('radius of curvature (mm)')
    ax.set_xlim(0, 1)


# plot radius of curvature vs normalized position
# ax = matplotlib axes
# norm_positions[f][i] = position s along path with 0 <= s <= 1
# curvature[f][i] = centerline curvature at s = s[i] in frame f
# frames and colors must be same length
def plot_curvature_vs_norm_position(ax, norm_positions, curvature, pressure, frames = [0, -1]):
    # plot display settings
    linewidth = 1.0

    colors = sns.color_palette(palette="viridis", n_colors=np.shape(frames)[0])

    i = 0
    for f in frames:
        color = colors[i]
        label = 'P = ' + str(int(pressure[f])) + ' kPa'
        ax.plot(norm_positions[f, :], curvature[f, :], color=color, linewidth=linewidth,
                label=label)  # plot curvature
        i += 1

    # format plot
    ax.set_xlabel('normalized position')
    ax.set_ylabel('FREE curvature (mm' + r'$^{-1}$' + ')')
    ax.set_xlim(0, 1)
    legend_without_duplicate_labels(ax)


# plot FREE radius vs. normlaized position
def plot_radius_vs_norm_position(ax, norm_positions, radii, pressure, frames = [0, -1]):
    # plot display settings
    linewidth = 1.0

    colors = sns.color_palette(palette="viridis", n_colors=np.shape(frames)[0])

    i = 0
    for f in frames:
        color = colors[i]
        label = 'P = ' + str(int(pressure[f])) + ' kPa'
        ax.plot(norm_positions[f, :], radii[f, :], color=color, linewidth=linewidth, label=label)
        i += 1

    # format plot
    ax.set_xlabel('normalized position')
    ax.set_ylabel('FREE radius (mm)')
    ax.set_xlim(0, 1)
    legend_without_duplicate_labels(ax)


# plot axial stretch vs pressure
def plot_lambda_1_vs_pressure(ax, pressure, lambda_1, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    ax.plot(pressure, lambda_1,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('axial stretch ratio ' + r'$\lambda_1$')


# plot radial stretch vs pressure
def plot_lambda_2_vs_pressure(ax, pressure, lambda_2, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.0

    ax.plot(pressure, lambda_2,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('radial stretch ratio ' + r'$\lambda_2$')


# plot axial rotation vs pressure
# same as plot_twist_vs_pressure, but more arguments for customization
def plot_delta_vs_pressure(ax, pressure, delta, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.0

    ax.plot(pressure, delta,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('net axial rotation ' + r'$\delta$ ($\degree$)')


# plot axial stretch vs fiber angle (assumes symmetric fibers)
def plot_lambda_1_vs_alpha(ax, alpha, lambda_1, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    ax.plot(alpha, lambda_1,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('fiber angle ' + r'$\alpha$ ($\degree$)')
    ax.set_ylabel('axial stretch ratio ' + r'$\lambda_1$')


# plot radial stretch vs fiber angle (assumes symmetric fibers)
def plot_lambda_2_vs_alpha(ax, alpha, lambda_2, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.0

    ax.plot(alpha, lambda_2,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('fiber angle ' + r'$\alpha$ ($\degree$)')
    ax.set_ylabel('radial stretch ratio ' + r'$\lambda_2$')


# plot axial rotation vs fiber angle (assumes symmetric fibers)
def plot_delta_vs_alpha(ax, alpha, delta, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.0

    ax.plot(alpha, delta,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            label=label)

    # axis labels
    ax.set_xlabel('fiber angle ' + r'$\alpha$ ($\degree$)')
    ax.set_ylabel('net axial rotation ' + r'$\delta$ ($\degree$)')


# plot radius of curvature vs pressure
# pressure[f] = pressure at frame f
# RoC[f] = mean centerline radius of curvature in frame f
def plot_mean_RoC_vs_pressure(ax, pressure, RoCs, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    # plot radius of curvature
    ax.plot(pressure, RoCs, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('radius of curvature (mm)')


# plot mean curvature vs pressure
# pressure[f] = pressure at frame f
# curvature[f] = mean centerline curvature in frame f
def plot_mean_curvature_vs_pressure(ax, pressure, curvature, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    # plot radius of curvature
    ax.plot(pressure, curvature, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('curvature (1/mm)')


# plot normalized curvature vs pressure
# pressure[f] = pressure at frame f
# norm_curvature[f] = normalized mean centerline curvature in frame f
def plot_norm_curvature_vs_pressure(ax, pressure, norm_curvature, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    # plot radius of curvature
    ax.plot(pressure, norm_curvature, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('normalized curvature')


# plot mean fiber angle vs pressure
# pressure[f] = pressure at frame f
# alpha[f] = mean fiber angle in frame f
def plot_mean_alpha_vs_pressure(ax, pressure, alpha, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    # plot radius of curvature
    ax.plot(pressure, alpha, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('alpha (deg)')

# plot radius * curvature vs pressure
def plot_mean_rk_vs_pressure(ax, pressure, curvature, radii, label, color='blue', linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    rk = curvature * radii

    # plot radius of curvature
    ax.plot(pressure, rk, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('r' + r'$\kappa$')


def plot_total_dist_force_vs_pressure(ax, pressure, force, color, label, linestyle='solid'):
    # plot display settings
    linewidth = 1.5

    # plot radius of curvature
    ax.plot(pressure, force, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # axis labels
    ax.set_xlabel('pressure (kPa)')
    ax.set_ylabel('net total distributed force,' + r'$F_d$' + ' (N)')


def plot_dist_force_vs_norm_position(ax, norm_positions, force, pressure, frames=[0, -1]):
    # plot display settings
    linewidth = 1.0

    colors = sns.color_palette(palette="viridis", n_colors=np.shape(frames)[0])

    i = 0
    for f in frames:
        color = colors[i]
        label = 'P = ' + str(int(pressure[f])) + ' kPa'
        ax.plot(norm_positions[f, :], force[f, :], color=color, linewidth=linewidth, label=label)  # plot radius / radius of curvature
        # ax.plot(norm_positions[f, :], r_ratio[f, :], linewidth=linewidth, label=label)  # plot radius / radius of curvature
        i += 1

    # format plot
    ax.set_xlabel('normalized position')
    ax.set_ylabel('net distributed force, d' + r'$F_d$' + ' (N/mm)')
    ax.set_xlim(0, 1)
    legend_without_duplicate_labels(ax)


# fix duplicate labels
# https://stackoverflow.com/questions/19385639/duplicate-items-in-legend-in-matplotlib
def legend_without_duplicate_labels(ax, ncol=1):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    leg = ax.legend(
        *zip(*unique),
        fontsize='small',
        frameon=False,
        ncol=ncol,
        handlelength=3.0)
    leg.set_draggable(state=True)

    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)


def format_plot(myFREE, label_suffix='', include_r=False, include_plot_method=False):
    alpha = int(myFREE.fiber_a_angle_deg)
    type = myFREE.path_type
    radius = int(myFREE.outer_radius)

    if include_plot_method:
        if myFREE.has_def_fib_angle:
            label_suffix = ', CFA'
        if not myFREE.has_def_fib_angle and type != 'straight':
            label_suffix = ', CP'
            radius = 8

    if type == 'straight':
        R = 1e6
        if include_r:
            plot_label = ('STR, ' + 'r' + str(radius) + ', ' +
                          r'$\alpha = $' + str(alpha) + r'$\degree$' + label_suffix)
        else:
            plot_label = ('STR, ' +
                          r'$\alpha = $' + str(alpha) + r'$\degree$' + label_suffix)

    elif type == 'arc':
        R = int(myFREE.path_parameters[0])
        if include_r:
            plot_label = (
                          'R' + str(R) + ', ' +
                          'r' + str(radius) + ', ' +
                          r'$\alpha = $' + str(alpha) + r'$\degree$' + label_suffix)
        else:
            plot_label = ('R' + str(R) + ', ' +
                          r'$\alpha = $' + str(alpha) + r'$\degree$' + label_suffix)
    else:
        R = 0
        plot_label = myFREE.label

    if alpha == 25:
        if radius == 5:  # blue
            c1 = '084594'
            c2 = '6baed6'
        elif radius == 8:  # purple
            c1 = '710bd6'
            c2 = 'bd62e8'
    elif alpha == 30:
        if radius == 5:  # green
            c1 = '0aa305'
            c2 = '7fcb6b'
        elif radius == 8:  # teal
            c1 = '078280'
            c2 = '05dede'
    elif alpha == 40:
        if radius == 5:  # orange
            c1 = '8c2d04'
            c2 = 'fe9929'
        elif radius == 8:  # red
            c1 = '820720'
            c2 = 'c1565e'
    else:  # bright pink
        color = 'magenta'
        linestyle = 'solid'

    if R == 1e6:  # gray
        # c1 = '252525'
        # c2 = '969696'
        color = create_color_list(4, color1=c1, color2=c2)[0]
        linestyle = 'solid'
    elif R == 20 or R == 32:
        color = create_color_list(4, color1=c1, color2=c2)[1]
        linestyle = 'dotted'
    elif R == 30 or R == 48:
        color = create_color_list(4, color1=c1, color2=c2)[2]
        linestyle = 'dashed'
    elif R == 60:  # red
        color = create_color_list(4, color1=c1, color2=c2)[3]
        linestyle = 'solid'
    else:  # bright pink
        color = 'magenta'
        linestyle = 'solid'

####
    #
    # if R == 1e6:  # black
    #     c1 = '000000'
    #     c2 = '969696'
    # elif R == 20:  # blue
    #     c1 = '08306b'
    #     c2 = '6baed6'
    # elif R == 30:  # green
    #     c1 = '00441b'
    #     c2 = '74c476'
    # elif R == 32:  # purple
    #     c1 = '422063'
    #     c2 = '807dba'
    # elif R == 40:  # teal
    #     c1 = '014636'
    #     c2 = '27b894'
    # elif R == 48:  # orange
    #     c1 = '7f2704'
    #     c2 = 'fd8d3c'
    # elif R == 60:  # red
    #     c1 = '7f2704'
    #     c2 = 'fd8d3c'
    # else:  # bright green
    #     c1 = '00ff00'
    #     c2 = '00ff00'
    #
    # if alpha == 25:
    #     color = create_color_list(3, color1=c1, color2=c2)[0]
    #     linestyle = 'dashed'
    # elif alpha == 30:
    #     color = create_color_list(3, color1=c1, color2=c2)[1]
    #     linestyle = 'dashdot'
    # elif alpha == 40:
    #     color = create_color_list(3, color1=c1, color2=c2)[2]
    #     linestyle = 'solid'
    # else:
    #     color = '000000'
    #     linestyle = 'dotted'

    return plot_label, color, linestyle


# create custom colormap from two hex values
def create_color_list(n_colors, color1='081d58', color2='c7e9b4'):
    color1_rgb = tuple(int(color1[i:i+2], 16) for i in (0, 2, 4))
    color2_rgb = tuple(int(color2[i:i+2], 16) for i in (0, 2, 4))
    colors = np.array([color1_rgb, color2_rgb]) / 256
    cmap = LinearSegmentedColormap.from_list('cmap', colors, N=n_colors)
    color_list = cmap(np.arange(n_colors))
    return color_list

