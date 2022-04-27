# process output data after saving output data with 2_Abaqus_FREE_output.py
# This file should be run OUTSIDE of Abaqus using Python 3.
# outputs CSV to Abaqus Working Directory

import csv

import numpy as np

from free_fea.initialization import *
from free_fea.post_processing import *
from free_fea.output_visualization import *

########## USER INPUT ##########
# initialize FREE object
# label='user_input' for a new FREE from free_fea.user_input.
# label='[input file name]' to load a FREE from an input file
myFREE = initialize_FREE(label='user_input')
# myFREE = initialize_FREE(label='cbez_0-0-0_20-60-20_40--40-20_80-20-0_r5_a30_b-30')
################################

# main code block
def main(myFREE):

    current_directory = os.getcwd()  # get current work directory
    output_directory = os.path.join(current_directory, 'output', myFREE.label)  # navigate to curved_FREEs\output\(myFREE.label)
    output = np.load(os.path.join(output_directory, 'output.npz'))  # load output file

    # load FEA output data
    pressure = output['pressure'] * 1e3  # convert pressure from MPa to kPa

    # fiber_a/b_coordinates[f][i][j] = np.array([x, y, z]) for fiber j at s = s[i] in frame f (group A/B)
    raw_fiber_a_coordinates = output['fiber_a_coordinates']
    raw_fiber_b_coordinates = output['fiber_b_coordinates']

    # smooth fiber paths
    fiber_a_coordinates = smooth_path_savgol(raw_fiber_a_coordinates, window_length=11, polyorder=3)  # smooth fiber paths
    fiber_b_coordinates = smooth_path_savgol(raw_fiber_b_coordinates, window_length=11, polyorder=3)  # smooth fiber paths
    fiber_a_coordinates = smooth_path_spline(fiber_a_coordinates, k=5, num_points=25)  # smooth fiber paths
    fiber_b_coordinates = smooth_path_spline(fiber_b_coordinates, k=5, num_points=25)  # smooth fiber paths
    fiber_a_coordinates[0] = raw_fiber_a_coordinates[0]
    fiber_b_coordinates[0] = raw_fiber_b_coordinates[0]

    num_frames = np.shape(fiber_a_coordinates)[0]  # retrieve number of FEA frames
    num_path_points = np.shape(fiber_a_coordinates)[1]  # retrieve number of path points

    # trim the endcap bulging effect from all output arrays
    # bulge_index = find_endcap_bulge_index(alpha, beta, tol=0.1)  # Default tolerance is 0.3. Adjust based on plots.
    bulge_index = 15
    bulge_free_zone = list(range(bulge_index, num_path_points - bulge_index))

    # endcap_coordinates[f] = position in frame f of point on free endcap offset from center by e1 in first frame
    if 'endcap_coordinates' in output:
        endcap_coordinates = output['endcap_coordinates']
    else:
        endcap_coordinates = np.zeros([num_frames, 2, 3])

    # bladder_coordinates[f][i][j] = point at bladder circumferential index j at s = s[i] in frame f
    raw_inner_bladder_coordinates = output['inner_bladder_coordinates']

    # find centerlines from inner bladder nodes
    inner_bladder_coordinates = np.zeros(np.shape(raw_inner_bladder_coordinates))
    for j in range(np.shape(raw_inner_bladder_coordinates)[2]):
        savgol_blad_coords = smooth_path_savgol(
            raw_inner_bladder_coordinates[:, :, j, :], window_length=21, polyorder=3)
        inner_bladder_coordinates[:, :, j, :] = smooth_path_spline(
            savgol_blad_coords, k=5, num_points=10)
    centerlines = get_centerlines_from_bladder(inner_bladder_coordinates)
    centerlines[0] = myFREE.path  # set centerline in initial frame to original path

    # extract FREE output measurements
    lengths = get_lengths(centerlines)  # lengths[f] = length of centerline in frame f
    radii = get_all_radii(centerlines, fiber_a_coordinates)  # radii[f][i] = radius at s = s[i] in frame f
    lambda_1, lambda_2 = get_stretch_ratios(lengths, radii)
    delta = get_twist(centerlines, endcap_coordinates)
    curvatures, RoCs = get_all_curves(centerlines)  # RoCs[f][i] = centerline radius of curvature at s = s[i] in frame f
    norm_positions = get_normalized_positions(centerlines)  # converts index range to length normalized positions (0 <= t <= 1)
    rk = radii * curvatures
    mean_curvatures = np.apply_along_axis(np.mean, 1, curvatures[:, bulge_free_zone])
    mean_radii = np.apply_along_axis(np.mean, 1, curvatures[:, bulge_free_zone])

    # fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0), midline (t = 1) and outside (t = 2) of the curvature
    alpha = get_fiber_angles(centerlines, fiber_a_coordinates)  # fiber angles for group A (deg)
    beta = get_fiber_angles(centerlines, fiber_b_coordinates)  # fiber angles for group B (deg)

    d_dist_force, total_dist_force = calculate_distributed_force(centerlines, raw_inner_bladder_coordinates, pressure)

    # write output data to csv file

    fields = define_output_fields()

    # data rows of csv file
    data = [{}] * num_frames * num_path_points
    for f in range(num_frames):
        for i in range(num_path_points):

            # write data to csv file
            row = f * num_path_points + i

            data[row] = {'Frame': f,
                         'Pressure': pressure[f],
                         'Index': i,
                         'Smooth CL X': centerlines[f][i][0],
                         'Smooth CL Y': centerlines[f][i][1],
                         'Smooth CL Z': centerlines[f][i][2],
                         'Radius': radii[f][i],
                         'Length': lengths[f],
                         'Alpha (inside)': alpha[f][i][0],
                         'Alpha (midline)': alpha[f][i][1],
                         'Alpha (outside)': alpha[f][i][2],
                         'Beta (inside)': beta[f][i][0],
                         'Beta (midline)': beta[f][i][1],
                         'Beta (outside)': beta[f][i][2],
                         'Curvature': curvatures[f][i]}

    # name of csv file
    output_csv = os.path.join(output_directory, myFREE.label + '.csv')
    
    # write to csv file
    with open(output_csv, 'w') as csvfile:
        # create a csv dict writer object
        csvwriter = csv.DictWriter(csvfile, fieldnames=fields, lineterminator='\n')

        # writing the header (field names)
        csvwriter.writeheader()

        # writing the data rows
        csvwriter.writerows(data)


    # plot output
    fig = plt.figure(figsize=(10, 8), dpi=100.0)  # initialize figure
    fig.suptitle(myFREE.label)

    # format subplot layout
    # plt.tight_layout()
    plt.subplots_adjust(top=0.93,
                        bottom=0.08,
                        left=0.09,
                        right=0.98,
                        hspace=0.25,
                        wspace=0.35)

    # sample pressures (frames) for multi-line plots
    plot_frames = np.concatenate((np.arange(0, num_frames, int(num_frames/6)), np.array([num_frames - 1])))
    if plot_frames[-2] == plot_frames[-1]:
        plot_frames = plot_frames[:-1]  # trim last entry off of plot_frames if it is a repeat

    # subplot 1
    # ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    plot_3D_FREE(ax, centerlines, fiber_a_coordinates, fiber_b_coordinates, frames=[0, -1])
    # ax = fig.add_subplot(2, 2, 1)
    # ax = fig.add_subplot(1, 1, 1)
    # plot_2D_FREE(ax, centerlines, fiber_a_coordinates, fiber_b_coordinates, frames=[0, -1])

    '''
    # subplot 2
    ax = fig.add_subplot(2, 2, 2)
    # plot_fiber_angles_vs_pressure(ax, pressure, alpha, group=1)
    # plot_fiber_angles_vs_pressure(ax, pressure, beta, group=2)
    # plot_angle_vs_norm_position(ax, norm_positions, alpha, group=1, frames=[0, -1])
    # plot_radius_vs_norm_position(ax, norm_positions, radii, pressure, frames=[0, -1])
    # plot_curvature_vs_norm_position(ax, norm_positions, curvatures, pressure, frames=plot_frames)
    # plot_angle_vs_norm_position(ax, norm_positions, beta, group=2, frames=[0, -1])
    # plt.vlines(x=norm_positions[-1][[bulge_index, -bulge_index - 1]],
    #            ymin=ax.get_ylim()[0],
    #            ymax=ax.get_ylim()[1],
    #            color='gray',
    #            linestyles='dashed',
    #            label='endcap bulge boundary')
    plot_mean_curvature_vs_pressure(
        ax=ax,
        pressure=pressure,
        curvature=mean_curvatures,
        color='red',
        label=myFREE.label)
    # legend_without_duplicate_labels(ax)
    # plot_total_dist_force_vs_pressure(
    #     ax=ax,
    #     pressure=pressure,
    #     force=np.apply_along_axis(np.min, 1, d_dist_force),
    #     color='blue',
    #     label=myFREE.label,
    #     linestyle='solid')

    # subplot 3
    ax = fig.add_subplot(2, 2, 3)
    # plot_RoC_vs_norm_position(ax, norm_positions, RoCs, frames=[0, -1], colors=['blue', 'red'])
    # plot_curvature_vs_norm_position(ax, norm_positions, curvatures, pressure, frames=plot_frames)
    # plt.vlines(x=norm_positions[-1][[bulge_index, -bulge_index - 1]],
    #            ymin=ax.get_ylim()[0],
    #            ymax=ax.get_ylim()[1],
    #            color='gray',
    #            linestyles='dashed',
    #            label='endcap bulge boundary')
    # plot_radius_vs_fiber_angle(ax=ax, radii=radii, fiber_angles=alpha, bulge_index=bulge_index, circ_position=1)  # plot mean midline fiber angle
    # plot_stretch_ratios_vs_pressure(ax, pressure, lambda_1, lambda_2)
    # plot_dist_force_vs_norm_position(ax, norm_positions, d_dist_force, pressure, frames=plot_frames)
    # plot_mean_curvature_vs_pressure(
    #     ax=ax,
    #     pressure=pressure,
    #     curvature=mean_curvatures,
    #     color='red',
    #     label=myFREE.label)
    plot_mean_rk_vs_pressure(ax, pressure, mean_curvatures, mean_radii, myFREE.label, 'blue', linestyle='solid')
    legend_without_duplicate_labels(ax)

    # plot twist angle if FREE has asymmetric fiber angles
    # if myFREE.fiber_a_angle_deg != -myFREE.fiber_b_angle_deg and not myFREE.has_def_fib_angle:
        # ax2 = ax.twinx()
        # plot_delta_vs_pressure(ax2, pressure, delta)
        # legend_without_duplicate_labels(ax2)

    # # subplot 4
    ax = fig.add_subplot(2, 2, 4)
    # plot_rk_vs_norm_position(ax, norm_positions, rk, frames=[0, -1], colors=['blue', 'red'])
    # plot_rk_vs_norm_position(ax, norm_positions, rk, pressure, frames=plot_frames)
    # plt.vlines(x=norm_positions[-1][[bulge_index, -bulge_index - 1]],
    #            ymin=ax.get_ylim()[0],
    #            ymax=ax.get_ylim()[1],
    #            color='gray',
    #            linestyles='dashed',
    #            label='endcap bulge boundary')
    # plot_total_dist_force_vs_pressure(
    #     ax=ax,
    #     pressure=pressure,
    #     force=total_dist_force,
    #     color='blue',
    #     label=myFREE.label,
    #     linestyle='solid')
    # plot_total_dist_force_vs_pressure(
    #     ax=ax,
    #     pressure=pressure,
    #     force=np.apply_along_axis(np.min, 1, d_dist_force) / pressure,
    #     color='blue',
    #     label=myFREE.label,
    #     linestyle='solid')

    # legend_without_duplicate_labels(ax)
    '''

    plt.savefig(fname=os.path.join(output_directory, myFREE.label + '.png'), dpi='figure', format='png', metadata=None, bbox_inches=None,
                pad_inches=0.1, facecolor='auto', edgecolor='auto')
    plt.savefig(fname=os.path.join(output_directory, myFREE.label + '.svg'), dpi='figure', format='svg', metadata=None, bbox_inches=None,
                pad_inches=0.1, facecolor='auto', edgecolor='auto')


    plt.show()

# specify field headers for output file
def define_output_fields():
    fields = ['Frame',
              'Pressure',
              'Index',
              'Smooth CL X',
              'Smooth CL Y',
              'Smooth CL Z',
              'Radius',
              'Length',
              'Alpha (inside)',
              'Alpha (midline)',
              'Alpha (outside)',
              'Beta (inside)',
              'Beta (midline)',
              'Beta (outside)',
              'Curvature']

    return fields


# run main()
if __name__ == "__main__":
    # batch process FREEs
    # all_labels = [
    #     'arc_20-330_r5_a25_b0_DFA',
    #     'arc_20-330_r5_a30_b0_DFA',
    #     'arc_20-330_r5_a40_b0_DFA',
    #
    #     'arc_30-220_r5_a25_b0_DFA',
    #     'arc_30-220_r5_a30_b0_DFA',
    #     'arc_30-220_r5_a40_b0_DFA',
    #
    #     'arc_60-110_r5_a25_b0_DFA',
    #     'arc_60-110_r5_a30_b0_DFA',
    #     'arc_60-110_r5_a40_b0_DFA',
    #
    #     'arc_20-330_r8_a25_b0_DFA',
    #     'arc_20-330_r8_a30_b0_DFA',
    #     'arc_20-330_r8_a40_b0_DFA',
    #
    #     'arc_30-220_r8_a25_b0_DFA',
    #     'arc_30-220_r8_a30_b0_DFA',
    #     'arc_30-220_r8_a40_b0_DFA',
    #
    #     'arc_60-110_r8_a25_b0_DFA',
    #     'arc_60-110_r8_a30_b0_DFA',
    #     'arc_60-110_r8_a40_b0_DFA',
    #
    #     'str_115_r5_a25_b-25',
    #     'str_115_r5_a30_b-30',
    #     'str_115_r5_a40_b-40',
    #
    #     'arc_20-330_r5_a25_b-25',
    #     'arc_20-330_r5_a30_b-30',
    #     'arc_20-330_r5_a40_b-40',
    #
    #     'arc_30-220_r5_a25_b-25',
    #     'arc_30-220_r5_a30_b-30',
    #     'arc_30-220_r5_a40_b-40',
    #
    #     'arc_60-110_r5_a25_b-25',
    #     'arc_60-110_r5_a30_b-30',
    #     'arc_60-110_r5_a40_b-40'
    # ]
    #
    # for label in all_labels:
    #     myFREE = initialize_FREE(label=label)
    #     main(myFREE)

    main(myFREE)
