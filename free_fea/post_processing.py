# contains functions for post-processing curved FREE FEA output in Python

# import necessary packages
import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import Delaunay
from scipy.spatial import cKDTree
from scipy.interpolate import InterpolatedUnivariateSpline
import pandas as pd
import networkx as nx
from free_fea.geometry_setup import *


# extract FREE centerline coordinates for all frames from bladder nodes
def get_centerlines_from_bladder(bladder_coordinates):
    num_frames = np.shape(bladder_coordinates)[0]
    num_path_points = np.shape(bladder_coordinates)[1]

    centerlines = np.zeros([num_frames, num_path_points, 3])

    for f in range(num_frames):

        for i in range(num_path_points):
            points_x = bladder_coordinates[f, i, :, 0]
            points_y = bladder_coordinates[f, i, :, 1]
            points_z = bladder_coordinates[f, i, :, 2]

            centerline_x = np.mean(points_x)
            centerline_y = np.mean(points_y)
            centerline_z = np.mean(points_z)

            centerlines[f][i] = np.array([centerline_x, centerline_y, centerline_z])

    return centerlines


# find all radii
def get_all_radii(centerlines, fiber_coordinates):
    num_frames = np.shape(fiber_coordinates)[0]
    num_path_points = np.shape(fiber_coordinates)[1]
    num_fibers = np.shape(fiber_coordinates)[2]

    radii = np.zeros([num_frames, num_path_points])
    radius_to_fibers = np.zeros(num_fibers)

    for f in range(num_frames):
        for i in range(num_path_points):
            for j in range(num_fibers):
                radius_to_fibers[j] = np.linalg.norm(fiber_coordinates[f][i][j] - centerlines[f][i])

            radii[f][i] = np.mean(radius_to_fibers)

    return radii


# find all centerline lengths
def get_lengths(centerlines):
    num_frames = np.shape(centerlines)[0]

    lengths = np.zeros([num_frames])

    for f in range(num_frames):
        lengths[f] = find_path_length(centerlines[f])

    return lengths


# find all centerline radii of curvature
# centerlines[f][i] = centerline point at s = s[i] in frame f
def get_all_curves(centerlines):
    num_frames = np.shape(centerlines)[0]
    num_path_points = np.shape(centerlines)[1]

    RoCs = np.zeros([num_frames, num_path_points])
    curvatures = np.zeros([num_frames, num_path_points])

    for f in range(num_frames):
        # curvatures[f][i] = centerline curvature at s = s[i] in frame f
        # RoCs[f][i] = centerline radius of curvature at s = s[i] in frame f
        curvatures[f], RoCs[f], x = find_local_curves(centerlines[f])

    return curvatures, RoCs


# calculate fiber angles using all fibers' coordinates
# fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0),
# midline (t = 1) and outside (t = 2) of the curvature
def get_fiber_angles(centerlines, fiber_coordinates):
    num_frames = np.shape(fiber_coordinates)[0]
    num_path_points = np.shape(fiber_coordinates)[1]
    num_fibers = np.shape(fiber_coordinates)[2]

    fiber_angles = np.zeros([num_frames, num_path_points, 3])

    for f in range(num_frames):

        parallel_T, parallel_U, parallel_V = \
            find_parallel_frame(centerlines[f])  # calculate parallel frame for FREE path in frame f

        fiber_T = np.zeros([num_path_points, num_fibers, 3])

        for j in range(num_fibers):
            # calculate parallel frame for FREE path in frame f
            fiber_T[:, j, :] = find_parallel_frame(fiber_coordinates[f, :, j])[0]

        for i in range(num_path_points):
            # local frame, based on parallel transport frame to minimize twist
            e1 = parallel_V[i]
            e2 = parallel_U[i]
            e3 = parallel_T[i]

            R_T_parallel = np.array(
                [e1, e2, e3])  # rotation matrix that maps global coordinates to local parallel frame

            local_fiber_angles = np.zeros([num_fibers])

            for j in range(num_fibers):
                local_fiber_angles[j] = math.acos(np.dot(e3, fiber_T[i][j]))  # fiber angle at fiber_coordinates[f,i,j]

                # adjust sign of fiber angle depending on fiber direction
                radial_vector = fiber_coordinates[f][i][j] - centerlines[f][i]
                local_radial_vector = apply_rotation(radial_vector, R_T_parallel)
                local_fiber_T = apply_rotation(fiber_T[i][j], R_T_parallel)
                proj_fiber_T = np.array([local_fiber_T[0], local_fiber_T[1], 0])
                cross_prod_z = np.cross(local_radial_vector, proj_fiber_T)[2]

                if cross_prod_z < 0:
                    local_fiber_angles[j] = -local_fiber_angles[j]

            fiber_angles[f, i, 0] = np.min(local_fiber_angles)
            fiber_angles[f, i, 1] = np.median(local_fiber_angles)
            fiber_angles[f, i, 2] = np.max(local_fiber_angles)

    fiber_angles = fiber_angles * 180 / math.pi

    return fiber_angles


# take in a path and smooth for more consistent curvature using a Savitzky-Golay filter
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
# path[f][i] = np.array([x, y, z]) at s = s[i] in frame f
# window_length = length of the filter window (i.e., the number of coefficients)
# polyorder = order of the polynomial used to fit the samples
# mode = type of extension to use for the padded signal
def smooth_path_savgol(path, window_length=21, polyorder=3, mode='interp'):
    smooth_path = np.zeros(np.shape(path))

    if len(np.shape(path)) == 3:
        for f in range(np.shape(path)[0]):
            x = path[f, :, 0]
            y = path[f, :, 1]
            z = path[f, :, 2]

            smooth_x = savgol_filter(x, window_length=window_length, polyorder=polyorder, mode=mode)
            smooth_y = savgol_filter(y, window_length=window_length, polyorder=polyorder, mode=mode)
            smooth_z = savgol_filter(z, window_length=window_length, polyorder=polyorder, mode=mode)

            smooth_path[f, :, 0] = smooth_x
            smooth_path[f, :, 1] = smooth_y
            smooth_path[f, :, 2] = smooth_z

    elif len(np.shape(path)) == 4:
        for f in range(np.shape(path)[0]):
            for j in range(np.shape(path)[2]):
                x = path[f, :, j, 0]
                y = path[f, :, j, 1]
                z = path[f, :, j, 2]

                smooth_x = savgol_filter(x, window_length=window_length, polyorder=polyorder, mode=mode)
                smooth_y = savgol_filter(y, window_length=window_length, polyorder=polyorder, mode=mode)
                smooth_z = savgol_filter(z, window_length=window_length, polyorder=polyorder, mode=mode)

                smooth_path[f, :, j, 0] = smooth_x
                smooth_path[f, :, j, 1] = smooth_y
                smooth_path[f, :, j, 2] = smooth_z

    return smooth_path


# return smoothing spline fit to input path
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
# path[f][i] = np.array([x, y, z]) at s = s[i] in frame f
# or path[f][i][j] = np.array([x, y, z]) at s = s[i] in frame f on path (fiber) j
# k = degree of smoothing spline, must be 1 <= k <= 5
def smooth_path_spline(path, k=5, num_points=20):
    smooth_path = np.zeros(np.shape(path))

    if len(np.shape(path)) == 3:  # this loop is for a set of single paths (e.g. centerlines)
        for f in range(np.shape(path)[0]):
            num_path_points = np.shape(path[f])[0]
            s = get_s_array(path[f])

            x = path[f, :, 0]
            y = path[f, :, 1]
            z = path[f, :, 2]

            indices = np.concatenate(
                (np.round(np.arange(0, num_path_points, num_path_points / (num_points - 1))),
                np.array([num_path_points - 1])))
            indices = indices.astype(int)

            x_spline = InterpolatedUnivariateSpline(s[indices], x[indices], k=k)
            y_spline = InterpolatedUnivariateSpline(s[indices], y[indices], k=k)
            z_spline = InterpolatedUnivariateSpline(s[indices], z[indices], k=k)

            smooth_path[f, :, 0] = x_spline(s)
            smooth_path[f, :, 1] = y_spline(s)
            smooth_path[f, :, 2] = z_spline(s)

    elif len(np.shape(path)) == 4:  # this loop for a set of multiple paths (e.g. fibers)
        for f in range(np.shape(path)[0]):
            for j in range(np.shape(path)[2]):
                num_path_points = np.shape(path[f, :, j])[0]
                s = get_s_array(path[f, :, j])

                x = path[f, :, j, 0]
                y = path[f, :, j, 1]
                z = path[f, :, j, 2]

                indices = np.concatenate(
                    (np.round(np.arange(0, num_path_points, num_path_points / (num_points - 1))),
                     np.array([num_path_points - 1])))
                indices = indices.astype(int)

                x_spline = InterpolatedUnivariateSpline(s[indices], x[indices], k=k)
                y_spline = InterpolatedUnivariateSpline(s[indices], y[indices], k=k)
                z_spline = InterpolatedUnivariateSpline(s[indices], z[indices], k=k)

                smooth_path[f, :, j, 0] = x_spline(s)
                smooth_path[f, :, j, 1] = y_spline(s)
                smooth_path[f, :, j, 2] = z_spline(s)

    return smooth_path


# returns array of normalized positions along a path
# path[f][i] = np.array([x, y, z]) at s = s[i] in frame f
# norm_positions[f][i] = position s along path with 0 <= s <= 1
def get_normalized_positions(path):
    norm_positions = np.zeros(np.shape(path)[0:2])

    for f in range(np.shape(path)[0]):  # loop through frames
        path_length = find_path_length(path[f])
        s = get_s_array(path[f])
        norm_positions[f] = s / path_length

    return norm_positions


# find the index for trimming off the FREE endcap bulge using the angle derivative w.r.t. path index
# fiber_angles[f, i, t] = fiber angle in frame f at index i on the inside (t = 0),
# midline (t = 1) and outside (t = 2) of the curvature
# tol = tolerance for moving average of the gradient of the midline fiber angles
def find_endcap_bulge_index(fiber_a_angles, fiber_b_angles, tol=0.3):
    num_path_points = np.shape(fiber_a_angles)[1]

    # gradient of the fiber angle w.r.t. index in the last (deformed) frame
    gradient_a = np.gradient(fiber_a_angles[-1, :, 1])
    gradient_b = np.gradient(fiber_b_angles[-1, :, 1])

    bulge_index = 0

    mid_index = int(num_path_points / 2)
    for i in list(range(2, mid_index)):  # trim 2 end nodes, then loop through indices toward the middle of the path
        # apply moving average to both ends of both fiber sets
        i_start = [i, i + 1, i + 2]
        i_end = [-i - 3, -i - 2, -i - 1]
        a_mean = np.array([np.mean(gradient_a[i_start]), np.mean(gradient_a[i_end])])
        b_mean = np.array([np.mean(gradient_b[i_start]), np.mean(gradient_b[i_end])])

        # find steepest slope between the two ends for each fiber group
        a_max = abs(a_mean).max()
        b_max = abs(b_mean).max()

        if a_max < tol and b_max < tol:  # if both ends' fiber angles have stopped changing rapidly
            bulge_index = i + 1  # trimmed_array = array[bulge_index : len(array)-bulge_index]
            break

    if i == mid_index - 1:
        warnings.warn('No end bulging detected.')

    return bulge_index


# takes in lengths and radii for all frames and returns stretch ratios
def get_stretch_ratios(lengths, radii):
    num_frames = np.shape(lengths)[0]

    mean_radii = np.zeros(num_frames)
    lambda_1 = np.zeros(num_frames)
    lambda_2 = np.zeros(num_frames)

    for f in range(num_frames):
        mean_radii[f] = np.mean(radii[f, :])
        lambda_1[f] = lengths[f] / lengths[0]
        lambda_2[f] = mean_radii[f] / mean_radii[0]

    return lambda_1, lambda_2


def get_twist(centerlines, endcap_coordinates):
    num_frames = np.shape(centerlines)[0]
    delta = np.zeros(num_frames)

    theta = np.zeros(num_frames)
    local_x = np.zeros(num_frames)
    local_y = np.zeros(num_frames)
    for f in range(1, num_frames):
        # calculate parallel frame and endpoint of FREE path
        parallel_T, parallel_U, parallel_V = find_parallel_frame(centerlines[f][:])

        # local frame at centerline end point, based on parallel transport frame to minimize twist
        e1 = parallel_V[-1]
        e2 = parallel_U[-1]
        e3 = parallel_T[-1]

        # rotation matrix that maps global coordinates to local parallel frame
        R_T_parallel = np.array([e1, e2, e3])

        twist_vector = endcap_coordinates[f][1] - endcap_coordinates[f][0]
        twist_vector = twist_vector / np.linalg.norm(twist_vector)  # normalize twist vector in case length > 1
        theta[f] = math.acos(np.dot(e1, twist_vector))

        local_twist_vector = apply_rotation(twist_vector, R_T_parallel)
        local_x[f] = local_twist_vector[0]
        local_y[f] = local_twist_vector[1]

        # ensure theta is in the correct quadrant for fiber angle calculation
        if local_y[f] < 0:
            theta[f] = 2 * math.pi - theta[f]

        # account for crossing the local d1 (x) axis
        new_theta = theta[f]
        old_theta = theta[f - 1]
        # if we cross the e1 axis moving counterclockwise...
        if local_y[f] >= 0 and local_y[f - 1] < 0 and local_x[f] > 0:
            old_theta = old_theta - 2 * math.pi
        # if we cross the e1 axis moving clockwise...
        elif local_y[f] < 0 and local_y[f - 1] >= 0 and local_x[f] > 0:
            old_theta = old_theta + 2 * math.pi

        d_theta = new_theta - old_theta

        delta[f] = delta[f - 1] + d_theta

    delta = delta * 180 / math.pi  # convert to degrees

    return delta


# calculate projected surface area on inside and outside of bladder curvature
# centerlines[f][i] = np.array([x, y, z]) on centerline c(s) at s = s[i] in frame f
# bladder_coordinates[f][i][j] = point at INNER bladder circumferential index j at s = s[i] in frame f
# pressure[f] = pressure in frame f
def calculate_distributed_force(centerlines, bladder_coordinates, pressure):
    num_frames = np.shape(centerlines)[0]
    num_path_points = np.shape(centerlines)[1]
    num_circ_points = np.shape(bladder_coordinates)[2]

    d_dist_force = np.zeros([num_frames, num_path_points])
    total_net_dist_force = np.zeros(num_frames)
    inward_dist_force = np.zeros(num_frames)
    outward_dist_force = np.zeros(num_frames)

    for f in range(num_frames):
        P = pressure[f] / 1e3  # convert pressure from kPa to MPa
        T, N, B = find_TNB_frame(centerlines[f])

        for i in range(num_path_points):

            ring_net_proj_area = 0  # net projected area for the ring surrounding c(s)
            for j in range(num_circ_points):  # loop through circumferential points
                # element axial length
                if i == 0:  # if at the start of the path, only use points ahead
                    i_f = i + 1
                    i_b = i
                elif i == num_path_points - 1:  # if at the end of the path, only use points behind
                    i_f = i
                    i_b = i - 1
                else:  # otherwise use points ahead and behind
                    i_f = i + 1
                    i_b = i - 1
                elem_ax_vec = bladder_coordinates[f][i_f][j] - bladder_coordinates[f][i_b][j]
                elem_ax_len = np.linalg.norm(elem_ax_vec) / 2

                # element circumferential length
                if j == 0:
                    j_f = 1
                    j_b = -1
                elif j == num_circ_points - 1:
                    j_f = 0
                    j_b = -2
                else:
                    j_f = j + 1
                    j_b = j - 1
                elem_circ_vec = bladder_coordinates[f][i][j_f] - bladder_coordinates[f][i][j_b]
                elem_circ_len = np.linalg.norm(elem_circ_vec) / 2

                elem_area = elem_ax_len * elem_circ_len

                # reverse sign of element surface normal if pointing toward centerline
                x_prod = np.cross(elem_ax_vec, elem_circ_vec)  # local element surface normal
                elem_surf_norm = x_prod / np.linalg.norm(x_prod)  # normalized surface normal
                radial_vector = bladder_coordinates[f, i, j] - centerlines[f, i]
                if np.dot(radial_vector, elem_surf_norm) < 0:
                    elem_surf_norm = -elem_surf_norm

                # projected area (+ = toward inside of curvature, - = toward outside)
                elem_proj_area = elem_area * np.dot(elem_surf_norm, N[i])

                ring_net_proj_area += elem_proj_area

            d_dist_force[f][i] = P * ring_net_proj_area

        total_net_dist_force[f] = sum(d_dist_force[f])

    return d_dist_force, total_net_dist_force,

