# This file contains all functions needed to set up generalized curved FREEs in Python

# import necessary packages
import math
import numpy as np
import warnings
from scipy.interpolate import interp1d


# define FREE path - FREE center axis will follow this discretized path
def generate_FREE_path(FREE):
    t = np.linspace(0, 1, FREE.num_path_points)  # normalized path position (always from 0 to 1)
    path = np.zeros((FREE.num_path_points, 3))  # initialize FREE_path (each row is [x, y, z] for one point)

    # FREE.path_parameters format:
    #     'straight': [length]
    #     'arc': [arc_radius(mm), arc_angle(deg)]
    #     'quadratic_Bezier': [[P0], [P1], [P2]] where Pi are Bezier control points in the form np.array([x, y, z)
    #     'cubic_Bezier': [[P0], [P1], [P2], [P3]] where Pi are Bezier control points in the form np.array([x, y, z)
    #     'from_experiment': [[exp_points_x], [exp_points_x], [exp_points_x]] where exp_points are experimental marker data (mm)

    if FREE.path_type == 'straight':
        # user-defined path parameters
        FREE_length = float(FREE.path_parameters[0])

        for i in range(FREE.num_path_points):
            x = 0.
            y = 0
            z = i * FREE_length / (FREE.num_path_points - 1)

            path[i] = np.array([x, y, z])

    elif FREE.path_type == 'arc':
        # user-defined arc parameters
        arc_rad = FREE.path_parameters[0]  # arc radius (mm)
        arc_ang_deg = FREE.path_parameters[1]  # arc angle (deg)

        # calculated arc path
        arc_ang = arc_ang_deg * math.pi / 180  # arc angle (rad)
        for i in range(FREE.num_path_points):
            x = arc_rad * (1 - np.cos(t[i] * arc_ang))
            y = 0.
            z = arc_rad * np.sin(t[i] * arc_ang)
            path[i] = np.array([x, y, z])

    elif FREE.path_type == 'quadratic_Bezier':
        # user-defined Bezier control points
        P0 = FREE.path_parameters[0]
        P1 = FREE.path_parameters[1]
        P2 = FREE.path_parameters[2]
        # bez_all_p = np.column_stack((P0, P1, P2))  # used for checking control points in a 3D plot

        # calculated Bezier curve path
        for i in range(FREE.num_path_points):
            path[i] = P0 * (1 - t[i]) ** 2 \
                      + P1 * (2 * (1 - t[i]) * t[i]) \
                      + P2 * t[i] ** 2

    elif FREE.path_type == 'cubic_Bezier':
        # user-defined Bezier control points
        P0 = FREE.path_parameters[0]
        P1 = FREE.path_parameters[1]
        P2 = FREE.path_parameters[2]
        P3 = FREE.path_parameters[3]
        # bez_all_p = np.column_stack((bez_p0, bez_p1, bez_p2, bez_p3))  # used for checking control points in a 3D plot

        # calculated Bezier curve path
        for i in range(FREE.num_path_points):
            path[i] = P0 * (1 - t[i]) ** 3 \
                      + P1 * (3 * (1 - t[i]) ** 2 * t[i]) \
                      + P2 * (3 * (1 - t[i]) * t[i] ** 2) \
                      + P3 * t[i] ** 3

    elif FREE.path_type == 'spiral':
        # user-defined arc parameters
        spiral_rad = FREE.path_parameters[0]  # spiral radius (mm)
        n_turns = FREE.path_parameters[1]  # number of spiral turns
        spiral_pitch = FREE.path_parameters[2]  # spiral pitch (dist. between coils) (mm)

        # calculated arc path
        spiral_ang = n_turns * 2 * math.pi  # swept angle (rad)
        for i in range(FREE.num_path_points):
            x = spiral_rad * (1 - np.cos(t[i] * spiral_ang))
            y = spiral_pitch * n_turns * i / (FREE.num_path_points - 1)
            z = spiral_rad * np.sin(t[i] * spiral_ang)
            path[i] = np.array([x, y, z])

    elif FREE.path_type == 'from_2D_experiment':
        # user-defined key points from experimental results
        exp_points_x = FREE.path_parameters[0]
        exp_points_y = FREE.path_parameters[1]

        f = interp1d(exp_points_x, exp_points_y, kind='cubic')  # interpolate experimental data points

        x_new = np.linspace(0, 1, num=FREE.num_path_points) * exp_points_x[-1]
        y_new = f(x_new)
        z_new = np.zeros(FREE.num_path_points)

        path = np.transpose(np.array([x_new, y_new, z_new]))

    else:
        raise ValueError('FREE path type not recognized.')

    return path


# manually modify path - for testing code
def modify_path(path):
    # add straight sections - used for testing local frame code
    path = np.append([[0, 0, -20], [0, 0, -10]], path, axis=0)
    path = np.append(path, [[60, 0, -10], [60, 0, -20]], axis=0)
    path = np.append(path, [[70, 0, -20]], axis=0)
    path = np.append(path, path + [80, 0, 0], axis=0)

    # rotate path - used for testing reorientation code
    R = rotation_matrix([1, 1, 1], -math.pi / 4)
    path = apply_rotation(path, R)

    return path


# calculate approximate path length (point-to-point)
def find_path_length(path):
    num_path_points = np.shape(path)[0]

    path_length = 0.

    for i in range(num_path_points - 1):
        point_1 = path[i]
        point_2 = path[i + 1]
        path_length += np.linalg.norm(point_2 - point_1)

    return path_length


# return position of points along a path (0 <= s <= L)
def get_s_array(path):
    num_path_points = np.shape(path)[0]

    s = np.zeros([num_path_points])

    for i in range(1, num_path_points):
        point_1 = path[i - 1]
        point_2 = path[i]
        s[i] = s[i - 1] + np.linalg.norm(point_2 - point_1)

    return s


# takes in a path and returns the local radius of curvature and curvature center for each point
def find_local_curves(path):
    num_path_points = np.shape(path)[0]

    T, N, B = find_TNB_frame(path)

    # get normalized positions along curve for all points
    s = get_s_array(path)

    dT = np.apply_along_axis(np.gradient, 0, T, s)

    curvature = np.apply_along_axis(np.linalg.norm, 1, dT)
    if np.amin(curvature) < 1e-6:  # if there are straight segments, cap radius of curvature at 1e6
        radius_of_curvature = 1e6 * np.ones(num_path_points)
        for i in range(num_path_points):
            if curvature[i] > 1e-6:
                radius_of_curvature[i] = 1 / curvature[i]
    else:
        radius_of_curvature = 1 / curvature
    radius_vectors = np.transpose(np.multiply(np.transpose(N), radius_of_curvature))
    curve_centers = np.add(path, radius_vectors)

    return curvature, radius_of_curvature, curve_centers


# check for curvature that is too tight for FREE radius
def check_FREE_curvature(path, FREE_outer_radius):
    path_radius = find_local_curves(path)[1]  # find all local radii
    if FREE_outer_radius > np.amin(path_radius):
        raise ValueError('Self-intersecting geometry: FREE outer radius exceeds minimum path radius.')


# takes in a 3D path and returns a Frenet (TNB) frame, i.e. local reference frame, for each point
# reference: https://janakiev.com/blog/framing-parametric-curves/
def find_TNB_frame(path):
    num_path_points = np.shape(path)[0]

    # get normalized positions along curve for all points
    s = get_s_array(path)

    # Calculate the first and second derivative of the points, taking uneven spacing into consideration
    dX = np.apply_along_axis(np.gradient, 0, path, s)
    ddX = np.apply_along_axis(np.gradient, 0, dX, s)

    # Normalize all tangents
    normalize = lambda m: m / np.linalg.norm(m)
    T = np.apply_along_axis(normalize, axis=1, arr=dX)

    # Calculate and normalize all binormals
    B = np.cross(dX, ddX)
    # if there are any straight sections, use previous binormal
    for i in range(num_path_points):
        if np.linalg.norm(B[i]) < 1e-3:
            if i == 0:
                B[i] = np.array([0, 1, 0])
            else:
                B[i] = B[i - 1]
    B = np.apply_along_axis(normalize, axis=1, arr=B)

    # Calculate all normals
    N = np.cross(B, T)

    # adjust starting tangent based on starting curvature
    dT = np.apply_along_axis(np.gradient, 0, T, s)
    if np.linalg.norm(dT[0]) > 1e-3:  # if the beginning of the FREE is not straight
        start_curvature = np.linalg.norm(dT[2])
        start_radius_of_curvature = 1 / start_curvature
        start_curve_center = np.add(path[2], start_radius_of_curvature * N[2])
        N[0] = normalize(start_curve_center - path[0])
        B[0] = normalize(np.cross((path[2] - path[0]), N[0]))
        T[0] = normalize(np.cross(N[0], B[0]))

    # adjust ending tangent based on ending curvature
    if np.linalg.norm(dT[-1]) > 1e-3:  # if the end of the FREE is not straight
        end_curvature = np.linalg.norm(dT[-3])
        end_radius_of_curvature = 1 / end_curvature
        end_curve_center = np.add(path[-3], end_radius_of_curvature * N[-3])
        N[-1] = normalize(end_curve_center - path[-1])
        B[-1] = normalize(np.cross((path[-1] - path[-3]), N[-1]))
        T[-1] = normalize(np.cross(N[-1], B[-1]))

    return T, N, B


# takes in a 3D path and returns normalized tangent vectors
# accounts for unevenly spaced points
def find_normalized_tangents(path):
    # get normalized positions along curve for all points
    s = get_s_array(path)

    # Calculate the first and second derivative of the points, taking uneven spacing into consideration
    dX = np.apply_along_axis(np.gradient, 0, path, s)

    # Normalize all tangents
    normalize = lambda m: m / np.linalg.norm(m)
    T = np.apply_along_axis(normalize, axis=1, arr=dX)

    return T


# takes in a 3D path and returns a parallel transport frame, i.e. local reference frame, for each point
# reference: https://janakiev.com/blog/framing-parametric-curves/
def find_parallel_frame(path):
    num_path_points = np.shape(path)[0]
    T, N, B = find_TNB_frame(path)

    # Initialize the first parallel-transported normal vector V
    V = np.zeros(np.shape(path))
    V[0] = N[0]
    V[0] = V[0] / np.linalg.norm(V[0])

    # Compute the values for V for each tangential vector from T
    for i in range(num_path_points - 1):
        b = np.cross(T[i], T[i + 1])
        if np.linalg.norm(b) < 1e-3:  # if the point is on a straight segment...
            V[i + 1] = V[i]
        else:
            b = b / np.linalg.norm(b)
            phi = np.arccos(np.dot(T[i], T[i + 1]))
            R = rotation_matrix(b, phi)
            V[i + 1] = np.dot(R, V[i])

    # Calculate the second parallel-transported normal vector U
    U = np.array([np.cross(t, v) for (t, v) in zip(T, V)])

    return T, U, V


# takes in a path and the TNB frame of the path's first point and reorients the
# path such that it initially points along the z-axis and curves about the y-axis
def reorient_path(path, T0, N0, B0):
    reorient_tol = 1e-6  # tolerance for reorienting path
    new_path = path  # initialize new path

    # shift path so path[0] is at [0, 0, 0] if not already
    if np.linalg.norm(path[0]) > reorient_tol:
        new_path = path - path[0]

    # rotate path to align N0, B0, T0 with x, y, z axes, unless already aligned
    if (np.linalg.norm(np.cross(T0, [0, 0, 1])) > reorient_tol or
            np.linalg.norm(np.cross(N0, [1, 0, 0])) > reorient_tol):
        R = np.array([N0, B0, T0])
        new_path = apply_rotation(new_path, R)

    return new_path


# Return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.
# reference: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


# takes in a path of one or more points and rotates them using rotation matrix R
def apply_rotation(points, R):
    new_points = np.zeros(np.shape(points))

    if np.shape(np.shape(points))[0] == 1:  # if applying to a single point...
        new_points = np.dot(R, points)
    else:  # if applying to a series of points...
        for i in range(np.shape(points)[0]):
            new_points[i] = np.dot(R, points[i])

    return new_points
