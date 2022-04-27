# user-defined FREE and FEA parameters
# this is the input file for curved_FREEs FEA scripts

import numpy as np


# FREE.path_type options:
#     'straight'
#     'arc'
#     'quadratic_Bezier'
#     'cubic_Bezier'
# Modify path parameters or add new paths in geometry_setup.define_FREE_path.

# User-defined FREE parameters. This is where the FREE is defined by the user before analysis.
def apply_user_input(newFREE):
    # FREE parameters
    newFREE.has_def_fib_angle = True  # FREE has definite fiber angle, including constant alpha
    newFREE.path_type = 'cubic_Bezier'  # FREE path type. See comment above.
    newFREE.label_addon = ''  # customizable label addon (default is '')
    newFREE.material = 'DS10 (NeoHk, Polygerinos 2019)'  # see abaqus_setup.define_Abaqus_materials() for options
    # newFREE.material = 'DS10 (Yeoh, Caasenbrood 2020)'
    # newFREE.material = 'DS30 (Yeoh, Chen 2018)'
    # newFREE.material = 'DS30 (M-R, Elsayed 2014)'
    newFREE.outer_radius = 5.0  # FREE outside radius (mm)
    newFREE.thickness = 2.0  # FREE thickness (mm)
    newFREE.fiber_a_angle_deg = 30.0  # initial fiber angle for fiber direction A, i.e. alpha OR for DFA, midline angle (deg)
    newFREE.fiber_b_angle_deg = 0.0  # initial fiber angle for fiber direction B, i.e. beta OR for DFA, offset angle (deg)
    newFREE.num_fibers = 8  # number of fibers in each direction
    newFREE.endcap_thickness = 1.0  # endcap thickness (mm)
    newFREE.final_pressure = 0.150  # FREE pressure (MPa)

    # Path parameters
    newFREE.num_path_points = 101  # Number of points along FREE path. Affects path fidelity and output resolution.

    # Abaqus mesh parameters
    newFREE.bladder_seed_size = 0.8  # bladder seed size (default = 0.8)
    newFREE.fiber_seed_size = 0.75 * newFREE.bladder_seed_size  # fiber seed size (default = 0.75 * newFREE.bladder_seed_size)
    newFREE.endcap_seed_size = 1.25 * newFREE.bladder_seed_size  # endcap seed size (default = 1.2 * newFREE.bladder_seed_size)
    newFREE.num_output_nodes_circ = 36  # number of output nodes in bladder mesh circumference
    newFREE.fiber_element = 'truss'  # Fiber element type. Use 'beam' when alpha = -beta. Use 'truss' for alpha != -beta.

    # Set path parameters. Change local variables to modify path. Do not change FREE.path_parameters definition.
    if newFREE.path_type == 'straight':  # format: [length]
        FREE_length = 115  # Free length (mm)

        newFREE.path_parameters = np.array([FREE_length])

    elif newFREE.path_type == 'arc':  # 'arc': [arc_radius(mm), arc_angle(deg)]
        arc_rad = 30      # arc radius (mm) | 20   30   32   48   60
        arc_ang_deg = 180  # arc angle (deg) | 330  220  206  138  110

        newFREE.path_parameters = np.array([arc_rad, arc_ang_deg])

    elif newFREE.path_type == 'quadratic_Bezier':  # 'quadratic_Bezier': [[P0], [P1], [P2]] where Pi are Bezier control points in the form np.array([x, y, z)
        # Bezier control points
        P0 = np.array([0, 0, 0])
        P1 = np.array([40, 40, 40])
        P2 = np.array([0, 60, 0])

        newFREE.path_parameters = np.array([P0, P1, P2])

    elif newFREE.path_type == 'cubic_Bezier':  #  'cubic_Bezier': [[P0], [P1], [P2], [P3]] where Pi are Bezier control points in the form np.array([x, y, z)
        # user-defined Bezier control points
        P0 = np.array([0, 0, 0])
        P1 = np.array([20, 60, 20])
        P2 = np.array([40, -40, 20])
        P3 = np.array([80, 20, 0])

        newFREE.path_parameters = np.array([P0, P1, P2, P3])

        # sample code for randomized control points
        # P0 = np.array([0, 0, 0])
        # P1 = np.array([random.uniform(-20, 20), random.uniform(-20, 20), random.uniform(-20, 20)])
        # P2 = np.array([random.uniform(-40, 40), random.uniform(-40, 40), random.uniform(-40, 40)])
        # P3 = np.array([random.uniform(-60, 60), random.uniform(-60, 60), random.uniform(-60, 60)])

    elif newFREE.path_type == 'spiral':  # 'from_2D_experiment': [[exp_points_x], [exp_points_y]] where exp_points are experimental marker data (mm)
        # spiral parameters
        spiral_rad = 10.0
        n_turns = 2.0
        spiral_pitch = 50.0  # must be > 2 * FREE.radius, should account for radial expansion

        newFREE.path_parameters = np.array([spiral_rad, n_turns, spiral_pitch])

    elif newFREE.path_type == 'from_2D_experiment':  # 'from_2D_experiment': [[exp_points_x], [exp_points_y]] where exp_points are experimental marker data (mm)
        # user-defined key points from experimental results
        exp_points_x = np.array([0, 4.8081, 12.2897, 22.0594, 33.0648, 43.1325, 51.5768, 57.6450, 60.7221])
        exp_points_y = np.array([0, 12.1111, 21.7737, 27.5226, 29.0152, 26.6876, 21.0902, 12.3598, 0.1530])

        newFREE.path_parameters = np.array([exp_points_x, exp_points_y])

    else:
        raise ValueError('FREE path type not recognized.')

    return newFREE