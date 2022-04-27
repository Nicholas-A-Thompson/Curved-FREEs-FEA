# Build curved FREE FEA model in Abaqus.
# This file can only be run in the Abaqus Python environment (Python 2.7).

# import the necessary curved FREE modules
import sys
import os
import time

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

# import custom modules
# navigate to output folder. Abaqus Work Directory and curved_FREEs must be in the same folder
current_directory = os.getcwd()  # get current work directory
parent_directory = os.path.split(current_directory)[0]  # go up one folder
new_path = os.path.join(parent_directory, 'curved_FREEs')  # navigate to curved_FREEs to import custom packages
if all(sys_path != new_path for sys_path in sys.path):
    sys.path.insert(1, new_path)

# reload modules to allow for changes without restarting Abaqus
import free_fea.initialization
reload(free_fea.initialization)
from free_fea.initialization import *

import free_fea.abaqus_setup
reload(free_fea.abaqus_setup)
from free_fea.abaqus_setup import *


########## USER INPUT ##########
# initialize FREE object
# label='user_input' for a new FREE from free_fea.user_input.
# label='[input file name]' to load a FREE from an input file
myFREE = initialize_FREE(label='user_input')
# myFREE = initialize_FREE(label='arc_60-110_r5_a40_b0_DFA')
################################

# main code block
def main():
    TNB_T, TNB_N, TNB_B = find_TNB_frame(myFREE.path)  # calculate TNB frame for FREE path after reorientation
    parallel_T, parallel_U, parallel_V = find_parallel_frame(myFREE.path)  # calculate parallel frame for FREE path

    # calculate FREE parameters for part creation
    output_nodes_separation_angle = 2 * math.pi / myFREE.num_output_nodes_circ  # angle difference between bladder mesh nodes (rad)
    fiber_a_angle = myFREE.fiber_a_angle_deg * math.pi / 180  # initial fiber angle for fiber direction A, i.e. alpha OR for CFA, midline angle (rad)
    fiber_b_angle = myFREE.fiber_b_angle_deg * math.pi / 180  # initial fiber angle for fiber direction B, i.e. beta OR for CFA, offset angle (rad)
    outer_radius = myFREE.outer_radius  # fiber helix radius (mm)
    fiber_separation_angle = 2 * math.pi / myFREE.num_fibers  # fiber separation angle (rad)

    # initialize fiber angular positions w.r.t. parallel frame
    fiber_a_theta_p = np.linspace(0, 2 * math.pi - fiber_separation_angle, myFREE.num_fibers)
    fiber_b_theta_p = np.linspace(0, 2 * math.pi - fiber_separation_angle, myFREE.num_fibers)

    # initialize fiber angular position w.r.t. Frenet / TNB frame
    fiber_a_theta_f = np.array(fiber_a_theta_p)
    fiber_b_theta_f = np.array(fiber_b_theta_p)

    # initialize fiber and output node arrays

    # Points on endcap used to track rotation. Starts as unit vector along local e1_p axis.
    endpoint_output_nodes = np.zeros([2, 3])

    # bladder_output_nodes[i][j] = point on inner bladder at s = s[i] and theta = j * d_theta
    bladder_output_nodes_in = np.zeros(
        [np.shape(myFREE.path)[0], myFREE.num_output_nodes_circ, np.shape(myFREE.path)[1]])

    # fiber_x_points[i][j] = point on fiber f at s = s[i] (direction x)
    fiber_a_points = np.zeros([np.shape(myFREE.path)[0], myFREE.num_fibers, np.shape(myFREE.path)[1]])
    fiber_b_points = np.zeros([np.shape(myFREE.path)[0], myFREE.num_fibers, np.shape(myFREE.path)[1]])

    curvature = find_local_curves(myFREE.path)[0]  # curvature[i] = curvature of the FREE path at index i

    # define fiber splines and output node positions
    for i in range(np.shape(myFREE.path)[0]):

        # local frame, based on parallel transport frame to minimize twist
        e1_p = parallel_V[i]
        e2_p = parallel_U[i]
        e3_p = parallel_T[i]

        R_p = np.array([e1_p, e2_p, e3_p]).T  # rotation matrix that maps local parallel frame coordinates to global

        # local frame, based on TNB / Frenet frame to keep track of curvature orientation
        e1_f = TNB_N[i]
        e2_f = TNB_B[i]
        e3_f = TNB_T[i]

        R_f = np.array([e1_f, e2_f, e3_f]).T  # rotation matrix that maps local Frenet / TNB frame coordinates to global

        # save last e1_p for endcap output positioning
        if i == np.shape(myFREE.path)[0] - 1:
            endpoint_output_nodes[0] = myFREE.path[i]
            endpoint_output_nodes[1] = e1_p + myFREE.path[i]

        # define bladder output nodes
        for j in range(myFREE.num_output_nodes_circ):
            point_in = np.array([math.cos(j * output_nodes_separation_angle),
                              math.sin(j * output_nodes_separation_angle),
                              0]) * (myFREE.outer_radius - myFREE.thickness)
            bladder_output_nodes_in[i][j] = myFREE.path[i] + apply_rotation(point_in, R_p)

        # define fiber points
        for j in range(myFREE.num_fibers):
            if i != 0:  # if not at the start of the path...
                dL = np.linalg.norm(myFREE.path[i] - myFREE.path[i - 1])  # length increment of previous path segment

                # determine the in-plane angular position of points on each fiber.
                if myFREE.has_def_fib_angle:  # definite fiber angle FREEs
                    # standard fiber angle assignment
                    # alpha = fiber_a_angle
                    # beta = fiber_b_angle

                    # experimental: vary fiber angle such that alpha = -beta = a - b * cos(theta)
                    alpha = fiber_a_angle - fiber_b_angle * (math.cos(fiber_a_theta_f[j]))
                    beta = -(fiber_a_angle - fiber_b_angle * (math.cos(fiber_b_theta_f[j])))

                    # theta_p for fiber group A
                    fiber_a_theta_p[j] += (math.tan(alpha)
                                           * (1.0 - outer_radius * curvature[i] * math.cos(fiber_a_theta_f[j]))
                                           * dL / outer_radius)  # increment fiber angular position, i.e. theta
                    # theta_p for fiber group B
                    fiber_b_theta_p[j] += (math.tan(beta)
                                           * (1.0 - outer_radius * curvature[i] * math.cos(fiber_b_theta_f[j]))
                                           * dL / outer_radius)  # increment fiber angular position, i.e. theta

                else:  # constant pitch FREEs
                    fiber_a_theta_p[j] += dL * math.tan(fiber_a_angle) / outer_radius  # increment fiber angular position, i.e. theta
                    fiber_b_theta_p[j] += dL * math.tan(fiber_b_angle) / outer_radius  # increment fiber angular position, i.e. theta

            # fiber_a_points[i][j] = point on fiber j at s = s[i] (direction A)
            point_p = np.array([math.cos(fiber_a_theta_p[j]),
                                math.sin(fiber_a_theta_p[j]),
                                0]) * outer_radius  # point in local frame w.r.t. parallel frame
            fiber_a_points[i][j] = myFREE.path[i] + apply_rotation(point_p, R_p)

            # calculate new theta_f (direction A)
            point_f = apply_rotation(point_p, np.dot(R_f.T, R_p))  # point in local frame w.r.t. Frenet / TNB frame
            rad = np.linalg.norm(point_f)
            if rad > outer_radius:
                point_f = point_f * outer_radius / rad  # eliminate small scaling errors
            fiber_a_theta_f[j] = math.acos(point_f[0] / outer_radius)
            # Note: quadrant correction not necessary since only cos(theta_f) is used.

            # fiber_b_points[i][j] = point on fiber j at s = s[i] (direction B)
            point_p = np.array([math.cos(fiber_b_theta_p[j]),
                                math.sin(fiber_b_theta_p[j]),
                                0]) * outer_radius
            fiber_b_points[i][j] = myFREE.path[i] + apply_rotation(point_p, R_p)

            # calculate new theta_f (direction B)
            point_f = apply_rotation(point_p, np.dot(R_f.T, R_p))  # point in local frame w.r.t. Frenet / TNB frame
            rad = np.linalg.norm(point_f)
            if rad > outer_radius:
                point_f = point_f * outer_radius / rad  # eliminate small scaling errors
            fiber_b_theta_f[j] = math.acos(point_f[0] / outer_radius)
            # Note: quadrant correction not necessary since only cos(theta_f) is used.

    # create Abaqus model
    create_Abaqus_model()

    # create Abaqus parts
    create_Abaqus_fibers(FREE=myFREE, fiber_coords_1=fiber_a_points, fiber_coords_2=fiber_b_points)
    create_Abaqus_bladder(path=myFREE.path, outer_radius=myFREE.outer_radius,
                          thickness=myFREE.thickness, bladder_seed_size=myFREE.bladder_seed_size)
    create_Abaqus_endcap(outer_radius=myFREE.outer_radius, FREE_thickness=myFREE.thickness,
                         endcap_thickness=myFREE.endcap_thickness, endcap_seed_size=myFREE.endcap_seed_size)
    create_Abaqus_output_nodes(endpoint_output_nodes, 'ENDCAP_OUTPUT')
    create_Abaqus_output_nodes(bladder_output_nodes_in, 'INNER_BLADDER_OUTPUT')
    create_Abaqus_output_nodes(fiber_a_points, 'FIBER_GROUP_A_OUTPUT')
    create_Abaqus_output_nodes(fiber_b_points, 'FIBER_GROUP_B_OUTPUT')

    # set up Abaqus assembly
    define_Abaqus_materials()
    define_Abaqus_sections(FREE=myFREE)
    assign_Abaqus_sections(FREE=myFREE)
    create_Abaqus_assembly()
    position_endcaps(path_end_point=myFREE.path[-1])
    setup_Abaqus_job(FREE=myFREE)

# run main()
if __name__ == "__main__":

    main()
