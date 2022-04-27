# This file contains all functions for setting up an Abaqus job for curved FREEs

# import necessary packages
import numpy as np

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


def create_Abaqus_model():
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-1')

# takes in a series of 3D fiber points, then builds and meshes Abaqus parts from them
# fiber_coords[i][j] = point on fiber j at s = s[i] (direction A)
def create_Abaqus_fibers(FREE, fiber_coords_1, fiber_coords_2):
    num_fibers = np.shape(fiber_coords_1)[1]

    # initialize parts
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='FIBER_GROUP_A', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='FIBER_GROUP_B', type=DEFORMABLE_BODY)

    # initialize fiber points
    fibers_start_points = list()

    # transfer fiber coordinates to fiber_points arrays
    for j in range(num_fibers):

        fiber_a_points = list()
        fiber_b_points = list()

        # define fiber nodes
        for i in range(np.shape(fiber_coords_1)[0]):
            x1 = fiber_coords_1[i][j][0]
            y1 = fiber_coords_1[i][j][1]
            z1 = fiber_coords_1[i][j][2]

            x2 = fiber_coords_2[i][j][0]
            y2 = fiber_coords_2[i][j][1]
            z2 = fiber_coords_2[i][j][2]

            fiber_a_points.append((x1, y1, z1))
            fiber_b_points.append((x2, y2, z2))

            if i == 0:
                fibers_start_points.append((x1, y1, z1))

        # add fiber to FIBERS1 part
        mdb.models['Model-1'].parts['FIBER_GROUP_A'].WireSpline(meshable=ON, points=fiber_a_points,
                                                                smoothClosedSpline=ON)

        # add fiber to fibers2 part
        mdb.models['Model-1'].parts['FIBER_GROUP_B'].WireSpline(meshable=ON, points=fiber_b_points,
                                                                smoothClosedSpline=ON)

    # define fiber surfaces and sets
    mdb.models['Model-1'].parts['FIBER_GROUP_A'].Set(
        edges=mdb.models['Model-1'].parts['FIBER_GROUP_A'].edges.findAt(coordinates=fibers_start_points),
        name='all')
    mdb.models['Model-1'].parts['FIBER_GROUP_A'].Surface(
        circumEdges=mdb.models['Model-1'].parts['FIBER_GROUP_A'].edges.getByBoundingSphere(
            center=(0.0, 0.0, 0.0),
            radius=1e6),
        name='circumferential')
    mdb.models['Model-1'].parts['FIBER_GROUP_B'].Set(
        edges=mdb.models['Model-1'].parts['FIBER_GROUP_B'].edges.findAt(coordinates=fibers_start_points),
        name='all')
    mdb.models['Model-1'].parts['FIBER_GROUP_B'].Surface(
        circumEdges=mdb.models['Model-1'].parts['FIBER_GROUP_B'].edges.getByBoundingSphere(
            center=(0.0, 0.0, 0.0),radius=1e6),
        name='circumferential')

    # seed, assign mesh controls, and mesh parts
    mdb.models['Model-1'].parts['FIBER_GROUP_A'].seedPart(
        deviationFactor=0.1, minSizeFactor=0.1, size=FREE.fiber_seed_size)
    mdb.models['Model-1'].parts['FIBER_GROUP_B'].seedPart(
        deviationFactor=0.1, minSizeFactor=0.1, size=FREE.fiber_seed_size)

    if FREE.fiber_element == 'beam':

        mdb.models['Model-1'].parts['FIBER_GROUP_A'].setElementType(
            elemTypes=(ElemType(elemCode=B32H, elemLibrary=STANDARD),),
            regions=(mdb.models['Model-1'].parts['FIBER_GROUP_A'].edges.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))

        mdb.models['Model-1'].parts['FIBER_GROUP_B'].setElementType(
            elemTypes=(ElemType(elemCode=B32H, elemLibrary=STANDARD),),
            regions=(mdb.models['Model-1'].parts['FIBER_GROUP_B'].edges.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))

    elif FREE.fiber_element == 'truss':

        mdb.models['Model-1'].parts['FIBER_GROUP_A'].setElementType(
            elemTypes=(ElemType(elemCode=T3D2, elemLibrary=STANDARD),),
            regions=(mdb.models['Model-1'].parts['FIBER_GROUP_A'].edges.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))

        mdb.models['Model-1'].parts['FIBER_GROUP_B'].setElementType(
            elemTypes=(ElemType(elemCode=T3D2, elemLibrary=STANDARD),),
            regions=(mdb.models['Model-1'].parts['FIBER_GROUP_B'].edges.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))

    else:

        raise Exception("Fiber element type not recognized. Check user input.")

    mdb.models['Model-1'].parts['FIBER_GROUP_A'].generateMesh()
    mdb.models['Model-1'].parts['FIBER_GROUP_B'].generateMesh()

# takes in FREE path, then builds and meshes Abaqus bladder part
def create_Abaqus_bladder(path, outer_radius, thickness, bladder_seed_size):
    FREE_inner_radius = outer_radius - thickness

    # initialize part
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='BLADDER_DEFINITION', type=DEFORMABLE_BODY)

    # define bladder path
    path_points = list()
    bladder_fixed_point = list()
    bladder_outside_point = list()
    bladder_inside_point = list()

    for i in range(np.shape(path)[0]):
        x = path[i][0]
        y = path[i][1]
        z = path[i][2]

        path_points.append((x, y, z))

        mid_thickness_radius = (FREE_inner_radius + outer_radius) / 2.
        if i == 0:
            bladder_fixed_point.append(
                (x, y + mid_thickness_radius, z))  # coordinates for selecting start endcap surface of bladder
        elif i == 1:
            bladder_outside_point.append(
                (x, y + outer_radius, z))  # coordinates for selecting outside surface of bladder
            bladder_inside_point.append(
                (x, y + FREE_inner_radius, z))  # coordinates for selecting inside surface of bladder

    # create a wire to use as the sweep path
    mdb.models['Model-1'].parts['BLADDER_DEFINITION'].WireSpline(
        meshable=ON,
        points=path_points,
        smoothClosedSpline=ON)

    # sketch bladder section
    mdb.models['Model-1'].parts['BLADDER_DEFINITION'].DatumPlaneByPrincipalPlane(
        offset=0.0,
        principalPlane=XYPLANE)
    mdb.models['Model-1'].parts['BLADDER_DEFINITION'].DatumAxisByPrincipalAxis(
        principalAxis=YAXIS)

    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__',
        sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0),
        point1=(outer_radius, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0),
        point1=(FREE_inner_radius, 0.0))

    mdb.models['Model-1'].parts['BLADDER_DEFINITION'].SolidSweep(
        path=mdb.models['Model-1'].parts['BLADDER_DEFINITION'].edges.getSequenceFromMask(('[#1 ]',), ),
        profile=mdb.models['Model-1'].sketches['__profile__'],
        sketchOrientation=RIGHT,
        sketchPlane=mdb.models['Model-1'].parts['BLADDER_DEFINITION'].datums[2],
        sketchUpEdge=mdb.models['Model-1'].parts['BLADDER_DEFINITION'].datums[3])

    del mdb.models['Model-1'].sketches['__profile__']

    # separate wire feature from solid bladder feature (only 'BLADDER' is used for FEA)
    mdb.models['Model-1'].Part(
        compressFeatureList=ON,
        name='BLADDER-Copy',
        objectToCopy=mdb.models['Model-1'].parts['BLADDER_DEFINITION'],
        separate=ON)
    mdb.models['Model-1'].parts.changeKey(
        fromName='BLADDER-Copy-1',
        toName='PATH')
    mdb.models['Model-1'].parts.changeKey(
        fromName='BLADDER-Copy-2',
        toName='BLADDER')

    # define bladder surfaces and sets
    mdb.models['Model-1'].parts['BLADDER'].Surface(
        name='inner',
        side1Faces=mdb.models['Model-1'].parts['BLADDER'].faces.findAt(coordinates=bladder_inside_point))
    mdb.models['Model-1'].parts['BLADDER'].Surface(
        name='outer',
        side1Faces=mdb.models['Model-1'].parts['BLADDER'].faces.findAt(coordinates=bladder_outside_point))
    mdb.models['Model-1'].parts['BLADDER'].Surface(
        name='fixed end',
        side1Faces=mdb.models['Model-1'].parts['BLADDER'].faces.findAt(coordinates=bladder_fixed_point))
    mdb.models['Model-1'].parts['BLADDER'].Surface(
        name='free end',
        side1Faces=mdb.models['Model-1'].parts['BLADDER'].faces.getByBoundingSphere(center=(x, y, z),
                                                                                    radius=outer_radius * 1.1))
    mdb.models['Model-1'].parts['BLADDER'].Set(
        name='all',
        cells=mdb.models['Model-1'].parts['BLADDER'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6))

    # mesh bladder
    mdb.models['Model-1'].parts['BLADDER'].seedPart(deviationFactor=0.1, minSizeFactor=0.1,
                                                    size=bladder_seed_size)
    mdb.models['Model-1'].parts['BLADDER'].setMeshControls(
        algorithm=MEDIAL_AXIS,
        regions=mdb.models['Model-1'].parts['BLADDER'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6))
    mdb.models['Model-1'].parts['BLADDER'].setElementType(
        elemTypes=(
            ElemType(
                elemCode=C3D8RH,
                elemLibrary=STANDARD,
                kinematicSplit=AVERAGE_STRAIN,
                hourglassControl=DEFAULT),),
        regions=(mdb.models['Model-1'].parts['BLADDER'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))
    mdb.models['Model-1'].parts['BLADDER'].generateMesh()


# takes in FREE path, then builds and meshes Abaqus bladder part
def create_Abaqus_endcap(outer_radius, FREE_thickness, endcap_thickness, endcap_seed_size):
    FREE_inner_radius = outer_radius - FREE_thickness

    # sketch endcap path
    mdb.models['Model-1'].ConstrainedSketch(
        name='__sweep__',
        sheetSize=200.0)
    mdb.models['Model-1'].sketches['__sweep__'].Line(point1=(0.0, 0.0), point2=(0.0, -endcap_thickness))

    # sketch endcap section
    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__',
        sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0),
                                                                          point1=(outer_radius, 0.0))

    # initialize endcap part
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='ENDCAP', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['ENDCAP'].BaseSolidSweep(
        path=mdb.models['Model-1'].sketches['__sweep__'],
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    del mdb.models['Model-1'].sketches['__sweep__']

    # partition endcap
    mdb.models['Model-1'].parts['ENDCAP'].DatumAxisByPrincipalAxis(
        principalAxis=ZAXIS)
    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__',
        sheetSize=200.0,
        transform=mdb.models['Model-1'].parts['ENDCAP'].MakeSketchTransform(
            sketchPlane=mdb.models['Model-1'].parts['ENDCAP'].faces[2],
            sketchPlaneSide=SIDE1,
            sketchUpEdge=mdb.models['Model-1'].parts['ENDCAP'].datums[2],
            sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models['Model-1'].parts['ENDCAP'].projectReferencesOntoSketch(
        filter=COPLANAR_EDGES,
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0),
        point1=(FREE_inner_radius, 0.0))
    mdb.models['Model-1'].parts['ENDCAP'].PartitionFaceBySketch(
        faces=mdb.models['Model-1'].parts['ENDCAP'].faces.findAt(coordinates=((0.0, 0.0, 0.0),)),
        sketch=mdb.models['Model-1'].sketches['__profile__'],
        sketchUpEdge=mdb.models['Model-1'].parts['ENDCAP'].datums[2])
    del mdb.models['Model-1'].sketches['__profile__']

    # define endcap surfaces and sets
    mdb.models['Model-1'].parts['ENDCAP'].Surface(
        name='inner',
        side1Faces=mdb.models['Model-1'].parts['ENDCAP'].faces.findAt(coordinates=((0.0, 0.0, 0.0),)))
    mdb.models['Model-1'].parts['ENDCAP'].Surface(
        name='ring',
        side1Faces=mdb.models['Model-1'].parts['ENDCAP'].faces.findAt(
            coordinates=(((outer_radius + FREE_inner_radius) / 2, 0.0, 0.0),)))
    mdb.models['Model-1'].parts['ENDCAP'].Set(
        cells=mdb.models['Model-1'].parts['ENDCAP'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),
        name='all')

    # mesh endcap
    mdb.models['Model-1'].parts['ENDCAP'].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=endcap_seed_size)
    mdb.models['Model-1'].parts['ENDCAP'].setMeshControls(
        algorithm=MEDIAL_AXIS,
        regions=mdb.models['Model-1'].parts['ENDCAP'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6), )
    mdb.models['Model-1'].parts['ENDCAP'].setElementType(
        elemTypes=(
            ElemType(
                elemCode=C3D8RH,
                elemLibrary=STANDARD,
                kinematicSplit=AVERAGE_STRAIN,
                hourglassControl=DEFAULT),),
        regions=(mdb.models['Model-1'].parts['ENDCAP'].cells.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6),))
    mdb.models['Model-1'].parts['ENDCAP'].generateMesh()


# converts output node coordinates into an Abaqus part with indexed nodes
def create_Abaqus_output_nodes(nodes, name):
    # initialize output part
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name=name, type=DEFORMABLE_BODY)

    # add specified output nodes to part
    # endpoint_output_nodes[n] = point n on free endcap
    # inner_bladder_nodes[i][j] = point on inner bladder at s = s[i] and theta = j * d_theta
    # fiber_a_points[i][j] = point on fiber f at s = s[i] (direction A)
    # fiber_b_points[i][j] = point on fiber f at s = s[i] (direction B)
    if len(np.shape(nodes)) == 1:  # single node
        x = nodes[0]
        y = nodes[1]
        z = nodes[2]

        mdb.models['Model-1'].parts[name].Node(coordinates=(x, y, z))

    elif len(np.shape(nodes)) == 2:  # multiple nodes, no grouping
        for i in range(np.shape(nodes)[0]):
            x = nodes[i][0]
            y = nodes[i][1]
            z = nodes[i][2]

            mdb.models['Model-1'].parts[name].Node(coordinates=(x, y, z))

    else:  # multiple nodes
        for j in range(np.shape(nodes)[1]):  # theta index for bladder, fiber number for fibers
            set_labels = []  # list of nodes for node sets in the axial direction

            for i in range(np.shape(nodes)[0]):  # index along path
                x = nodes[i][j][0]
                y = nodes[i][j][1]
                z = nodes[i][j][2]

                node_label = int(j * np.shape(nodes)[0] + i + 1)
                set_labels.append(node_label)

                mdb.models['Model-1'].parts[name].Node(coordinates=(x, y, z))

            set_name = 'Set-' + str(j + 1)

            mdb.models['Model-1'].parts[name].Set(
                name=set_name,
                nodes=mdb.models['Model-1'].parts[name].nodes.sequenceFromLabels(set_labels))

    mdb.models['Model-1'].parts[name].Set(
        name='all',
        nodes=mdb.models['Model-1'].parts[name].nodes.getByBoundingSphere(center=(0.0, 0.0, 0.0), radius=1e6))


# defines Abaqus materials
def define_Abaqus_materials():
    mdb.models['Model-1'].Material(name='steel (elastic)')
    mdb.models['Model-1'].materials['steel (elastic)'].Elastic(table=((200000.0, 0.3),))

    mdb.models['Model-1'].Material(name='Kevlar (elastic, Connolly 2015)')
    mdb.models['Model-1'].materials['Kevlar (elastic, Connolly 2015)'].Elastic(table=((31067, 0.36),))

    mdb.models['Model-1'].Material(name='Ecoflex 00-30 (M-R, Steck 2019)')
    mdb.models['Model-1'].materials['Ecoflex 00-30 (M-R, Steck 2019)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((0.048, -0.152, 0.0),),
        testData=OFF,
        type=MOONEY_RIVLIN,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='Ecoflex 00-30 (Yeoh, Steck 2019)')
    mdb.models['Model-1'].materials['Ecoflex 00-30 (Yeoh, Steck 2019)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((0.017, -0.0002, 2.3e-05, 0.0, 0.0, 0.0),),
        testData=OFF,
        type=YEOH,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='Ecoflex 00-50 (Yeoh, Kulkarni 2015)')
    mdb.models['Model-1'].materials['Ecoflex 00-50 (Yeoh, Kulkarni 2015)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((1.9e-02, 9.0e-04, -4.75e-06, 0.0, 0.0, 0.0),),
        testData=OFF,
        type=YEOH,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='Ecoflex 00-50 (Yeoh, Xue 2015)')
    mdb.models['Model-1'].materials['Ecoflex 00-50 (Yeoh, Xue 2015)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((0.1, 0.02, 0., 0., 0., 0.),),
        testData=OFF,
        type=YEOH,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='DS10 (Yeoh, Caasenbrood 2020)')
    mdb.models['Model-1'].materials['DS10 (Yeoh, Caasenbrood 2020)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((36.e-3, 0.25e-3, 0.023e-3, 0., 0., 0.),),
        testData=OFF,
        type=YEOH,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='DS10 (NeoHk, Polygerinos 2019)')
    mdb.models['Model-1'].materials['DS10 (NeoHk, Polygerinos 2019)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((0.0425, 0.0),),
        testData=OFF,
        type=NEO_HOOKE,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='DS30 (M-R, Elsayed 2014)')
    mdb.models['Model-1'].materials['DS30 (M-R, Elsayed 2014)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((1.19e-3, 23.028e-3, 0.0),),
        testData=OFF,
        type=MOONEY_RIVLIN,
        volumetricResponse=VOLUMETRIC_DATA)

    mdb.models['Model-1'].Material(name='DS30 (Yeoh, Chen 2018)')
    mdb.models['Model-1'].materials['DS30 (Yeoh, Chen 2018)'].Hyperelastic(
        materialType=ISOTROPIC,
        table=((96e-3, 9.5e-3, 0., 0., 0., 0.),),
        testData=OFF,
        type=YEOH,
        volumetricResponse=VOLUMETRIC_DATA)


# define Abaqus sections
def define_Abaqus_sections(FREE):
    mdb.models['Model-1'].HomogeneousSolidSection(
        material=FREE.material,
        name='bladder',
        thickness=None)

    mdb.models['Model-1'].HomogeneousSolidSection(
        material='steel (elastic)',
        name='steel solid',
        thickness=None)

    if FREE.fiber_element == 'beam':

        mdb.models['Model-1'].CircularProfile(
            name='Kevlar profile',
            r=0.0889)  # fiber radius (mm) from Connolly 2015

        mdb.models['Model-1'].BeamSection(
            consistentMassMatrix=False,
            integration=DURING_ANALYSIS,
            material='Kevlar (elastic, Connolly 2015)',
            name='Kevlar beam',
            poissonRatio=0.36,
            profile='Kevlar profile',
            temperatureVar=LINEAR)

    elif FREE.fiber_element == 'truss':

        mdb.models['Model-1'].TrussSection(
            area=1.0,
            material='steel (elastic)',
            name='steel truss')

    else:

        raise Exception("Fiber element type not recognized. Check user input.")


# assign Abaqus sections
def assign_Abaqus_sections(FREE):
    mdb.models['Model-1'].parts['BLADDER'].SectionAssignment(
        offset=0.0,
        offsetField='',
        offsetType=MIDDLE_SURFACE,
        region=mdb.models['Model-1'].parts['BLADDER'].sets['all'],
        sectionName='bladder',
        thicknessAssignment=FROM_SECTION)

    mdb.models['Model-1'].parts['ENDCAP'].SectionAssignment(
        offset=0.0,
        offsetField='',
        offsetType=MIDDLE_SURFACE,
        region=mdb.models['Model-1'].parts['ENDCAP'].sets['all'],
        sectionName='steel solid', thicknessAssignment=FROM_SECTION)

    if FREE.fiber_element == 'beam':

        mdb.models['Model-1'].parts['FIBER_GROUP_A'].SectionAssignment(
            offset=0.0,
            offsetField='',
            offsetType=MIDDLE_SURFACE,
            region=mdb.models['Model-1'].parts['FIBER_GROUP_A'].sets['all'],
            sectionName='Kevlar beam', thicknessAssignment=FROM_SECTION)

        mdb.models['Model-1'].parts['FIBER_GROUP_A'].assignBeamSectionOrientation(
            method=N1_COSINES,
            n1=(0.0, 0.0, -1.0),
            region=mdb.models['Model-1'].parts['FIBER_GROUP_A'].sets['all'])

        mdb.models['Model-1'].parts['FIBER_GROUP_B'].SectionAssignment(
            offset=0.0,
            offsetField='',
            offsetType=MIDDLE_SURFACE,
            region=mdb.models['Model-1'].parts['FIBER_GROUP_B'].sets['all'],
            sectionName='Kevlar beam', thicknessAssignment=FROM_SECTION)

        mdb.models['Model-1'].parts['FIBER_GROUP_B'].assignBeamSectionOrientation(
            method=N1_COSINES,
            n1=(0.0, 0.0, -1.0),
            region=mdb.models['Model-1'].parts['FIBER_GROUP_B'].sets['all'])

    elif FREE.fiber_element == 'truss':

        mdb.models['Model-1'].parts['FIBER_GROUP_A'].SectionAssignment(
            offset=0.0,
            offsetField='',
            offsetType=MIDDLE_SURFACE,
            region=mdb.models['Model-1'].parts['FIBER_GROUP_A'].sets['all'],
            sectionName='steel truss', thicknessAssignment=FROM_SECTION)

        mdb.models['Model-1'].parts['FIBER_GROUP_B'].SectionAssignment(
            offset=0.0,
            offsetField='',
            offsetType=MIDDLE_SURFACE,
            region=mdb.models['Model-1'].parts['FIBER_GROUP_B'].sets['all'],
            sectionName='steel truss', thicknessAssignment=FROM_SECTION)

    else:

        raise Exception("Fiber element type not recognized. Check user input.")


# assemble Abaqus model
def create_Abaqus_assembly():
    # add part instances to assembly
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='BLADDER-1',
        part=mdb.models['Model-1'].parts['BLADDER'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='ENDCAP-1',
        part=mdb.models['Model-1'].parts['ENDCAP'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='ENDCAP-2',
        part=mdb.models['Model-1'].parts['ENDCAP'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='FIBER_GROUP_A-1',
        part=mdb.models['Model-1'].parts['FIBER_GROUP_A'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='FIBER_GROUP_B-1',
        part=mdb.models['Model-1'].parts['FIBER_GROUP_B'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='ENDCAP_OUTPUT-1',
        part=mdb.models['Model-1'].parts['ENDCAP_OUTPUT'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='INNER_BLADDER_OUTPUT-1',
        part=mdb.models['Model-1'].parts['INNER_BLADDER_OUTPUT'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='FIBER_GROUP_A_OUTPUT-1',
        part=mdb.models['Model-1'].parts['FIBER_GROUP_A_OUTPUT'])
    mdb.models['Model-1'].rootAssembly.Instance(
        dependent=ON,
        name='FIBER_GROUP_B_OUTPUT-1',
        part=mdb.models['Model-1'].parts['FIBER_GROUP_B_OUTPUT'])


# move endcap 2 to end of FREE in Abaqus
def position_endcaps(path_end_point):
    # rotate fixed endcap to align with fixed end of FREE
    mdb.models['Model-1'].rootAssembly.rotate(
        instanceList=('ENDCAP-1',),
        angle=90.0,
        axisDirection=(1.0, 0.0, 0.0),
        axisPoint=(0.0, 0.0, 0.0))

    # rotate free endcap to align with free end of FREE
    mdb.models['Model-1'].rootAssembly.ParallelFace(
        fixedPlane=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].faces[2],
        flip=ON,
        movablePlane=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-2'].faces[0])

    # translate free endcap to free end of FREE
    x = path_end_point[0]
    y = path_end_point[1]
    z = path_end_point[2]
    mdb.models['Model-1'].rootAssembly.translate(
        instanceList=('ENDCAP-2',),
        vector=(x, y, z))


# set up Abaqus constraints, loading, and initialize job
def setup_Abaqus_job(FREE):
    # initialize FEM step
    mdb.models['Model-1'].StaticStep(
        initialInc=0.01,
        maxInc=0.01,
        maxNumInc=1000,
        name='Step-1',
        nlgeom=ON,
        previous='Initial')

    # apply constraints
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-1'].surfaces['ring'],
        name='ENDCAP-1 to BLADDER',
        positionToleranceMethod=COMPUTED,
        slave=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['fixed end'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-2'].surfaces['ring'],
        name='ENDCAP-2 to BLADDER',
        positionToleranceMethod=COMPUTED,
        slave=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['free end'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['outer'],
        name='FIBER_A to BLADDER',
        positionToleranceMethod=COMPUTED,
        slave=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_A-1'].surfaces['circumferential'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['outer'],
        name='FIBER_B to BLADDER',
        positionToleranceMethod=COMPUTED,
        slave=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_B-1'].surfaces['circumferential'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['inner'],
        name='INNER_BLADDER_OUTPUT to BLADDER',
        positionTolerance=1.0,
        positionToleranceMethod=SPECIFIED,
        slave=mdb.models['Model-1'].rootAssembly.instances['INNER_BLADDER_OUTPUT-1'].sets['all'],
        thickness=ON,
        tieRotations=ON,
        constraintEnforcement=NODE_TO_SURFACE)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_A-1'].surfaces['circumferential'],
        name='FIBER_A_OUTPUT to FIBER_A',
        positionTolerance=1.0,
        positionToleranceMethod=SPECIFIED,
        slave=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_A_OUTPUT-1'].sets['all'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_B-1'].surfaces['circumferential'],
        name='FIBER_B_OUTPUT to FIBER_B',
        positionTolerance=1.0,
        positionToleranceMethod=SPECIFIED,
        slave=mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_B_OUTPUT-1'].sets['all'],
        thickness=ON,
        tieRotations=ON)
    mdb.models['Model-1'].Tie(
        adjust=ON,
        master=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-2'].surfaces['inner'],
        name='ENDCAP_OUTPUT to ENDCAP-2',
        positionTolerance=1.0,
        positionToleranceMethod=SPECIFIED,
        slave=mdb.models['Model-1'].rootAssembly.instances['ENDCAP_OUTPUT-1'].sets['all'],
        thickness=ON,
        tieRotations=ON,
        constraintEnforcement=NODE_TO_SURFACE)

    # apply boundary conditions
    mdb.models['Model-1'].EncastreBC(
        createStepName='Step-1',
        localCsys=None,
        name='fixed ENDCAP-1',
        region=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-1'].sets['all'])

    # apply loads
    mdb.models['Model-1'].Pressure(
        amplitude=UNSET, createStepName='Step-1',
        distributionType=UNIFORM,
        field='',
        magnitude=FREE.final_pressure,
        name='BLADDER pressure',
        region=mdb.models['Model-1'].rootAssembly.instances['BLADDER-1'].surfaces['inner'])
    mdb.models['Model-1'].Pressure(
        amplitude=UNSET,
        createStepName='Step-1',
        distributionType=UNIFORM,
        field='',
        magnitude=FREE.final_pressure,
        name='ENDCAP-1 pressure',
        region=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-1'].surfaces['inner'])
    mdb.models['Model-1'].Pressure(
        amplitude=UNSET,
        createStepName='Step-1',
        distributionType=UNIFORM,
        field='',
        magnitude=FREE.final_pressure,
        name='ENDCAP-2 pressure',
        region=mdb.models['Model-1'].rootAssembly.instances['ENDCAP-2'].surfaces['inner'])

    # create job
    job = mdb.Job(
        atTime=None,
        contactPrint=OFF,
        description='',
        echoPrint=OFF,
        explicitPrecision=SINGLE,
        getMemoryFromAnalysis=True,
        historyPrint=OFF,
        memory=90,
        memoryUnits=PERCENTAGE,
        model='Model-1',
        modelPrint=OFF,
        multiprocessingMode=DEFAULT,
        name=FREE.label,
        nodalOutputPrecision=SINGLE,
        numCpus=4,
        numDomains=4,
        numGPUs=1,
        queue=None,
        resultsFormat=ODB,
        scratch='',
        type=ANALYSIS,
        userSubroutine='',
        waitHours=0,
        waitMinutes=0)
