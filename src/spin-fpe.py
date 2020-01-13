# Solves the FPE formulation of monodomain LLG on the surface of a unit sphere
#
# With only uniaxial unisotropy and initial uniform distribution of PDF over the
# surface of the unit sphere, the result of the simulation should give a
# distribution with two peaks, one at each poles of the sphere defined by the
# axis of the anisotropy
#
from fipy import FaceVariable, CellVariable, Gmsh2DIn3DSpace, Viewer, TransientTerm, DiffusionTerm, DefaultSolver
from fipy.variables.variable import Variable
from fipy.tools import numerix

# define mesh
mesh = Gmsh2DIn3DSpace('''
    radius = 1.0;
    cellSize = 0.01;

    // create inner 1/8 shell
    Point(1) = {0, 0, 0, cellSize};
    Point(2) = {-radius, 0, 0, cellSize};
    Point(3) = {0, radius, 0, cellSize};
    Point(4) = {0, 0, radius, cellSize};
    Circle(1) = {2, 1, 3};
    Circle(2) = {4, 1, 2};
    Circle(3) = {4, 1, 3};
    Line Loop(1) = {1, -3, 2} ;
    Ruled Surface(1) = {1};

    // create remaining 7/8 inner shells
    t1[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{1};}};
    t2[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{1};}};
    t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{1};}};
    t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2} {Duplicata{Surface{1};}};
    t5[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{t4[0]};}};
    t6[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{t4[0]};}};
    t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{t4[0]};}};

    // create entire inner and outer shell
    Surface Loop(100)={1,t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
''', order=2).extrude(extrudeFunc=lambda r: 1.05 * r) # doctest: +GMSH
#

mmag = FaceVariable(name=r"$mmag$", mesh=mesh) # doctest: +GMSH
gridCoor = mesh.faceCenters

## Constants
kBoltzmann = 1.38064852e-23
mu0 = numerix.pi * 4.0e-7

## LLG parameters
##gamFac = 1.7608e11 * pi * 4.0e-7
gamFac = 2.2128e5
alphaDamping = 0.01
Temperature = 300
Msat = 1050e3
magVolume = 2.0e-9 * (25e-9 * 25e-9) * numerix.pi
D = alphaDamping * gamFac * kBoltzmann * Temperature / ((1 + alphaDamping) * Msat * magVolume)

# Uniaxial anisotropy
Ku2 = 800e3
uAxis = numerix.array((0., 0., 1.))

# Define arrays storing the torque terms in LLG
TeffBase = numerix.zeros((3,len(gridCoor[0])), 'd')
TuniaxBase = numerix.zeros((3,len(gridCoor[0])), 'd')

# Define array of HeffBase that does into LLG
# HeffBase contains Heff terms that only vary with m (independent of t)
HeffBase = numerix.zeros((3,len(gridCoor[0])), 'd')
HuniaxBase = numerix.zeros((3,len(gridCoor[0])), 'd')

# Normalize uniaxial anisotropy axis vector
uAxisNorm2 = (uAxis[0]*uAxis[0]) + (uAxis[1]*uAxis[1]) + (uAxis[2]*uAxis[2])
uAxisNorm = numerix.sqrt(uAxisNorm2)
uAxisUnit = uAxis / uAxisNorm
HuniScaleFac = -2.0 * Ku2 / Msat

print('Looping over')
print(len(gridCoor[0]))
for idx in range(len(gridCoor[0]) - 1):
    m0x = gridCoor[0][idx]
    m0y = gridCoor[1][idx]
    m0z = gridCoor[2][idx]
    mag2 = (m0x*m0x) + (m0y*m0y) + (m0z*m0z)
    mag = numerix.sqrt(mag2)
    mx = m0x / mag
    my = m0y / mag
    mz = m0z / mag
    mdotu = mx*uAxisUnit[0] + my*uAxisUnit[1] + mz*uAxisUnit[2]
    HuniaxBase[0][idx] = HuniScaleFac * mdotu * uAxisUnit[0]
    HuniaxBase[1][idx] = HuniScaleFac * mdotu * uAxisUnit[1]
    HuniaxBase[2][idx] = HuniScaleFac * mdotu * uAxisUnit[2]
    mxHx = my*HuniaxBase[2][idx] - mz*HuniaxBase[1][idx]
    mxHy = mz*HuniaxBase[0][idx] - mx*HuniaxBase[2][idx]
    mxHz = mx*HuniaxBase[1][idx] - my*HuniaxBase[0][idx]

    mxmxHx = my*mxHz - mz*mxHy
    mxmxHy = mz*mxHx - mx*mxHz
    mxmxHz = mx*mxHy - my*mxHx
    TuniaxBase[0][idx] = -gamFac * (mxHx + alphaDamping * mxmxHx)
    TuniaxBase[1][idx] = -gamFac * (mxHy + alphaDamping * mxmxHy)
    TuniaxBase[2][idx] = -gamFac * (mxHz + alphaDamping * mxmxHz)
    print('Progress')
    print(100.0 * idx/len(gridCoor[0]))
