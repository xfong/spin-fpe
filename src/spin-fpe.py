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
from joblib import Parallel, delayed
import multiprocessing

def NormDotProduct(mx, my, mz, uax, uay, uaz, scaleFac):
    mag2 = (mx*mx) + (my*my) + (mz*mz)
    mag = numerix.sqrt(mag2)
    unorm2 = (uax*uax) + (uay*uay) + (uaz*uaz)
    unorm = numerix.sqrt(unorm2)
    dotProd = scaleFac * (((mx*uax) + (my*uay) + (mz*uaz)) / (mag * unorm))
    return dotProd

def ScaleVectorComp(ax, scaleFac):
    return ax*scaleFac

def UniaxialXComp(mx, my, mz, uax, uay, uaz, scaleFac):
    mag2 = (mx*mx) + (my*my) + (mz*mz)
    mag = numerix.sqrt(mag2)
    unorm2 = (uax*uax) + (uay*uay) + (uaz*uaz)
    unorm = numerix.sqrt(unorm2)
    dotProd = scaleFac * (((mx*uax) + (my*uay) + (mz*uaz)) / (mag * unorm))
    return dotProd * uax

def UniaxialYComp(mx, my, mz, uax, uay, uaz, scaleFac):
    mag2 = (mx*mx) + (my*my) + (mz*mz)
    mag = numerix.sqrt(mag2)
    unorm2 = (uax*uax) + (uay*uay) + (uaz*uaz)
    unorm = numerix.sqrt(unorm2)
    dotProd = scaleFac * (((mx*uax) + (my*uay) + (mz*uaz)) / (mag * unorm))
    return dotProd * uay

def UniaxialZComp(mx, my, mz, uax, uay, uaz, scaleFac):
    mag2 = (mx*mx) + (my*my) + (mz*mz)
    mag = numerix.sqrt(mag2)
    unorm2 = (uax*uax) + (uay*uay) + (uaz*uaz)
    unorm = numerix.sqrt(unorm2)
    dotProd = scaleFac * (((mx*uax) + (my*uay) + (mz*uaz)) / (mag * unorm))
    return dotProd * uaz

def VectorCrossProduct(ax, ay, az, bx, by, bz, scaleFac):
    outputVector = numerix.zeros(3,'d')
    outputVector[0] = scaleFac*(ay*bz - az*by)
    outputVector[1] = scaleFac*(az*bx - ax*bz)
    outputVector[2] = scaleFac*(ax*by - ay*bx)
    return outputVector

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

num_cores = multiprocessing.cpu_count()

HuniaxBase[0] = Parallel(n_jobs=num_cores)(delayed(UniaxialXComp)(mesh.faceCenters[0][idx], mesh.faceCenters[1][idx], mesh.faceCenters[2][idx], uAxisUnit[0], uAxisUnit[1], uAxisUnit[2], HuniScaleFac) for idx in range(len(mesh.faceCenters[0])-1))
#HuniaxBase[1] = HuniScaleFac * mdotu * uAxisUnit[1]
#HuniaxBase[2] = HuniScaleFac * mdotu * uAxisUnit[2]
