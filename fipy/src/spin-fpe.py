# Solves the FPE formulation of monodomain LLG on the surface of a unit sphere
#
# With only uniaxial unisotropy and initial uniform distribution of PDF over the
# surface of the unit sphere, the result of the simulation should give a
# distribution with two peaks, one at each poles of the sphere defined by the
# axis of the anisotropy
#
from fipy import FaceVariable, CellVariable, Gmsh2DIn3DSpace, Viewer, TransientTerm, ExplicitDiffusionTerm, DiffusionTerm, ExponentialConvectionTerm, DefaultSolver
from fipy.variables.variable import Variable
from fipy.tools import numerix

def HeffUniaxialAnisotropyKM(mUnit, uAxis, Ku2, Msat):
    uAxisNorm = numerix.linalg.norm(uAxis)
    uAxisUnit = uAxis / uAxisNorm
    mNorm = numerix.linalg.norm(mUnit)
    mArray = mUnit / mNorm
    uAxisArr = numerix.tile(uAxisUnit, (len(mUnit[0]), 1))
    uAxisArr = numerix.transpose(uAxisArr)
    mdotu = numerix.dot(mArray, uAxisArr)
    scaleFac = numerix.multiply(mdotu, (-2.0 * Ku2 / Msat))
    Heff = numerix.zeros((3, len(scaleFac)), 'd')
    Heff[0] = numerix.multiply(scaleFac, uAxisArr[0])
    Heff[1] = numerix.multiply(scaleFac, uAxisArr[1])
    Heff[2] = numerix.multiply(scaleFac, uAxisArr[2])
    return Heff

def HeffUniaxialAnisotropyFac(mUnit, uAxis, Hfac):
    uAxisNorm = numerix.linalg.norm(uAxis)
    uAxisUnit = uAxis / uAxisNorm
    mNorm = numerix.linalg.norm(mUnit)
    mArray = mUnit / mNorm
    uAxisArr = numerix.tile(uAxisUnit, (len(mUnit[0]), 1))
    uAxisArr = numerix.transpose(uAxisArr)
    mdotu = numerix.dot(mArray, uAxisArr)
    scaleFac = numerix.multiply(mdotu, Hfac)
    Heff = numerix.zeros((3, len(scaleFac)), 'd')
    Heff[0] = numerix.multiply(scaleFac, uAxisArr[0])
    Heff[1] = numerix.multiply(scaleFac, uAxisArr[1])
    Heff[2] = numerix.multiply(scaleFac, uAxisArr[2])
    return Heff

def calculateFieldLikeTorque(mUnit, uAxis):
    uAxisNorm = numerix.linalg.norm(uAxis)
    uAxisUnit = uAxis / uAxisNorm
    mNorm = numerix.linalg.norm(mUnit)
    mArray = mUnit / mNorm
    uAxisArr = numerix.tile(uAxisUnit, (len(mUnit[0]), 1))
    uAxisArr = numerix.transpose(uAxisArr)
    m_x_u = numerix.cross(mArray, uAxisArr)
    return numerix.transpose(m_x_u)

def calculateDampingLikeTorque(mUnit, uAxis):
    m_x_u = calculateFieldLikeTorque(mUnit, uAxis)
    mNorm = numerix.linalg.norm(mUnit)
    mArray = mUnit / mNorm
    m_x_m_x_u = numerix.cross(mArray, m_x_u)
    return numerix.transpose(m_x_m_x_u)

def calculatePrecessionTerm(mUnit, Heff):
    PrecessionBase = numerix.cross(mUnit, Heff, axisa=0, axisb=0)
    return numerix.transpose(PrecessionBase)

def calculateDampingTerm(mUnit, Heff):
    Precession = calculatePrecessionTerm(mUnit, Heff)
    DampingBase = numerix.cross(mUnit, Precession, axisa=0, axisb=0)
    return numerix.transpose(DampingBase)


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
gridCoor = mesh.cellCenters

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
uAxis = numerix.array([[0., 0., 1.]])

# Calculate H-field for uniaxial anisotropy
HuniaxBase = HeffUniaxialAnisotropyKM(gridCoor, uAxis, Ku2, Msat)

HeffBase = HuniaxBase

TuniaxBase = calculateDampingTerm(gridCoor / numerix.linalg.norm(gridCoor), HeffBase)
TuniaxBase = numerix.multiply(TuniaxBase, alphaDamping)
TuniaxBase = TuniaxBase + calculatePrecessionTerm(gridCoor / numerix.linalg.norm(gridCoor), HeffBase)
TuniaxBase = numerix.multiply(TuniaxBase, (-1.0 * gamFac))

TeffBase = TuniaxBase
Teff = CellVariable(mesh=mesh, value=TeffBase)

phi = CellVariable(name=r"$\phi", mesh=mesh)
phi.setValue(0.25 / numerix.pi)

eqI = (TransientTerm()
      == DiffusionTerm(coeff=D)
      - ExponentialConvectionTerm(coeff=Teff)) # doctest: +GMSH

eqX = (TransientTerm()
      == ExplicitDiffusionTerm(coeff=D)
      - ExponentialConvectionTerm(coeff=Teff)) # doctest: +GMSH

steps = 10
timeStepDuration = 1e-18
time = 0
for step in range(steps):
    print("Time point is")
    print(time)
    if step < 0 :
        eqX.solve(var=phi,
                 dt=timeStepDuration)
    else:
        eqI.solve(var=phi,
                 dt=timeStepDuration)
    time = time + timeStepDuration
