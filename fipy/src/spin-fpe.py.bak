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

# define function to calculate uniaxial anisotropy
def uniaxialHeff(m0, uax, coeff):
    uaxx = uax[0]
    uaxy = uax[1]
    uaxz = uax[2]
    unorm2 = (uaxx * uaxx) + (uaxy * uaxy) + (uaxz * uaxz)
    unorm = numerix.sqrt(unorm2)
    heffx = uaxx / unorm
    heffy = uaxy / unorm
    heffz = uaxz / unorm
    mdotu = (m0[0] * heffx) + (m0[1] * heffy) + (m0[2] * heffz)
    HeffMag = [[ (mdotu * heffx), (mdotu * heffy), (mdotu * heffz) ]]

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
Xgrid = mesh.faceCenters[0]
Ygrid = mesh.faceCenters[1]
Zgrid = mesh.faceCenters[2]
mmag.setValue(uniaxialHeff([[ Xgrid, Ygrid, Zgrid ]], numerix.array((0., 0., 1.0)), 1.0))
