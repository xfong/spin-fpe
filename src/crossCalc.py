# Solves the FPE formulation of monodomain LLG on the surface of a unit sphere
#
# With only uniaxial unisotropy and initial uniform distribution of PDF over the
# surface of the unit sphere, the result of the simulation should give a
# distribution with two peaks, one at each poles of the sphere defined by the
# axis of the anisotropy
#
from fipy import FaceVariable, CellVariable, Grid3D
from fipy.tools import numerix

mesh = Grid3D(dx = 0.1, dy = 0.1, dz = 0.1,
              nx = 10, ny = 10, nz = 10)

uax = 0.
uay = 0.
uaz = 1.
coeff = 10.

xVector = FaceVariable(mesh=mesh, value=[0., 0., 0.])
mVal = mesh.faceCenters
baseMask = numerix.zeros(len(mVal[0]), 'l')
chkVector = numerix.zeros((3,len(mVal[0])), 'd')

for idx in range(len(mVal[0]) - 1):
    mx0 = mVal[0][idx]
    my0 = mVal[1][idx]
    mz0 = mVal[2][idx]
    mnorm2 = (mx0 * mx0) + (my0 * my0) + (mz0 * mz0)
    mnorm = numerix.sqrt(mnorm2)
    mx = mx0 / mnorm
    my = my0 / mnorm
    mz = mz0 / mnorm
    mxu_x = my*uaz - mz*uay
    mxu_y = mz*uax - mx*uaz
    mxu_z = mx*uay - my*uax
    chkVector[0][idx] = mxu_x
    chkVector[1][idx] = mxu_y
    chkVector[2][idx] = mxu_z
    print('XVal')
    print(chkVector[0][idx])
    print('YVal')
    print(chkVector[1][idx])
    print('ZVal')
    print(chkVector[2][idx])

xVector.setValue(chkVector)
print(xVector)
