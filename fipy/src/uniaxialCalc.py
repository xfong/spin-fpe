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

udot2 = (uax*uax) + (uay*uay) + (uaz*uaz)
udot = numerix.sqrt(udot2)

uniAxis = FaceVariable(mesh=mesh, value=[uax/udot, uay/udot, uaz/udot])
##Heff = FaceVariable(mesh=mesh, value=[uax/udot, uay/udot, uaz/udot])

scaleHeff = numerix.dot(mesh.faceCenters, uniAxis)
Heff = coeff * scaleHeff * uniAxis

#for idx in range(len(Heff[0]) - 1):
#    numerix.put(Heff[0], idx, Heff[0][idx] * scaleHeff[idx])
#    numerix.put(Heff[1], idx, Heff[1][idx] * scaleHeff[idx])
#    numerix.put(Heff[2], idx, Heff[2][idx] * scaleHeff[idx])
print(scaleHeff)
print(Heff)
