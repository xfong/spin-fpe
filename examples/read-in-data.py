print('Import starts')
from fipy import FaceVariable, CellVariable, Gmsh2DIn3DSpace, VTKViewer, TransientTerm, ExplicitDiffusionTerm, DiffusionTerm, ExponentialConvectionTerm, DefaultSolver
from fipy.variables.variable import Variable
from fipy.tools import numerix, dump
from scipy.integrate import dblquad
import time
from shutil import copyfile
print('Import complete')

##### Internal functions
def getCellVariableDatapoint(coord):
  # expect coord to be a nested list
  #   coord[0][0] = x-coordinate
  #   coord[1][0] = y-coordinate
  #   coord[2][0] = z-coordinate
  global phi
  return phi(coord, order=1)

def sphericalDatapoint(theta, rho):
  # Expect theta and rho to be in radians
  #   theta is the angle from the +z-axis to the coordinate vector
  #   rho is the angle from the +x-axis to projection of the coordinate
  #     vector onto the xy-plane
  sin_theta = numerix.sin(theta)
  cos_theta = numerix.cos(theta)
  sin_rho = numerix.sin(rho)
  cos_rho = numerix.cos(rho)
  return sin_theta*getCellVariableDatapoint([[sin_theta*cos_rho], [sin_theta*sin_rho], [cos_theta]])

phi = dump.read(filename="phi.dat")
temp_location = [[0.], [1.], [0.]]

print(getCellVariableDatapoint(temp_location))

f = lambda y, x: sphericalDatapoint(y, x)

print(dblquad(f, 0, 2.0*numerix.pi, lambda x:  0, lambda x: 1.0*numerix.pi, epsabs=1e-13, epsrel=1e-5))
print(dblquad(f, 0, 2.0*numerix.pi, lambda x:  0, lambda x: 0.5*numerix.pi, epsabs=1e-13, epsrel=1e-5))
print(dblquad(f, 0, 2.0*numerix.pi, lambda x: 0.5*numerix.pi, lambda x: 1.0*numerix.pi, epsabs=1e-13, epsrel=1e-8))
