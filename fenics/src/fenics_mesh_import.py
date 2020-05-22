import petsc4py
import fenics
from dolfin import *

msh = Mesh()
with XDMFFile(MPI.comm_world,
              "sphere_mf.xdmf") as xdmf_infile:
    xdmf_infile.read(msh)
