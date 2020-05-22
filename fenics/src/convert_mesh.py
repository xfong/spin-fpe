import meshio
msh = meshio.read("sphere_msh20.msh")

for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data
    elif  cell.type == "tetra":
        tetra_cells = cell.data

for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]
        triangle_mesh =meshio.Mesh(points=msh.points,
                                   cells=[("triangle", triangle_cells)],
                                   cell_data={"name_to_read":[triangle_data]})
        meshio.write("sphere_mf.xdmf", triangle_mesh)
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]
        tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})
        meshio.write("sphere_mesh.xdmf", tetra_mesh)

