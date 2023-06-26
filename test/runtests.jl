using Mesh, MeshDataset, Test

mname = "Vejiga_11cm"
mpath = meshfile(mname)
mesh = loadmesh(mpath, LightMesh)

eimesh = convert(MeshEdgeIter, mesh)


# Convertir dos veces devuelve la misma mesh
@test mesh == convert(LightMesh, eimesh)

# Que deduzca la extensión
@test loadmesh(mpath) == loadmesh(mpath, LightMesh)


