# # Tutorial

# Como siempre, con el paquete instalado, para usarlo hay que correr `using Mesh`.

using Mesh


# Esto trae varios tipos de dato y funciones que se detallan en la documentación.

# Un vértice es de tipo `Vertex`, y consiste de 3 `Float32` (precisión simple).


v1 = Vertex(0, 0, 0)
v2 = Vertex(1, 0, 0)
v3 = Vertex(0, 1, 0)

v1, v2, v3

# Una cara es de tipo `Face`, y se arma con 3 vértices.

f = Face(v1, v2, v3)


# Le puedo extraer la tupla de vértices con `vertices`

vertices(f)

# La normal a una cara o a un grupo de vértices con `normal`

normal(f)

normal(vertices(f))

normal(v1, v2, v3)




# ## Mesh de ejemplo
# ​
# Los meshes que usamos están en el paquete `MeshDataset`.

using MeshDataset

# Se pueden listar los meshes que hay con `listmeshnames`

listmeshnames()[40:50]

path = getmeshpath("Esferoide_5k.lmesh")


# y conseguir el path para cargarlos con 

lm = loadMesh(path, LightMesh);



# # Mesh

# Las meshes son del tipo abstracto `AbstractMesh`.

# Uno de esos tipos es `LightMesh`, que son las que guardan internamente la data con las estructuras `selv` y `vertex`. Es la representación más comprimida: una matriz de vértices, y otra con los índices de las caras.

# El segundo argumento de `loadMesh` es el tipo de dato. La extensión `.lmesh` hace referencia a `LightMesh`.

# No debería haber necesidad de esto: podría deducir el segundo agrumento de `loadMesh` a través de la extensión. Algún día...

# # Tipos de dato