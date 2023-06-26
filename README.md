# Mesh


Package to handle meshes.

STL can be loaded and saved both in binary (`loadmesh_bin`/`savemesh_bin`) and ASCII (`loadmesh_ASCII`/`savemesh_ASCII`) formats.

Meshes are transformed into one of multiple types: `LightMesh`, `SimpleMesh`, `MeshAdj` and `MeshEdgeIter`,
which differ in internal representation, efficiency, and operations allowed.
For example, `LightMesh` is the most compact memory-wise, and `MeshEdgeIter` allows for iterating over the edges.
These can be freely converted into each other with `Base.convert`, and saved as a binary retaining the representation
for more efficient future importing with functions `loadmesh` and `savemesh`.

All mesh types are subtypes of `AbstractMesh`, which itself is an `AbstractVector` whose elements are the mesh faces.
Faces are represented by the type `Face` and are tuples of vertices, which are `StaticVector` of length 3, each coordinate of type `Float32`. 

Edges are of type `Edge` and obtained from the function `edges(::AbstractMesh)`. They can be indexed to get the edge
endpoints, and they support some natural functions like `norm(::Edge)`.

Function `vertices` returns a vector of the vertices of the mesh, face, or edge.

## Integration

Function `integrate_triangle` allows for integrating an arbitrary function defined on a face or a full mesh. This
is done using a default quadrature rule that can be customized with the keyword argument `quad`.
Quadrature rules are represented by the type `QuadratureRule` which contain points and their respective weights under
a canonical triangle. Macros `@quad_gquts(n::Int)` and `@quad_gqutm(n::Int)` are provided for easy and efficient (compile-time) generation of quadrature rules of $n^2$ and $n * (n + 1) / 2$ points respectively.

## Other

Other functions provided include

- `area` computes the area of a face or a mesh. 
- `edgelenstats` provides an overview of statistics of the edges of the mesh.
- `boundingbox`
- `faceIter`/`vertexIter` gives an iterator over the faces/vertices
- `refineFlat` refines the mesh into smaller triangles while retaining the geometry.

There's also a primitive incomplete implementation of curved meshes.