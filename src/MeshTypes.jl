using StaticArrays
using Printf: @printf
import Base: size, getindex, setindex!, similar, IndexStyle, convert, promote_rule, (>>)

export loadmesh, LightMesh, MeshAdj, MeshEdgeIter, Vertex, FaceIndex, FloatType, Face
export Edge, AbstractMesh, savemesh, SimpleMesh, normal, vertices, edges
export vertexIndices, FaceAdj, faceIndex, doIfNextFace, AbstractFace, cleanupVertices
export meshMegaMatrix, CurvedMesh, CMesh, CMeshAdj, CMeshEdgeIter, CFace, CFaceAdj
export controlPoints, AbstractCurvedFace, AbstractCurvedMesh, makeMegaSelvVerticesMatrix
export adjacent
export emptyMesh, emptyMeshAdj, emptyMeshEdgeIter

#---------- Self Comments - prov
#=
Cada type tiene
    Si es grupo
        * Specs de lo que requiere para ser parte del grupo
    * Specs de lo que implementa propio

    * implementación de cosas que hereda
    * implementación de cosas propias
    * implementaciones propias generales quizá, que dependen de lo requerido

=#


#-------------------------------------
# Type definitions -------------------
#-------------------------------------


# Basic
const FloatType = Float32
const Vertex = SVector{3, FloatType} 
const FaceIndex = UInt32


"""
`AbstractFace`

Requires:

- `vertices`

Implements:

- `normal`

"""
abstract type AbstractFace <: AbstractVector{Vertex} end

Base.size(f::AbstractFace) = (3,)
getindex(f::AbstractFace, i::Integer) = vertices(f)[i]



normal(f::AbstractFace) = normal(vertices(f))
"Normal a triángulo"
normal(triang :: NTuple{3}) = cross(
    triang[2] .- triang[1], triang[3] .- triang[1]) |> normalize 
normal(v1, v2, v3) = normal((v1, v2, v3))


#-----
""" 
`AbstractMesh{FACE}` 

Requires:

- Indexing returns a `FACE`.
- `vertices.

Optional:

- Indexar un range: recortar la mesh, junto con sus conexiones
- `vertexIndices`
"""
abstract type AbstractMesh{FACE} <: AbstractVector{FACE}  end

getindex(am::AbstractMesh, i::Integer, j::Integer) = am[i][j]



#-------
"Versión más simple de un `AbstractFace`"
struct Face <: AbstractFace 
    vertices :: NTuple{3, Vertex}
end

Face(v1, v2, v3) = Face((v1, v2, v3))
vertices(f::Face) = f.vertices



#-------
" `LightMesh` es la versión más comprimida de AbstractMesh, con `selv` y `vertex`"
struct LightMesh <: AbstractMesh{Face}
    vertices :: Vector{Vertex}
    selv     :: Vector{NTuple{3, FaceIndex}}
end

    # Req
Base.size(m::LightMesh) = (length(m.selv),)
@inbounds function getindex(m::LightMesh, i::Integer)
    vertexIndices = m.selv[i]
    vertices = m.vertices

    return Face( ( 
            vertices[vertexIndices[1]], 
            vertices[vertexIndices[2]], 
            vertices[vertexIndices[3]] ))
end
vertices(m::LightMesh) = m.vertices
vertexIndices(m::LightMesh) = m.selv


    # Opt
@inbounds getindex(m::LightMesh, range) = LightMesh(m.vertices, m.selv[range])

#------------
" `SimpleMesh` es la versión más manejable, con vector de triángulos (v1, v2, v3)"
struct SimpleMesh <: AbstractMesh{Face}
    vec::Vector{NTuple{3, Vertex}}
end

    # Req
Base.size(m::SimpleMesh) = size(m.vec)
getindex(m::SimpleMesh, i::Integer) = Face(m.vec[i])
vertices(m::SimpleMesh) = unique(Vertex[ v for f in m.vec for v in f ])
vertexIndices(m::SimpleMesh) = error("Not implemented")

    # Opt
getindex(m::SimpleMesh, range) = SimpleMesh(m.vec[range])


#-------------
"""
`FaceAdjAux{MESH}`

A face that knows about adjacencies.

* `adjacent`: recorrer caras adyacentes (alias `>>`)
* `setindex!` vértices

El parámetro `MESH` es el tipo de mesh al que te da acceso, para type stability

Implementa:
* `faceIndex`
* `vertexIndices`
"""
struct FaceAdjAux{MESH} <: AbstractFace
    parentMesh :: MESH
    faceIndex  :: FaceIndex
end

vertices(f::FaceAdjAux) = begin
    vs = vertices(f.parentMesh)
    vi1, vi2, vi3 = vertexIndices(f)

    return vs[vi1], vs[vi2], vs[vi3]
end

vertexIndices(f::FaceAdjAux) = vertexIndices(f.parentMesh)[faceIndex]

faceIndex(f::FaceAdjAux) = f.faceIndex





#------------
" `MeshAdj` devuelve caras que permiten ir a una mesh `MeshAdj`"
struct MeshAdj <: AbstractMesh{FaceAdjAux{MeshAdj}}
    mesh :: LightMesh
    faceConnections :: Vector{NTuple{3, FaceIndex}}
end

    # Req
Base.size(am::MeshAdj) = size(am.mesh)
@inbounds getindex(ma::MeshAdj, i::Integer) = FaceAdj(ma, i)
vertices(am::MeshAdj) = am.mesh.vertices
vertexIndices(am::MeshAdj) = vertexIndices(am.mesh)

    # Opt
@inbounds getindex(ma::MeshAdj, range) = begin
    mesh = ma.mesh[range]

	newFaceConnections = ma.faceConnections[range]
	# Las nuevas conexiones entre caras deben reajustar sus índices
    rule = Dict{Int, Int}(zip(collect(range), 1:length(range)))
	for i in eachindex(newFaceConnections)
        newFaceConnections[i] = map(fc -> get(rule, fc, 0x000000), newFaceConnections[i])
	end

	MeshAdj(mesh, newFaceConnections)
end

#------------
"Cara que va a un `MeshAdj`"
const FaceAdj = FaceAdjAux{MeshAdj}

@inbounds adjacent(f::FaceAdj) = begin
    amesh = f.parentMesh
    indices = amesh.faceConnections[f.faceIndex]

    return FaceAdj[ amesh[index] for index in indices if index!=0]
end

# ESTO ESTA PARA EL OJETE
@inbounds adjacentNothing(f::FaceAdj, i::Integer) = begin
    amesh = f.parentMesh
    index = amesh.faceConnections[f.faceIndex][i]

    if index == 0
        return nothing
    end
    
    return amesh[index]
end

#--- FIN DE OJETE
@inbounds adjacent(f::FaceAdj, i::Integer) = begin
    amesh = f.parentMesh
    index = amesh.faceConnections[f.faceIndex][i]

    if index == 0
        error("No había cara ahí")
    end
    
    return amesh[index]
end

(>>)(f::FaceAdj, i::Integer) = adjacent(f, i)


vertexIndices(f::FaceAdj) = f.parentMesh.mesh.selv[f.faceIndex]

@inbounds Base.setindex!(f::FaceAdj, x::Vertex, i::Integer) = begin
    mesh = f.parentMesh.mesh
    mesh.vertices[mesh.selv[f.faceIndex][i]] = x
end




#= Falta
edges(f::FaceAdj) = begin
end
=#




#-----------
"""
MeshAdj que también permite seguir aristas
O sea, navegando se vuelve a un `MeshAdj`
"""
struct MeshEdgeIter <: AbstractMesh{FaceAdj}
    amesh      :: MeshAdj
    edgeVector :: Vector{Tuple{FaceIndex, UInt8}} #tupla con cara, y arista de cara 
end

    # Req
vertices(eim::MeshEdgeIter) = eim.amesh.mesh.vertices
vertexIndices(eim::MeshEdgeIter) = vertexIndices(eim.amesh)
size(mei::MeshEdgeIter) = size(mei.amesh)
@inbounds getindex(mei::MeshEdgeIter, i::Integer) = mei.amesh[i]

    # Opt
@inbounds getindex(mei::MeshEdgeIter, range) = MeshEdgeIter(mei.amesh[range], mei.edgeVector)


edges(m::MeshEdgeIter)::EdgeVector{Edge} = EdgeVector{Edge}(m)
edgeLength(eim::MeshEdgeIter) = length(eim.edgeVector)
edge(mesh::MeshEdgeIter, i::Integer)::Edge = begin
    
    fIdx, eIdx = mesh.edgeVector[i]
    face = mesh[fIdx]
    otherFace = adjacentNothing(face, eIdx)
    vsFace = vertices(face)
    vsEdge = vsFace[eIdx], vsFace[eIdx%3+1]
    
    return Edge(vsEdge, (face, otherFace != nothing ? otherFace : face ::FaceAdjAux))
end



#-----------
"""
Requiere:

* `adjacent` permite acceder a las caras adyacentes
* `vertices` indexa vértices
"""
abstract type AbstractEdge <: AbstractVector{Vertex} end

size(e::AbstractEdge) = (2,)
@inbounds getindex(e::AbstractEdge, i::Integer) = vertices(e)[i]




#-----------
"`AbstractEdge` simple"
struct Edge <: AbstractEdge
    vertices :: NTuple{2, Vertex}
    faces    :: NTuple{2, FaceAdj}
end

vertices(e::Edge) = e.vertices
adjacent(e::Edge) = e.faces
adjacent(e::Edge, i) = e.faces[i]




#---------
"""
Requiere
* Indexar un edge

Debe guardar un tipo de mesh que incorpore 
* `edges(m::AbstractMesh)::EdgeVector{E<:AbstractEdge} = EdgeVector{E}(m)`
* `edge(::AbstractMesh, i)::EDGE`
* `edgeLength(::AbstractMesh)`
"""
struct EdgeVector{EDGE} <: AbstractVector{EDGE}
    mesh ::AbstractMesh
end



size(ev::EdgeVector) = (edgeLength(ev.mesh),)
getindex(ev::EdgeVector{EDGE}, i::Integer) where EDGE = edge(ev.mesh, i)::EDGE
#IndexStyle(::Type{<:EdgeVector}) = IndexLinear()

    



#------------
struct CurvedEdge <: AbstractEdge
    edge :: Edge
    control  :: Vertex
end



vertices(e::CurvedEdge) = vertices(e.edge)
adjacent(e::Edge, r...) = adjacent(e.edge, r...)





#-------------
"Implementa `controlPoints` y `vertices`"
abstract type AbstractCurvedFace <: AbstractFace end



#------------
struct CFaceAdj <: AbstractCurvedFace
    face :: FaceAdj
    controls :: Vector{NTuple{3, Vertex}}
end


vertices(cf::CFaceAdj) = vertices(cf.face)
controlPoints(cf::CFaceAdj) = cf.controls[face.faceIndex]


#------------
struct CFace <: AbstractCurvedFace
    face :: Face
    control :: NTuple{3, Vertex}
end

vertices(cf::CFace) = vertices(cf.face)
controlPoints(cf::CFace) = cf.control



#------------
"""
`AbstractCurvedMesh` representa mesh curva, 
* `controlPoints`
* Indexa una `FACE`.
* `vertices`
"""
abstract type AbstractCurvedMesh{FACE} <: AbstractMesh{FACE} end




#------------
"""
* Indexa una `FACE`.
"""
struct CurvedMesh{MESH,FACE} <: AbstractCurvedMesh{FACE}
    mesh :: MESH
    controls :: Vector{NTuple{3, Vertex}}
end

    # Req
vertices(cm::CurvedMesh) = vertices(cm.mesh)
controlPoints(cm::CurvedMesh) = cm.controls
size(c::CurvedMesh) = size(c.mesh)

    # Opt
getindex(cm::CurvedMesh, range) = begin
    typeof(cm)(
        cm.mesh[range],
        controls[range]
    )
end



edges(m::CurvedMesh) = EdgeVector{CurvedEdge}(m)
edge(m::CurvedMesh, i::Integer) = CurvedEdge(edge(m.mesh, i), m.controls[i])
edgeLength(m::CurvedMesh) = edgeLength(m.mesh)

#----------------
const CMeshAdj = CurvedMesh{MeshAdj, CFaceAdj}
const CMesh    = CurvedMesh{LightMesh, CFace}
const CMeshEdgeIter = CurvedMesh{MeshEdgeIter, CFaceAdj}
const CSimpleMesh = CurvedMesh{SimpleMesh, CFace}



getindex(cm::CurvedMesh{MESH, CFace} where MESH, i::Integer) = CFace(cm.mesh[i], cm.controls[i]) 
getindex(cm::CurvedMesh{MESH, CFaceAdj} where MESH, i::Integer) = CFaceAdj(cm.mesh[i], cm.controls)


#----------------------------------------------------------------
# FIN DE TYPES --------------------------------------------------
#----------------------------------------------------------------


    # Empty meshes
const emptyMesh = LightMesh([],[])
const emptyMeshAdj = MeshAdj(emptyMesh, [])
const emptyMeshEdgeIter = MeshEdgeIter(emptyMeshAdj, [])




extractMesh(::Type{MeshEdgeIter}, e::EdgeVector{Edge}) = e.emesh
extractMesh(::Type{MeshEdgeIter}, e::EdgeVector{CurvedEdge}) = e.emesh.mesh



@inbounds doIfNextFace(fun, f::FaceAdj, i::Integer) = begin
    amesh = f.parentMesh
    index = amesh.faceConnections[f.faceIndex][i]

    if index > 0
        fun(amesh[index]) 
    end

    return nothing
end





edges(::CMesh) = error("undefined")
edges(m::AbstractMesh) = begin
    # @warn "Converting to MeshEdgeIter"
    return edges(convert(MeshEdgeIter, m))
end




@inbounds Base.setindex!(m::LightMesh, x::Vertex, i::Integer, j::Integer) = begin
	m.vertices[m.selv[i][j]] = x
end
Base.setindex!(ma::MeshAdj, x::Vertex, i::Integer, j::Integer) = (ma.mesh[i,j]=x)
Base.setindex!(mei::MeshEdgeIter, x::Vertex, i::Integer, j::Integer) = (mei.amesh[i,j]=x)




# toma 3 tupla, devuelve 2 tupla de "eje" según `eIdx`
function edgeVertexIndices(fTup :: NTuple{3}, eIdx :: Integer)
    v1, v2 = fTup[eIdx], fTup[eIdx%3+1]
    v1 < v2 ? (v1, v2) : (v2, v1)
end
edgeVertexIndices(fTup) = map(i->edgeVertexIndices(fTup, i), (1, 2, 3))





"""
Borra los vértices que no pertenecen a ninguna cara.
"""
cleanupVertices(cmesh::CMesh) = CMesh(
    cleanupVertices(cmesh.mesh), 
    cmesh.controls)
cleanupVertices(mesh::LightMesh) = convert(LightMesh, convert(SimpleMesh, mesh))

                            
                            

savemesh(s::String, am::AbstractCurvedMesh) =  open(s, "w") do io      
    
    selv, vertices = makeMegaSelvVerticesMatrix(am)
                                
    write(io, length(vertices))
    for v in vertices, el in v
        write(io, el)
    end
    write(io, length(selv))
    for f in selv, i in f
        write(io, i)
    end
end
 
                                
#=                   
savemesh_ASCII(s::String, am::AbstractCurvedMesh) =                                
    open(s, "w") do io           
        @printf io "%i %i\n" length(vertices) length(selv)

        for i in 1:length(vertices)
            v = vertices[i]
            @printf io "%e %e %e %e\n" i v[1] v[2] v[3]
        end

        for i in 1:length(selv)
            f = selv[i]
            @printf io "%i %i %i %i %i %i %i %i\n" i 206 f[1] f[2] f[3] f[4] f[5] f[6]
        end
    end 
end
=#

meshMegaMatrix(m::SimpleMesh; normals=false) = begin
    if !normals
        return reinterpret(Float32, m.vec, (3, 3, length(m.vec)))
    else
        mnvec = map(tup -> (tup..., normal(tup)), m.vec)
        return reinterpret(Float32, mnvec, (3, 4, length(m.vec)))
    end
end

makeMegaSelvVerticesMatrix(am::AbstractMesh) = begin
    m = convert(LightMesh, am)
    
    selv = reinterpret(UInt32, m.selv, (length(m.selv), 3))
    vertices = reinterpret(Float32, m.vertices, (length(m.vertices), 3))
end

makeMegaSelvVerticesMatrix(am::AbstractCurvedMesh) = begin
    m = convert(LightMesh, am.mesh)

    vlen = length(m.vertices)
    selvBase = m.selv
    
    vertices = append!(copy(m.vertices), Base.Iterators.flatten(am.controls))
    selv = [ riffle(selvBase[i], vlen+3*(i-1).+(1:3)) for i in eachindex(selvBase)]
    
    return selv, vertices
end





# CONVERSION --------------------------
convert(::Type{Face}, fa::FaceAdj) = fa.parentMesh.mesh[fa.faceIndex]
convert(::Type{CurvedMesh{MESH, FACE}}, cm::CurvedMesh{MESHFROM, F}) where {MESH, FACE, MESHFROM, F} = begin
    CurvedMesh{MESH, FACE}(convert(LightMesh, cm.mesh), cm.controls)
end

convert(::Type{LightMesh}, m :: SimpleMesh) = begin
    mesh = m.vec
    vs   = unique(v for f in mesh for v in f)    
    selv = map.(nestedReplaceByIndex(vs, UInt32), mesh)

    return LightMesh(vs, selv)
end
convert(::Type{MeshEdgeIter}, m::SimpleMesh) = convert(MeshEdgeIter, convert(LightMesh, m))
convert(::Type{SimpleMesh}, x::SimpleMesh) = x
convert(::Type{SimpleMesh}, m::AbstractMesh) = SimpleMesh(map(vertices, m))
convert(::Type{MeshEdgeIter}, m::LightMesh) = begin
    selv = m.selv
    
    edgeDict = group(
        feidxs -> edgeVertexIndices(selv[feidxs[1]], feidxs[2]), 
        (fIdx::FaceIndex, eIdx::UInt8) 
        for fIdx in UInt32(1):UInt32(length(selv)) 
            for eIdx in 0x01:0x03
    )

    edgeVector = collect(Tuple{FaceIndex, UInt8}, map(first, values(edgeDict)))
    
    faceConnections = map(UInt32(1):UInt32(length(selv))) do fIdx
        faceVertexIndices = selv[fIdx]
        
        map(edgeVertexIndices(faceVertexIndices)) do evi
            # Busco los índices (fIdx, eIdx) de las caras con ese arista
            facesWithEdge = edgeDict[evi]
            out ::UInt32 = zero(UInt32)
            for (fi, _) in facesWithEdge
                if fi != fIdx
                    out = fi
                end
            end
            out
        end
    end
                
    meshAdj = MeshAdj(m, faceConnections)
    return MeshEdgeIter(meshAdj, edgeVector)
end

convert(::Type{MeshEdgeIter}, m::MeshAdj) = convert(MeshEdgeIter, convert(LightMesh, m))
convert(::Type{LightMesh}, m::MeshAdj) = m.mesh
convert(::Type{LightMesh}, m::MeshEdgeIter) = m.amesh.mesh
convert(::Type{MeshAdj}, m::MeshEdgeIter) = m.amesh
convert(::Type{T}, x::T) where T <: AbstractMesh = x
convert(::Type{T}, x::F) where F <: AbstractMesh where T <: AbstractMesh = convert(T, convert(promote_type(T, F), x))
promote_rule(::Type{<:AbstractMesh}, ::Type{<:AbstractMesh}) = MeshEdgeIter


const exts = Dict(
    "smesh"=>SimpleMesh, 
    "lmesh"=>LightMesh, 
    "ma"=>MeshAdj, 
    "mei"=>MeshEdgeIter)

extension(fn) = last(splitext(fn))[2:end]
loadmesh(s::String) = loadmesh(s, exts[extension(s)])


# Load and save
savemesh(s::String, m) = open(io->savemesh(io, m), s, "w")
loadmesh(s::String, t) = open(io->loadmesh(io, t), s)

function savemesh(io::IOStream, m::SimpleMesh)
    write(io, length(m))
    for i in m.vec, v in i, p in v
        write(io, p)
    end
end
function savemesh(io::IOStream, m::LightMesh) 
    write(io, length(m.vertices))
    for v in m.vertices, el in v
        write(io, el)
    end
    write(io, length(m.selv))
    for f in m.selv, i in f
        write(io, i)
    end
end

function savemesh(io::IOStream, ma::MeshAdj)
    savemesh(io, ma.mesh)
    faceConnections = ma.faceConnections
    fcc = length(faceConnections)
    write(io, fcc)
    for fc in faceConnections, i in fc
        write(io, i)
    end
end

function savemesh(io::IOStream, mei::MeshEdgeIter)
    savemesh(io, mei.amesh)
    edgeVector = mei.edgeVector
    ecc = length(edgeVector)
    write(io, ecc)
    for e in edgeVector, i in e
        write(io, i)
    end
end

function loadmesh(io::IOStream, ::Type{SimpleMesh})
    fn = read(io, Int64)
    vec = Vector{NTuple{3, Vertex}}(fn)
    for i in 1:fn
        vec[i] = (
                 SVector(read(io, FloatType), read(io, FloatType), read(io, FloatType)),
                 SVector(read(io, FloatType), read(io, FloatType), read(io, FloatType)),
                 SVector(read(io, FloatType), read(io, FloatType), read(io, FloatType))
                )
    end

    return SimpleMesh(vec)
end

function loadmesh(io::IOStream, ::Type{LightMesh})
    vn = read(io, Int64)
    vertices = Vector{Vertex}(undef, vn)
    for i in 1:vn
        vertices[i] = SVector(read(io, FloatType), read(io, FloatType), read(io, FloatType))
    end
    fn = read(io, Int64)
    selv = Vector{NTuple{3, FaceIndex}}(undef, fn)
    for i in 1:fn
        selv[i] = (read(io, FaceIndex), read(io, FaceIndex), read(io, FaceIndex))
    end
    
    LightMesh(vertices, selv)
end

function loadmesh(io::IOStream, ::Type{MeshAdj})
    mesh = loadmesh(io, Mesh)
    fcs = read(io, Int64)
    faceConnections = Vector{NTuple{3, FaceIndex}}(undef, fcs)
    
    for i in 1:fcs
        faceConnections[i] = (read(io, FaceIndex), read(io, FaceIndex), read(io, FaceIndex))
    end
    
    MeshAdj(mesh, faceConnections)
end

function loadmesh(io::IOStream, ::Type{MeshEdgeIter})
    amesh = loadmesh(io, MeshAdj)
    ecs = read(io, Int64)
    edgeVector = Vector{Tuple{FaceIndex, UInt8}}(ecs)
    
    for i in 1:ecs
        edgeVector[i] = (read(io, FaceIndex), read(io, UInt8))
    end
    
    MeshEdgeIter(amesh, edgeVector)
end


include("import-aux-ed.jl")


#= BROKEN
loadmesh(s::IOStream, T::Type{<:AbstractCurvedMesh}) = convert(T, loadmesh(s, CMesh))
loadmesh(s::IOStream, ::Type{CMesh}) = begin


    nv::Int, nf::Int = readline(io) |> split .|> parse


    vertices = Vertex[ 
        readline(io) |> rest ∘ split .|> parse |> SVector{3} 
        for i in 1:nv ]

    selv     = Vector{NTuple{3, FaceIndex}}(undef, nf)
    controls = Vector{NTuple{3, Vertex}}(undef, nf)

    for i in 1:nf
        # Leo la próxima línea correspondiente a triángulo curvo
        sline = readline(io) |> split
        while sline[2] != "206" # triángulo curvo
            sline = readline(io) |> split
        end
        line = parse.(sline)

        selv[i]     = FaceIndex.((line[3], line[5], line[7]))
        controls[i] = convert.(Vertex, (
            vertices[line[4]],
            vertices[line[6]],
            vertices[line[8]]))
    end

    cmesh = CMesh(Mesh(vertices, selv), controls)

    return cleanupVertices(cmesh)
    end
                                                        
                                                        

end

                                                                                                     =#
####### Algunas medio generales

function nestedReplaceByIndex(v :: Vector, t=Int)
    d = Dict(zip(v, UnitRange{t}(1, length(v))))
    c -> d[c]
end


"Agrupar comunes"
function group(f, itr)

    T = eltype(itr)
    K = Base.promote_op(f, T)
    d = Dict{K, Vector{T}}()
    
    for x in itr
        key = f(x)
        if !haskey(d, key)
            d[key] = [x]
        else
            push!(d[key], x)
        end
    end
    d
end                                                        

riffle(v, c) = zip(v, c) |> Base.Iterators.flatten |> collect
rest(l) = l[2:end]
