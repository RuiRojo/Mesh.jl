module Mesh

using Reexport

using Statistics: quantile, mean
using Requires

@reexport using StaticArrays

#using DistributedArrays

export loadmesh_ASCII, loadmesh_bin, isclosed, savemesh_ASCII, savemesh_bin
export faceIter, vertexIter, vertexIndexIter, tangentPlane
export eachcoordinate
export edgelenstats
export area, boundingbox

include("MeshTypes.jl")

#DMeshType = DArray{Face}

function __init__()
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" include("cuda.jl") #=begin
        @require FoldsCUDA="6cd66ae4-5932-4b96-926d-e73e578e42cc" include("cuda.jl")
    end=#

    @require GeometryBasics="5c1252a2-5f33-56bf-86c9-59e7332b4326" begin
        
        function Base.convert(::Type{GeometryBasics.Mesh}, m::AbstractMesh)
            lm = convert(LightMesh, m)

            mv = GeometryBasics.Point.(vertices(lm))
            mf = GeometryBasics.TriangleFace.(lm.selv)
            GeometryBasics.Mesh(mv, mf)
        end
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin

            import GeometryBasics
            for fun in (:plot, :plot!, :mesh, :mesh!, :wireframe, :wireframe!)
                @eval function Makie.$fun(fig::Union{Makie.Figure, Makie.Block}, m::AbstractMesh, r...; kwargs...)
                    Makie.$fun(fig, convert(GeometryBasics.Mesh, m), r...; kwargs...)
                end
                @eval function Makie.$fun(m::AbstractMesh, r...; kwargs...)
                    Makie.$fun(convert(GeometryBasics.Mesh, m), r...; kwargs...)
                end
            end
        end
    end

end

boundingbox(m::AbstractMesh) = Tuple(extrema(c) for c in eachcoordinate(vertices(m)))

"""
    area(v1, v2, v3)
    area(face)

Return the area of a triangle.

# Example

```jldoctest
julia> using Mesh

julia> area((1, 2, 3), (4, 5, 6.2), (7, 8, 9))
0.8485281374238061
```
"""
area((v1, v2, v3) :: Face) = area(v1, v2, v3)
function area(v1, v2, v3)

    v12  = v2 .- v1 
    v13  = v3 .- v1
    ev21 = dot(v12, v12)
    ev31 = dot(v13, v13)
    dp   = dot(v12, v13)

    return sqrt(ev21 * ev31 - dp * dp) / 2
    
end


"Saves the mesh `m` to file `s` as binary STL"
function savemesh_bin(s::String, m::AbstractMesh)
    open(s, "w") do io

        nF = length(m)
        
        for i in 1:80 # write empty header
            write(io,0x00)
        end

        write(io, UInt32(nF)) # write triangle count
        for i = 1:nF
            f = m[i]
            vs = vertices(f)
            n = normal(vs)
            
            for j=1:3; write(io, n[j]); end # write normal

            for v in vs
                for j = 1:3
                    write(io, v[j]) # write vertex coordinates
                end
            end
            write(io,0x0000) # write 16bit empty bit
        end        
    end
end

"Saves the mesh `m` to file `s` as an ASCII STL"
function savemesh_ASCII(s::String, m::AbstractMesh)
    open(s, "w") do io
        
        nF = length(m)

        # write the header
        write(io,"solid vcg\n")

        # write the data
        for i = 1:nF
            f = m[i]
            v1, v2, v3 = vs = vertices(f)
            n = normal(vs)
            
            @printf io "  facet normal %e %e %e\n" n[1] n[2] n[3]
            write(io,"    outer loop\n")

            @printf io "      vertex  %e %e %e\n" v1[1] v1[2] v1[3]
            @printf io "      vertex  %e %e %e\n" v2[1] v2[2] v2[3]
            @printf io "      vertex  %e %e %e\n" v3[1] v3[2] v3[3]

            write(io,"    endloop\n")
            write(io,"  endfacet\n")
        end
        write(io,"endsolid vcg\n")
    end
    nothing
end

isFacingTowards(x, y) = sign(dot(x, y)) > zero(eltype(x))

#edgeLengthQuantile(mesh::DMeshType, q) = edgeLengthQuantile(convert(Array, mesh), q)
edgeLengthQuantile(mesh::MeshEdgeIter, q) = quantile([ norm(e[1].-e[2]) for e in edges(mesh)], q)


"""
    eachcoordinate(v)
    eachcoordinat(v, i)

Given an iterator of SVectors{N}, return the iterators for each column, 
or just for column i.
"""
function eachcoordinate(v::AbstractVector)
    return Tuple(eachcoordinate(v, i) for i in eachindex(first(v)))
end
function eachcoordinate(v::AbstractVector, i)
    return imap(pt -> pt[i], v)
end


# No recuerdo esto si sirve para algo
#=
function loadMeshBinVector(fs::String)
    open(fs) do io
        read(io, 80) # throw out header

        triangle_count = read(io, UInt32)

        mesh       = Vector{NTuple{3, Vector{Float32}}}(triangle_count)
        
        i = 0
        while !eof(io)
            read(io, Float32), read(io, Float32), read(io, Float32)
            mesh[i+1] = (
                         SVector(read(io, Float32), read(io, Float32), read(io, Float32)),
                         [read(io, Float32), read(io, Float32), read(io, Float32)],
                 [read(io, Float32), read(io, Float32), read(io, Float32)]
               )

            skip(io, 2) # throwout 16bit attribute
            i += 1
        end
        
        return mesh
    end
end=#

isclosed(m::AbstractMesh) = prod(prod.(convert(MeshAdj, m).faceConnections))==0


"Carga STL binario. Por Rui"
function loadmesh_bin(fs::String)::SimpleMesh

    open(fs) do io
        read(io, 80) # throw out header

        triangle_count = read(io, UInt32)

        mesh = Vector{NTuple{3, SVector{3, Float32}}}(undef, triangle_count)
        
        i = 0
        while !eof(io)
            read(io, Float32), read(io, Float32), read(io, Float32)
            mesh[i+1] = (
                (read(io, Float32), read(io, Float32), read(io, Float32)),
                (read(io, Float32), read(io, Float32), read(io, Float32)),
                (read(io, Float32), read(io, Float32), read(io, Float32))
                )

            skip(io, 2) # throwout 16bit attribute
            i += 1
        end
        
        return SimpleMesh(mesh)
    end
end

function loadmesh_ASCII(fs::String)::SimpleMesh
    open(fs) do io
        mesh = Vector{NTuple{3, SVector{3, Float32}}}()
        
        while !eof(io)
            line = split(lowercase(readline(io)))
            if !isempty(line) && line[1] == "facet"
                readline(io) # Throw away outerloop
                face = (
                        parse.(Float32, split(readline(io))[2:4]),
                        parse.(Float32, split(readline(io))[2:4]),
                        parse.(Float32, split(readline(io))[2:4])
                       ) .|> SVector{3, Float32}
                push!(mesh, face)
                readline(io) # throwout endloop
                readline(io) # throwout endfacet
            end
        end
        
        SimpleMesh(mesh)
    end
end


# Iteradores
struct BreadthFirstIter{N, T}
    startData :: T
    nmax :: Int
end
struct BreadthFirstIterState{N, T}
    level :: Int
    childNum :: Int
    data :: T
    parentState :: BreadthFirstIterState{N, T}

    BreadthFirstIterState{N, T}(l, cn, d::T) where {N, T} = new(l, cn, d)
    BreadthFirstIterState{N, T}(a, b, c::T, d) where {N, T} = new(a, b, c, d)
end


# Probablemente esté mal
function Base.iterate(
    it::BreadthFirstIter{N, T},
    state = BreadthFirstIterState{N, T}(1, 3, it.startData)
   ) where {N, T}

    return state.level > it.nmax ? nothing : (state.data,  brotherState(st)) 

end

#= Pre 1.0 implementation 
start(it::BreadthFirstIter{N, T}) where {N, T} = 
    BreadthFirstIterState{N, T}(0, 3, it.startData)
done(it::BreadthFirstIter, st::BreadthFirstIterState) = st.level > it.nmax
next(it::BreadthFirstIter{N, T}, st::BreadthFirstIterState{N, T}) where {N, T} = 
    st.data, brotherState(st)
    =#  

brotherState(st::BreadthFirstIterState{N, T}) where {N, T} = begin

    lv = st.level
    if lv > 0
        ps = st.parentState
        cn = st.childNum        
        chdata = st.data # Al pedo, pero si no no la reconoce como necesariamente definida

        valid = false
        while !valid && cn<N
            cn = cn + 1
            vch = childData(ps.data, cn)
            chdata = vch[2]
            valid = vch[1]
        end

        if valid
            return BreadthFirstIterState{N, T}(
                lv, 
                cn, 
                chdata,
                ps)
        else
            bro = brotherState(ps)
            return firstSonState(bro)
        end
    else        
        # caso próxima generación 
        return firstSonState(st)
    end
end

firstSonState(st::BreadthFirstIterState{N}) where N = begin
    valid = false
    cn = 0

    chdata = st.data # Al pedo, pero si no no la reconoce como necesariamente definidaa

    while !valid && cn<N
        cn = cn + 1
        valid, chdata = childData(st.data, cn)
    end

    if !valid
        error("Unhandled case, weird")
    end
    
    return typeof(st)(
        st.level + 1,
        cn,
        chdata,
        st
    )
end

faceIter(f, n) = BreadthFirstIter{3, FaceAdj}(f, n) |> unique
vertexIter(f, n) = map(vertices, faceIter(f, n)) |> Base.Iterators.Flatten |> unique
vertexIndexIter(f, n) = map(f -> f.parentMesh[f.faceIndex], faceIter(f, n)) |> Base.Iterators.Flatten |> unique

# Esto debería ser parte del tipo de dato de los iteradores supongo
# Horrible, para cubrir el caso donde no hay hijos. Diseño desagradable
childData(f::FaceAdj, n::Integer) = begin
    valid = false
    cf = f

    doIfNextFace(f, n) do child
        valid = true
        cf = child
    end

    return valid, cf
end

tangentplane(f::AbstractFace) = tangentplane(vertices(f))
tangentplane(t::NTuple{3, Vertex}) = tangentplane(t[1], t[2], t[3])
tangentplane(v1, v2, v3) = function (x, y)
    @. v1 + x * (v2-v1) + y * (v3-v1)
end

changeBasis(fFrom::AbstractFace, fTo::AbstractFace) = changeBasis(vertices(fFrom), vertices(fTo))
changeBasis(ptsFrom::NTuple{3}, ptsTo::NTuple{3}) = begin
    v0 = ptsFrom[1]
    v1 = ptsFrom[2] .- ptsFrom[1]
    v2 = ptsFrom[3] .- ptsFrom[1]
    w0 = ptsTo[1]
    w1 = ptsTo[2] .- ptsTo[1]
    w2 = ptsTo[3] .- ptsTo[1]
    
    vMat = hcat(v1, v2, cross(v1, v2))
    wMat = hcat(w1, w2, cross(w1, w2))
    
    
    wInv = inv(wMat)
    wv0  = v0 .- w0
    
    function cb(αv::SVector{3})::SVector{3}
        wInv*(vMat*αv .+ wv0)
    end    
    
    return cb
end

tangentplane(fplane::AbstractFace, fbasis::AbstractFace) = function (x, y)
    xn, yn, zn = changeBasis(fbasis, fplane)(Vertex(x, y, zero(typeof(x))))
    tangentplane(fplane)(xn, yn)
end


include("integration.jl")
using .Integration
export integrate_triangle, QuadratureRule
export Quad7, Quad21, Quad44, Quad65, @quad_gquts, @quad_gqutm, quad_gquts, quad_gqutm

include("Refine.jl")
include("convenience-tools.jl")


end # module
