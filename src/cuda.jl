using .CUDA

export CuLightMesh

struct CuLightMesh <: AbstractMesh{Face}
    vertices :: CuVector{Vertex}
    selv     :: CuVector{NTuple{3, FaceIndex}}
end

CuLightMesh(lm:: LightMesh) = CuLightMesh(CuArray(lm.vertices), CuArray(lm.selv))
LightMesh(clm:: CuLightMesh) = LightMesh(Array(clm.vertices), Array(lm.selv))

# Req
Base.size(m::CuLightMesh) = (length(m.selv),)
@inbounds function getindex(m::CuLightMesh, i::Integer)
    #vertexIndices = m.selv[i]
    vi1, vi2, vi3 = m.selv[i]
    vertices = m.vertices

    return Face( ( 
            vertices[vi1], 
            vertices[vi2], 
            vertices[vi3] ))
end
vertices(m::CuLightMesh) = m.vertices
vertexIndices(m::CuLightMesh) = m.selv