# Como la mesh viene importada para el orto, sin reutilizar vértices y duplicando, 
# esto mejor es llamarlo de a pares. Papelón
using IterTools

export refineFlat, refineSmooth, makeCurvedMesh

REFINE_VERBOSE :: Bool = false;

refineFlat(mesh ::LightMesh, target ::Integer, maxfrac :: Real = 0.2) = begin
    currFacesLog = log(length(mesh))
    targetLog = log(target)
    maxStepLog = log(1+maxfrac)
    numSteps = 1 + ceil((targetLog-currFacesLog)/maxStepLog)|>Int
    ksLog = range(currFacesLog; stop=targetLog, length=numSteps)
    ksTot = Int.(round.(exp.(ksLog)))
    ksTot[end] = target
    ks = diff(ksTot)
    refineFlat(mesh, ks)
end
refineFlat(mesh::AbstractMesh, rest...) = begin
    T = typeof(mesh)
    m = convert(LightMesh, mesh)
    return convert(T, refineFlat(m, rest...))
end
refineFlat(mesh ::LightMesh) = refineFlat(mesh, 1, 1)
refineFlat(mesh ::LightMesh, i :: Integer, k :: Integer) = refineFlat(mesh, fill(i, k))
refineFlat(mesh ::LightMesh, ks :: Vector{<:Integer}) = begin
    vertices, faces = copy(mesh.vertices), copy(mesh.selv)
    for i in eachindex(ks)
        #@info "Iteración $i de $(length(ks))"
        REFINE_VERBOSE && println("Iteración $i de $(length(ks))")
        refineFlat!(vertices, faces, ks[i])
    end
    LightMesh(vertices, faces)
end

function refineFlat!(vertices, faces, k :: Integer) 
    # Buscar aristas y catalogar caras
    edges = getEdges(faces)
    
    
    # Los ejes de máxima longitud        
    maxEdges = partialsort(keys(edges)|>collect, 1:k, by=(-) ∘ edgelength(vertices))
    
    
    # Para cada eje largo, modifico
    for e in maxEdges
        divideEdge!(edges, vertices, faces, e)
    end
    
    vertices, faces
end

function getEdges(faces)

    edges = Dict{Tuple{UInt32, UInt32}, Vector{UInt32}}()
    for i in eachindex(faces), j in [ (1, 2), (2, 3), (1, 3) ] 
        appendToDictKey!(edges, tsort(faces[i][j[1]], faces[i][j[2]]), i)
    end

    return edges
end

tsort(x, y) = x > y ? (y, x) : (x, y)


appendToDictKey!(d::Dict, key, val) = begin
    
    if haskey(d, key)
        push!(d[key], val)
    else
        d[key] = [val]
    end
    d
end

#global buf = []

function divideEdge!(edges, vertices, faces, e)
#    push!(buf, map(deepcopy, (edges, vertices, faces)))

    # Adaptá `vertices`
        # Agregá el vértice nuevo al final
    push!(vertices, (vertices[e[1]]+vertices[e[2]])/2)
    newVertexIdx = length(vertices)
#@debug "New vertex = $newVertexIdx"
#@debug "Edge = $(Int.(e))"
    # Adaptá `faces`
        # Borrá las 2 caras viejas que contenian a la arista, y armo 2 por cada una
    for fIdx in edges[e]
        oldFace = faces[fIdx]
        faces[fIdx] = replaceVal.(oldFace, e[2], newVertexIdx)
        push!(faces, replaceVal.(oldFace, e[1], newVertexIdx))

        faceExtraVertexIdx = tupleRemaining(oldFace, e)
        edges[tsort(faceExtraVertexIdx, newVertexIdx)] = [ fIdx, length(faces) ]
        replaceElement!(edges[tsort(faceExtraVertexIdx, e[2])], fIdx, length(faces))
    end
    newFacesId = length(faces).- (1:length(edges[e])) .+ 1 |> collect

    
    # Adaptá `edges`
        # Borra la arista y agrego los dos nuevos, cada uno con sus 2 caras asociadas
    edges[tsort(e[1], newVertexIdx)] = [ edges[e]... ]
    edges[tsort(e[2], newVertexIdx)] = newFacesId
    delete!(edges, e)
end

replaceElement!(a :: Array, from, to ) = begin
    for i in eachindex(a)
        if a[i] == from
            a[i] = to
        end
    end
    a
end
tupleRemaining(tup :: NTuple{3, T}, els :: NTuple{2, T}) where T = begin
    for i in tup
        if i != els[1] && i != els[2]
            return i
        end
    end
    error()
    return tup[1]
end

replaceVal(trip, oVal, nVal) = trip == oVal ? nVal : trip 
            
edgelength(vertices) = edge -> norm(vertices[edge[1]]-vertices[edge[2]])
getVertices(mesh::LightMesh) = map(Vec, mesh.vertices)
getFaces(mesh::LightMesh) = faceToTuple.(mesh.faces)





#--------------------------------------
#------ otra versión
meanSide(f, i) = begin
    vs = vertices(f)
    if i==1
        (vs[1]+vs[2])/2
    elseif i==2
        (vs[2]+vs[3])/2
    else
        (vs[3]+vs[1])/2
    end
end


refineFacesUnstitched(aggregationFunction, s)::Matrix = begin
    coords = ( (0.5, 0.), (0.5, 0.5), (0., 0.5) )
    [ begin

            u, v = coords[i]

            # caras adyacente a los lados no-`i`
            a1 = adjacentNothing(f, i%3+1)
            a2 = adjacentNothing(f, (i+1)%3+1)
            pl1 = tangentPlane(a1 != nothing ? a1 : a2::FaceAdjAux, f)
            pl2 = tangentPlane(a2 != nothing ? a2 : a1::FaceAdjAux, f)


            orig = meanSide(f, i)
            int1 = pl1(u, v)
            int2 = pl2(u, v)

            aggregationFunction(orig, int1, int2)
        end
        for i in 1:3, f in s ]
end

# asume que existe 1 vértice no compartido entre caras adyacentes
commonEdge(fs) = commonEdge(fs[1], fs[2])
commonEdge(f1, f2) = begin
    vs1 = vertices(f1)
    vs2 = vertices(f2)
    
    @inbounds index(v, vs) = v == vs[1] ? 1 : v == vs[2] ? 2 : 3
    

    v1Not2Idx = index(setdiff(vs1, vs2)[1]::Vertex, vs1)
    v2Not1Idx = index(setdiff(vs2, vs1)[1]::Vertex, vs2)
    
    
    return (f1.faceIndex, v1Not2Idx%3+1), (f2.faceIndex, v2Not1Idx%3+1)
end

partnerFaces(edge) = edge.faces |> commonEdge

stitch!(stitchingFunction, controls::Matrix, edgePairs) = begin
    for (fe1, fe2) in edgePairs
        f1, e1 = fe1
        f2, e2 = fe2
        agg = stitchingFunction(controls[e1, f1], controls[e2, f2])
        controls[e1, f1] = controls[e2, f2] = agg
    end
end
mean3(x, y, z) = (x+y+z)/3
mean2(x, y) = (x+y)/2
makeCurvedFaces(s; aggregationFunction=mean3, stitchingFunction=mean2) = 
    makeCurvedFaces(
        convert(MeshEdgeIter, s), 
        aggregationFunction=aggregationFunction, 
        stitchingFunction=stitchingFunction)
makeCurvedFaces(s::MeshEdgeIter; aggregationFunction=mean3, stitchingFunction=mean2) = begin
    cs = refineFacesUnstitched(aggregationFunction, s)

    es = edges(s)
    stitch!(stitchingFunction, cs, [partnerFaces(e) for e in es if !(e.faces[1] === e.faces[2])])
    
    controls = reinterpret(NTuple{3, Vertex}, cs, (size(cs, 2), ))
    
    return CMeshEdgeIter(s, controls)
end

refineSmooth(s; aggregationFunction=mean3, stitchingFunction=mean2) = begin
    cm =  makeCurvedFaces(s, 
        aggregationFunction=aggregationFunction, 
        stitchingFunction=stitchingFunction) |> curvedMeshToPlane

    convert(typeof(s), cm)
end



fixOriginalFaces!(of, cn) = begin
    for i in eachindex(of)
        ofi = of[i]
        cni = cn[i]
        of[i] = (ofi[1], cni[1], cni[3])
        push!(of, (cni[1], ofi[2], cni[2]))
        push!(of, (cni[1], cni[2], cni[3]))
        push!(of, (cni[3], cni[2], ofi[3]))
    end
end


curvedMeshToPlane(cm::CurvedMesh) = begin    
    s = convert(SimpleMesh, cm.mesh)
    fixOriginalFaces!(s.vec, cm.controls)
    return s
end

