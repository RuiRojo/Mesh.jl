function edgelenstats(mesh) 
    edgelens = [ norm(v1 .- v2) for face in mesh for (v1, v2) in subsets(face, 2) ]

    return (
        deciles = map(q -> quantile(edgelens, q), 0.1:0.1:0.9),
        mean = mean(edgelens),
        minimum = minimum(edgelens),
        maximum = maximum(edgelens)
    )        
end