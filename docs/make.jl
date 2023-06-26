using Mesh
using Documenter

DocMeta.setdocmeta!(Mesh, :DocTestSetup, :(using Mesh); recursive=true)

makedocs(;
    modules=[Mesh],
    authors="Rui Rojo <rui.rojo@gmail.com>",
    repo="https://gitlab.com/das-ara/priv/julia/Mesh.jl/blob/{commit}{path}#{line}",
    sitename="Mesh.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://das-ara.gitlab.io/priv/julia/Mesh.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

