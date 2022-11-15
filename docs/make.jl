using DiffusionFlux
using Documenter

DocMeta.setdocmeta!(DiffusionFlux, :DocTestSetup, :(using DiffusionFlux); recursive=true)

makedocs(;
    modules=[DiffusionFlux],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/DiffusionFlux.jl/blob/{commit}{path}#{line}",
    sitename="DiffusionFlux.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/DiffusionFlux.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/DiffusionFlux.jl",
    devbranch="main",
)
