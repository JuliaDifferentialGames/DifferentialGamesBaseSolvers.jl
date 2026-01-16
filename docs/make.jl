using DifferentialGamesBaseSolvers
using Documenter

DocMeta.setdocmeta!(DifferentialGamesBaseSolvers, :DocTestSetup, :(using DifferentialGamesBaseSolvers); recursive=true)

makedocs(;
    modules=[DifferentialGamesBaseSolvers],
    authors="BennetOutland <bennet.outland@pm.me> and contributors",
    sitename="DifferentialGamesBaseSolvers.jl",
    format=Documenter.HTML(;
        canonical="https://BennetOutland.github.io/DifferentialGamesBaseSolvers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/BennetOutland/DifferentialGamesBaseSolvers.jl",
    devbranch="main",
)
