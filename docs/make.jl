using Documenter, QuixoticSimulating

makedocs(;
    modules=[QuixoticSimulating],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/BradenDKelly/QuixoticSimulating.jl/blob/{commit}{path}#L{line}",
    sitename="QuixoticSimulating.jl",
    authors="Braden D. Kelly",
    assets=String[],
)
