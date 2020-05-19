using Documenter, QuixoticSimulating

makedocs(;
    modules=[QuixoticSimulating],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/BradenDKelly/QuixoticSimulating/blob/{commit}{path}#L{line}",
    sitename="QuixoticSimulating",
    authors="Braden D. Kelly",
    assets=String[],
)
