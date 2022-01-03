using LoopPoly
using Documenter

DocMeta.setdocmeta!(LoopPoly, :DocTestSetup, :(using LoopPoly); recursive=true)

makedocs(;
    modules=[LoopPoly],
    authors="Yingbo Ma <mayingbo5@gmail.com> and contributors",
    repo="https://github.com/YingboMa/LoopPoly.jl/blob/{commit}{path}#{line}",
    sitename="LoopPoly.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://YingboMa.github.io/LoopPoly.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YingboMa/LoopPoly.jl",
    devbranch="master",
)
