using SIMDPolynomials
using Documenter

DocMeta.setdocmeta!(SIMDPolynomials, :DocTestSetup, :(using SIMDPolynomials); recursive=true)

makedocs(;
    modules=[SIMDPolynomials],
    authors="Yingbo Ma <mayingbo5@gmail.com> and contributors",
    repo="https://github.com/YingboMa/SIMDPolynomials.jl/blob/{commit}{path}#{line}",
    sitename="SIMDPolynomials.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://YingboMa.github.io/SIMDPolynomials.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YingboMa/SIMDPolynomials.jl",
    devbranch="master",
)
