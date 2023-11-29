using QRWithoutQ
using Documenter

DocMeta.setdocmeta!(QRWithoutQ, :DocTestSetup, :(using QRWithoutQ); recursive = true)

makedocs(;
    modules = [QRWithoutQ],
    authors = "Filip Tronarp <filip.tronarp@matstat.lu.se> and contributors",
    repo = "https://github.com/filtron/QRWithoutQ.jl/blob/{commit}{path}#{line}",
    sitename = "QRWithoutQ.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://filtron.github.io/QRWithoutQ.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/filtron/QRWithoutQ.jl", devbranch = "main")
