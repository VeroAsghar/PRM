using Documenter, PRM

makedocs(
    modules = [PRM],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "gralvarek",
    sitename = "PRM.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/Vero/PRM.jl.git",
    push_preview = true
)
