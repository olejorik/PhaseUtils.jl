# using Pkg
# Pkg.develop(PackageSpec(; path=pwd()))
# Pkg.instantiate()
using Documenter, Literate
using PhaseUtils
# varinfo(PhaseUtils)

DocMeta.setdocmeta!(PhaseUtils, :DocTestSetup, :(using PhaseUtils); recursive=true)

# TODO Generate tutorials from example files
@info "current dir =$(@__DIR__)"
tutorials_folder = (@__DIR__) * "/../examples"
docs_tutorials_folder = (@__DIR__) * "/src/examples"
@info tutorials_folder
for f in readdir(tutorials_folder; join=true)
    Literate.markdown(f, docs_tutorials_folder)
end

makedocs(;
    sitename="PhaseUtils",
    modules=[PhaseUtils],
    authors="Oleg Soloviev",
    # repo="https://github.com/olejorik/PhaseUtils.jl/blob/{commit}{path}#L{line}",
    checkdocs=:exports,
    # doctest=:fix,
    format=Documenter.HTML(;
        # Use clean URLs, unless built as a "local" build
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseUtils.jl/stable/",
        assets=["assets/favicon.ico"],
        highlights=["yaml"],
    ),
    clean=false,
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        "Examples" => [
            "examples/Poisson.md",
            "examples/Unwrapping.md",
            "examples/Unwrapping_method.md",
        ],
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
deploydocs(; repo="github.com/olejorik/PhaseUtils.jl.git", target="build")
