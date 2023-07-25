# using Pkg
# Pkg.develop(PackageSpec(; path=pwd()))
# Pkg.instantiate()
using Documenter
using PhaseUtils
# varinfo(PhaseUtils)

DocMeta.setdocmeta!(PhaseUtils, :DocTestSetup, :(using PhaseUtils); recursive=true)

makedocs(;
    sitename="PhaseUtils",
    format=Documenter.HTML(),
    modules=[PhaseUtils],
    pages=["Home" => "index.md"],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
deploydocs(; repo="github.com/olejorik/PhaseUtils.jl.git", target="build")
