# Simple utilities to read and save files
#


"""
    natural_sort(filenames)

Sort filenames by the first integer found in each name, then lexicographically.
Useful for files named `1.tif`, `2.tif`, …, `10.tif` which would otherwise sort
as `1, 10, 2, …` under plain alphabetical ordering.
"""
function natural_sort(filenames)
    return sort(filenames; by=f -> begin
        m = match(r"(\d+)", f)
        m === nothing ? (0, f) : (parse(Int, m[1]), f)
    end)
end

"""
    readfiles(filedir, filenamemask; loadfunction=load, sorter=sort)

Load all files in `filedir` whose names match `filenamemask` (a `String`, `Regex`,
or array thereof).  Returns `(filenames, fileimages)` where `filenames` is the
filtered list of base names and `fileimages` is the corresponding vector of loaded
objects.

# Keyword arguments
- `loadfunction`: callable used to load each file (default: `FileIO.load`).
- `sorter`: function applied to the filtered filename list before loading
  (default: `sort`, i.e. lexicographic order).  Pass `natural_sort` to sort by
  the first integer found in each filename — useful when files are named
  `1.tif`, `2.tif`, …, `10.tif` which would otherwise sort as `1, 10, 2, …`.

# Examples
```julia
## Lexicographic order (default)
names, imgs = readfiles("/data/run42", r".*\\.tif\$")

## Numeric order
names, imgs = readfiles("/data/run42", r".*\\.tif\$"; sorter=natural_sort)
```
"""
function readfiles(filedir, filenamemask::Array; loadfunction=load, sorter=sort)
    filenames = readdir(filedir)
    filenames = filter(x -> any(occursin.(filenamemask, Ref(x))), filenames)
    filenames = sorter(filenames)
    fileimages = [loadfunction(joinpath(filedir, filename)) for filename in filenames]

    return filenames, fileimages
end

readfiles(filedir, filenamemask::Union{String,Regex}; kwargs...) =
    readfiles(filedir, [filenamemask]; kwargs...)
