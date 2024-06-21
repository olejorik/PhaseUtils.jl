# Simple utilitilities to read and save files
#


function readfiles(filedir, filenamemask::Array; loadfunction=load)
    filenames = readdir(filedir)
    filenames = filter(x -> any(occursin.(filenamemask, Ref(x))), filenames)
    fileimages = [loadfunction(joinpath(filedir, filename)) for filename in filenames]

    return filenames, fileimages
end

readfiles(filedir, filenamemask::Union{String,Regex}; kwargs...) =
    readfiles(filedir, [filenamemask]; kwargs...)
