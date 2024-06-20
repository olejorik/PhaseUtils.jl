# Simple utilitilities to read and save files
#


function readfiles(filedir, filenamemask::Array)
    filenames = readdir(filedir)
    filenames = filter(x -> any(occursin.(filenamemask, Ref(x))), filenames)
    fileimages = [load(joinpath(filedir, filename)) for filename in filenames]

    return filenames, fileimages
end

readfiles(filedir, filenamemask::Union{String,Regex}) = readfiles(filedir, [filenamemask])
