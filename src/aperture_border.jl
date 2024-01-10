
alldirs = (:left, :right, :up, :down)

cw = (; :up => :left, :left => :down, :down => :right, :right => :up)

ccw = (; :up => :right, :left => :up, :down => :left, :right => :down)

# Using these symbols, we can define also the steps in the needed direction
step = (
    up=CartesianIndex(-1, 0),
    left=CartesianIndex(0, -1),
    down=CartesianIndex(1, 0),
    right=CartesianIndex(0, 1),
)


"""
    _find_edges(ap)

Return arrays of the left, right, top, and bottom edges of the aperture
"""
function _find_edges(ap)
    LE = findall(ap - circshift(ap, (0, 1)) .== 1)
    RE = findall(ap - circshift(ap, (0, -1)) .== 1)
    TE = findall(ap - circshift(ap, (1, 0)) .== 1)
    BE = findall(ap - circshift(ap, (-1, 0)) .== 1)
    return LE, RE, TE, BE
end  # function _find_edges



"""
    _find_set_edges(ap)

Return sets of the left, right, top, and bottom edges of the aperture
"""
function _find_set_edges(ap)
    return NamedTuple{alldirs}(_find_edges(ap))
end  # function _find_edges


function _find_cw_border_push(ap)
    edge = _find_set_edges(ap)

    # set any direction #TODO rewrtie as afunction
    direction = :left
    # We begin wiht the position in the left edge, so we cannot move to the left
    # Let's check the other directions

    startpos = first(edge[direction])
    for d in setdiff(alldirs, [direction])
        if startpos + step[d] in edge[ccw[d]]
            direction = d
            break
        end
    end


    # Now we can move along this direction as long as we haven't come to the perpedicular edge

    # This is approach without preallocation
    cw_cont_ap = CartesianIndex{2}[]
    push!(cw_cont_ap, startpos)
    pos = startpos + step[direction]
    while pos != startpos
        push!(cw_cont_ap, pos)
        # @show pos
        # println([pos in edge[d] for d in alldirs])
        if pos in edge[ccw[direction]] # we are still on the border
            # @info "we are on the $(ccw[direction]) edge going $direction"
            while pos in edge[direction] # then change direction
                direction = cw[direction]
                # @info "Direction changed to $direction"
            end
        else # we're out of the border, turn back
            direction = ccw[direction]
        end
        # @show direction # debug
        pos = pos + step[direction]
    end
    return cw_cont_ap, edge

end


function _find_cw_border_alloc(ap; outside=false)
    if outside
        (cwlike, ccwlike) = (ccw, cw)
        ap .= 1 .- ap
    else
        (cwlike, ccwlike) = (cw, ccw)
    end

    edge = _find_set_edges(ap)

    # set any direction #TODO rewrtie as afunction
    direction = :left
    # We begin wiht the position in the left edge, so we cannot move to the left
    # Let's check the other directions

    startpos = first(edge[direction])
    for d in setdiff(alldirs, [direction])
        if startpos + step[d] in edge[ccwlike[d]]
            direction = d
            break
        end
    end


    # Approach with preallocation
    cw_cont_ap = zeros(CartesianIndex{2}, sum([length(e) for e in edge]) - 4) # total amount of edges - left turns (belongs to two edges) + right turns (no edge). Should be improved if we have points belonging to three edges


    pos = startpos
    for i in eachindex(cw_cont_ap)
        cw_cont_ap[i] = pos
        # @show pos
        # println([pos in edge[d] for d in alldirs])
        if pos in edge[ccwlike[direction]] # we are still on the border
            # @info "we are on the $(ccwlike[direction]) edge going $direction"
            while pos in edge[direction] # then change direction
                direction = cwlike[direction]
                # @info "Direction changed to $direction"
            end
        else # we're out of the border, turn back
            direction = ccwlike[direction]
        end
        # @show direction # debug
        pos = pos + step[direction]
    end
    return cw_cont_ap, edge
end

# Both functions seem to be equally fast. Expport one of them
find_cw_border = _find_cw_border_alloc
export find_cw_border
