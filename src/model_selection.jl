function step_through_folds(flds, evaluate_fun, geotable,save_distances=false, save_coords=false)
    dif_blocks     = [Float64[] for i in flds]
    if save_distances
        NNdist_blocks  = [Float64[] for i in flds]
    end
    if save_coords
        xcoord_blocks  = [Float64[] for i in flds]
        ycoord_blocks  = [Float64[] for i in flds]
    end
    for (j,fs) in enumerate(flds)
        # find the neighbors that the folds routine (https://github.com/JuliaEarth/GeoStatsBase.jl/blob/master/src/folding/block.jl) leaves out
        # there might be a mistake in the partitioning routine in Meshes.jl, the neighbors don't make sense (also not tested well)
        neighbors = Vector(1:length(geotable.geometry))
        deleteat!(neighbors, unique(sort([fs[1];fs[2]])))
        append!(fs[1],neighbors)

        sdat  = view(geotable, fs[1])
        stest = view(domain(geotable), fs[2])
        @assert length(sdat.Z) > length(stest)

        y_pred = evaluate_fun(fs[1],fs[2])
        dif_blocks[j] = y_pred .- geotable.Z[fs[2]]
        # Plots.scatter!(x_[fs[1]], y_[fs[1]], markersize=2, color="grey",markerstrokewidth=0)
        # Plots.scatter!(x_[fs[2]], y_[fs[2]], markersize=2, markerstrokewidth=0)
        # if j == 1
        #     break
        # end
        if save_distances
            # save distances
            ids_nn  = zeros(Int,length(stest))
            dist_nn = zeros(length(stest))
            for (i,st) in enumerate(stest)
                id_st, dist = searchdists(st, KNearestSearch(domain(sdat),1))
                ids_nn[i]  = id_st[1]
                dist_nn[i] = dist[1]
            end
            @assert length(dist_nn) == length(test_data)
            NNdist_blocks[j] = dist_nn
        end
        if save_coords
            # save coordinates
            crds = coordinates.(stest)
            xcoord_blocks[j] = first.(crds)
            ycoord_blocks[j] = last.(crds)
        end
    end
    rt = (vcat(dif_blocks...))
    if save_distances
        rt = (rt..., vcat(NNdist_blocks))
    end
    if save_coords
        rt = (rt..., vcat(xcoord_blocks), vcat(ycoord_blocks))
    end
    return rt
end
