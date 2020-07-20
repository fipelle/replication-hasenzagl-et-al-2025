function rw_benchmark(data_v, nM, nQ, data_order, last_not_na)
    #=
    Name: rw_benchmark.jl
    Description: Random walk benchmark
    =#

    # Mixed-frequency
    if nM > 0 && nQ > 0
        drift_monthly   = hcat([mean(diff(data_v[:, i])[diff(data_v[:, i]).!==missing], dims=1) for i=1:nM]...);
        drift_quarterly = hcat([mean(diff(data_v[3:3:end, i])[diff(data_v[3:3:end, i]).!==missing], dims=1) for i=nM+1:nM+nQ]...);
        drift           = [drift_monthly drift_quarterly];

    # Mono frequency
    else
        drift = hcat([mean(diff(data_v[:, i])[diff(data_v[:, i]).!==missing], dims=1) for i=1:nM+nQ]...);
    end

    output = Array{Float64}(undef, h, nM+nQ);
    for i=1:nM+nQ
        output[:,i] = data_v[last_not_na[i], i] .+ (collect(1:h).*drift[i]);
    end
    output = output[:, data_order];

    return output;
end
