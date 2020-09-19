function standardize_data(data, nM, nQ, h, data_order, MNEMONIC=[], σʸ=[])

    is_σʸ_empty = isempty(σʸ) == true;

    if is_σʸ_empty

        # Mixed-frequency
        if nM > 0 && nQ > 0
             σʸ_monthly   = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM]...);
             σʸ_quarterly = hcat([std(collect(skipmissing(diff(data[3:3:end, i]))), dims=1) for i=nM+1:nM+nQ]...);
             σʸ           = [σʸ_monthly σʸ_quarterly];

        # Mono frequency
        else
             σʸ = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM+nQ]...);
        end

        data = data./σʸ;
    end

    # Re-order
    quarterly_position = [zeros(nM); ones(nQ)];

    if isempty(data_order) == false
        data               = data[:, data_order];
        quarterly_position = quarterly_position[data_order];

        if is_σʸ_empty
            σʸ = σʸ[data_order];
        end

        if isempty(MNEMONIC) == false
            MNEMONIC = MNEMONIC[data_order];
        end
    end

    if ~is_σʸ_empty
        data = data./σʸ;
    end

    return data, MNEMONIC, quarterly_position, σʸ;
end
