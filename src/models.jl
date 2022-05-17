# Function for IP weighting
function cal_ipw(data, adjustments, second_degree_terms, treatment)

    for term in second_degree_terms
        data[!, term*"^2"] = data[!, term] .^ 2
    end

    fit = glm(Term(Symbol(treatment)) ~ sum(Term.(Symbol.(adjustments))) + sum(Term.(Symbol.(second_degree_terms .* "^2"))), data, Binomial(), ProbitLink())
    pred = predict(fit, data)
    data[!, treatment*"_w"] = 1 .- pred
    data[data[!, treatment] .== 1, treatment*"_w"] = pred[data[!, treatment] .== 1]
    #println(data[!, :low_quality])
    
    return fit
end

function estimate_ipw(data, adjustments, second_degree_terms, treatment, target, bootstrap)

    effects = zeros(bootstrap)

    for i in 1:bootstrap

        sample_rows = sample(1:nrow(data), nrow(data), replace=true)

        tmp_data = copy(data[sample_rows, :])

        fit_ip_weight = cal_ipw(tmp_data, adjustments, second_degree_terms, treatment)

        tmp_data[!, target*"_confound"] = tmp_data[!, target] .* (1 ./ tmp_data[!, treatment*"_w"])
        tmp_data[tmp_data[!, treatment] .== 0, target*"_confound"] = tmp_data[tmp_data[!, treatment] .== 0, target] .* (1 ./ (1 .- tmp_data[tmp_data[!, treatment] .== 0, treatment*"_w"]))

        fit_effect = lm(Term(Symbol(target*"_confound")) ~ Term(Symbol(treatment)), tmp_data)

        effect = GLM.coef(fit_effect)[2]

        effects[i] = effect
    
    end

    return effects, mean(effects), [percentile(effects, 2.5), percentile(effects, 97.5)]

end