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

    return effects

end

function standardization(data, adjustments, second_degree_terms, treatment, target)

    for term in second_degree_terms
        data[!, term*"^2"] = data[!, term] .^ 2
    end

    block1 = copy(data)
    block1[!, treatment] .= 1

    block2 = copy(data)
    block2[!, treatment] .= 0

    fit = lm(Term(Symbol(target)) ~ Term(Symbol(treatment)) + sum(Term.(Symbol.(adjustments))) + sum(Term.(Symbol.(second_degree_terms.*"^2"))) + sum([Term(Symbol(treatment))&Term(Symbol(var)) for var in adjustments]), data) 
    #fit = lm(Term(Symbol(target)) ~ Term(Symbol(treatment)) + sum(Term.(Symbol.(adjustments))) + sum(Term.(Symbol.(second_degree_terms*"^2"))), data) 

    pred_1 = predict(fit, block1)
    pred_2 = predict(fit, block2)

    effect = mean(pred_1) - mean(pred_2)

    return fit, effect

end

function estimate_standardization(data, adjustments, second_degree_terms, treatment, target, bootstrap)
    effects = zeros(bootstrap)

    for i in 1:bootstrap

        sample_rows = sample(1:nrow(data), nrow(data), replace=true)

        tmp_data = copy(data[sample_rows, :])

        fit, effect = standardization(tmp_data, adjustments, second_degree_terms, treatment, target)

        effects[i] = effect
    
    end

    return effects
end