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
    if bootstrap == 1

        fit_ip_weight = cal_ipw(data, adjustments, second_degree_terms, treatment)

        data[!, target*"_confound"] = data[!, target] .* (1 ./ data[!, treatment*"_w"])
        data[data[!, treatment] .== 0, target*"_confound"] = data[data[!, treatment] .== 0, target] .* (1 ./ (1 .- data[data[!, treatment] .== 0, treatment*"_w"]))

        fit_effect = lm(Term(Symbol(target*"_confound")) ~ Term(Symbol(treatment)), data)

        effect = GLM.coef(fit_effect)[2]

        return fit_ip_weight, effect
    
    else

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

    pred_1 = predict(fit, block1)
    pred_2 = predict(fit, block2)

    effect = mean(pred_1) - mean(pred_2)

    return fit, effect

end

function estimate_standardization(data, adjustments, second_degree_terms, treatment, target, bootstrap)
    effects = zeros(bootstrap)
    if bootstrap == 1
        
        return standardization(data, adjustments, second_degree_terms, treatment, target)

    else

        for i in 1:bootstrap

            sample_rows = sample(1:nrow(data), nrow(data), replace=true)

            tmp_data = copy(data[sample_rows, :])

            fit, effect = standardization(tmp_data, adjustments, second_degree_terms, treatment, target)

            effects[i] = effect
        
        end

        return effects
    end
end

function estimate_doubly_robust(data, adjustments, second_degree_terms, treatment, target, bootstrap)
    
    effects = zeros(bootstrap)
    
    if bootstrap == 1
        fit_ipw = estimate_ipw(data, adjustments, second_degree_terms, treatment, target, 1)
        ws = 1 ./ data[!, treatment*"_w"]
        data[!, treatment*"_R"] = ws
        data[data[!, treatment] .== 0, treatment*"_R"] = ws[data[!, treatment] .== 0] .* (-1)

        adjustments = [adjustments; treatment*"_R"]

        return estimate_standardization(data, adjustments, second_degree_terms, treatment, target, 1)
    else

        for i in 1:bootstrap

            sample_rows = sample(1:nrow(data), nrow(data), replace=true)

            tmp_data = copy(data[sample_rows, :])
            tmp_adjustments = copy(adjustments)

            fit_ipw = estimate_ipw(tmp_data, tmp_adjustments, second_degree_terms, treatment, target, 1)
            ws = 1 ./ tmp_data[!, treatment*"_w"]
            tmp_data[!, treatment*"_R"] = ws
            tmp_data[tmp_data[!, treatment] .== 0, treatment*"_R"] = ws[tmp_data[!, treatment] .== 0] .* (-1)
    
            tmp_adjustments = [tmp_adjustments; treatment*"_R"]
    
            fit, effect = estimate_standardization(tmp_data, tmp_adjustments, second_degree_terms, treatment, target, 1)

            effects[i] = effect
        
        end

        return effects

    end
end

function CI(effects)
    Non_NA_vec = effects[.!isnan.(effects)]
    println("Original runs: "*string(length(effects)))
    println("Non NA runs: "*string(length(Non_NA_vec)))

    mean_value = mean(Non_NA_vec)
    CI = [percentile(Non_NA_vec, 2.5), percentile(Non_NA_vec, 97.5)]
    return mean_value, CI
end