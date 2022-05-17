# Function for IP weighting
function cal_ipw(data, adjustments, second_degree_terms, treatment)

    for term in second_degree_terms
        data[!, term*"^2"] = data[!, term] .^ 2
    end

    fit = glm(Term(Symbol(treatment)) ~ sum(Term.(Symbol.(adjustments))) + sum(Term.(Symbol.(second_degree_terms .* "^2"))), data, Binomial(), ProbitLink())
    pred = predict(fit, data)
    data[!, treatment*"_w"] = 1 .- pred
    data[data[!, treatment] .= 1, treatment*"_w"] = pred
    
    return fit
end