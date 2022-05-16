# Function for IP weighting
function ip_weighting(data, adjustments, treatment, target)
    fit = glm(Term(Symbol(target)) ~ sum(Term.(Symbol.(adjustments))), data, binomial(), ProbitLink())

end