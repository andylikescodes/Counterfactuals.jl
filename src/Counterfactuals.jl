module Counterfactuals

using Revise
using CSV
using DataFrames
using GLM
using StatsBase

# Write your package code here.
include("models.jl")

export cal_ipw, estimate_ipw

end
