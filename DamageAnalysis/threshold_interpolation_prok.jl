using DataFrames
using JLD2
using Pipe
using GLM

ftax = "Bacteria"; pip = "amaw"; frank="species";
ftax = "Archaea"; pip = "amaw"; frank="species";

before = load("./DmgMixtureModel/output/$(pip)_$(ftax)_$(frank)_filter_thres_before.csv") |> DataFrame # this are manually set thresholds
before = load("./DmgMixtureModel/output/$(pip)_$(ftax)_$(frank)_filter_thres_before_ANI95_strict.csv") |> DataFrame # this are manually set thresholds

reg = lm(@formula(lower ~ yrBP), before)

test = filter(:lower => l -> ismissing(l), before)
test.lower = predict(reg, DataFrame(yrBP=test.yrBP))

after = vcat(dropmissing(before), test)
sort!(after, :yrBP)

save("./DmgMixtureModel/output/$(pip)_$(ftax)_$(frank)_filter_thres_after.csv", after)
save("./DmgMixtureModel/output/$(pip)_$(ftax)_$(frank)_filter_thres_after_ANI95_strict.csv", after)
