#main_DecayFunctions_test.jl

using DataFrames, CSV, PlotlyJS

include("Test_module_main.jl")

using .Test_module

# neut dist. parameters (lognormal(mu, sig))
mu = 1.0
sig = 0.5

c50_in = -0.5
k_in = 5.0
delta_t = 300.0
n = 100

(neuts, eff) = generate_efficacy_distribution(mu, sig, c50_in, k_in, delta_t, n)

df = DataFrame(log_titers = log.(neuts), Efficacy = eff)
println("finished making dataframe")


a = plot(df, x=:log_titers, kind = "histogram")
b = plot(df, x=:Efficacy, kind = "histogram")

display(a)
display(b)

println("finished plotting")