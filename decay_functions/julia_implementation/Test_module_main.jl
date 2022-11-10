module Test_module

include("DecayFunctions.jl")

#__precompile__(false)

#if !(pwd() in LOAD_PATH)
#    push!(LOAD_PATH, pwd())
#end

#using .DecayFunctions
using .DecayFunctions 


using StatsBase
using Distributions 
using Random

export generate_efficacy_distribution
function generate_efficacy_distribution(sig::AbstractFloat,
                                        mu::AbstractFloat,
                                        c50_in::AbstractFloat, 
                                        k_in::AbstractFloat, 
                                        delta_t::AbstractFloat,
                                        n::Int)::Tuple{Vector{Float64}, Vector{Float64}} 
    # testing sigmoid mapping:
    x = 1.0
    c50 = c50_in
    k = k_in

    println("c50: $c50")

    y = sigmoid_1(c50, k, log(x))
    println("sigmoid test complete:")
    println("y = $y")

    # testing double-exponential mapping: 
    w1 = 0.5 #weight of decay component 1
    w2 = 1.0 - w1 #weight of decay componenet 2 

    v0 = 1.0

    v0_1 = v0 * w1
    v0_2 = v0 * w2
    tau_1 = 1.0 / 0.0085
    tau_2 = 1.0 / 0.00085

    #delta_t = 100.0 

    (y1, v1, v2) = double_exponential(v0_1, v0_2, tau_1, tau_2, delta_t)
    println("double exponential test complete:")
    println("y1 = $y1")



    # generating a heterogeneous population: 
    #n = 10
    neuts = zeros(n)
    efficacy = zeros(n)

    seed_neuts = 1
    rng_neuts = MersenneTwister(seed_neuts)

    #mu = 1.0
    #sig = 0.5

    neut_dist = LogNormal{Float64}(mu, sig)
    rand!(rng_neuts, neut_dist, neuts)

    #println(neuts)

    #TODO next: map each value in neuts through the sigmoidal function sigmoid_1()
    iterator = 0
    for neuts_i in neuts 
        iterator += 1
        #apply decay function: 
        (neuts_dt, v1, v2) = double_exponential(neuts_i * w1, neuts_i * w2, tau_1, tau_2, delta_t)
        neuts[iterator] = neuts_dt
        efficacy[iterator] = sigmoid_1(c50, k, log(neuts_dt))
    end

    println("finished computing efficacy and neut values")

    return (neuts, efficacy)

end

end 