#main_DecayFunctions_test.jl

#include("DecayFunctions.jl")

#__precompile__(false)

if !(pwd() in LOAD_PATH)
    push!(LOAD_PATH, pwd())
end

using DecayFunctions

# testing sigmoid mapping:
x = 1.0 
c50 = 0.5 
k = 5.0
y = sigmoid_1(c50, k, x)
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

delta_t = 100.0 

(y1, v1, v2) = double_exponential(v0_1, v0_2, tau_1, tau_2, delta_t)
println("double exponential test complete:")
println("y1 = $y1")