#double_exponential.jl

export double_exponential 
function double_exponential(v0_1::AbstractFloat,
                            v0_2::AbstractFloat,
                            tau_1::AbstractFloat,
                            tau_2::AbstractFloat,
                            delta_t::AbstractFloat)::Tuple{Float64, Float64, Float64}

    v_1 = v0_1 * exp( -(1/tau_1) * delta_t)

    v_2 = v0_2 * exp( -(1/tau_2) * delta_t) 
    
    v_t = v_1 + v_2

    return (v_t, v_1, v_2)

end
