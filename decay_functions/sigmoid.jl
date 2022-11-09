
function sigmoid_1(c50::AbstractFloat, k::AbstractFloat, x::AbstractFloat)::Float64

    y =  1 / (1 + exp(-k * (x - c50)))

    return y

end

export sigmoid_1 #test if this can be written in "sigmoid.jl" instead. 