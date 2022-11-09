%double_exponential

function [v_t, v_1, v_2] = double_exponential(v0_1, v0_2, tau_1, tau_2, delta_t) 


v_1 = v0_1 .* exp( -1 * (1/tau_1) .* delta_t);

v_2 = v0_2 .* exp ( -1 * (1/tau_2) .* delta_t) ;

v_t = v_1 + v_2;

