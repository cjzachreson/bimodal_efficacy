function Eff = logit_normal_mean_quasi_MC(sig, mu, K)
 
x_dist = makedist('Normal', 'mu', mu, 'sigma', sig);

p_vals = [1:K-1] ./ K;

x_sample = icdf(x_dist, p_vals);

logit_normal_sample = standard_logistic(x_sample);

Eff =  sum(logit_normal_sample) ./ (K-1);


end

function y = standard_logistic(x)

    y = 1 ./ (1 + exp(-x));

end