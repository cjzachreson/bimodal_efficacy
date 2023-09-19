% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 
% author: Cameron Zachreson


function [mu_est, Eff_out] = fit_mu_quasi_MC(LogitNormal_mean, sig, n_iter, tolerance, K)

% a bespoke search function, estimating a value of mu
% that corresponds to a desired Logit-Normal distribution mean
% for fixed sigma. 

mu_est = NaN;
Eff_out = NaN;

target = LogitNormal_mean;

mu_i = rand(); 

x_i = logit_normal_mean_quasi_MC(sig, mu_i, K);
d0 = target - x_i;
d_last = d0;
sz_i = 0.1;

for i = 1:n_iter 

    step_i = sign(d_last) * sz_i;

    mu_i =  mu_i + step_i;

    x_i = logit_normal_mean_quasi_MC(sig, mu_i, K);

    d_i = target - x_i;

    if abs(d_i) < tolerance
        mu_est = mu_i;
        Eff_out = x_i;
        disp(['mu estimate converged to: ' num2str(mu_est)])
        break
    end

    % step size gets smaller if the 
    % target is overshot. 
    if sign(d_i) ~= sign(d_last)
        sz_i = sz_i / 2;
    end

    d_last = d_i;

end

if isnan(mu_est)
    disp('*warning: mu estimator did not converge*')
end

if abs(x_i - target) > tolerance
    disp('*warning: mu estimator is allowing unacceptable errors*')
end
    
    
end



