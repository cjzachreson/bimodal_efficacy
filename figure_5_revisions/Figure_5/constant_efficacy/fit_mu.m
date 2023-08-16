function [mu_est, Eff_out] = fit_mu(LN_mean, sig, n_iter, tolerance, n)

mu_est = NaN;
Eff_out = NaN;

target = LN_mean;

mu_i = rand(); 

x_i = logit_normal_mean(sig, mu_i, n);
d0 = target - x_i;
d_last = d0;
sz_i = 0.1;

for i = 1:n_iter 

    step_i = sign(d_last) * sz_i;

    mu_i =  mu_i + step_i;

    x_i = logit_normal_mean(sig, mu_i, n);

    d_i = target - x_i;

    if abs(d_i) < tolerance
        mu_est = mu_i;
        Eff_out = x_i;
        disp(['mu estimate converged to: ' num2str(mu_est)])
        break
    end

    if sign(d_i) ~= sign(d_last)
        sz_i = sz_i / 3;
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



