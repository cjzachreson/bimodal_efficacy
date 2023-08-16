function [mu_est] = fit_mu_descent(LN_mean, sig, n_iter, tolerance, n)

% TODO: this needs more testing - not robust. 

mu_est = NaN;

x_target = LN_mean;


eps = 0.1;
mu_0 = 0;
mu_0_eps = mu_0 + eps * (rand() - 0.5);
x_0 = logit_normal_mean(sig, mu_0, n);
x_0_eps = logit_normal_mean(sig, mu_0_eps, n);

dx_0_loc = x_0_eps - x_0;
dmu_0_loc = mu_0_eps - mu_0;

dmu_dx_0_loc = dmu_0_loc / dx_0_loc;

dx_0 = x_0 - x_target;

dmu_i = dx_0 * dmu_dx_0_loc;

mu_i = mu_0 + dmu_i;

x_i = logit_normal_mean(sig, mu_i, n);

dmu_dx_i_loc = dmu_dx_0_loc;

for i = 1:n_iter 
    
    %while abs(dmu_dx_i_loc) > 5
    
        mu_i_eps =  mu_i + eps * sign(rand() - 0.5);
        x_i_eps = logit_normal_mean(sig, mu_i_eps, n);

        dmu_i_loc = mu_i_eps - mu_i;
        dx_i_loc = x_i_eps - x_i;

        dmu_dx_i_loc = dmu_i_loc / dx_i_loc;
    %end
    
    dx_i = x_target - x_i;
    
    
    dmu_i = dx_i * dmu_dx_i_loc;
    
    %try to avoid divergence: 
    dmu_i = sign(dmu_i) * min(abs(dmu_i), 1);
    
    mu_i = mu_i + dmu_i;
    
    x_i = logit_normal_mean(sig, mu_i, n);
    
    dx_i = x_i - x_target;
    
    if abs(dx_i) < tolerance
        mu_est = mu_i;
        disp(['mu estimate converged to: ' num2str(mu_est) ])
        disp(['   corresponding Efficacy: ' num2str(x_i) ])
        break
    end
    
    


end

if isnan(mu_est)
    disp('*warning: mu estimator did not converge*')
end

if abs(x_i - x_target) > tolerance
    disp('*warning: mu estimator is allowing unacceptable errors*')
end
    
    
end



