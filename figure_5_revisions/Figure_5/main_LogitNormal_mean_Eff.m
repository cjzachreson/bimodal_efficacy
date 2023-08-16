% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 


% examining efficacy over logit-normal parameter space. 

clear all
close all 

n = 100000; % n = 10^5 samples from each dist. 


sig_1 = 0.05;
d_sig = 0.05;
sig_f = 6;

sig_vals = sig_1:d_sig:sig_f;

    mu_1 = -3;
    d_mu = 1;
    mu_f = 0;

mu_vals = mu_1:d_mu:mu_f;

n_sig_vals = numel(sig_vals);
n_mu_vals = numel(mu_vals);

Eff_mean = zeros(n_sig_vals, n_mu_vals);


sig_i = 0;

for sig = sig_vals   
    
    sig_i = sig_i + 1;
    mu_i = 0;
    
    for mu = mu_vals
        
        mu_i = mu_i + 1;

        x_dist = makedist('Normal', 'mu', mu, 'sigma', sig);
        
        x_sample = random(x_dist, [n ,1]);
        
        logit_normal_sample = standard_logistic(x_sample);
        
        Eff_mean(sig_i, mu_i) =  mean(logit_normal_sample);

    end
    
end

figure(2)
imagesc(Eff_mean)


dlmwrite('Eff_mean_mu_vs_sig_R1_s6.csv', Eff_mean);


function y = standard_logistic(x)

    y = 1 ./ (1 + exp(-x));

end
