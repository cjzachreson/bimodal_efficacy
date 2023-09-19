% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 
% author: Cameron Zachreson

% examining efficacy over logit-normal parameter space. 

clear all
close all 

n = 100000; % n = 10^5 samples from each dist. 


sig_1 = 0.05;
d_sig = 0.05;
sig_f = 4;

sig_vals = sig_1:d_sig:sig_f;

mu_1 = -5;
d_mu = 0.1;
mu_f = 5;

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

        disp(['sigma = ' num2str(sig) ', mu = ' num2str(mu)])

        Eff_mean(sig_i, mu_i) =  logit_normal_mean_quasi_MC(sig, mu, n);

    end
    
end

figure(2)
imagesc(Eff_mean)


dlmwrite('Eff_mean_mu_vs_sig_QMC.csv', Eff_mean);

