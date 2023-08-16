% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 


% examining efficacy over logit-normal parameter space. 
% this script produces a vector of mu as a function of sigma, holding 
% efficacy constant (to a specified tolerance) 

clear all
close all 



Efficacy = 0.5 %the distribution mean, held constant

Eff_label = strrep('.', 'p', num2str(Efficacy));

n = 100000; % n = 10^5 samples from each dist. 

max_iter = 1000%for the estimator of mu

tolerance = 0.0001

% systematically vary sigma, and find mu for each to hold Eff const. 
% 
sig_1 = 0.05;
d_sig = 0.05;
sig_f = 6;

sig_vals = sig_1:d_sig:sig_f;
mu_vals = NaN(size(sig_vals));
Eff_vals = NaN(size(sig_vals));

mu_1 = -1; %initial guess for mu. 


n_sig_vals = numel(sig_vals);



for i = 1:numel(sig_vals)   
    
    sig_i = sig_vals(i);
   
    disp(['sigma = ' num2str(sig_i)])
    
    [mu_i, eff_i] = fit_mu(Efficacy, sig_i, max_iter, tolerance, n) ;
    
    mu_vals(i) = mu_i;
    Eff_vals(i) = eff_i;

    
end

figure(1)
plot(sig_vals, mu_vals)

output_table = table(sig_vals, mu_vals, Eff_vals); 

writetable(output_table, ['Eff_' Eff_label '_const_mu_vs_sig_R1.csv']);





