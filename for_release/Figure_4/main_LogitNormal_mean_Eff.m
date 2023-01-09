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

        x_dist = makedist('Normal', 'mu', mu, 'sigma', sig);
        
        x_sample = random(x_dist, [n ,1]);
        
        logit_normal_sample = standard_logistic(x_sample);
        
        Eff_mean(sig_i, mu_i) =  mean(logit_normal_sample);
        
        %do a ton of samples and compute the mean. 

    end
    
end

figure(2)
imagesc(Eff_mean)


dlmwrite('Eff_mean_mu_vs_sig.csv', Eff_mean);




function p_LN = logit_normal_pdf(mu, sig, x)

    term_1 = 1./(sig * sqrt(2*pi));

    term_2 = 1./(x.*(1-x));

    logit_x = log(x./(1 - x));

    term_3 = exp(- (logit_x - mu).^2 ./ (2*sig.^2));

    p_LN = term_1 .* term_2 .* term_3;

end

function y = standard_logistic(x)

    y = 1 ./ (1 + exp(-x));

end

function y = general_logistic(x, L, k, x0)

    y = L ./ (1 + exp(-k * (x - x0)));

end