%simulate a phase III trial, i.e., simulate individual infection outcomes for a large
%cohort of vaccinated individuals. Then see what the distribution of the
%mean efficacy is for differently sized subsamples. See if the logit-normal
%model produces different stats than the 'leaky' or 'all-or-nothing'
%models. 

% parameters from
%khoury et al: 
%https://www.nature.com/articles/s41591-021-01377-8#Sec20
%supp: 
%https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01377-8/MediaObjects/41591_2021_1377_MOESM1_ESM.pdf

clear all
close all


% define neut distribution:

% global offset does not match the curves
% this means an offset is needed, but is a function of mu(t)

% logistic transform parameters
%from khoury et al. 
k = 3 
c50 = log10(0.2)
L = 1;

% neut distribution:
% exponential decay (linear decay of log neuts.)
% khoury et al:

t_half = 108;
lambda = log(2)/t_half;
%decay: 
t1 = 0;
t2 = 10 * 7; % decay period in days


% neutralisation level: 

%moderna
%mu_neuts = log10(654/94)

%sinovac:
%mu_neuts = log10(28/164)

%pfizer
%mu_neuts = log10(223/94)

%AstraZeneca
mu_neuts = log10(32/59)

%Johnson
%mu_neuts = log10(246/522)

%convalescent
%mu_neuts = 0;

%from Cromer et al. 
sd_neuts = 0.465;

%change of base. 
k = k / log(10);
c50 = base10_to_base_e(c50);
mu_neuts = base10_to_base_e(mu_neuts);
sd_neuts = base10_to_base_e(sd_neuts);



% simulate the phase III trial: 
p_base = 0.1; % p(infected | unvaccinated)

n = 100000;

mu_neuts = mu_neuts - mu_neuts * lambda * (t2 - t1);
% approximate efficacy dist with finite samples:
log_neut_dist = makedist('Normal', 'mu', mu_neuts, 'sigma', sd_neuts);

log_neuts = random(log_neut_dist, [n, 1]);
eff_ind_LN = general_logistic(log_neuts, L, k, c50);

% this should give the value from Khoury et al. 
eff_mean = mean(eff_ind_LN);

%all or nothing: 
eff_ind_AoN = NaN([n, 1]);
n_all = round(n * eff_mean);
eff_ind_AoN(1:n_all) = 1.0;
eff_ind_AoN((n_all + 1):end) = 0.0;
  

% Leaky: 
eff_ind_Leaky = ones(n, 1) * eff_mean;

k = 1000; %number of bootstraps
m = 10; %number of samples in each bootstrap


bstrap_AoN = bootstrap_outcomes(k, m, eff_ind_AoN, p_base);

bstrap_Leaky = bootstrap_outcomes(k, m, eff_ind_Leaky, p_base);

bstrap_LN = bootstrap_outcomes(k, m, eff_ind_LN, p_base);

figure(1)
histogram(bstrap_AoN, 'numbins', 10)
figure(2)
histogram(bstrap_Leaky, 'numbins', 10)
figure(3)
histogram(bstrap_LN, 'numbins', 10)
   
%TODO: test this function. 
function eff_bstrap_dist = bootstrap_outcomes(k, m, eff_dist, p0)

% k is number of bootstrap samples to examine
% m is number of samples in each bootstrap
% eff_dist is the underlying 'ground truth' efficacy distribution 

    eff_bstrap_dist = NaN([k, 1]);

    for i = 1:k

        bs_i = randsample(eff_dist, m);
        
        p_true = p0 .* (1 - bs_i);
        
        n_infections = sum(rand([m, 1]) < p_true); 

        eff_bstrap_dist(i) = 1 - (n_infections/m) / p0;
        
    end

end

function y = general_logistic(x, L, k, x0)

y = L ./ (1 + exp(-k * (x - x0)));

end

function n_base_e = base10_to_base_e(n_base_10)
    n_base_e = n_base_10 / log10(exp(1));
end
