%correlated output example MAIN

%khoury et al: 
%https://www.nature.com/articles/s41591-021-01377-8#Sec20
%supp: 
%https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01377-8/MediaObjects/41591_2021_1377_MOESM1_ESM.pdf

clear all


% define neut distribution:

% global offset does not match the curves
% this means an offset is needed, but is a function of mu(t)

% logistic transform parameters
%from khoury et al. 
k = 3 
c50 = log10(0.2)
L = 1;

% neut distribution:
% don't care about decay for now... 
% khoury et al:
%tau = 0.0064

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

% try change of base. 
k = k / log(10);
c50 = base10_to_base_e(c50);
mu_neuts = base10_to_base_e(mu_neuts);
sd_neuts = base10_to_base_e(sd_neuts);


% approximate efficacy dist with finite samples:
log_neut_dist = makedist('Normal', 'mu', mu_neuts, 'sigma', sd_neuts);
n = 100000;
log_neuts = random(log_neut_dist, [n, 1]);
eff_acquisition_ind = general_logistic(log_neuts, L, k, c50);

% this should give the value from Khoury et al. 
eff_ind_mean = mean(eff_acquisition_ind)


%eff_pop_mean = general_logistic(mu_neuts, L, k_pop, c50_pop)











function y = general_logistic(x, L, k, x0)

y = L ./ (1 + exp(-k * (x - x0)));

end

function n_base_e = base10_to_base_e(n_base_10)
    n_base_e = n_base_10 / log10(exp(1));
end
