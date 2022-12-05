% testing variance of distributions derived from the data listed in Khoury et al. 

%khoury et al: 
%https://www.nature.com/articles/s41591-021-01377-8#Sec20
%supp: 
%https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01377-8/MediaObjects/41591_2021_1377_MOESM1_ESM.pdf

%pfizer: 

vac_mean = 223
vac_sd = 0.38 %standard deviation of the base10 logarithm. 
vac_sd = vac_sd / log10(exp(1));

mu = log(vac_mean) - ((vac_sd^2)/2)

conv_mean = 94 %normalisation factor

n = 10000

vac_dist = makedist('lognormal', 'mu', mu, 'sigma', vac_sd);

vac_sample = random(vac_dist, [n, 1]);

histogram(vac_sample)

check = fitdist(vac_sample, 'lognormal')

mean_test = exp(check.mu + check.sigma^2/2)
std_test = check.sigma / log(10)

% now try normalising: 

vac_sample_norm = vac_sample ./ conv_mean

figure(2)
histogram(vac_sample_norm)

norm_dist_check = fitdist(vac_sample_norm, 'lognormal')

norm_std_test = norm_dist_check.sigma / log(10)