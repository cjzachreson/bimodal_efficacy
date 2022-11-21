% testing parameterisation of logit-normal distribution 

clear all

n = 1000000;

dx = 0.001;
x_logit_normal = (0:dx:1)';


%% approximating logit-normal moments at week 0 (peak) and week 12. 

% parameters of underlying neut dist. for logit-normal
mu_neuts_t12 = -1.4; %approx. 10 wks waning of PF3 (peak is -0.8)
mu_neuts_t0 = -0.8; %these were derived from waning model (converted to base e)
sig_neuts = 1;

% parameters of logistic transform from neuts -> efficacy for logit-normal
L = 1;
k =  2.4;
x0 = -1;

general_normal_t0 = makedist('Normal', 'mu', mu_neuts_t0, 'sigma', sig_neuts);
general_normal_t10 = makedist('Normal', 'mu', mu_neuts_t12, 'sigma', sig_neuts);

% take n smaples from normal and estimate mean and variance. 
%% for peak (day 0)
neut_samples_t0 = random(general_normal_t0, [n, 1]);
efficacy_samples_logit_normal_t0 = general_logistic(neut_samples_t0, L, k, x0);
mean_logit_normal_t0 = mean(efficacy_samples_logit_normal_t0)
var_logit_normal_t0 = var(efficacy_samples_logit_normal_t0)

%% for week 10: 
neut_samples_t10 = random(general_normal_t10, [n, 1]);
efficacy_samples_logit_normal_t10 = general_logistic(neut_samples_t10, L, k, x0);
mean_logit_normal_t10 = mean(efficacy_samples_logit_normal_t10)
var_logit_normal_t10 = var(efficacy_samples_logit_normal_t10)


% make a comparable uniform distribution: 
% uniform distribution: 
%% for t0 (peak)
mean_t0 = mean_logit_normal_t0;

[a_t0, b_t0] =  Uni_ab_from_mean(mean_t0);
Uni_dist_t0 = makedist('uniform', 'lower', a_t0, 'upper', b_t0);
efficacy_samples_Uni_t0 = random(Uni_dist_t0, [n, 1]);

mean_Uni_t0 = mean(efficacy_samples_Uni_t0)
var_Uni_t0 = var(efficacy_samples_Uni_t0)
% maximum variance without being bi-modal. 

%%
mean_t10 = mean_logit_normal_t10;
[a_t10, b_t10] =  Uni_ab_from_mean(mean_t10);
Uni_dist_t10 = makedist('uniform', 'lower', a_t10, 'upper', b_t10);
efficacy_samples_Uni_t10 = random(Uni_dist_t10, [n, 1]);

mean_Uni_t10 = mean(efficacy_samples_Uni_t10)
var_Uni_t10 = var(efficacy_samples_Uni_t10)


function [a,b] = Uni_ab_from_mean(mu)

r = min(mu, 1-mu);

a = mu - r;
b = mu + r;

end





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

%% notes 
% sig = k * sig_1;
% mu = k * (mu_1 - x0);
% pdf_logit_normal = logit_normal_pdf(mu, sig, x_logit_normal);

