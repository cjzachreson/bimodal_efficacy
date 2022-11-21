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
general_normal_t12 = makedist('Normal', 'mu', mu_neuts_t12, 'sigma', sig_neuts);

% take n smaples from normal and estimate mean and variance. 
%% for peak (day 0)
neut_samples_t0 = random(general_normal_t0, [n, 1]);
efficacy_samples_logit_normal_t0 = general_logistic(neut_samples_t0, L, k, x0);
mean_logit_normal_t0 = mean(efficacy_samples_logit_normal_t0)
var_logit_normal_t0 = var(efficacy_samples_logit_normal_t0)

%% for week 12: 
neut_samples_t12 = random(general_normal_t12, [n, 1]);
efficacy_samples_logit_normal_t12 = general_logistic(neut_samples_t12, L, k, x0);
mean_logit_normal_t12 = mean(efficacy_samples_logit_normal_t12)
var_logit_normal_t12 = var(efficacy_samples_logit_normal_t12)


% make some comparable single-mode distributions: 
% beta distribution is defined [0, 1], with flexible shape so should be appropriate
%% for t0 (peak)
mean_t0 = mean_logit_normal_t0;
var_t0 = var_logit_normal_t0;
[Beta_a_t0, Beta_b_t0] =  Beta_ab_from_mean_and_var(mean_t0, var_t0);
Beta_dist_t0 = makedist('Beta', 'a', Beta_a_t0, 'b', Beta_b_t0);
efficacy_samples_Beta_t0 = random(Beta_dist_t0, [n, 1]);

function [a, b] =  Beta_ab_from_mean_and_var(mu_1, v_1)

    v_t0 = (mu_1*(1 - mu_1) / v_1) - 1;
    a = mu_1 * v_t0;
    b = (1 - mu_1) * v_t0;

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

