%correlated output example MAIN

clear all
close all


% define neut distribution: 

% t1 

% t10

% set some parameters: 

% neut distribution: 
tau = 0.0085
mu_0 = -0.8
sd_neuts = 1.0
%decay of mean: 
t1 = 0;
t2 = 10; % time in weeks

t2 = t2 * 7; %time in days. 
mu_neuts = mu_0 - tau * (t2 - t1);
log_neut_dist = makedist('Normal', 'mu', mu_neuts, 'sigma', sd_neuts);

% logistic transform parameters
k = 2.4;
L = 1;
c50_acquisition = -1;
c50_death = -2.7;

% approximate efficacy dist with finite samples: 

n = 100000;
log_neuts = random(log_neut_dist, [n, 1]);

eff_acquisition = general_logistic(log_neuts, L, k, c50_acquisition);
eff_death = general_logistic(log_neuts, L, k, c50_death);

% what is efficacy against death, given infection? 
% ratio of relative risks: 
rr_aq = 1 - eff_acquisition;
rr_death = 1 - eff_death; 
% compute risk of death, given infection. 
rr_death_given_infection = rr_death ./ rr_aq;
eff_death_given_infection = 1 - rr_death_given_infection;

CFR = 0.1; % base cfr p(death | infection, not vaccinated)

% evaluate deaths as bernoulli trials from p(death | infection, vaccinated) 
p_death_given_infection = CFR * (1 - eff_death_given_infection);
%verbose, but explicit. 

% sample cases from cdf of rr_aq
n_cases = 10000;
cases = NaN(n_cases, 1);
cases_eff_death_given_infection = NaN(n_cases, 1);
cases_eff_death = NaN(n_cases, 1);
cases_eff_acquisition = NaN(n_cases, 1);

% finite distribution function
cdf_rr_aq = cumsum(rr_aq) ./ sum(rr_aq);

for i = 1:n_cases
    j = eval_cdf(cdf_rr_aq, [1:n]);
    cases_eff_death_given_infection(i) = eff_death_given_infection(j);
    cases_eff_death(i) = eff_death(j);
    cases_eff_acquisition(i) = eff_acquisition(j);
end


%     j = randsample(n, n_cases, true, rr_aq.* 0.001);
%     cases_eff_death_given_infection = eff_death_given_infection(j);
%     cases_eff_death = eff_death(j);
%     cases_eff_acquisition = eff_acquisition(j);


% evaluate p_death_given_infection for sampled cases. 
r = rand([n_cases, 1]);
p_death = CFR * (1 - cases_eff_death_given_infection);

death = r < p_death; 

CFR_vac = sum(death) / n_cases; 

eff_obs = 1 - CFR_vac / CFR

eff_mean = mean(eff_death)



% figure(1)
% scatter(log_neuts, rr_death_given_infection, '.')
% hold on 
% scatter(log_neuts, rr_death, 'r.')
% hold off
% 
% figure(2)
% histogram(rr_death)
% hold on 
% histogram(rr_death_given_infection)
% histogram(rr_aq)
% hold off
% 
% figure(3)
% scatter(eff_acquisition, rr_death_given_infection, '.')




    



    

        
        

    
   
    







function y = general_logistic(x, L, k, x0)

    y = L ./ (1 + exp(-k * (x - x0)));

end


function [a,b] = Uni_ab_from_mean(mu)
r = min(mu, 1-mu);
a = mu - r;
b = mu + r;
end

function n_base_e = base10_to_base_e(n_base_10)
    n_base_e = n_base_10 / log10(exp(1));
end