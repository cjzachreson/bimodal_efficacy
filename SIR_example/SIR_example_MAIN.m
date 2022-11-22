%SIR example MAIN

clear all
close all

% define vaccine efficacy distribution
 eff_dist_label = 'constant'; %354
 %eff_dist_label = 'uniform';
% eff_dist_label = 'logit-normal'; %274


%waning_period_label = '0 weeks'; %
waning_period_label = '10 weeks';

eff_dist_label
waning_period_label


% uniform population, frequency-dependent transmission.
k_runs = 1000; %number of stochastic runs for SIR simulation. 

final_size = zeros(k_runs, 1);

%% define system parameters:

% reproductive ratio
R0 = 2;

% recovery rate
gamma = 0.2; %average infectious period of 5 days

%transmission rate per contact per day
beta = R0 * gamma;

% time step.
dt = 0.1;

% population size.
N = 1000;

%raw transmission probability per susceptible contact.
%note this is altered by vaccination.
p_trans = 1 - exp(-1 * dt * beta/N);

p_recover = 1 - exp(-1 * dt * gamma);



for i = 1:k_runs
    %% set up initial conditions
    %initialise population
    S = ones(N, 1);
    I = zeros(N, 1);
    R = zeros(N, 1);
    
    %efficacy of vaccination
    Eff = zeros(N, 1); %initialise to 0
    
    % vaccination: 
    Eff = sample_efficacy_dist(N, eff_dist_label, waning_period_label);
    
%     figure(1)
%     histogram(Eff)
    
    FoI_S = zeros(N, 1); %force of infection on each susceptible
    
    %infect an index case:
    i_index = randperm(N, 1);
    S(i_index) = 0;
    I(i_index) = 1;
    
    t = 0;
    day = 0;
    
    cases_per_day = [];
    total_recovered = [];
    total_infected = [];
    total_susceptible = [];
    
    new_cases_today = sum(I);
    
    while sum(I) > 0
        
        t = t + dt;
        % some timekeeping.
        day = floor(t);
        day_last = floor((t - dt));
        if day > day_last
            %day
            cases_per_day = [cases_per_day; new_cases_today];
            total_infected = [total_infected; sum(I)];
            total_recovered = [total_recovered; sum(R)];
            total_susceptible = [total_susceptible; sum(S)];
            
            new_cases_today = 0;
        end
        
        
        %aggregate force of infection:
        FoI = beta * dt * sum(I) / N;
        
        %compute reduced force of infection for each susceptible individual
        % based on vaccination status:
        FoI_S = FoI .* ((1 - Eff) .* S);
        
        p_infect = 1 - exp(-FoI_S);
        
        %noting this is probably a faster way in matlab because iteration is slow
        newly_infected = rand(N, 1) < p_infect;
        
        new_cases_today = new_cases_today + sum(newly_infected);
        
        %infect people
        I(newly_infected) = 1;
        S(newly_infected) = 0;
        
        % compute probability of recovery
        p_recovery = p_recover .* I;
        
        newly_recovered = rand(N, 1) < p_recovery;
        
        I(newly_recovered) = 0;
        R(newly_recovered) = 1;
        
        
    end
    
    final_size(i) = sum(R);
    
    % do some checks here...
    figure(1)
    plot(total_infected, 'Color', [1, 0, 0, 0.2])
    hold on
    figure(2)
    plot(total_recovered,  'Color', [1, 0, 0, 0.2])
    hold on
    figure(3)
    plot(total_susceptible,  'Color', [1, 0, 0, 0.2])
    hold on
    
    
end


figure(4)
histogram(final_size)

q90 = quantile(final_size, 0.9)

p_outbreak = nnz(final_size > 1) / k_runs













%% vaccine efficacy
function Eff = sample_efficacy_dist(N, eff_dist_label, waning_period_label)

%mean vaccine efficacy
%flag defining distribution of vaccine efficacy.
% parameters of logit-normal efficacy distribution

% sigmoidal mapping (does not change)
% parameters of logistic transform from neuts -> efficacy for logit-normal
L = 1;
k =  2.4;
x0 = -1;

% parameters of underlying neut dist. for logit-normal
mu_neuts_t10 = -1.4; %approx. 10 wks waning of PF3 (peak is -0.8)
mu_neuts_t0 = -0.8; %these were derived from waning model (converted to base e)
sig_neuts = 1;
%t0 (peak protection)
%note: average efficacy estimated for large n = 10^6
eff_mu_t0 = 0.56;
eff_mu_t10 = 0.37;

%placeholder
average_efficacy = 0;
mu_neuts = 0;

if strcmp(waning_period_label, '0 weeks')
    average_efficacy = eff_mu_t0;
    mu_neuts = mu_neuts_t0;
    
elseif(strcmp(waning_period_label, '10 weeks'))
    average_efficacy = eff_mu_t10;
    mu_neuts = mu_neuts_t10;
end


if strcmp(eff_dist_label, 'constant')
    % everyone gets the same efficacy
    Eff = ones(N, 1) * average_efficacy;
elseif strcmp(eff_dist_label, 'uniform')
    % uniformly-distributed for defined mean efficacy
    [a, b] = Uni_ab_from_mean(average_efficacy)
    eff_dist = makedist('Uniform', 'upper', b, 'lower', a);
    Eff = random(eff_dist, [N, 1]);
elseif strcmp(eff_dist_label, 'logit-normal')
    neut_dist = makedist('Normal', 'mu', mu_neuts, 'sigma', sig_neuts);
    neuts = random(neut_dist, [N, 1]);
    Eff = general_logistic(neuts, L, k, x0);
    % logit-normal distribution of efficacy.
end


%% functions defining efficacy distribution

    function y = general_logistic(x, L, k, x0)
    y = L ./ (1 + exp(-k * (x - x0)));
    end


    function [a,b] = Uni_ab_from_mean(mu)
    r = min(mu, 1-mu);
    a = mu - r;
    b = mu + r;
    end

end