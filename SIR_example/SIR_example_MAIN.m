%SIR example MAIN

clear all
close all

model_types = {'logit-normal'}%{'constant', 'all-or-nothing', 'logit-normal'}

for m = 1:numel(model_types)
    
    % uniform population, frequency-dependent transmission.
    k_runs = 100; %number of stochastic runs for each SIR simulation.
    
    % define vaccine efficacy distribution
    %eff_dist_label = 'constant';
    %eff_dist_label = 'all-or-nothing';
    % eff_dist_label = 'logit-normal';
    
    eff_dist_label = model_types{m};
    
    sig_1 = 0.05;
    d_sig = 0.05;
    sig_f = 4;
    
    sig_vals = sig_1:d_sig:sig_f;
    
    mu_1 = -5;
    d_mu = 0.1;
    mu_f = 5;
    
    mu_vals = mu_1:d_mu:mu_f;
    
    % lookup table for the same mu and sig vals:
    Eff_mat_fname = ['Eff_mean_mu_vs_sig.csv'];
    mean_eff_mat = dlmread(Eff_mat_fname);
    
    p_outbreak = NaN(size(mean_eff_mat));
    final_size = NaN(size(mean_eff_mat));
    final_size_q90 = NaN(size(mean_eff_mat));
    
    sig_i = 0;
    
    mu_indices = 1:numel(mu_vals);
    
    parfor sig_i = 1:numel(sig_vals)
        
        sig = sig_vals(sig_i);
        
        for mu_i = mu_indices%1:numel(mu_vals)
            
            mu = mu_vals(mu_i);
            
            disp(['mu = ' num2str(mu), ' ; sig = ' num2str(sig)]);
            
            % pre-computed from 10^5 samples.
            average_efficacy = mean_eff_mat(sig_i, mu_i);
            
            final_size_sig_mu = zeros(k_runs, 1);
            
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
            %p_trans = 1 - exp(-1 * dt * beta/N);
            
            p_recover = 1 - exp(-1 * dt * gamma);
            
            for i = 1:k_runs
                %% set up initial conditions
                %initialise population
                S = ones(N, 1);
                I = zeros(N, 1);
                R = zeros(N, 1);
                
                %efficacy of vaccination
                % vaccination:
                
                
                Eff = sample_efficacy_dist(N, eff_dist_label, mu, sig, average_efficacy);
                
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
                
                final_size_sig_mu(i) = sum(R);
                
                % do some checks here...
                %             figure(1)
                %             plot(total_infected, 'Color', [1, 0, 0, 0.2])
                %             hold on
                %             figure(2)
                %             plot(total_recovered,  'Color', [1, 0, 0, 0.2])
                %             hold on
                %             figure(3)
                %             plot(total_susceptible,  'Color', [1, 0, 0, 0.2])
                %             hold on
                
                
            end
            
            final_size(sig_i, mu_i) = mean(final_size_sig_mu);
            
            final_size_q90(sig_i, mu_i) = quantile(final_size_sig_mu, 0.9);
            
            p_outbreak(sig_i, mu_i) = nnz(final_size_sig_mu > 1) / k_runs;
            
            
        end
        
    end
    
    flabel_final_size = ['final_size_k_' num2str(k_runs) '_' eff_dist_label '.csv' ];
    flabel_p_outbreak = ['p_outbreak_k_' num2str(k_runs) '_' eff_dist_label '.csv'];
    flabel_final_size_q90 = ['final_size_q90_k_' num2str(k_runs) '_' eff_dist_label '.csv'];
    
    dlmwrite(flabel_final_size, final_size);
    dlmwrite(flabel_p_outbreak, p_outbreak);
    dlmwrite(flabel_final_size_q90, final_size_q90);
    
end






%% vaccine efficacy
function Eff = sample_efficacy_dist(N, eff_dist_label, mu_neuts, sig_neuts, average_efficacy)

%mean vaccine efficacy
%flag defining distribution of vaccine efficacy.
% parameters of logit-normal efficacy distribution

% standard sigmoidal mapping (does not change)
L = 1;
k =  1;
x0 = 0;

%placeholder
Eff = NaN(N, 1);

% parameters of underlying neut dist. for logit-normal
if strcmp(eff_dist_label, 'constant')
    % everyone gets the average efficacy
    Eff = ones(N, 1) * average_efficacy;
    
elseif strcmp(eff_dist_label, 'all-or-nothing')
    
    n_all = round(N * average_efficacy);
    Eff(1:n_all) = 1.0;
    Eff((n_all + 1):end) = 0.0;
    
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