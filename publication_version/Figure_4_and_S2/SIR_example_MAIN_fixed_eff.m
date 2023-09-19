% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 
% author: Cameron Zachreson


%SIR example MAIN


clear all
close all

 model_types = {'constant','all_or_nothing', 'logit_normal'}
%model_types = {'all_or_nothing', 'logit_normal'}
%model_types = {'logit_normal'}
%WARNING: these model type labels link to the efficacy
%distribution function, should make this more robust (globals?) 

%NOTE: should output variance or IQR of final size over k samples. 
% adding 10th, 90th, 


%first, pre-compute mu for each value of sigma and efficacy in the 
% logit-normal model. 
Efficacy_vals = [0.5, 0.3, 0.1]

data_dirname_pre_processing = [pwd(), '\mu_vs_sigma_fixed_efficacy\'];
if ~isfolder(data_dirname_pre_processing)
    mkdir(data_dirname_pre_processing)
end
% toggle overwrite_flag = true to pre-compute new values of mu for sigma-efficacy combos
% already computed previously. Any parameter combos not already processed
% will be computed.
overwrite_flag = false;

% pre-computation runs for values of sigma:
% [0.05, 0.25:0.25:20], these can be changed within the function
pre_process_mu_vs_sigma_fixed_efficacy(Efficacy_vals,... 
                                       data_dirname_pre_processing,...
                                       overwrite_flag)

OB_threshold = 10; %used for selecting a subset of outbreak stats not subject to stochastic die-out. 

% set seed for stochastic simulations
seed = 1;

for m = 1:numel(model_types)
    
    eff_dist_label = model_types{m};
    
    output_dirname = [pwd(), '\SIR_results\', eff_dist_label '\' ];

    if ~isfolder(output_dirname)
        mkdir(output_dirname)
    end

    % reseeding rng before each set 
    rng(seed)
    
    % uniform population, frequency-dependent transmission.
    k_runs = 1000; %number of stochastic runs for each SIR simulation.
    
    % For revised version, trace mu and sigma for fixed efficacy

    % lookup table for the same mu and sig vals:
    
    for eff_i = 1:numel(Efficacy_vals)
    
        efficacy = Efficacy_vals(eff_i);
        
        eff_label = strrep(num2str(efficacy), '.', 'p');
        % note this file name specification must be matched to the one in
        % the preprocessing function that generates these. 
        Eff_vs_mu_sig_fname = ['Eff_' eff_label '_const_mu_vs_sigma.csv'];
        mean_eff_table = readtable([data_dirname_pre_processing, Eff_vs_mu_sig_fname]);

        output_fname = [output_dirname, 'summary_output_Eff_' eff_label '.csv'];
        

        % summary stats: 
        p_outbreak = NaN(size(mean_eff_table, 1), 1);
        final_size_mean = NaN(size(mean_eff_table, 1), 1);
        final_size_q90 = NaN(size(mean_eff_table, 1), 1);

        final_size_q25 = NaN(size(mean_eff_table, 1), 1);
        final_size_q50 = NaN(size(mean_eff_table, 1), 1);
        final_size_q75 = NaN(size(mean_eff_table, 1), 1);

        %summary stats with case threshold applied:
        p_outbreak_th = NaN(size(mean_eff_table, 1), 1);
        final_size_mean_th = NaN(size(mean_eff_table, 1), 1);
        final_size_q90_th = NaN(size(mean_eff_table, 1), 1);

        final_size_q25_th = NaN(size(mean_eff_table, 1), 1);
        final_size_q50_th = NaN(size(mean_eff_table, 1), 1);
        final_size_q75_th = NaN(size(mean_eff_table, 1), 1);


        sig_vals = mean_eff_table.sig_vals;
        mu_vals = mean_eff_table.mu_vals;
        eff_vals = mean_eff_table.Eff_vals;

        

        parfor sig_mu_i = 1:numel(sig_vals)

            sig = sig_vals(sig_mu_i);
            

                mu = mu_vals(sig_mu_i);

                disp(['mu = ' num2str(mu), ' ; sig = ' num2str(sig)]);

                % pre-computed from estimator
                average_efficacy = eff_vals(sig_mu_i)

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



                end

                final_size_mean(sig_mu_i) = mean(final_size_sig_mu);

                final_size_q90(sig_mu_i) = quantile(final_size_sig_mu, 0.9);
                final_size_q25(sig_mu_i) = quantile(final_size_sig_mu, 0.25);
                final_size_q50(sig_mu_i) = quantile(final_size_sig_mu, 0.5);
                final_size_q75(sig_mu_i) = quantile(final_size_sig_mu, 0.75);

                p_outbreak(sig_mu_i) = nnz(final_size_sig_mu > 1) / k_runs;

                final_size_sig_mu_th = final_size_sig_mu(final_size_sig_mu >= OB_threshold);

                %apply case threshold and compute summary stats: 
                p_outbreak_th(sig_mu_i) = nnz(final_size_sig_mu >= OB_threshold) / k_runs;

                final_size_mean_th(sig_mu_i) = mean(final_size_sig_mu_th);

                final_size_q90_th(sig_mu_i) = quantile(final_size_sig_mu_th, 0.9);
                final_size_q25_th(sig_mu_i) = quantile(final_size_sig_mu_th, 0.25);
                final_size_q50_th(sig_mu_i) = quantile(final_size_sig_mu_th, 0.5);
                final_size_q75_th(sig_mu_i) = quantile(final_size_sig_mu_th, 0.75);



            

        end

        %summary stats: 
        output_table = table(sig_vals, mu_vals, eff_vals, p_outbreak_th,...
                        final_size_mean_th, final_size_q90_th,final_size_q75_th,...
                        final_size_q50_th, final_size_q25_th, p_outbreak,...
                        final_size_mean, final_size_q90,final_size_q75,...
                        final_size_q50, final_size_q25);
                    
        writetable(output_table, output_fname)
                    
         
        
    end
    
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
    
elseif strcmp(eff_dist_label, 'all_or_nothing')
    
    n_all = round(N * average_efficacy);
    Eff(1:n_all) = 1.0;
    Eff((n_all + 1):end) = 0.0;
    
elseif strcmp(eff_dist_label, 'logit_normal')
    neut_dist = makedist('Normal', 'mu', mu_neuts, 'sigma', sig_neuts);
    neuts = random(neut_dist, [N, 1]);
    Eff = general_logistic(neuts, L, k, x0);
    % logit-normal distribution of efficacy.
end


%% functions defining efficacy distribution

    function y = general_logistic(x, L, k, x0)
        y = L ./ (1 + exp(-k * (x - x0)));
    end


end