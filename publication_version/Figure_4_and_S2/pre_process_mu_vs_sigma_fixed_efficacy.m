% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 
% author: Cameron Zachreson


% examining efficacy over logit-normal parameter space. 
% this script produces a vector of mu as a function of sigma, holding 
% efficacy constant (to a specified tolerance) 
function pre_process_mu_vs_sigma_fixed_efficacy(Efficacy_vals, output_dirname, overwrite_flag)

    disp(['pre-processing mu vs. sigma for fixed efficacy values...'])

    for eff_i = 1:numel(Efficacy_vals)

        Efficacy = Efficacy_vals(eff_i); %the distribution mean, held constant (must be between 0 and 1)

        Eff_label = strrep(num2str(Efficacy), '.', 'p');

        output_fname = ['Eff_' Eff_label '_const_mu_vs_sigma.csv'];
        
        if isfile([output_dirname, output_fname])
            if ~overwrite_flag
                disp(['mu vs. sigma is already computed for efficacy ' Eff_label])
                disp(['toggle overwrite_flag = true to replace existing data in directory: \n' output_dirname])
                continue
            else
                disp(['mu vs. sigma is already computed for efficacy ' Eff_label])
                disp(['overwriting existing data in directory: ' output_dirname])
            end
        end

        max_iter = 1000; %number of iterations to use for mu estimate

        K = 100000; % number of slices to use in discrete estimate

        tolerance = 10^(-7); % desired precision of estimated mean

        % systematically vary sigma, and find mu for each to hold Eff const. 
        % 
        sig_0 = 0.05;
        d_sig = 0.25;
        sig_1 = d_sig;
        sig_f = 20;
        sig_vals = [sig_0, sig_1:d_sig:sig_f]';

        mu_vals = NaN(size(sig_vals));
        Eff_vals = NaN(size(sig_vals));

        for i = 1:numel(sig_vals)   

            sig_i = sig_vals(i);

            disp(['sigma = ' num2str(sig_i)])

            [mu_i, mean_i] = fit_mu_quasi_MC(Efficacy, sig_i, max_iter, tolerance, K) ;

            mu_vals(i) = mu_i;
            Eff_vals(i) = mean_i;


        end

        figure(1)
        plot(sig_vals, mu_vals)
        hold on

        output_table = table(sig_vals, mu_vals, Eff_vals); 

        writetable(output_table, [output_dirname, output_fname]);

    end
    
    disp(['pre-processing mu vs. sigma for fixed efficacy values... complete'])

end





