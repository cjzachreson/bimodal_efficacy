%simulate decay rates that vary for each individual.

% base parameters from
%khoury et al:
%https://www.nature.com/articles/s41591-021-01377-8#Sec20
%supp:
%https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01377-8/MediaObjects/41591_2021_1377_MOESM1_ESM.pdf

%pooled standard deviation sd_neuts is from Cromer et al., 2021
% https://doi.org/10.1016/S2666-5247(21)00267-6 table S3

clear all
close all


    
    
    lambda_sig_label = ['const_lambda'];   
    
    % define neut distribution
    
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
    t1 = 0;
    
    
    % neutralisation level:
    
    %moderna
    %mu_neuts = log10(654/94)
    
    %sinovac:
    %mu_neuts = log10(28/164)
    
    %pfizer
    %mu_neuts = log10(223/94)
    
    %AstraZeneca
    %mu_neuts = log10(32/59)
    
    %Johnson
    %mu_neuts = log10(246/522)
    
    %convalescent
    mu_neuts = 0;
    neut_label = 'conv';
    
    
    %from Cromer et al. (need link)
    sd_neuts = 0.465;
    
    
    %change of base.
    k = k / log(10);
    c50 = base10_to_base_e(c50);
    mu_neuts = base10_to_base_e(mu_neuts);
    sd_neuts = base10_to_base_e(sd_neuts);
    
    
    
    n = 100000;
    neut_dist = makedist('normal', 'mu', mu_neuts, 'sigma', sd_neuts);
    
    neuts_i = random(neut_dist, [n, 1]);

   
    %output for fancy plotting:
    
    output_LN_pdf = {};
    output_varnames = {};
    
    %plot pdf (Efficacy) at t0
    eff_t0 = general_logistic(neuts_i, L, k, c50);
    
    
    
    x_neut_dist_plot = [mu_neuts - 4 * sd_neuts: sd_neuts/100: mu_neuts + 4 * sd_neuts]';
    
    neut_dist_pdf_t0 = pdf(neut_dist, x_neut_dist_plot);
    
    x_logistic_mapping  = x_neut_dist_plot;
    
    Logistic_mapping = general_logistic(x_logistic_mapping, L, k, c50);
    
    % record neut dist at time t0
    output_table_neut_dist = table(x_neut_dist_plot, neut_dist_pdf_t0, x_logistic_mapping, Logistic_mapping);
    
    
    
    
    dx = 0.001;
    x_logit_normal = (0:dx:1)';
    
    x_edges = [0; x_logit_normal(2:end)-dx/2; 1+dx/2];
    
    sig_t0 = k * sd_neuts;
    mu_t0 = k * (mu_neuts - c50);
    pdf_logit_normal_t0 = logit_normal_pdf(mu_t0, sig_t0, x_logit_normal);
    
    figure(2)
    plot(x_logit_normal, pdf_logit_normal_t0, 'k--')
    hold on
    
    h_t0 = histcounts(eff_t0, 'BinEdges', x_edges, 'normalization', 'pdf');
    figure(2)
    plot(x_logit_normal, h_t0, 'k')
    hold on
    
    
    t0_varname = ['dt_0wk_' neut_label];
    
    output_varnames{1} = 'Eff_vals_pdf';
    output_LN_pdf{1} = x_logit_normal;
    
    output_finite_pdf{1} = x_logit_normal;
    
    output_varnames{2} = t0_varname;
    output_LN_pdf{2} = pdf_logit_normal_t0;
    output_finite_pdf{2} = h_t0';
    
    % plot for different periods of decay:
    
    %decay:
    
    t2_vals = [10, 30, 50]%5:5:50 %decay time in weeks
    
    n_t = numel(t2_vals)
    
    for t_i = 1:n_t
        
        t2 = t2_vals(t_i) * 7 % decay period in days
        
        cR = 1 - (1 - t_i/n_t);
        cG = 0;
        cB = 0;
        c = [cR, cG, cB];
        
        
        
        neuts_t = neuts_i - lambda*(t2 - t1);
        mu_neuts_t = mu_neuts - lambda*(t2 - t1);
        
        x_neut_dist_plot_t = [mu_neuts_t - 4 * sd_neuts: sd_neuts/100: mu_neuts_t + 4 * sd_neuts]';
        neut_dist_t = makedist('normal', 'mu', mu_neuts_t, 'sigma', sd_neuts);
        neut_dist_pdf_t = pdf(neut_dist_t, x_neut_dist_plot_t);
        
        varname_x = ['x_neut_dist_plot_t' num2str(t2_vals(t_i)) 'wk'];
        varname_pdf = ['neut_dist_pdf_t' num2str(t2_vals(t_i)) 'wk'];
        
        % record neut dist at time t2
        output_table_neut_dist.(varname_x) = x_neut_dist_plot_t;
        output_table_neut_dist.(varname_pdf) = neut_dist_pdf_t;
        
        
        eff_t2 = general_logistic(neuts_t, L, k, c50);
        
        % std(neuts_i)
        % std(neuts_t)
        
        sig_t2 = k * sd_neuts;
        mu_t2 = k * (mu_neuts_t - c50);
        pdf_logit_normal_t2 = logit_normal_pdf(mu_t2, sig_t2, x_logit_normal);
        
        figure(2)
        plot(x_logit_normal, pdf_logit_normal_t2, 'color', c, 'LineStyle', '--')
        
        
        h_t2 = histcounts(eff_t2, 'BinEdges', x_edges, 'normalization', 'pdf');
        figure(2)
        plot(x_logit_normal, h_t2,  'color', c, 'LineStyle', '-')
        hold on
        
        
        t2_varname = ['dt_' num2str(t2_vals(t_i)) 'wk_' neut_label];
        
        output_varnames{t_i + 2} = t2_varname;
        output_LN_pdf{t_i + 2} = pdf_logit_normal_t2;
        output_finite_pdf{t_i + 2} = h_t2';
        
    end
    
    figure(2)
    hold off
    
%     figure(3)
%     hold off
    
    output_table_LN = array2table([output_LN_pdf{:}], 'VariableNames', output_varnames);
    output_table_finite = array2table([output_finite_pdf{:}], 'VariableNames', output_varnames);
    
    output_dirname = [pwd(), '/' neut_label '/' lambda_sig_label];
    output_LN_pdf_fname = [output_dirname, '/', 'pdf_Eff_vs_t_' neut_label '_' lambda_sig_label '.csv'];
    output_finite_pdf_fname = [output_dirname, '/', 'pdf_finite_Eff_vs_t_' neut_label '_' lambda_sig_label '.csv'];
    
    output_neut_dist_fname = [output_dirname, '/', 'pdf_neut_dist_and_logistic.csv'];
    
    if ~isfolder(output_dirname)
        mkdir(output_dirname)
    end
    
    writetable(output_table_LN, output_LN_pdf_fname);
    writetable(output_table_finite, output_finite_pdf_fname);
    
    writetable(output_table_neut_dist, output_neut_dist_fname);
    

function p_LN = logit_normal_pdf(mu, sig, x)

term_1 = 1./(sig * sqrt(2*pi));

term_2 = 1./(x.*(1-x));

logit_x = log(x./(1 - x));

term_3 = exp(- (logit_x - mu).^2 ./ (2*sig.^2));

p_LN = term_1 .* term_2 .* term_3;

end



function y = general_logistic(x, L, k, x0)

y = L ./ (1 + exp(-k * (x - x0)));

end

function n_base_e = base10_to_base_e(n_base_10)
n_base_e = n_base_10 / log10(exp(1));
end
