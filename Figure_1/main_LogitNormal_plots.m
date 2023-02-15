% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 

%This script reproduces the data used in Figure 1 of the above manuscript. 

% plotting logit-normal distributions:

%NOTE: this script implements the general parameterisation, but sets the logistic
%function to the standard logistic. 


% x values to evaluate pdf
dx = 0.001; % resolution for normal dist x values (for plotting)

x_logit_normal = (0:dx:1)'; % x vlaues for plotting the logit-normal dist.


% control variable dist. parameters:
mu_vals = [-1, 0, 1];

n = numel(mu_vals);

%output sturctures: 
outtable_Fig1a = table();
outtable_Fig1b = table(); 

for i = 1:n

    mu_1 = mu_vals(i);
    
    sig_1 = 2;
    
    x_normal = ((mu_1 - 4*sig_1):dx:(mu_1 + 4*sig_1))';
    outtable_Fig1a.x_normal = x_normal;
    
    general_normal = makedist('Normal', 'mu', mu_1, 'sigma', sig_1);
    pdf_normal = pdf(general_normal, x_normal);
    
    mu_label = ['mu_' num2str(mu_1)];
    N_varname = ['Normal_pdf_' mu_label];
    outtable_Fig1a.(N_varname) = pdf_normal;
    
    figure(1)
    yyaxis left
    plot(x_normal, pdf_normal, '-', 'color', [i/n, i/n/2, i/n/3])
    
    %logistic parameters:
    L = 1;
    k =  1;
    x0 = 0;
    y2_general_logistic = general_logistic(x_normal, L, k, x0);
    
    % plot of the normal dist. control variable
    % overlayed with logistic mapping.
    figure(1)
    yyaxis right
    plot(x_normal, y2_general_logistic, 'g-')
    hold on
    
    % this sets up the logit-normal distribution using the control variable
    % dist. and logistic mapping specified above.
    sig = k * sig_1;
    mu = k * (mu_1 - x0);
    pdf_logit_normal = logit_normal_pdf(mu, sig, x_logit_normal);
    
    figure(3)
    plot(x_logit_normal, pdf_logit_normal, 'color', [i/n, i/n/2, i/n/3])
    hold on
    
    outtable_Fig1b.x_logit_normal = x_logit_normal;
    
    LN_label = ['pdf_logit_normal_' mu_label];
    
    outtable_Fig1b.(LN_label) = pdf_logit_normal;
    
end

writetable(outtable_Fig1a, 'data_Figure_1a.csv');
writetable(outtable_Fig1b, 'data_Figure_1b.csv');


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