% examining bimodal parameter space. 

clear all
close all 

dx = 0.001;
x_logit_normal = (0:dx:1)';


% %% relaxing both the underlying dist. and the logistic mapping: 
% mu_1 = -1;
% sig_1 = 2;
% 
% x_normal = ((mu_1 - 4*sig_1):dx:(mu_1 + 4*sig_1))';
% 
% general_normal = makedist('Normal', 'mu', mu_1, 'sigma', sig_1); %call this N(mu_1, sig_1)
% pdf_normal = pdf(general_normal, x_normal);
% 
% figure(1)
% plot(x_normal, pdf_normal, 'k-')
% 
% L = 1
% k =  1
% x0 = 0
% y2_general_logistic = general_logistic(x_normal, L, k, x0);
% 
% figure(1)
% hold on
% plot(x_normal, y2_general_logistic, 'r-')
% hold off


% sig = k * sig_1;
% mu = k * (mu_1 - x0);

sig_1 = 0.5;
d_sig = 0.1;
sig_f = 4;

sig_vals = sig_1:d_sig:sig_f;

mu_1 = -5;
d_mu = 0.1;
mu_f = 5;

mu_vals = mu_1:d_mu:mu_f;

n_sig_vals = numel(sig_vals);
n_mu_vals = numel(mu_vals);

test = zeros(n_sig_vals, n_mu_vals);

sig_i = 0;

for sig = sig_vals   
    
    sig_i = sig_i + 1;
    mu_i = 0;
    
    for mu = mu_vals
        
        mu_i = mu_i + 1;

        pdf_logit_normal = logit_normal_pdf(mu, sig, x_logit_normal);
        
        pdf_logit_normal(1) = 0;
        pdf_logit_normal(end) = 0;
        
        d_pdf = diff(pdf_logit_normal);
        
        sign_diff = sign(d_pdf);
        
        d_sign_diff = diff(sign_diff);
        
        d_sign_diff = d_sign_diff(~isnan(d_sign_diff));
        
        n_modes = sum(d_sign_diff < 0);
        
        basins = d_sign_diff > 0;
       
        if n_modes == 2
            basin_loc = find(basins);
            cum_dens_0_to_b = nansum(pdf_logit_normal(1:basin_loc));
            cum_dens_b_to_end = nansum(pdf_logit_normal(basin_loc:end));
            
            ratio = min([cum_dens_0_to_b, cum_dens_b_to_end]) / max([cum_dens_0_to_b, cum_dens_b_to_end]) ;
            
            test(sig_i, mu_i) = ratio;
            
        end
        
%         if  test(sig_i, mu_i) > 0
%             figure(sig_i)
%             plot(x_logit_normal, pdf_logit_normal)
%             hold on
%         end
            
            
            
    
    end
    
end

figure(2)
imagesc(test)



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