% examining logit-normal parameter space. 

clear all
close all 

dx = 0.001;
x_logit_normal = (0:dx:1)';

% specifying the range and step size for sigma parameter scan
sig_1 = 0.05;
d_sig = 0.05;
sig_f = 4;

sig_vals = sig_1:d_sig:sig_f;

% specifying the range and step size for mu parameter scan 
mu_1 = -5;
d_mu = 0.1;
mu_f = 5;

mu_vals = mu_1:d_mu:mu_f;

n_sig_vals = numel(sig_vals);
n_mu_vals = numel(mu_vals);

mode_density_ratio = zeros(n_sig_vals, n_mu_vals);

mode_distance = zeros(n_sig_vals, n_mu_vals);

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
        
        maxima = d_sign_diff < 0; 
        
        n_modes = sum(maxima);
        
        minima = d_sign_diff > 0;
       
        if n_modes == 2
            % analysis of mode densities
            min_loc = find(minima);
            cum_dens_0_to_min = nansum(pdf_logit_normal(1:min_loc));
            cum_dens_min_to_end = nansum(pdf_logit_normal(min_loc:end));
            
            ratio = min([cum_dens_0_to_min, cum_dens_min_to_end]) / max([cum_dens_0_to_min, cum_dens_min_to_end]) ;
            
            mode_density_ratio(sig_i, mu_i) = ratio;
            
            %analysis of mode separation:
            max_loc = find(maxima);
            
            x_max = x_logit_normal(max_loc);
            
            mode_dist = x_max(2) - x_max(1);
            
            mode_distance(sig_i, mu_i) = mode_dist;
            
            
            
        end  
    
    end
    
end

figure(2)
imagesc(mode_density_ratio)

figure(3)
imagesc(mode_distance)

dlmwrite('mode_dist_mu_vs_sig.csv', mode_distance);
dlmwrite('mode_ratio_mu_vs_sig.csv', mode_density_ratio);



function p_LN = logit_normal_pdf(mu, sig, x)

    term_1 = 1./(sig * sqrt(2*pi));

    term_2 = 1./(x.*(1-x));

    logit_x = log(x./(1 - x));

    term_3 = exp(- (logit_x - mu).^2 ./ (2*sig.^2));

    p_LN = term_1 .* term_2 .* term_3;

end
