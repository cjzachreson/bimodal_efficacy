function [eff_acquisition,...
          eff_transmission,... 
          eff_symptoms,... 
          eff_hospitalisation,... 
          eff_death,...   
          neuts_t] = compute_neuts_base_e_v3(vac_type, t_last_vac, t)

% NOTE: 2022 09 24, modified to include neuts in output. 

% map vaccine type, time since vaccination and (possibly) subsequent
% infection to NAT concentration, and efficacy at time t = now.

% model is from here: https://github.com/goldingn/neuts2efficacy/blob/master/methods.md
% parameters sent by Ruarai Tobin on 2022 07 20.

neuts_t = 0;

if ~strcmp(vac_type, 'naive')
    
    % convert to base e
    c50_hosp = base10_to_base_e(-1.206);
    c50_death = base10_to_base_e(-1.184);
    c50_acquisition = base10_to_base_e(-0.4717);
    c50_transmission = base10_to_base_e(0.01846);
    c50_symptoms = base10_to_base_e(-0.6349);
    
    omicron_log_neut_fold = base10_to_base_e(-0.8006);
    
    sd_log_neut_titres = base10_to_base_e(0.4647);
    
    % this is the one I'm not sure about 
    log_k = 1.707;
    k = exp(log_k);
    k =  k / log(10);
    
    neut_decay_rate = 0.0085;
    
    %NOTE: the 4th dose gives a stronger response than the 3rd dose
    %(log10(1.33))
    
    vac_levels =           {'naive';'IAI';'PF1';   'PF2';  'PF3';  'PF4';  'PF4IAI'};
    log10_mean_neut_vals = {-Inf;    0;   -0.2296; 0.1706; 0.4644; 0.4644 + log10(1.33); 0.4644 + log10(1.33); };
    vac_2_neut = containers.Map(vac_levels, log10_mean_neut_vals);
    %note that infection acquired immunity iterates an individual forward to the next.
    % so there are no additional steps needed to account for permutations via
    % this simplifying assumption.
    
    if ~any(strcmp(vac_levels, vac_type))
        disp(['WARNING: vaccine type not found. Specify as one of:'])
        vac_levels
        return
    else
        log10_neut = vac_2_neut(vac_type);
    end
    
    log_neut = base10_to_base_e(log10_neut);
    
    log_neuts_dist = makedist('normal', 'mu', log_neut, 'sigma', sd_log_neut_titres);
    
    log_neuts_0_i = random(log_neuts_dist);
    
    % for omicron
    % NOTE: the neut_fold factor should only be used if the person has NOT
    % been infected with Omicron before. Unclear precisely how to implement
    % multistrain immune evasion. 
    log_neuts_0_i = log_neuts_0_i + omicron_log_neut_fold;
    
    dt_vac = (t - t_last_vac) * 7; % convert weeks to days.
    
    neuts_0 = exp(log_neuts_0_i);
    
    neuts_t = neut_decay(neuts_0, dt_vac, neut_decay_rate);
    
    log_neuts_t = log(neuts_t);
    
    eff_acquisition = log_neuts_to_efficacy(log_neuts_t, c50_acquisition, k);
    
    eff_transmission = log_neuts_to_efficacy(log_neuts_t, c50_transmission, k);
    
    eff_symptoms = log_neuts_to_efficacy(log_neuts_t, c50_symptoms, k);
    
    eff_hospitalisation = log_neuts_to_efficacy(log_neuts_t, c50_hosp, k);
    
    eff_death = log_neuts_to_efficacy(log_neuts_t, c50_death, k);
    
    
end
end

function eff = log_neuts_to_efficacy(log_neuts_t, c50, k)
eff = 1 / (1 + exp(-k*(log_neuts_t - c50)));
end

function neuts_t = neut_decay(neut_0, t, decay_rate)
neuts_t = neut_0 * exp(-decay_rate * t);
end


function n_base_e = base10_to_base_e(n_base_10)

    n_base_e = n_base_10 / log10(exp(1));
    

end

% if strcmp(vac_type, 'PF1')
% log10_mean_neut_vac = -0.2296;
% elseif strcmp(vac_type, 'PF2')
%log10_mean_neut_vac = 0.1706;
% elseif strcmp(vac_type, 'PF3')
%log10_mean_neut_vac = 0.4644;
% elseif strcmp(vac_type, 'PF4')
% log10_mean_neut_vac = 0.4644;
% end