function [eff_acquisition,...
    eff_transmission,...
    eff_symptoms,...
    eff_hospitalisation,...
    eff_death,...
    neuts_t] = ...
    compute_vaccine_efficacy_Omicron_v2(vac_type, t_last_vac, t)

% NOTE: 2022 09 24, modified to include neuts in output. 

% map vaccine type, time since vaccination and (possibly) subsequent
% infection to NAT concentration, and efficacy at time t = now.

% model is from here: https://github.com/goldingn/neuts2efficacy/blob/master/methods.md
% parameters sent by Ruarai Tobin on 2022 07 20.

eff_acquisition = 0;
eff_transmission = 0;
eff_symptoms = 0;
eff_hospitalisation = 0;
eff_death = 0;
neuts_t = 0;

if ~strcmp(vac_type, 'naive')
    
    c50_hosp = -1.206;
    c50_death = -1.184;
    c50_acquisition = -0.4717;
    c50_transmission = 0.01846;
    c50_symptoms = -0.6349;
    
    omicron_log10_neut_fold = -0.8006;
    
    sd_log10_neut_titres = 0.4647;
    
    log_k = 1.707;
    k = exp(log_k);
    
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
    
    
    log10_neuts_dist = makedist('normal', 'mu', log10_neut, 'sigma', sd_log10_neut_titres);
    
    log10_neuts_0_i = random(log10_neuts_dist);
    
    % for omicron
    % NOTE: the neut_fold factor should only be used if the person has NOT
    % been infected with Omicron before. Unclear precisely how to implement
    % multistrain immune evasion. 
    log10_neuts_0_i = log10_neuts_0_i + omicron_log10_neut_fold;
    
    dt_vac = (t - t_last_vac) * 7; % convert weeks to days.
    
    neuts_0 = 10^(log10_neuts_0_i);
    
    neuts_t = neut_decay(neuts_0, dt_vac, neut_decay_rate);
    
    log10_neuts_t = log10(neuts_t);
    
    eff_acquisition = log10_neuts_to_efficacy(log10_neuts_t, c50_acquisition, k);
    
    eff_transmission = log10_neuts_to_efficacy(log10_neuts_t, c50_transmission, k);
    
    eff_symptoms = log10_neuts_to_efficacy(log10_neuts_t, c50_symptoms, k);
    
    eff_hospitalisation = log10_neuts_to_efficacy(log10_neuts_t, c50_hosp, k);
    
    eff_death = log10_neuts_to_efficacy(log10_neuts_t, c50_death, k);
    
    
end
end

function eff = log10_neuts_to_efficacy(log_neuts_t, c50, k)
eff = 1 / (1 + exp(-k*(log_neuts_t - c50)));
end

function neuts_t = neut_decay(neut_0, t, decay_rate)
neuts_t = neut_0 * exp(-decay_rate * t);
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
