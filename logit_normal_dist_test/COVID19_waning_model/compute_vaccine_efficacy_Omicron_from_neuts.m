function [eff_acquisition,...
    eff_transmission,...
    eff_symptoms,...
    eff_hospitalisation,...
    eff_death] = ...
    compute_vaccine_efficacy_Omicron_from_neuts(neuts)

% model is from here: https://github.com/goldingn/neuts2efficacy/blob/master/methods.md
% parameters sent by Ruarai Tobin on 2022 07 20.

eff_acquisition = 0;
eff_transmission = 0;
eff_symptoms = 0;
eff_hospitalisation = 0;
eff_death = 0;

if neuts > 0
    
    c50_hosp = -1.206;
    c50_death = -1.184;
    c50_acquisition = -0.4717;
    c50_transmission = 0.01846;
    c50_symptoms = -0.6349;
    
    log_k = 1.707;
    k = exp(log_k);
    
    log10_neuts = log10(neuts);
    
    eff_acquisition = log10_neuts_to_efficacy(log10_neuts, c50_acquisition, k);
    
    eff_transmission = log10_neuts_to_efficacy(log10_neuts, c50_transmission, k);
    
    eff_symptoms = log10_neuts_to_efficacy(log10_neuts, c50_symptoms, k);
    
    eff_hospitalisation = log10_neuts_to_efficacy(log10_neuts, c50_hosp, k);
    
    eff_death = log10_neuts_to_efficacy(log10_neuts, c50_death, k);
    
    
end
end

function eff = log10_neuts_to_efficacy(log_neuts, c50, k)
eff = 1 / (1 + exp(-k*(log_neuts - c50)));
end

