clear all 
close all

n = 10000;


%% residents

mu_residents = -1.47;
sigma_residents = 1.34;

neut_dist_residents = makedist('Lognormal', mu_residents, sigma_residents);

eff_aq_residents = zeros(n, 1); 
eff_t_residents = zeros(n, 1); 
eff_s_residents = zeros(n, 1); 
eff_h_residents = zeros(n, 1); 
eff_d_residents = zeros(n, 1); 

neuts_residents = random(neut_dist_residents, [n, 1]); 


for i = 1:n
    [eff_aq_residents(i, 1),...
    eff_t_residents(i, 1),...
    eff_s_residents(i, 1),...
    eff_h_residents(i, 1),...
    eff_d_residents(i, 1)] = compute_vaccine_efficacy_Omicron_from_neuts(neuts_residents(i));
end

%% general staff
mu_staff_g = -1.79;
sigma_staff_g = 1.43;

neut_dist_staff_g = makedist('Lognormal', mu_staff_g, sigma_staff_g);

eff_aq_staff_g = zeros(n, 1); 
eff_t_staff_g = zeros(n, 1); 
eff_s_staff_g = zeros(n, 1); 
eff_h_staff_g = zeros(n, 1); 
eff_d_staff_g = zeros(n, 1); 

neuts_staff_g = random(neut_dist_staff_g, [n, 1]); 


for i = 1:n
    [eff_aq_staff_g(i, 1),...
    eff_t_staff_g(i, 1),...
    eff_s_staff_g(i, 1),...
    eff_h_staff_g(i, 1),...
    eff_d_staff_g(i, 1)] = compute_vaccine_efficacy_Omicron_from_neuts(neuts_staff_g(i));
end

%% medical staff
mu_staff_m = -1.64;
sigma_staff_m = 1.51;

neut_dist_staff_m = makedist('Lognormal', mu_staff_m, sigma_staff_m);

eff_aq_staff_m = zeros(n, 1); 
eff_t_staff_m = zeros(n, 1); 
eff_s_staff_m = zeros(n, 1); 
eff_h_staff_m = zeros(n, 1); 
eff_d_staff_m = zeros(n, 1); 

neuts_staff_m = random(neut_dist_staff_m, [n, 1]); 


for i = 1:n
    [eff_aq_staff_m(i, 1),...
    eff_t_staff_m(i, 1),...
    eff_s_staff_m(i, 1),...
    eff_h_staff_m(i, 1),...
    eff_d_staff_m(i, 1)] = compute_vaccine_efficacy_Omicron_from_neuts(neuts_staff_m(i));
end














%pd_0 = makedist('lognormal', 'mu', mu_0, 'sigma', sig_0);

% figure(1)
% h = histogram(log_s_0, 'normalization', 'pdf');
% hold on
% plot(x2, y0_2)
% 
% figure(2)
% histogram(eff_0, 'normalization', 'pdf')
% hold on


% pd_normal = fitdist(log_s_1, 'Normal');
% x2 = log(x1);
% y1_2 = pdf(pd_normal, x2);







