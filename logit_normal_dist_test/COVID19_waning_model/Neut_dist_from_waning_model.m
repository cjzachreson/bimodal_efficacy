clear all 
close all

n = 10000;

t1 = 0

% time in weeks. 
t2 = 0

%t2 = 34

%t2 = 1

%t2 = 16

vac_type = 'PF3'

eff_aq = zeros(n, 1); 
eff_t = zeros(n, 1); 
eff_s = zeros(n, 1); 
eff_h = zeros(n, 1); 
eff_d = zeros(n, 1); 
neuts = zeros(n, 1); 


for i = 1:n
    [eff_aq(i, 1),...
    eff_t(i, 1),...
    eff_s(i, 1),...
    eff_h(i, 1),...
    eff_d(i, 1),...
    neuts(i, 1)] = compute_vaccine_efficacy_Omicron_v2(vac_type,t1, t2);
end



pd_lognormal = fitdist(neuts, 'Lognormal')














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







