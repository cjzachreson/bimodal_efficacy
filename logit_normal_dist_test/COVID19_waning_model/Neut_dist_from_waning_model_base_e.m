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

eff_aq_b10 = zeros(n, 1); 
eff_t_b10 = zeros(n, 1); 
eff_s_b10 = zeros(n, 1); 
eff_h_b10 = zeros(n, 1); 
eff_d_b10 = zeros(n, 1); 

neuts_10 = zeros(n, 1); 


eff_aq_be = zeros(n, 1); 
eff_t_be = zeros(n, 1); 
eff_s_be = zeros(n, 1); 
eff_h_be = zeros(n, 1); 
eff_d_be = zeros(n, 1); 

neuts_e =  zeros(n, 1); 


for i = 1:n
    [eff_aq_b10(i, 1),...
    eff_t_b10(i, 1),...
    eff_s_b10(i, 1),...
    eff_h_b10(i, 1),...
    eff_d_b10(i, 1),...
    neuts_10(i, 1)] = compute_vaccine_efficacy_Omicron_v2(vac_type,t1, t2);

    [eff_aq_be(i, 1),...
     eff_t_be(i, 1),... 
     eff_s_be(i, 1),... 
     eff_h_be(i, 1),... 
     eff_d_be(i, 1),... 
     neuts_e(i, 1)] = compute_neuts_base_e_v3(vac_type, t1, t2);

end


% these must be the same (if my change-of-base alterations are correct)
pd_lognormal_e = fitdist(neuts_e, 'Lognormal')

pd_lognormal_10 = fitdist(neuts_10, 'Lognormal')

figure(1)
histogram(eff_aq_b10, 'BinEdges', 0:0.01:1)
hold on
histogram(eff_aq_be, 'BinEdges', 0:0.01:1)
hold off














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






