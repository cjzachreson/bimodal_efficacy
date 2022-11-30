% testing parameterisation of logit-normal distribution 

clear all

%n = 1000000;

dx = 0.001;
x_logit_normal = (0:dx:1)';


%% relaxing both the underlying dist. and the logistic mapping: 
%mu_1 =  0; % convalescent serum - Khoury et al. 

%AstraZeneca
mu_1 = -0.61

t1 = 0
t2 = 10 %waning period in weeks
t2 = t2 * 7; %waning period in days. 

t_half = 108;
lambda = log(2)/t_half;

%linear decay of log neuts: 
mu_1 = mu_1 - lambda * (t2 - t1)

sig_1 = 1.07; % Cromer et al. (changed to base e)

x_normal = ((mu_1 - 4*sig_1):dx:(mu_1 + 4*sig_1))';

general_normal = makedist('Normal', 'mu', mu_1, 'sigma', sig_1); %call this N(mu_1, sig_1)
pdf_normal = pdf(general_normal, x_normal);

figure(2)
plot(x_normal, pdf_normal, 'k-')

L = 1
k =  1.3 %Khoury et al (converted to base e)
x0 = -1.6 %Khoury et al. (converted to base e)

y2_general_logistic = general_logistic(x_normal, L, k, x0);

figure(2)
hold on
plot(x_normal, y2_general_logistic, 'r-')
hold off


sig = k * sig_1;
mu = k * (mu_1 - x0);
pdf_logit_normal = logit_normal_pdf(mu, sig, x_logit_normal);

figure(3)
plot(x_logit_normal, pdf_logit_normal)

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