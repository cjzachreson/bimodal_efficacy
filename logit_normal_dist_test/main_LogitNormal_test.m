% testing parameterisation of logit-normal distribution 

n = 1000000;

dx = 0.001;
x = dx-(dx/2):dx:1-(dx/2);

%% standard parameterisation: general normal -> standard logistic
mu = -1;
sig = 2;

y = logit_normal_pdf(mu, sig, x);

custom_normal = makedist('Normal', 'mu', mu, 'sigma', sig);

y0_samples = random(custom_normal, [n,1]);
y1_standard_logistic = standard_logistic(y0_samples);

% figure(1)
% h = histogram(y1_standard_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
% hold on
% plot(x, y, 'r')
% hold off

figure(2)
% plot(x, h.Values)
% hold on
plot(x, y, 'r--')
%hold off

%% alternate parameterisation
%standard normal -> general logistic (reparameterisation). 
standard_normal = makedist('Normal', 'mu', 0, 'sigma', 1);

y0_standard_normal = random(standard_normal, [n, 1]);

L = 1; 
k = sig;
x0 = -(mu/sig);
y1_general_logistic = general_logistic(y0_standard_normal, L, k, x0);
%y = histcounts(y1_general_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
% figure(2)
% hold on 
% %plot(x, h1.Values)
% plot(x, y, 'k-')
%hold off

%% relaxing both the underlying dist. and the logistic mapping: 
mu_1 = -5;
sig_1 = 0.5;
general_normal = makedist('Normal', 'mu', mu_1, 'sigma', sig_1); %call this N(mu_1, sig_1)

y2_general_normal = random(general_normal, [n, 1]);



L = 1; 

% re-parameterisation

sig_scaler = 1/sig_1 ; %we need to multiply N(mu_1, sig_1) by this to get the variance of a standard normal

y2_test = y2_general_normal * sig_scaler;

mu_offset = mu_1 / sig_1; %we need to subtract this from the stretched dist. with sig = 1 to get to standard normal.

y2_test = y2_test - mu_offset;

figure(5)
histogram(y2_test, 'normalization', 'pdf')

% the formulae below were re-arranged until the transform matched, followed
% by simplification. 
k =  sig * sig_scaler
x0 = -(mu/k) + mu_1


y2_general_logistic = general_logistic(y2_general_normal, L, k, x0);

% h1 = histogram(y1_general_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
y = histcounts(y2_general_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
figure(2)
hold on 
plot(x, y, 'k-')
hold off




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