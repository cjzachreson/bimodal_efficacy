% testing parameterisation of logit-normal distribution 

n = 1000000;

dx = 0.001;
x = dx-(dx/2):dx:1-(dx/2);

mu = 3;
sig = 3;

y = logit_normal_pdf(mu, sig, x);

custom_normal = makedist('Normal', 'mu', mu, 'sigma', sig);
%custom_normal_pdf = pdf(custom_normal, x);

y0_samples = random(custom_normal, [n,1]);
y1_standard_logistic = standard_logistic(y0_samples);

figure(1)
h = histogram(y1_standard_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
hold on
plot(x, y, 'r')
hold off

figure(2)
plot(x, h.Values)
hold on
plot(x, y, 'r--')
hold off


standard_normal = makedist('Normal', 'mu', 0, 'sigma', 1);

y0_standard_normal = random(standard_normal, [n, 1]);

L = 1; 

fit_mu = [];
fit_sig = [];

for k = sig%1:0.1:3

x0 = -(mu/sig) ;

y1_general_logistic = general_logistic(y0_standard_normal, L, k, x0);

% h1 = histogram(y1_general_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
y = histcounts(y1_general_logistic, 'BinEdges', 0:dx:1, 'Normalization', 'pdf');
figure(2)
hold on 
%plot(x, h1.Values)
plot(x, y, 'k-')
end

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