




%% decaying control parameter (i.e., neut concentration and B-cell pop.)
w1 = 0.5;
w2 = 1 - w1;

v0 = 1

v0_1 = w1 * v0;  

v0_2 = w2 * v0;

tau_1 = (1/0.0085);
tau_2 = (1/0.00085);

t = 0:1:(2 * 365)

[v_t, v_1, v_2] = double_exponential(v0_1, v0_2, tau_1, tau_2, t);

figure(1)
semilogy(t, v_t)

figure(2)
plot([t', t'], [v_1', v_2'])

%% sigmoid 
% maps from log(neuts) to efficacy.

c50 = -0.5; % transition point for log-transformed control parameter v_t 
k = 5; % steepness of transition

y = sigmoid(c50, k, log(v_t)) ;

figure(3)
plot(log(v_t), y)

figure(4)
plot(t, y)