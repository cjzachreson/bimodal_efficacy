
% make a finite distribution of trajectories 
n = 100

% distribution of log neuts
mu = 1
sig = 0.5

v0_dist = makedist('lognormal', 'mu', mu, 'sigma', sig);

for i = 1:n 
    v0 = random(v0_dist);



%% decaying control parameter (i.e., neut concentration and B-cell pop.)
    w1 = 0.5;
    w2 = 1 - w1;



    v0_1 = w1 * v0;  

    v0_2 = w2 * v0;

    tau_1 = (1/0.0085);
    tau_2 = (1/0.00085);

    t = 0:1:(2 * 365)

    [v_t, v_1, v_2] = double_exponential(v0_1, v0_2, tau_1, tau_2, t);

    figure(1)
    hold on
    semilogy(t, v_t)

    figure(2)
    hold on
    plot([t', t'], [v_1', v_2'])

%% sigmoid 
% maps from log(neuts) to efficacy.

    c50 = -0.5; % transition point for log-transformed control parameter v_t 
    k = 5; % steepness of transition

    y = sigmoid(c50, k, log(v_t)) ;

%     figure(3)
%     hold on
%     plot(log(v_t), y)

    figure(4)
    hold on
    plot(t, y)

end