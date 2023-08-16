function Eff = logit_normal_mean(sig, mu, n)

Eff = NaN(numel(sig), 1);
if n< 1
    n = 1;
end
    for s_i = 1:numel(sig)

            x_dist = makedist('Normal', 'mu', mu, 'sigma', sig(s_i));

            x_sample = random(x_dist, [n ,1]);

            logit_normal_sample = standard_logistic(x_sample);

            Eff(s_i) =  mean(logit_normal_sample);
    end

end

