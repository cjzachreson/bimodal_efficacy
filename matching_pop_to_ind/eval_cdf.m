function i_out = eval_cdf(cdf_discrete, possible_vals)

% inverse transform sample of discrete distribution

x = rand();

%s = abs(cdf_discrete - x);

s = cdf_discrete - x;
s(s<0) = Inf;
[~, ind_s] = min(s);

i_out = possible_vals(ind_s);


