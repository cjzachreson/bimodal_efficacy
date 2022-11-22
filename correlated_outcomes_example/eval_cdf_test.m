% eval_cdf_test 

clear all
close all

n = 10

w = rand([1, n]);
%w = sort(w);

w = w./sum(w);

x = 1:n;

cdf_w = cumsum(w) ./ sum(w); 

k = 1000000

out = NaN(k, 1);

for i = 1:k
    
    out(i) = eval_cdf(cdf_w, x);
    
end

w_out = NaN(size(w));

for j = 1:numel(x)
    
    w_out(j) = (sum(out == x(j))/k);
    
end

w'

w_out'

scatter(w, w_out)