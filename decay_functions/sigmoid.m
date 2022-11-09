%sigmoid

function [y] = sigmoid(c50, k, x) 


    y = 1 ./ (1 + exp(-k .* (x - c50)));

end

