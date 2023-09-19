% Supplementary file to accompany the manuscript by Zachreson et al., 2023
% "Individual variation in vaccine immune response can produce bimodal distributions of protection" 
% author: Cameron Zachreson

function y = standard_logistic(x)

    y = 1 ./ (1 + exp(-x));

end

