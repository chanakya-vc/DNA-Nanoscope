% Compute a weighted error between embedded and measured distances for the distance geometry problem
% Usage:
%   J = DGPcostFunction(theta, M, weights)
% theta is an 2*n by 1 column vector contains x coordinates for all points
% followed by y coordinates for all points
% M is an n by n lower triangular matrix containing pairwise distances
% between points.
% weights is an n by n lower triangular matrix containing weights for the
% pairwise distances

function J = DGPcostFunction(theta, M, weights)

% restructure theta
temp_n = size(theta,1);
theta = [theta(1:temp_n/2) theta(temp_n/2 + 1:end)];

% initialize some useful values
n = size(theta,1); % number of points

% pairwise distances from parameters
D = zeros(size(M)); % initialize
for i=1:n
    for j=1:i-1
        D(i,j) = sqrt((theta(j,2) - theta(i,2))^2 + (theta(j,1) - theta(i,1))^2);
    end
end

pairwise_error = tril(D-M,-1); % 
pairwise_sqerr = pairwise_error.^2; % matrix of squared errors
J = sum(sum(weights.*pairwise_sqerr));

end