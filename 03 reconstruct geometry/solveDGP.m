% Compute the most likely geometry of a group of points from pairwise distance values
% Usage: 
% [theta, prune, score] = solveDGP(M, weights, threshold)
% where:
% M is the output from finddist_geometry and is an array of measured distances for each pair of points
% weights is the output from finddist_geometry and is the peak height (in bases) corresponding to the distances measured. It is a measure of the confidence in the measurement.
% threshold is a parameter used to prune unreliable measurements and to generate weights for measurements reflecting their reliability (see supplementary text for details on how opt_threshold is auto-set),
% theta is a list of coordinates, specifying the final embedding,
% prune is a logical-valued array indicating which target points were dropped from the final reconstruction, and
% score is a measure of the internal consistency of the embedding. 

function [theta, prune, score] = solveDGP(M, weights, threshold)

if nargin == 1 % if no weights supplied, set all non-filtered weights to be equal to 1
    weights = tril(ones(size(M,1)),-1);
    threshold = 0; % lower threshold to accept all distance measurements with equal weight
elseif nargin == 2
    threshold = 20; % choose a reasonable threshold for filtering by weight
end

%% Get an initial guess for the embedding through ExpVecEDM
M_th = M.*(weights > threshold); % only keep distances above threshold
Dpartial = tril(M_th,-1).^2 +triu(M_th',1).^2;
[initial_theta, prune, theta_refined] = ExpVecEDM(Dpartial,[],2);

if ~isreal(theta_refined)
    disp('Warning! ExpVecEDM gave us complex initialization points. Will do random initialization instead ....');
    theta_refined = [30*rand(size(M,1),1) 30*rand(size(M,1),1)];
    prune = true(1,size(M,1));
end

n = size(theta_refined,1);

%% plot initial answer (uncomment to plot)
% figure; 
% hold on; % plot
% scatter(initial_theta(:,1), initial_theta(:,2),50,[0.5 0.5 0.5],'filled');
% %text(1+initial_theta(:,1), 1+initial_theta(:,2),num2str([1:n]'),'FontSize',18);
% axis equal
% hold off;
% % 
% figure; 
% hold on; % plot
% scatter(theta_refined(:,1), theta_refined(:,2),50,[0.5 0.5 0.5],'filled');
% %text(1+theta_refined(:,1), 1+theta_refined(:,2),num2str([1:n]'),'FontSize',18);
% axis equal
% hold off;
%% recast theta for fminunc
initial_theta = [initial_theta(:,1); initial_theta(:,2)]; % recast initial answer
theta_refined = [theta_refined(:,1); theta_refined(:,2)];  % recast initial answer

%% Prune M and weights
num_points_pruned = size(M,1) - sum(prune);
% fprintf('Pruned %d points', num_points_pruned);
M = M(prune,prune);
weights = weights(prune,prune);

%% find embedding that minimizes sum of weighted squared-errors
weights = tril(max(sigmoid(weights-threshold,0.8)-0.5,0),-1); % rectified positive sigmoid, lower triangular
weights = weights./sum(sum(weights)); % Normalize weights

% Set options for fminunc
opt = optimoptions('fminunc','MaxFunEvals',1e6,'MaxIter',50000,'TolFun',1e-9,'TolX',1e-9,'FiniteDifferenceType','forward','UseParallel',false);%,'Display','none');
[theta, funcval] = fminunc(@(theta)(DGPcostFunction(theta, M, weights)),theta_refined,opt); % run minimization
theta = [theta(1:n) theta(n+1:end)]; % recast theta
prune_penalty = 0.5;
score = funcval + prune_penalty*(num_points_pruned);
fprintf('Minimal score obtained is: %4.2f\n',score);

%% Plot the answer
figure; 
hold on; % plot
scatter(theta(:,1), theta(:,2),50,[0.5 0.5 0.5],'filled');
%text(1+theta(:,1), 1+theta(:,2),num2str([1:n]'),'FontSize',18);

xmin = min(theta(:,1));
xmax = max(theta(:,1)); % X-axis and Y-axis range of the points
ymin = min(theta(:,2));
ymax = max(theta(:,2));
xlim([xmin-5 xmax+5]);
ylim([ymin-5 ymax+5]);
axis equal; hold off;

%toc;
end