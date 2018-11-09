function [xmean, C, sigma] = CMAES(varargin)
xmean = [];
C = [];
contextVector = [];
maxGen = [];
maxFEs = [];
lambda = [];
dims = [];
Global = [];
sigma = [];
proStr = {'sigma', 'xmean', 'C', 'contextVector', 'maxGen', 'lambda', 'maxFEs', 'dims', 'Global'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if isempty(sigma)
    sigma = 0.5;
end
if isempty(contextVector)
    if isempty(Global.bestIndividual)
        Global.evaluate(Global.problem.lowerbound + (Global.problem.upperbound - Global.problem.lowerbound) .* rand(1, Global.problem.dimension));
    end
    contextVector = Global.bestIndividual;
end
dimension = numel(dims);
if isempty(lambda)
    lambda = 4 + floor(3 * log(dimension));
end
if isempty(maxFEs)
    if ~isempty(maxGen)
        maxFEs = maxGen * lambda;
    else
        error('max evaluation consumption unassigned');
    end
end
% --------------------  Initialization --------------------------------
% User defined input parameters (need to be edited)
if isempty(xmean)
    xmean = rand(dimension, 1);    % objective variables initial point
end
% Strategy parameter setting: Selection
mu = lambda/2;               % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);     % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); %�?�代方差有效数�? variance-effectiveness of sum w_i x_i
% Strategy parameter setting: Adaptation
cc = (4 + mueff/dimension) / (dimension+4 + 2*mueff/dimension); % time constant for cumulation for C
cs = (mueff+2) / (dimension+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((dimension+1.3)^2+mueff);    %rank-one的学习率 learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((dimension+2)^2+mueff));  % rank-mu的学习率 for rank-mu update
damps = 1 + 2 * max(0, sqrt((mueff-1)/(dimension+1))-1) + cs; % damping for sigma usually close to 1
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(dimension,1); ps = zeros(dimension,1);   % evolution paths for C and sigma
if isempty(C)
    B = eye(dimension,dimension);                       % B defines the coordinate system
    D = ones(dimension,1);                      % diagonal D defines the scaling
    C = B * diag(D.^2) * B';            % covariance matrix C
else
    [B, D] = eig(C);
    D = sqrt(abs(diag(D)));        % D contains standard deviations now
end
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN = dimension^0.5*(1-1/(4*dimension)+1/(21*dimension^2));  % expectation of ||dimension(0,I)|| == norm(randn(dimension,1))
G = 0;
% -------------------- Generation Loop --------------------------------
startFEs = Global.evaluated;
ub = Global.problem.upperbound(dims)';
lb = Global.problem.lowerbound(dims)';
while Global.evaluated - startFEs < maxFEs %#ok<*BDSCI>
    G = G + 1;
    % Generate and evaluate lambda offspring
    arx = repmat(xmean, 1, lambda);
    for k = 1: lambda
        arx(:, k) = arx(:, k) + sigma * B * (D .* randn(dimension, 1)); % m + sig * Normal(0,C)
    end
    arx = trimPopulation(arx, lb, ub);
    arfitness = Global.evaluate(combine(arx', contextVector, dims));
    
    % Sort by fitness and compute weighted mean into xmean
    [~, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value（新的�?�值）
    xmean = trimPopulation(xmean, lb, ub);
    
    % Cumulation: Update evolution paths
    ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = sum(ps.^2)/(1-(1-cs)^(2*(Global.evaluated - startFEs)/lambda))/dimension < 2 + 4/(dimension+1);
    pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C + c1 * (pc * pc' + (1-hsig) * cc*(2-cc) * C) + cmu * artmp * diag(weights) * artmp';
    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    % Update B and D from C
    if (Global.evaluated - startFEs) - eigeneval > lambda/(c1+cmu)/dimension/10  % to achieve O(dimension^2)
        eigeneval = (Global.evaluated - startFEs);
        C = triu(C) + triu(C, 1)'; % enforce symmetry
        [B, D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
        D = sqrt(abs(diag(D)));        % D contains standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
end
end