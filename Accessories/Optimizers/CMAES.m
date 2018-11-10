% Nikolaus Hansen. "The CMA Evolution Strategy: A Tutorial"
% https://arxiv.org/abs/1604.00772
function [xmean, C, sigma] = CMAES(varargin)
[lb, ub, objfunc, N, sigma, xmean, C, B, D, lambda, maxFEs] = parseargs(varargin);
% Strategy parameter setting: Selection
mu = lambda / 2; % number of parents/points for recombination
weights = log(mu + 1 / 2) - log(1: mu)'; % muXone array for weighted recombination
mu = floor(mu);
weights = weights / sum(weights); % normalize recombination weights array
mueff = sum(weights) ^ 2 / sum(weights .^ 2); %variance-effectiveness of sum w_i x_i
% Strategy parameter setting: Adaptation
cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N); % time constant for cumulation for C
cs = (mueff + 2) / (N + mueff + 5); % t-const for cumulation for sigma control
c1 = 2 / ((N + 1.3) ^ 2 + mueff); % learning rate for rank-one update of C
cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2) ^ 2 + mueff));  % rank-mu的学习率 for rank-mu update
damps = 1 + 2 * max(0, sqrt((mueff - 1)/(N + 1)) - 1) + cs; % damping for sigma usually close to 1
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N, 1); % evolution paths for C and sigma
ps = zeros(N, 1); % evolution paths for sigma
invsqrtC = B * diag(D .^ -1) * B'; % C ^ -0.5
eigeneval = 0; % track update of B and D
chiN = N ^ 0.5 * (1 - 1 / (4 * N) + 1 / (21 * N ^ 2)); % expectation of ||N(0,I)|| == norm(randn(N,1))
evaluated = 0;
while evaluated < maxFEs
    % Generate and evaluate lambda offspring
    arx = repmat(xmean, 1, lambda);
    for k = 1: lambda
        arx(:, k) = arx(:, k) + sigma * B * (D .* randn(N, 1)); % m + sig * Normal(0,C)
    end
    arx = trimPopulation(arx, lb, ub);
    arfitness = feval(objfunc, arx');
    evaluated = evaluated + lambda;
    % Sort by fitness and compute weighted mean into xmean
    [~, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:, arindex(1: mu)) * weights;  % recombination, new mean value
    % Cumulation: Update evolution paths
    ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * invsqrtC * (xmean - xold) / sigma;
    hsig = sum(ps .^ 2) / (1 - (1 - cs) ^ (2 * evaluated / lambda)) / N < 2 + 4 / (N + 1);
    pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * (xmean - xold) / sigma;
    % Adapt covariance matrix C
    artmp = (1 / sigma) * (arx(:, arindex(1: mu)) - repmat(xold, 1, mu));  % mu difference vectors
    C = (1 - c1 - cmu) * C + c1 * (pc * pc' + (1 - hsig) * cc * (2 - cc) * C) + cmu * artmp * diag(weights) * artmp';
    % Adapt step size sigma
    sigma = sigma * exp((cs / damps) * (norm(ps) / chiN - 1));
    % Update B and D from C
    if evaluated - eigeneval > lambda / (c1 + cmu) / N /10  % to achieve O(N^2)
        eigeneval = evaluated;
        C = triu(C) + triu(C, 1)'; % enforce symmetry
        [B, D] = eig(C); % eigen decomposition, B==normalized eigenvectors
        D = sqrt(abs(diag(D))); % D contains standard deviations now
        invsqrtC = B * diag(D .^ -1) * B';
    end
end
end

function [lb, ub, objfunc, N, sigma, xmean, C, B, D, lambda, maxFEs] = parseargs(varargin)
varargin = varargin{:};
xmean = [];
C = [];
maxFEs = [];
lambda = [];
sigma = [];
objfunc = [];
lb = [];
ub = [];
N = [];
proStr = {'lb', 'ub', 'sigma', 'xmean', 'N', 'C', 'objfunc', 'lambda', 'maxFEs'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if isempty(sigma)
    sigma = 0.5;
end
if isempty(lambda)
    lambda = 4 + floor(3 * log(N));
end
if isempty(maxFEs)
    if ~isempty(maxGen)
        maxFEs = maxGen * lambda;
    else
        error('max evaluation consumption unassigned');
    end
end
if isempty(xmean)
    xmean = rand(N, 1); % objective variables initial point
end
if isempty(C)
    B = eye(N, N); % B defines the coordinate system
    D = ones(N, 1); % diagonal D defines the scaling
    C = B * diag(D .^ 2) * B'; % covariance matrix C
else
    [B, D] = eig(C);
    D = sqrt(abs(diag(D)));        % D contains standard deviations now
end
end