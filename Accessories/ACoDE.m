function varargout = ACoDE(varargin)
contextVector = [];
contextFitness = [];
initPop = [];
maxGen = [];
NP = [];
maxFEs = [];
dims = [];
Global = [];
proStr = {'contextVector', 'contextFitness', 'initPop', 'maxGen', 'NP', 'maxFEs', 'dims', 'Global'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~,Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if isempty(Global)
    error('no Global');
end
startFEs = Global.evaluated;
if isempty(dims)
    D = Global.problem.dimension;
	dims = 1: D;
else
	D = numel(dims);
end
if isempty(NP)
	error('NP unassigned');
end
if isempty(maxFEs)
   if isempty(maxGen)
       error('cannot allocate maxFEs since maxGen unassigned');
   else
       maxFEs = 3 * NP * maxGen;
   end
end
if isempty(contextVector)
	contextVector = Global.bestIndividual;
	contextFitness = Global.bestFitness;
elseif isempty(contextFitness)
	contextFitness = Global.evaluate(contextVector);
end
success = 0;
fail = 0;
LEP = 50;
p = repmat(Global.problem.lowerbound( dims), NP, 1) + rand(NP, D) .* (repmat(Global.problem.upperbound( dims) - Global.problem.lowerbound( dims), NP, 1));
PS = size(initPop, 1);
if PS > 0
    p(1: PS, :) = initPop(:, dims);
end
% Evaluate the objective function values
fit = Global.evaluate(combine(p, contextVector, dims));
for j = 1: 3
    % Success memory for each trial contextVector generation strategy
    sucMemo{j} = zeros(LEP, 3);
    % Failure memory for each trial contextVector generation strategy
    failMemo{j} = zeros(LEP, 3);
end
gen = 1;
while Global.evaluated - startFEs <= maxFEs && ~Global.terminated
    [best, ibest] = min(fit);
    if best < contextFitness
        contextVector(dims) = p(ibest, :);
        contextFitness = best;
        success = success + 1;
        fprintf('CYCLE %d: %.2e\n', gen, Global.bestFitness);
    else
        fail = fail + 1;
    end
    if gen < LEP
        % The control parameter settings have the same probability
        % to be selected for each trial contextVector generation strategy.
        paraProb = ones(3, 3) * 1/3;
    else
        for k = 1: 3
            % Compute the success rate of each control parameter
            % setting for each trial contextVector generation strategy
            paraSR = sum(sucMemo{k}) ./ (sum(sucMemo{k}) + sum(failMemo{k})) + 0.01;
            % Normalized the success rates
            paraProb(k, :) = paraSR ./ sum(paraSR);
        end
    end
    pTemp = p;
    fitTemp = fit;
    % uSet: the set of trial vectors
    uSet = zeros(3 * NP, D);
    % numP{i}: at each generation, record which control parameter
    % setting is used to produce the trial contextVector for each target
    % contextVector, please note that three trial vectors are produced for
    % each target contextVector
    numP = cell(NP, 1);
    for i = 1: NP
        numP{i} = zeros(3);
    end
    % numS: at each generation, record which control parameter
    % setting is used to produce the trial contextVector entering the next
    % population successfully
    numS = zeros(3);
    for i = 1: NP
        % Generate the trail vectors
        [u, tempPara] = reproduce(p, [Global.problem.lowerbound(dims); Global.problem.upperbound(dims)], i, NP, D, paraProb);
        uSet(i * 3-2: 3 * i, :) = u;
        for k = 1: 3
            % Judge which control parameter setting is used to
            % produce the trial contextVector for each target contextVector,
            % please note that three trial vectors are generated
            numP{i}(k, 1) = numP{i}(k, 1) + length(find(tempPara(k, 1) == 1));
            numP{i}(k, 2) = numP{i}(k, 2) + length(find(tempPara(k, 1) == 2));
            numP{i}(k, 3) = numP{i}(k, 3) + length(find(tempPara(k, 1) == 3));
        end
    end
    % Evaluate the trial vectors
    fitSet = Global.evaluate(combine(uSet, contextVector, dims));
    for i = 1: NP
        % Find the best trial contextVector of the three trial vectors
        % generated for each target contextVector
        % Minindex: denotes which trial contextVector generation strategy
        % is used to generate the best trial contextVector
        [minVal, minIndex] = min(fitSet(3 * i - 2: 3 * i, :));
        bestInd = uSet(3 * (i - 1) + minIndex, :);
        bestIndFit = fitSet(3 * (i - 1) + minIndex, :);
        % Selection
        if fit(i) >= bestIndFit
            pTemp(i, :) = bestInd;
            fitTemp(i, :) = bestIndFit;
            % Judge which control parameter setting is used to
            % generate the best trial entering the next population
            temp = find(numP{i}(minIndex, :) == 1);
            numS(minIndex, temp) = numS(minIndex, temp) + 1;
        end
    end
    p = pTemp;
    fit = fitTemp;
    % Update the success memory
    for k = 1: 3
        sucMemo{k}(1, :) = [];
        sucMemo{k}(LEP, :) = numS(k, :);
    end
    % Record the total number of each control parameter setting
    % used at each generation for each trial contextVector generation
    % strategy
    totalNum = zeros(3);
    for i = 1: NP
        totalNum = totalNum + numP{i};
    end
    % Update the failure memory
    for k = 1: 3
        failMemo{k}(1, :) = [];
        failMemo{k}(LEP, :) = totalNum(k, :) - numS(k, :);
    end
    gen = gen + 1;
end
if nargout == 2
    varargout = {success, fail};
elseif nargout == 3
    varargout = {success, fail, combine(p, contextVector, dims)};
end
end
function [u, tempPara] = reproduce(p, lu, i, NP, n, paraProb)
% The three control parameter settings
F    =  [1.0 1.0 0.8];
CR =  [0.1 0.9 0.2];
tempPara = zeros(3, 1);
%.... "rand/1/bin" strategy ....%
% Choose the indices for mutation
indexSet = 1: NP;
indexSet(i) = [];
% Choose the first index
temp = floor(rand * (NP - 1)) + 1;
nouse(1) = indexSet(temp);
indexSet(temp) = [];
% Choose the second index
temp = floor(rand * (NP - 2)) + 1;
nouse(2) = indexSet(temp);
indexSet(temp) = [];
% Choose the third index
temp = floor(rand * (NP - 3)) + 1;
nouse(3) = indexSet(temp);
% Use the roulette wheel selection to select one control parameter setting
paraIndex = length(find(rand > cumsum(paraProb(1, :)))) + 1;
% Record which control parameter setting is used
tempPara(1) = paraIndex;
% Mutation
v1 = p(nouse(1), :) + F(paraIndex) .* (p(nouse(2), :) - p(nouse(3), :));
% Handle the elements of the mutant contextVector which violate the boundary
vioLow = find(v1 < lu(1, :));
if ~isempty(vioLow)
    v1(1, vioLow) = 2 .* lu(1, vioLow) - v1(1, vioLow);
    vioLowUpper = find(v1(1, vioLow) > lu(2, vioLow));
    if ~isempty(vioLowUpper)
        v1(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
    end
end
vioUpper = find(v1 > lu(2, :));
if ~isempty(vioUpper)
    v1(1, vioUpper) = 2 .* lu(2, vioUpper) - v1(1, vioUpper);
    vioUpperLow = find(v1(1, vioUpper) < lu(1, vioUpper));
    if ~isempty(vioUpperLow)
        v1(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
    end
end
% Binomial crossover
j_rand = floor(rand * n) + 1;
t = rand(1, n) < CR(paraIndex);
t(1, j_rand) = 1;
t_ = 1 - t;
u(1, :) = t .* v1 + t_ .* p(i, :);
%... "current to rand/1" strategy ...%
% The mechanism to choose the indices for mutation is slightly different from that of the classic
% DE, we found that using the following mechanism to choose the indices for
% mutation can improve the performance to certain degree
nouse = floor(rand(1, 3) * NP) + 1;
% Use the roulette wheel selection to select one control parameter setting
paraIndex = length(find(rand > cumsum(paraProb(2, :)))) + 1;
% Record which control parameter setting is used
tempPara(2) = paraIndex;
% Mutation
v2 = p(i, :) + rand * (p(nouse(1), :) - p(i, :)) + F(paraIndex) .* (p(nouse(2), :) - p(nouse(3), :));
% Handle the elements of the mutant contextVector which violate the boundary
vioLow = find(v2 < lu(1, :));
if ~isempty(vioLow)
    v2(1, vioLow) = 2 .* lu(1, vioLow) - v2(1, vioLow);
    vioLowUpper = find(v2(1, vioLow) > lu(2, vioLow));
    if ~isempty(vioLowUpper)
        v2(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
    end
end
vioUpper = find(v2 > lu(2, :));
if ~isempty(vioUpper)
    v2(1, vioUpper) = 2 .* lu(2, vioUpper) - v2(1, vioUpper);
    vioUpperLow = find(v2(1, vioUpper) < lu(1, vioUpper));
    if ~isempty(vioUpperLow)
        v2(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
    end
end
% Binomial crossover is not used to generate the trial contextVector
u(2, :) = v2;
%... "rand/2/bin" strategy ...%
% Choose the indices for mutation
indexSet = 1: NP;
indexSet(i) = [];
% Choose the first index
temp = floor(rand * (NP - 1)) + 1;
nouse(1) = indexSet(temp);
indexSet(temp) = [];
% Choose the second index
temp = floor(rand * (NP - 2)) + 1;
nouse(2) = indexSet(temp);
indexSet(temp) = [];
% Choose the third index
temp = floor(rand * (NP - 3)) + 1;
nouse(3) = indexSet(temp);
indexSet(temp) = [];
% Choose the fourth index
temp = floor(rand * (NP - 4)) + 1;
nouse(4) = indexSet(temp);
indexSet(temp) = [];
% Choose the fifth index
temp = floor(rand * (NP - 5)) + 1;
nouse(5) = indexSet(temp);
% Use the roulette wheel selection to select one control parameter setting
paraIndex = length(find(rand > cumsum(paraProb(3, :)))) + 1;
% Record which control parameter setting is used
tempPara(3) = paraIndex;
% Mutation
% The first scaling factor is randomly chosen from 0 to 1
v3 = p(nouse(1), :) + rand .* (p(nouse(2), :) - p(nouse(3), :)) + F(paraIndex) .* (p(nouse(4), :) - p(nouse(5), :));
% Handle the elements of the mutant contextVector which violate the boundary
vioLow = find(v3 < lu(1, :));
if ~isempty(vioLow)
    v3(1, vioLow) = 2 .* lu(1, vioLow) - v3(1, vioLow);
    vioLowUpper = find(v3(1, vioLow) > lu(2, vioLow));
    if ~isempty(vioLowUpper)
        v3(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
    end
end
vioUpper = find(v3 > lu(2, :));
if ~isempty(vioUpper)
    v3(1, vioUpper) = 2 .* lu(2, vioUpper) - v3(1, vioUpper);
    vioUpperLow = find(v3(1, vioUpper) < lu(1, vioUpper));
    if ~isempty(vioUpperLow)
        v3(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
    end
end
% Binomial crossover
j_rand = floor(rand * n) + 1;
t = rand(1, n) < CR(paraIndex);
t(1, j_rand) = 1;
t_ = 1 - t;
u(3, :) = t .* v3 + t_ .* p(i, :);
end
