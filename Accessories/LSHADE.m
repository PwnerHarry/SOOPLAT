function varargout = LSHADE(varargin)
% stagnantCheck([]);
contextVector = [];
contextFitness = [];
initPop = [];
initFit = [];
maxGen = [];
initNP = [];
meanNP = [];
maxFEs = [];
NP = [];
dims = [];
minNP = [];
Global = [];
HM = [];
proStr = {'initFit', 'HM', 'contextVector', 'contextFitness', 'initPop', 'maxGen', 'NP', 'initNP', 'meanNP', 'maxFEs', 'dims', 'minNP', 'Global'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
%% set initNP, initPop, initialize Track, Historical Memory and Archive
startFEs = Global.evaluated;
if isempty(dims)
    D = Global.problem.dimension;
    dims = 1: D;
else
    D = numel(dims);
end
if isempty(contextVector)
    if isempty(Global.bestIndividual)
        Global.evaluate(Global.problem.lowerbound + (Global.problem.upperbound - Global.problem.lowerbound) .* rand(1, Global.problem.dimension));
    end
    contextVector = Global.bestIndividual;
    contextFitness = Global.bestFitness;
elseif isempty(contextFitness)
    contextFitness = Global.evaluate(contextVector);
end
if isempty(NP)
    if isempty(minNP)
        minNP = 4;
    end
    if isempty(initNP)
        if ~isempty(meanNP)
            initNP = 2 * meanNP - minNP;
        else
            error('NP unassigned');
        end
    end
else
    initNP = NP;
    minNP = NP;
end
if isempty(maxFEs)
    if ~isempty(maxGen)
        if isempty(meanNP)
            maxFEs = maxGen * (initNP + minNP) / 2;
        else
            maxFEs = maxGen * meanNP;
        end
    else
        error('max evaluation consumption unassigned');
    end
end
currentNP = initNP;
% X = repmat(Global.problem.lowerbound(dims), currentNP, 1) + rand(currentNP, D) .* repmat(Global.problem.upperbound(dims) - Global.problem.lowerbound(dims), currentNP, 1);
PS = size(initPop, 1);
if PS > currentNP
    error('initPop too large');
elseif PS == currentNP
    X = initPop(:, dims);
    if isempty(initFit)
        valParents = Global.evaluate(combine(X, contextVector, dims));
    else
        valParents = initFit;
    end
else
    X = NaN(currentNP, numel(dims));
    if isempty(initPop)
        X(1, :) = contextVector(dims);
        X(2: end, :) = normPopulation(dims, contextVector, currentNP - 1, Global);
        % X(2: end, :) = initPopulation(dims, currentNP - 1, Global);
        valParents = Global.evaluate(combine(X, contextVector, dims));
    else
        X(1: PS, :) = initPop(:, dims);
        X(PS + 1: end, :) = normPopulation(dims, contextVector, currentNP - PS, Global);
        if isempty(initFit)
            valParents = Global.evaluate(combine(X, contextVector, dims));
        else
            valParents = NaN(currentNP, 1);
            valParents(1: PS) = initFit;
            valParents(PS + 1: end) = Global.evaluate(combine(X(PS + 1: end, :), contextVector, dims));
        end
    end
end
CR = zeros(currentNP, 1);
F = zeros(currentNP, 1);
if isempty(HM)
    H = 6;% Size of HM = 6, according to Tanabe's paper
    HM = initHM(H, 0.5, 0.5);
end
Afactor = 2.6;
archive = Initialize_Archive(Afactor, initNP);
k = 1;
terminal_MCR = -1;
p = 0.11;
G = 1;
activity = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ascIndex] = sort(valParents, 'ascend');
ibest = ascIndex(1);
best = valParents(ibest);
if best < contextFitness
    contextFitness = best;
    contextVector(dims) = X(ibest, :);
    activity = [activity, 1];
else
    activity = [activity, 0];
end
% Main Loop
while Global.evaluated - startFEs < maxFEs && ~Global.terminated %#ok<*BDLGI>
    Variance = NaN(size(X));
    r = randi(numel(HM), currentNP, 1); % ri the index randomly generated from [1, H]
    % MUTATION
    for i = 1: currentNP
        % current-to-pbest/1/bin and archive
        [CR(i), F(i)] = generateCRF(HM(r(i)).MCR, HM(r(i)).MF, 0.1, 0.1, terminal_MCR);
        pbest = randpick(1: max(round(p * currentNP), 2));
        pbest_index = ascIndex(pbest);
        % r1 = randpick(setdiff(1: currentNP, [i, pbest_index]));
        % r2 = randpick(setdiff(1: currentNP + size(archive.X, 1), [i, pbest_index, r1]));
        r1 = randpick([1: i - 1, i + 1: currentNP]);
        r2 = randpick(setdiff(1: currentNP + size(archive.X, 1), [i, r1]));
        if r2 > size(X, 1)
            x2 = archive.X(r2 - size(X, 1), :);
        else
            x2 = X(r2, :);
        end
        Variance(i, :) = F(i)' .* (X(pbest_index, :) - X(i, :) + X(r1, :) - x2);
    end
    U = BoundConstraint(X + Variance, X, Global.problem.lowerbound(dims), Global.problem.upperbound(dims));
    % CROSSOVER
    CR_MAP = repmat(CR, 1, D);
    CR_MAP(:, randi(D)) = 1;
    mask = rand(currentNP, D) > CR_MAP; % mask is used to indicate which elements of U comes from the parent
    U(mask) = X(mask);
    % SELECTION
    valOffspring = Global.evaluate(combine(U, contextVector, dims));
    Good_Offspring_Indice = valOffspring <= valParents;
    Bad_Parent_Indice = find(valOffspring < valParents);
    % update Archive
    archive = updateArchive(archive, X(Bad_Parent_Indice, :), valParents(Bad_Parent_Indice));
    % update HM
    if ~isempty(Bad_Parent_Indice) %update HM(k).MCR and HM(k).MF based on SCR and SF
        SCR = CR(Bad_Parent_Indice); % CR used by successful individuals
        SF = F(Bad_Parent_Indice); % F used by successful individuals
        
        w = zeros(size(SCR));
        for j = 1: numel(Bad_Parent_Indice)
            w(j) = norm(U(Bad_Parent_Indice(j), :) - X(Bad_Parent_Indice(j), :));
        end
        w = w ./ sum(w); % weight contextVector
        if HM(k).MCR == terminal_MCR || max(SCR) == 0
            HM(k).MCR = terminal_MCR;
        else
            HM(k).MCR = w' * SCR; % Weighted Mean
        end
        HM(k).MF = sum(w .* (SF .^ 2)) / sum(w .* SF); % Weighted Lehmer Mean proposed by Peng et al.
        k = mod(k, numel(HM)) + 1;
    end
    valParents(Good_Offspring_Indice) = valOffspring(Good_Offspring_Indice);
    X(Good_Offspring_Indice, :) = U(Good_Offspring_Indice, :);
    [~, ascIndex] = sort(valParents, 'ascend');
    ibest = ascIndex(1);
    best = valParents(ibest);
    if best < contextFitness
        contextFitness = best;
        contextVector(dims) = U(ibest, :); %#ok<*AGROW>
        activity = [activity, 1];
    else
        activity = [activity, 0];
    end
    %%%%%%%%%%%%%%%%%%%%%
    newNP = round(initNP - (initNP - minNP) * (Global.evaluated - startFEs) / maxFEs);
    if newNP < currentNP
        currentNP = newNP;
        archive.NP = round(Afactor * currentNP);
        archive = updateArchive(archive, [], []);
        remain_index = ascIndex(1: currentNP);
        X = X(remain_index, :);
        valParents = valParents(remain_index);
        ascIndex = 1: currentNP;
        CR = zeros(currentNP, 1);
        F = zeros(currentNP, 1);
    end
%     if stagnantCheck(X)
%         HM = [];
%         break;
%     end
    G = G + 1;
end
if nargout == 1
    varargout = {activity};
elseif nargout == 2
    varargout = {activity, X};
elseif nargout == 3
    varargout = {activity, X, HM};
end
end
function flag = stagnantCheck(P)
global espace_counter;
persistent M STD counter;
if isempty(P)
    M = [];
    STD = [];
    counter = 0;
end
currentM = mean(P, 1);
currentSTD = std(P, 0, 1);
if isequal(M, currentM) && isequal(STD, currentSTD)
    counter = counter + 1;
    if counter >= size(P, 2)
        flag = true;
        espace_counter = espace_counter + 1;
        fprintf('stagnant escaped for the %dth time\n', espace_counter);
        return;
    end
else
    M = currentM;
    STD = currentSTD;
end
flag = false;
end

%% Generate CR and F with MCR, MF, SigmaCR, SigmaF
function [CR, F] = generateCRF(MCR, MF, sCR, sF, terminal_MCR)
if MCR == terminal_MCR
    CR = 0;
else
    CR = max(0, min(1, (normrnd(MCR, sCR))));
end
F = 0;
while F <= 0
    F = min(1, cauchyrnd(MF, sF, 1, 1));
end
end
%% Generate random with Cauchy Distribution
function result = cauchyrnd(mu, delta, m, n)
result = mu + delta .* tan(pi .* (rand(m, n) - 0.5));
end
%% Check the bounds, if bad, throw it back
function vi = BoundConstraint (vi, X, lb, ub)
% if the boundary constraint is violated
% set the value to be the middle of the previous value and the bound
% check the lower bound
xl = repmat(lb, size(X, 1), 1);
pos = vi < xl;
vi(pos) = (X(pos) + xl(pos)) / 2;
% check the upper bound
xu = repmat(ub, size(X, 1), 1);
pos = vi > xu;
vi(pos) = (X(pos) + xu(pos)) / 2;
end
%% Update the archive with bad vectors to maintain diversity
function archive = updateArchive(archive, X, funvalue)
archive.X = [archive.X; X];
archive.funvalues = [archive.funvalues; funvalue];
[~, IX]= unique(archive.X, 'rows');
if length(IX) < size(archive.X, 1)
    archive.X = archive.X(IX, :);
    archive.funvalues = archive.funvalues(IX, :);
end
if size(archive.X, 1) > archive.NP % add all new individuals
    r = randperm(size(archive.X, 1));
    rndpos = r(1: archive.NP);
    archive.X = archive.X(rndpos, :);
    archive.funvalues = archive.funvalues(rndpos, :);
end
end
%% Historical memory initialization
function HM = initHM(H, MCR, MF)
HM.MCR = MCR;
HM.MF = MF;
HM = repmat(HM, 1, H);
end
%% Archive initialization
function archive = Initialize_Archive(Afactor, NP)
archive.NP = round(Afactor * NP); % the maximum size of the archive
archive.X = []; % the soGlobal.boundarytions stored in te archive
archive.funvalues = []; % the function vaGlobal.boundarye of the archived soGlobal.boundarytions
end