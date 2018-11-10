function varargout = LSHADE2(varargin)
[evaluated, D, HM, X, valParents, lb, ub, objfunc, currentNP, initNP, minNP, maxFEs, fbest] = parseARGS(varargin);
CR = zeros(currentNP, 1);
F = zeros(currentNP, 1);
Afactor = 2.6;
archive = initArchive(Afactor, initNP);
k = 1;
terminal_MCR = -1;
p = 0.11;
G = 1;
activity = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ascIndex] = sort(valParents, 'ascend');
ibest = ascIndex(1);
best = valParents(ibest);
if best < fbest
    fbest = best;
    activity = [activity, 1];
else
    activity = [activity, 0];
end
while evaluated < maxFEs
    Variance = NaN(size(X));
    r = randi(numel(HM), currentNP, 1); % ri the index randomly generated from [1, H]
    % MUTATION
    for i = 1: currentNP
        % current-to-pbest/1/bin and archive
        [CR(i), F(i)] = generateCRF(HM(r(i)).MCR, HM(r(i)).MF, 0.1, 0.1, terminal_MCR);
        pbest = randpick(1: max(round(p * currentNP), 2));
        pbest_index = ascIndex(pbest);
        r1 = randpick([1: i - 1, i + 1: currentNP]);
        r2 = randpick(setdiff(1: currentNP + size(archive.X, 1), [i, r1]));
        if r2 > size(X, 1)
            x2 = archive.X(r2 - size(X, 1), :);
        else
            x2 = X(r2, :);
        end
        Variance(i, :) = F(i)' .* (X(pbest_index, :) - X(i, :) + X(r1, :) - x2);
    end
    U = BoundConstraint(X + Variance, X, lb, ub);
    % CROSSOVER
    CR_MAP = repmat(CR, 1, D);
    CR_MAP(:, randi(D)) = 1;
    mask = rand(currentNP, D) > CR_MAP; % mask is used to indicate which elements of U comes from the parent
    U(mask) = X(mask);
    % SELECTION
    valOffspring = feval(objfunc, U);
    evaluated = evaluated + size(U, 1); %Is this right?
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
        w = w ./ sum(w); % weight
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
    if best < fbest
        fbest = best;
        activity = [activity, 1];
    else
        activity = [activity, 0];
    end
    %%%%%%%%%%%%%%%%%%%%%
    newNP = round(initNP - (initNP - minNP) * evaluated / maxFEs);
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
function archive = initArchive(Afactor, NP)
archive.NP = round(Afactor * NP); % the maximum size of the archive
archive.X = []; % the soGlobal.boundarytions stored in te archive
archive.funvalues = []; % the function vaGlobal.boundarye of the archived soGlobal.boundarytions
end

function [evaluated, D, HM, X, valParents, lb, ub, objfunc, currentNP, initNP, minNP, maxFEs, xbest, fbest] = parseARGS(varargin)
varargin = varargin{:};
initPop = [];
initFit = [];
maxGen = [];
initNP = [];
meanNP = [];
maxFEs = [];
NP = [];
minNP = [];
HM = [];
xbest = [];
fbest = [];
lb = [];
ub = [];
objfunc = [];
proStr = {'lb', 'ub', 'objfunc', 'initPop', 'initFit', 'xbest', 'fbest', 'HM', 'maxGen', 'NP', 'initNP', 'meanNP', 'maxFEs', 'dims', 'minNP'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
%% set initNP, initPop, initialize Track, Historical Memory and Archive
D = length(lb);
if D ~= length(ub)
    error('lb and ub lengths inconsistent');
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
PS = size(initPop, 1);
if PS > currentNP
    error('initPop too large');
elseif PS == currentNP
    X = initPop;
    if isempty(initFit)
        valParents = feval(objfunc, X);
        evaluated = currentNP;
    else
        valParents = initFit;
        evaluated = 0;
    end
elseif PS == 0
    valParents = NaN(currentNP, 1);
    X = normPopulation2(currentNP, lb, ub, xbest, diag(0.5 * ub - lb));
    valParents(2: end) = feval(objfunc, X(2: end));
    evaluated = currentNP - 1;
else % 0 < PS < currentNP
    newX = normPopulation2(currentNP - PS + 1, lb, ub, xbest, diag(0.5 * ub - lb));
    X = [initPop; newX(2: end, :)];
    if isempty(initFit)
        valParents = feval(objfunc, X);
    else
        valParents = NaN(currentNP, 1);
        valParents(1: PS) = initFit;
        valParents(PS + 1: end) = feval(objfunc, X(PS + 1: end, :));
    end
    evaluated = currentNP - PS;
end
if isempty(HM)
    H = 6;% Size of HM = 6, according to Tanabe's paper
    HM = initHM(H, 0.5, 0.5);
end
end