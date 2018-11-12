function varargout = SaNSDE2(varargin)
[activity, D, evaluated, pop, val, xbest, fbest, ccm, objfunc, lb, ub, maxFEs, NP] = parseARGS(varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = zeros(NP, 1);
linkp = 0.5;
l1 = 1; l2 = 1;
nl1 = 1; nl2 = 1;
fp = 0.5;
ns1 = 1; nf1 = 1;
ns2 = 1; nf2 = 1;
ui = zeros(NP, D); % intermediate population of perturbed vectors
rot = (0: NP - 1); % rotating index array (size NP)
G = 0;
while evaluated < maxFEs || G == 0
    popold = pop; % save the old population
    ind = randperm(4); % index pointer array
    a1 = randperm(NP);                      pm1 = popold(a1, :);
    a2 = a1(rem(rot + ind(1), NP) + 1);     pm2 = popold(a2, :);  
    a3 = a2(rem(rot + ind(2), NP) + 1);     pm3 = popold(a3, :); 
    a4 = a3(rem(rot + ind(3), NP) + 1);     pm4 = popold(a4, :);
    bm = ones(NP, 1) * xbest;
    if rem(G, 24) == 0
        if G ~= 0 && ~isempty(cc_rec)
            ccm = sum(f_rec .* cc_rec) / sum(f_rec);
        end
        cc_rec = [];
        f_rec = [];
    end
    if rem(G, 5) == 0
        index = [];
        while numel(index) < NP
            cc = normrnd(ccm, 0.1, NP * 3, 1);
            index = find((cc < 1) & (cc > 0));
        end
        cc = cc(index(1: NP));
    end
    fst1 = (rand(NP, 1) <= fp); fst2 = ~fst1;
    tmp = normrnd(0.5, 0.3, NP, 1);
    F(fst1) = tmp(fst1);
    tmp = normrnd(0, 1, NP, 1) ./ normrnd(0, 1, NP, 1);
    F(fst2) = tmp(fst2);
    F = abs(F);
    aa = (rand(NP, D) < repmat(cc, 1, D));
    index = find(sum(aa, 2) == 0); % TODO: this should be written using any(...)
    tmpsize = length(index);
    for k = 1: tmpsize
        bb = ceil(D * rand());
        aa(index(k), bb) = 1;
    end
    mui = aa;
    mpo = mui < 0.5;                % inverse mask to mui
    aaa = (rand(NP, 1) <= linkp); bbb = 1 - aaa;
    aindex = find(aaa == 0);
    bindex = find(aaa ~= 0);
    if ~isempty(bindex) % mutation and crossover
        ui(bindex,:) = popold(bindex,:)+repmat(F(bindex,:),1,D).*(bm(bindex,:)-popold(bindex,:)) + repmat(F(bindex,:),1,D).*(pm1(bindex,:) - pm2(bindex,:) + pm3(bindex,:) - pm4(bindex,:));
        ui(bindex,:) = popold(bindex,:).*mpo(bindex,:) + ui(bindex,:).*mui(bindex,:);
    end
    if ~isempty(aindex)
        ui(aindex,:) = pm3(aindex,:) + repmat(F(aindex,:),1,D).*(pm1(aindex,:) - pm2(aindex,:));
        ui(aindex,:) = popold(aindex,:).*mpo(aindex,:) + ui(aindex,:).*mui(aindex,:);
    end
    %-----Select which vectors are allowed to enter the new population-------
    ui = trimPopulation(ui, lb, ub, 'mod');
    tempval = feval(objfunc, ui);   evaluated = evaluated + length(tempval);
    for i = 1: NP
        if tempval(i) <= val(i)
            if tempval(i) < val(i)
                cc_rec = [cc_rec, cc(i, 1)]; %#ok<*AGROW>
                f_rec = [f_rec, (val(i) - tempval(i))];
            end
            pop(i, :) = ui(i, :);   val(i) = tempval(i);
            l1 = l1 + aaa(i);       l2 = l2 + bbb(i);
            ns1 = ns1 + fst1(i);    ns2 = ns2 + fst2(i);
        else
            nl1 = nl1 + aaa(i);     nl2 = nl2 + bbb(i);
            nf1 = nf1 + fst1(i);    nf2 = nf2 + fst2(i);
        end
    end
    [best, ibest] = min(val);
    if best < fbest
        fbest = best;   xbest = pop(ibest, :);
        activity = [activity, 1];
    else
        activity = [activity, 0];
    end
    if rem(G, 24) == 0 && G ~= 0
        linkp = (l1 / (l1 + nl1)) / (l1 / (l1 + nl1) + l2 / (l2 + nl2));
        l1 = 1;     l2 = 1;
        nl1 = 1;    nl2 = 1;
        fp = (ns1 * (ns2 + nf2)) / (ns2 * (ns1 + nf1) + ns1 * (ns2 + nf2));
        ns1 = 1;    nf1 = 1;
        ns2 = 1;    nf2 = 1;
    end
    G = G + 1;
end
if nargout == 1
    varargout = {pop};
elseif nargout == 2
    varargout = {pop, ccm};
elseif nargout == 3
    varargout = {pop, ccm, activity};
elseif nargout == 4
    varargout = {pop, ccm, sum(activity), length(activity) - sum(activity)};
end
end

function [activity, D, evaluated, pop, val, xbest, fbest, ccm, objfunc, lb, ub, maxFEs, NP] = parseARGS(varargin)
varargin = varargin{:};
xbest = [];
fbest =[];
ccm = [];
objfunc = [];
lb = [];
ub = [];
maxFEs = [];
NP = [];
initPop = [];
initFit = [];
proStr = {'xbest', 'fbest', 'initPop', 'initFit', 'ccm', 'objfunc', 'lb', 'ub', 'maxGen', 'NP', 'maxFEs'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if isempty(maxFEs) && ~isempty(maxGen)
    maxFEs = NP * maxGen;
end
if isempty(ccm)
    ccm = 0.5;
end
val = NaN(NP, 1);
D = intersect(length(lb), length(ub));
if size(initPop, 1) ~= NP %#ok<*BDSCI>
    error('size of initPop is not NP');
end
pop = initPop;
[tf, xbest_index] = ismember(xbest, pop, 'rows');
if tf
    val(xbest_index) = fbest;
end
if isempty(initFit)
    non_best_index = setdiff(1: NP, xbest_index);
    val(non_best_index) = feval(objfunc, pop(non_best_index, :));
    evaluated = length(non_best_index);
else
    if size(initPop, 1) ~= length(initFit)
        error('sizes of initPop and initFit inconsistent');
    end
    val = initFit;
    evaluated = 0;
end
[best, ibest] = min(val);
if best < fbest
    fbest = best;
    xbest = pop(ibest, :);
    activity = 1;
else
    activity = 0;
end
end