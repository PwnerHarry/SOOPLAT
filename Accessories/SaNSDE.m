function varargout = SaNSDE(varargin)
ccm = [];
contextVector = [];
contextFitness = [];
initPop = [];
initFit = [];
maxGen = [];
NP = [];
maxFEs = [];
dims = [];
Global = [];
proStr = {'ccm', 'contextVector', 'contextFitness', 'initPop', 'initFit', 'maxGen', 'NP', 'maxFEs', 'dims', 'Global'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if isempty(Global)
    error('no Global');
end
if isempty(ccm)
    ccm = 0.5;
end
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
        maxFEs = NP * maxGen;
    end
end
if isempty(contextVector)
    contextVector = Global.bestIndividual;
    contextFitness = Global.bestFitness;
elseif isempty(contextFitness)
    contextFitness = Global.evaluate(contextVector);
end
Lbound = repmat(Global.problem.lowerbound(dims), NP, 1);
Ubound = repmat(Global.problem.upperbound(dims), NP, 1);
if isempty(initPop)
    pop = initPopulation(dims, NP, Global);
    val = Global.evaluate(combine(pop, contextVector, dims));
else
    if NP == size(initPop, 1)
        pop = initPop(:, dims);
        if isempty(initFit)
            val = Global.evaluate(combine(pop, contextVector, dims));
        else
            val = initFit;
        end
    elseif size(initPop, 1) < NP
        newpop = initPopulation(dims, NP - size(initPop, 1), Global);
        pop = [initPop(:, dims); newpop];
        if isempty(initFit)
            val = Global.evaluate(combine(pop, contextVector, dims));
        else
            newval = Global.evaluate(combine(newpop, contextVector, dims));
            val = [initFit, newval];
        end
    else
        error('NP does not equal to initPop size');
    end
end
[best, ibest] = min(val);
if best < contextFitness
    contextVector(dims) = pop(ibest, :);
    contextFitness = best;
    success = 1;
    fail = 0;
else
    success = 0;
    fail = 1;
end
startFEs = Global.evaluated;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = zeros(NP,1);
linkp = 0.5;
l1 = 1;
l2 = 1;
nl1 = 1;
nl2 = 1;
fp = 0.5;
ns1 = 1;
nf1 = 1;
ns2 = 1;
nf2 = 1;
ui  = zeros(NP, D);              % intermediate population of perturbed vectors
rot = (0: NP - 1);               %#ok<*BDSCI> % rotating index array (size NP)
gen = 0;
while Global.evaluated - startFEs < maxFEs && ~Global.terminated %#ok<*BDLGI>
    popold = pop;                   % save the old population
    ind = randperm(4);              % index pointer array
    a1 = randperm(NP);             % shuffle locations of vectors
    a2 = a1(rem(rot + ind(1), NP) + 1);                 % rotate vector locations
    a3 = a2(rem(rot + ind(2), NP) + 1);
    a4 = a3(rem(rot + ind(3), NP) + 1);
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    bm = ones(NP, 1) * contextVector(1, dims);
    if rem(gen,24)==0
        if gen~=0 && ~isempty(cc_rec)
            ccm = sum(f_rec.*cc_rec)/sum(f_rec);
        end
        cc_rec = [];
        f_rec = [];
    end
    if rem(gen,5)==0
        index = [];
        while numel(index) < NP
            cc = normrnd(ccm, 0.1, NP * 3, 1);
            index = find((cc < 1) & (cc > 0));
        end
        cc = cc(index(1:NP));
    end
    fst1 = (rand(NP,1) <= fp);
    fst2 = 1-fst1;
    fst1_index = find(fst1 ~= 0);
    fst2_index = find(fst1 == 0);
    tmp = normrnd(0.5, 0.3, NP, 1);
    F(fst1_index) = tmp(fst1_index);
    tmp = normrnd(0, 1, NP, 1) ./ normrnd(0, 1, NP, 1);
    F(fst2_index) = tmp(fst2_index);
    F = abs(F);
    aa = rand(NP,D) < repmat(cc,1,D);
    index = find(sum(aa') == 0);% TODO: this should be written using any(...)
    tmpsize = size(index, 2);
    for k=1:tmpsize
        bb = ceil(D*rand);
        aa(index(k), bb) = 1;
    end
    mui=aa;
    mpo = mui < 0.5;                % inverse mask to mui
    aaa = (rand(NP,1) <= linkp);
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
    bbb = 1 - aaa;
    %-----Select which vectors are allowed to enter the new population-------
    index = find(ui > Ubound);
    ui(index) = Ubound(index) - mod((ui(index)-Ubound(index)), (Ubound(index)-Lbound(index)));
    index = find(ui < Lbound);
    ui(index) = Lbound(index) + mod((Lbound(index)-ui(index)), (Ubound(index)-Lbound(index)));
    tempval = Global.evaluate(combine(ui, contextVector, dims));
    for i = 1: NP
        if tempval(i) <= val(i)
            if tempval(i) < val(i)
                cc_rec = [cc_rec, cc(i,1)];
                f_rec = [f_rec, (val(i) - tempval(i))];
            end
            pop(i, :) = ui(i, :);
            val(i) = tempval(i);
            l1 = l1 + aaa(i);
            l2 = l2 + bbb(i);
            ns1 = ns1 + fst1(i);
            ns2 = ns2 + fst2(i);
        else
            nl1 = nl1 + aaa(i);
            nl2 = nl2 + bbb(i);
            nf1 = nf1 + fst1(i);
            nf2 = nf2 + fst2(i);
        end
    end
    [best, ibest] = min(val);
    if best < contextFitness
        success = success + 1;
        contextFitness = best;
        contextVector(dims) = pop(ibest, :);
    else
        fail = fail + 1;
    end
    if rem(gen, 24) == 0 && gen ~= 0
        linkp = (l1/(l1+nl1))/(l1/(l1+nl1)+l2/(l2+nl2));
        l1 = 1;
        l2 = 1;
        nl1 = 1;
        nl2 = 1;
        fp = (ns1 * (ns2 + nf2))/(ns2 * (ns1 + nf1) + ns1 * (ns2 + nf2));
        ns1 = 1;
        nf1 = 1;
        ns2 = 1;
        nf2 = 1;
    end
    gen = gen + 1;
end
if nargout == 1
    varargout = {pop};
elseif nargout == 2
    varargout = {pop, ccm};
elseif nargout == 4
    varargout = {pop, ccm, success, fail};
end