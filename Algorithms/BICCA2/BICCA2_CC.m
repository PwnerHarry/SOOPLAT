function BICCA2_CC(Global)
global CC_RESTART_FLAG CC_OPTIMIZER;
recordGroupings('-Global', Global);
CC_RESTART_FLAG = true;
CC_OPTIMIZER = 'LSHADE';
GLOBAL_RESTART_COUNTER = 0;
cycles = 0;
NP = 50;
[~, SOA] = sort(rand(50, Global.problem.dimension));
X = SOA ./ NP .* repmat(Global.problem.upperbound-Global.problem.lowerbound, NP, 1) + repmat(Global.problem.lowerbound, NP, 1);
Global.evaluate(X);
patterns = [];
while ~Global.terminated
    [cratio, patterns] = CC(Global, patterns);
    if cratio < 0.05
        GLOBAL_RESTART_COUNTER = GLOBAL_RESTART_COUNTER + 1;
        if GLOBAL_RESTART_COUNTER >= 3
            GLOBAL_RESTART_COUNTER = 0;
            fprintf('BICCA2: GLOBAL RESTART\n');
            CC_RESTART_FLAG = true;
        end
    else
        GLOBAL_RESTART_COUNTER = 0;
    end
end
fprintf('cycles: %d\n', cycles);
end

function [cratio, patterns] = CC(Global, patterns)
global CC_RESTART_FLAG CC_OPTIMIZER;
persistent GPOP U maxGen NP H CC_RESTART_COUNTER;
if strcmp(CC_OPTIMIZER, 'SaNSDE')
    persistent ccm; %#ok<*TLEV>
elseif strcmp(CC_OPTIMIZER, 'LSHADE')
    persistent HM;
end
if CC_RESTART_FLAG || isempty(patterns)
    CC_RESTART_COUNTER = 0;
    H = 100;
    maxGen = 250;
    NP = 50;
    if strcmp(CC_OPTIMIZER, 'SaNSDE')
        ccm = 0.5;
    elseif strcmp(CC_OPTIMIZER, 'LSHADE')
        HM = [];
    end
    clear patterns;
    if isempty(U)
        groups = randomGroup(Global.problem.dimension, H);
    else
        groups = brandnewGroups(U, H, Global);
    end
    for i = 1: numel(groups)
        patterns(i) = PATTERN(groups{i});
    end
    U = zeros(Global.problem.dimension, Global.problem.dimension);
    CC_RESTART_FLAG = false;
    GPOP = initPopulation([], NP, Global);
end
fitness_before_CC = Global.bestFitness;
cvector_before_CC = Global.bestIndividual;
fprintf('CC: ');
% ASYNCRONOUS OPTIMIZATION
group = {};
for i = 1: numel(patterns)
    group{i} = pattern(i).core;
end
activeness = [patterns.success] ./ ([patterns.success] + [patterns.fail]);
recordGroupings('-activeness', mean(activeness), '-groups', group, '-Global', Global);
[~, I] = sort(activeness, 'descend');
for i = I
    if strcmp(CC_OPTIMIZER, 'SaNSDE')
        [GPOP(:, patterns(i).core), ccm, success, fail] = SaNSDE('-ccm', ccm, '-contextFitness', Global.bestFitness, '-contextVector', Global.bestIndividual, '-initPop', GPOP, '-maxGen', maxGen, '-NP', NP, '-dims', patterns(i).core, '-Global', Global);
    elseif strcmp(CC_OPTIMIZER, 'LSHADE')
        [activity, GPOP(:, patterns(i).core), HM] = LSHADE('-HM', HM, '-contextFitness', Global.bestFitness, '-contextVector', Global.bestIndividual, '-initPop', GPOP, '-maxGen', maxGen, '-NP', NP, '-dims', patterns(i).core, '-Global', Global);
        success = sum(activity);
        fail = length(activity) - success;
    end
    patterns(i).success = patterns(i).success + success;
    patterns(i).fail = patterns(i).fail + fail;
    U(patterns(i).core, patterns(i).core) = 1;
    renderCurve(Global);
end
Global.evaluate(GPOP);
cvector_after_CC = Global.bestIndividual;
cvector_difference = cvector_before_CC - cvector_after_CC;
fitness_after_CC = Global.bestFitness;
cratio = min(1, abs((fitness_before_CC - fitness_after_CC) / fitness_before_CC));
if cratio < 0.05
    CC_RESTART_COUNTER = CC_RESTART_COUNTER + 1;
    if CC_RESTART_COUNTER >= 3
        CC_RESTART_COUNTER = 0;
        fprintf('BICCA2: CC RESTART\n');
        CC_RESTART_FLAG = true;
    end
else
    CC_RESTART_COUNTER = 0;
    CC_RESTART_FLAG = false;
end
if cratio > 0
    global BICCA2_LS_STAGNANT_FLAG;
    BICCA2_LS_STAGNANT_FLAG = false;
end
hold on;
fill([0, 1, 1, 0], [fitness_after_CC, fitness_after_CC, fitness_before_CC, fitness_before_CC], [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', '0.2');
hold off;
% cratio = abs(fitness_before_CC / fitness_after_CC);
% PATTERN EVOLUTION
activeness = [patterns.success] ./ ([patterns.success] + [patterns.fail]);
normalized_activeness = mapminmax(activeness, 0, 1);
for i = 1: numel(patterns)
    patterns(i).compactness = normalized_activeness(i);
    members = patterns(i).core;
    dim_difference = cvector_difference(members);
    [~, I] = sort(dim_difference, 'ascend');
    patterns(i).scatters = members(I(1: ceil((1 - patterns(i).compactness) * length(members))));
    patterns(i).core = setdiff(members, patterns(i).scatters);
end
for i = 1: numel(patterns)
    scatters_to_pick = H - length(patterns(i).core);
    if scatters_to_pick <= 0
        continue;
    else
        patterns(i).success = 0;
        patterns(i).fail = 0;
        other_scatters = setdiff([patterns.scatters], patterns(i).scatters);
        other_scatters = other_scatters(randperm(length(other_scatters)));
        if scatters_to_pick > length(other_scatters)
            gathered_scatters = [other_scatters, patterns(i).scatters];
            patterns(i).scatters = [];
        else
            gathered_scatters = other_scatters(1: scatters_to_pick);
        end
        patterns(i).core = [patterns(i).core, gathered_scatters];
        for j = setdiff(1: numel(patterns), i)
            patterns(j).scatters = setdiff(patterns(j).scatters, gathered_scatters);
        end
    end
end
fprintf('CRATIO = %.3f%%\n', 100 * cratio);
end

function groups = brandnewGroups(U, H, Global)
NP = 10000;
groups_population = {};
groups_fitness = NaN(NP, 1);
for i = 1: NP
    groups_population{i} = randomGroup(Global.problem.dimension, H);
    groups_fitness(i) = evaluateGrouping(groups_population{i}, U);
end
[~, I] = min(groups_fitness);
groups = groups_population{I};
end

function fit = evaluateGrouping(groups, U)
fit = 0;
for i = 1: numel(groups)
    dims = groups{i};
    fit = fit + sum(sum(U(dims, dims))) / 2;
end
fit = fit / numel(groups);
end