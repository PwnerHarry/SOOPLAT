function [cratio, patterns] = CC(Global, patterns)
global CC_RESTART_FLAG;
persistent r U maxGen NP H CC_RESTART_COUNTER HM;
if CC_RESTART_FLAG || isempty(patterns)
    CC_RESTART_COUNTER = 0;
    H = 100;
    maxGen = 500;
    NP = 15;
    HM = [];
    clear patterns;
    if isempty(U) || size(U, 2) ~= Global.problem.dimension
        groups = randomGroup(Global.problem.dimension, H);
    else
        groups = brandnewGroups(U, H, Global);
    end
    for i = 1: numel(groups)
        patterns(i) = PATTERN(groups{i});
    end
    r = 3;
    U.cube = zeros(r, Global.problem.dimension, Global.problem.dimension);
    U.pointer = 1;
    CC_RESTART_FLAG = false;
end
fitness_before_CC = Global.bestFitness;
cvector_before_CC = Global.bestIndividual;
fprintf('CC: ');
% ASYNCRONOUS OPTIMIZATION
activeness = [patterns.success] ./ ([patterns.success] + [patterns.fail]);
[~, I] = sort(activeness, 'descend');
for i = I
    dims = patterns(i).core;
    [activity, ~, ~] = LSHADE('-HM', [], '-contextFitness', Global.bestFitness, '-contextVector', Global.bestIndividual, '-initPop', Global.bestIndividual, '-initFit', Global.bestFitness, '-maxGen', maxGen, '-meanNP', NP, '-dims', dims, '-Global', Global);
    patterns(i).success = patterns(i).success + sum(activity);
    patterns(i).fail = patterns(i).fail + length(activity) - sum(activity);
    U.cube(U.pointer, patterns(i).core, patterns(i).core) = 1;
    renderCurve(Global);
end
group = {};
for i = 1: numel(patterns)
    group{i} = patterns(i).core;
end
activeness = [patterns.success] ./ ([patterns.success] + [patterns.fail]);
recordGroupings('-activeness', mean(activeness), '-groups', group, '-Global', Global);
U.pointer = mod(U.pointer, r) + 1;
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
drawnow;
hold off;
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
M = sum(U.cube, 1);
M(M > 0) = 1;
for i = 1: NP
    groups_population{i} = randomGroup(Global.problem.dimension, H);
    groups_fitness(i) = evaluateGrouping(groups_population{i}, M);
end
[~, I] = min(groups_fitness);
groups = groups_population{I};
end

function fit = evaluateGrouping(groups, M)
fit = 0;
for i = 1: numel(groups)
    dims = groups{i};
    fit = fit + sum(sum(sum(M(1, dims, dims))));
end
fit = fit / numel(groups);
end