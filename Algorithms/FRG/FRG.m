function FRG(Global)
% Frequent Random Grouping
Global.evaluate(initPopulation([], 200, Global));
groupSize = 50;
rounds = ceil(Global.problem.dimension / groupSize) * 30;
groups = 1: ceil(Global.problem.dimension / groupSize);
NP = 50;
gen = 0;
operator = 'ACoDE';
while ~Global.terminated
    gen = gen + 1;
    fprintf('FRG: CYCLE %d\n', gen);
    dims = randperm(Global.problem.dimension);
    for i = groups
        subDim = dims(1 + 50 * (i - 1): min(50 * i, Global.problem.dimension));
        if strcmp(operator, 'LSHADE')
            LSHADE('-initPop', Global.bestIndividual, '-dims', subDim, '-meanNP', NP, '-maxFEs', Global.evaluation / rounds, '-Global', Global);
        elseif strcmp(operator, 'ACoDE')
            ACoDE('-initPop', Global.bestIndividual, '-dims', subDim, '-NP', NP, '-maxFEs', Global.evaluation / rounds, '-Global', Global);
        end
        if Global.terminated
            return;
        end
        renderCurve(Global);
    end
end
end