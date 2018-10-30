function DECCR(Global)
recordGroupings('-Global', Global);
NP = 15;
groupSize = 50;
maxGen = 250;
pop = initPopulation([], NP, Global);
Global.evaluate(pop);
while ~Global.terminated
    group = randomGroup(Global.problem.dimension, groupSize);
    activeness = zeros(1, numel(group));
    for i = 1: numel(group)
        renderCurve(Global);
        dims = group{i};
        [activity, ~, ~] = LSHADE('-HM', [], '-contextFitness', Global.bestFitness, '-contextVector', Global.bestIndividual, '-initPop', Global.bestIndividual, '-initFit', Global.bestFitness, '-maxGen', maxGen, '-meanNP', NP, '-dims', dims, '-Global', Global);
        success = sum(activity);
        fail = length(activity) - success;
        activeness(i) = success / (success + fail);
    end
    recordGroupings('-activeness', mean(activeness), '-groups', group, '-Global', Global);
end
end