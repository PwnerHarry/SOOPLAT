function CCCMAES(Global)
Global.evaluate(initPopulation([], 200, Global));
m = 50; % dimensionality of the subspace
FEs = 10 * m; % FEs allocated for each subspace
lambda = 50; % sample size
adapter = ADAPTER({'MiVD', 'MaVD', 'RaVD'}, 5);
X_w = Global.bestIndividual;
C = eye(Global.problem.dimension, Global.problem.dimension);
sigma = 0.5 * (min(Global.problem.upperbound) - max(Global.problem.lowerbound));
while ~Global.terminated
    renderCurve(Global);
    old_fitness = Global.bestFitness;
    action = adapter.decide();
    groups = feval(action, Global.problem.dimension, m, C);
    for i = 1: numel(groups)
        dims = groups{i};
        [X_w(dims), C(dims, dims), sigma] = CMAES('-sigma', sigma, '-xmean', X_w(dims)', '-C', C(dims, dims), '-contextVector', Global.bestIndividual, '-lambda', lambda, '-maxFEs', FEs, '-dims', dims, '-Global', Global);
    end
    contribution_ratio = abs(old_fitness - Global.bestFitness) / old_fitness;
    adapter.update(action, contribution_ratio);
    fprintf('%s: %.2f%%\n', action, 100 * contribution_ratio);
end
end

