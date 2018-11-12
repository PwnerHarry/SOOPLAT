function DECCG(Global)
NP = 100;
groupSize = 100;
maxGen = 199;
ccm = 0.5; % the initial crossover rate for SaNSDE
pop = initPopulation([], NP, Global);
Global.evaluate(pop);
while ~Global.terminated
    group = randomGroup(Global.problem.dimension, groupSize);
	for i = 1: numel(group)
        dims = group{i};
        lb = Global.problem.lowerbound(dims);
        ub = Global.problem.upperbound(dims);
        OBJFUNC = @(X) Global.evaluate(combine(X, Global.bestIndividual, dims));
        [pop(:, dims), ccm] = SaNSDE2('-ub', ub, '-lb', lb, '-objfunc', OBJFUNC, '-ccm', ccm, '-NP', NP, '-initPop', pop(:, dims), '-xbest', Global.bestIndividual(dims), '-fbest', Global.bestFitness, '-maxGen', maxGen);
        renderCurve(Global);
    end
    val = Global.evaluate(pop);
    de_weight(Global.bestIndividual, Global.bestFitness, Global.problem.lowerbound, Global.problem.upperbound, NP, maxGen, group, Global);
    [mem_val, mem_id] = max(val);
    [newmember, newmem_val] = de_weight(pop(mem_id, :), Global.bestFitness, Global.problem.lowerbound, Global.problem.upperbound, NP, maxGen, group, Global);
    if newmem_val < mem_val
        val(mem_id) = newmem_val;
        pop(mem_id, :) = newmember;
    end
    mem_id = randi(NP);
    mem_val = val(mem_id);
    [newmember, newmem_val] = de_weight(pop(mem_id, :), Global.bestFitness, Global.problem.lowerbound, Global.problem.upperbound, NP, maxGen, group, Global);
    if newmem_val < mem_val
        val(mem_id) = newmem_val;
        pop(mem_id, :) = newmember;
    end
end
end