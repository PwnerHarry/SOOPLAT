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
        renderCurve(Global);
        dims = group{i};
        [pop(:, dims), ccm] = SaNSDE('-ccm', ccm, '-NP', NP, '-dims', dims, '-initPop', pop, '-contextVector', Global.bestIndividual, '-contextFitness', Global.bestFitness, '-maxGen', maxGen, '-Global', Global);
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