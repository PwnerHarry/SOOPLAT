function P = initPopulation(subdims, NP, Global)
if isempty(subdims)
    subdims = 1: Global.problem.dimension;
end
P = repmat(Global.problem.lowerbound(subdims), NP, 1) + repmat(Global.problem.upperbound(subdims) - Global.problem.lowerbound(subdims), NP, 1) .* rand(NP, numel(subdims));
end