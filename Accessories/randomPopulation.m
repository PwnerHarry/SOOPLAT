function P = randomPopulation(NP, evaluated_flag, Global)
X = repmat(Global.problem.lowerbound, NP, 1) + repmat(Global.problem.upperbound - Global.problem.lowerbound, NP, 1) .* rand(NP, numel);
for i = 1: NP
    P(i) = INDIVIDUAL(X(i, :));
end
if evaluated_flag
    Global.evaluateIndividuals(P);
end
end