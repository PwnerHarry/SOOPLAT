function P = normPopulation(dims, c, NP, Global)
P = NaN(NP - 1, numel(dims));
for i = 1: numel(dims)
    real_index = dims(i);
    P(:, i) = normrnd(Global.bestIndividual(real_index), 0.5 * (Global.problem.upperbound(real_index) - Global.problem.lowerbound(real_index)), [1, NP - 1]);
end
xl = repmat(Global.problem.lowerbound(dims), NP - 1, 1);% check the lower bound
less_pos = P < xl;
xu = repmat(Global.problem.upperbound(dims), NP - 1, 1);% check the upper bound
greater_pos = P > xu;
mask = logical(less_pos + greater_pos);
R = initPopulation(dims, NP - 1, Global);
P(mask) = R(mask);
P = [c(dims); P];
end