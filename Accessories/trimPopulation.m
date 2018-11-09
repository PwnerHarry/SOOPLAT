function X = trimPopulation(X, lb, ub)
if iscolumn(lb) && iscolumn(ub)
    NP = size(X, 2);
    MIN = repmat(lb, 1, NP);
    MAX = repmat(ub, 1, NP);
elseif isrow(lb) && isrow(ub)
    NP = size(X, 1);
    MIN = repmat(lb, NP, 1);
    MAX = repmat(ub, NP, 1);
end
overflow_ub_index = X > MAX;
overflow_lb_index = X < MIN;
X(overflow_ub_index) = MAX(overflow_ub_index);
X(overflow_lb_index) = MIN(overflow_lb_index);
end