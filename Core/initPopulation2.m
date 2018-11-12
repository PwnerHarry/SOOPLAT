function X = initPopulation2(D, NP, lb, ub)
if isrow(lb) && isrow(ub)
    X = ones(NP, 1) * lb + (ones(NP, 1) * (ub - lb)) .* rand(NP, D);
elseif iscolumn(lb) && iscolumn(ub)
    X = lb * ones(1, NP) + ((ub - lb) * ones(1, NP)) .* rand(D, NP);
end
end