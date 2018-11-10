function P = normPopulation2(NP, lb, ub, mu, C)
P = mvnrnd(mu, C, NP - 1);
P = trimPopulation([mu; P], lb, ub, 'random');
end