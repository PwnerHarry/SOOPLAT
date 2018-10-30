function [ranks, p, chi2] = easyFriedman(X)
NA = size(X, 2) / 2;
X(:, 2 * (1: NA)) = [];
[p, tbl, stats] = friedman(X, 1, 'off');
ranks = stats.meanranks;
chi2 = tbl{2, 5};
end