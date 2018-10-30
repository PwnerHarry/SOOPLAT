function subInfo = MiVD(D, m, C)
dC = diag(C);
s = ceil(D / m);
subInfo = {};
[~, sortedIndex] = sort(dC, 'ascend');
for i = 1: s
    subInfo = [subInfo, sort(sortedIndex((i - 1) * m + 1: min(i * m, D)), 'ascend')];
end
end
