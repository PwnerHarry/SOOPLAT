function subInfo = MaVD(D, m, C)
dC = diag(C);
s = ceil(D / m);
subInfo = {};
[~, sortedIndex] = sort(dC, 'ascend');
for i = 1: s
    subInfo = [subInfo, sort(sortedIndex(i: s: D), 'ascend')];
end
end