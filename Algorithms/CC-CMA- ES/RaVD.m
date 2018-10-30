function subInfo = RaVD(D, m, ~)
I = randperm(D);
s = ceil(D / m);
subInfo = {};
for i = 1: s
    subInfo = [subInfo, I((i - 1) * m + 1: min(i * m, D))];
end
end