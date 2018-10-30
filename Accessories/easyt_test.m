function strs = easyt_test(X)
NA = size(X, 2) / 2;
strs = {};
for i = 1: NA - 1
    strs{i} = t_test(X(:, 1), X(:, 2), X(:, 2 * i + 1), X(:, 2 * i + 2));
end
end