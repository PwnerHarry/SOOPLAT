function accuracy = groupingAccuracy(groups, problem)
% problem = feval(problem_str, problem.dimension);
if isempty(problem.idealgroups)
    accuracy = 1;
elseif numel(problem.idealgroups) == 1
    accuracy = 1;
else
    C = zeros(numel(problem.idealgroups), numel(groups));
    nonseparableproblem.dimension = 0;
    for i = 1: numel(problem.idealgroups)
        nonseparableproblem.dimension = nonseparableproblem.dimension + numel(problem.idealgroups{i});
        for j = 1: numel(groups)
            C(i, j) = numel(intersect(groups{j}, problem.idealgroups{i}));
        end
    end
    accuracy = sum(max(C, [], 2)) / nonseparableproblem.dimension;
end
end