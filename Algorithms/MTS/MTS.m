function MTS(Global)
M = 5;
ofForeground = 3;
ofLocalSearchTest = 3;
ofLocalSearch = 10;
ofLocalSearchBest = 150;
BONUS1 = 10;
BONUS2 = 1;
N = Global.problem.dimension;
[~, SOA] = sort(rand(M, N));
X = SOA ./ M .* repmat(Global.problem.upperbound-Global.problem.lowerbound, M, 1) + repmat(Global.problem.lowerbound, M, 1);
y = Global.evaluate(X);
Enable = true(M, 1);
Improve = true(M, 1);
GradeX = NaN(M, 1);
SR = ones(M, 1) * 0.5 * (mean(Global.problem.upperbound) - mean(Global.problem.lowerbound));
while ~Global.terminated
    for i = 1: M
        if Enable(i)
            GradeX(i) = 0;
            TestGrades = zeros(3, 1);
            for j = 1: ofLocalSearchTest
                [TestGrades(1), X(i, :), y(i), SR(i), Improve(i)] = LocalSearch1(X(i, :), y(i), SR(i), Improve(i), TestGrades(1), BONUS1, BONUS2, Global);
                [TestGrades(2), X(i, :), y(i), SR(i), Improve(i)] = LocalSearch2(X(i, :), y(i), SR(i), Improve(i), TestGrades(2), BONUS1, BONUS2, Global);
                [TestGrades(3), X(i, :), y(i), SR(i), Improve(i)] = LocalSearch3(X(i, :), y(i), SR(i), Improve(i), TestGrades(3), BONUS1, BONUS2, Global);
            end
            [~, K] = max(TestGrades);
            LocalSearchK = sprintf('LocalSearch%d', K);
            for j = 1: ofLocalSearch
                [GradeX(i), X(i, :), y(i), SR(i), Improve(i)] = feval(LocalSearchK, X(i, :), y(i), SR(i), Improve(i), GradeX(i), BONUS1, BONUS2, Global);
            end
        end
    end
    [~, best_solution_index] = min(y);
    for i = 1: ofLocalSearchBest
        [~, X(best_solution_index, :), y(best_solution_index), SR(best_solution_index), Improve(best_solution_index)] = LocalSearch1(X(best_solution_index, :), y(best_solution_index), SR(best_solution_index), Improve(best_solution_index), [], BONUS1, BONUS2, Global);
    end
    Enable = false(M, 1);
    [~, ind] = sort(GradeX, 'descend');
    Enable(ind(1: ofForeground)) = true;
end
renderCurve(Global);
end