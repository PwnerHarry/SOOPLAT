function [grade, x, y, SR, improve] = LocalSearch3(x, y, SR, improve, grade, BONUS1, BONUS2, Global)
% 'SR' and 'improve' will not be updated in local search 3
startFEs = Global.evaluated;
N = Global.problem.dimension;
cbf = Global.bestFitness;
for i = randperm(N)
    oldx = x;
    oldy = y;
    x1 = x;
    x1(i) = truncate2boundary(x1(i) + 0.1, i, Global);
    y1 = Global.evaluate(x1);
    x2 = x;
    x2(i) = truncate2boundary(x2(i) - 0.1, i, Global);
    y2 = Global.evaluate(x2);
    x3 = x;
    x3(i) = truncate2boundary(x3(i) - 0.2, i, Global);
    y3 = Global.evaluate(x3);
    if y1 < cbf
        grade = grade + BONUS1;
        cbf = y1;
    end
    if y2 < cbf
        grade = grade + BONUS1;
        cbf = y2;
    end
    if y3 < cbf
        grade = grade + BONUS1;
        cbf = y3;
    end
    D1 = y - y1;
    if D1 > 0
        grade = grade + BONUS2;
    end
    D2 = y - y2;
    if D2 > 0
        grade = grade + BONUS2;
    end
    D3 = y - y3;
    if D3 > 0
        grade = grade + BONUS2;
    end
    a = 0.3 + 0.1 * rand();
    b = 0.1 + 0.2 * rand();
    c = rand();
    x(i) = truncate2boundary(x(i) + a * (D1 - D2) + b * (D3 - 2 * D1) + c, i, Global);
	if isequal(x, oldx)
        y = oldy;
    else
        y = Global.evaluate(x);
    end
    if y >= oldy
        y = oldy;
        x = oldx;
    else
        grade = grade + BONUS2;
    end % TODO: I don't know if this is right, should be
end
renderCurve(Global);
end