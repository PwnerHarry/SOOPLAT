function [grade, xk, yk, SR, improve] = LocalSearch2(xk, yk, SR, improve, grade, BONUS1, BONUS2, Global)
startFEs = Global.evaluated;
if ~improve
    SR = SR / 2;
    if SR < 1e-15
        SR = 0.4 * (mean(Global.problem.upperbound) - mean(Global.problem.lowerbound));
    end
end
improve = false;
N = Global.problem.dimension;
cbf = Global.bestFitness;
for l = randperm(N)
    oldxk = xk;
    oldyk = yk;
    D = rand(1, N);
    D(D > 0.5) = 1;
    D(D <= 0.5) = -1;
    r = rand(1, N);
    for i = 1: N
        if r(i) < 1/4
            xk(i) = truncate2boundary(xk(i) - SR * D(i), i, Global);
        end
    end
    if isequal(xk, oldxk)
        yk = oldyk;
    else
        yk = Global.evaluate(xk);
    end
    if yk < cbf
        grade = grade + BONUS1;
        cbf = yk;
    end
    if yk == oldyk
        xk = oldxk;
    else
        if yk > oldyk
            xk = oldxk;
            % yk = oldyk;
            for i = 1: N
                if r(i) < 1/4
                    xk(i) = truncate2boundary(xk(i) + 0.5 * SR * D(i), i, Global);
                end
            end
            if isequal(xk, oldxk)
                yk = oldyk;
            else
                yk = Global.evaluate(xk);
            end
            if yk < cbf
                grade = grade + BONUS1;
                cbf = yk;
            end
            if yk >= oldyk
                xk = oldxk;
                yk = oldyk;
            else
                grade = grade + BONUS2;
                improve = true;
            end
        else
            grade = grade + BONUS2;
            improve = true;
        end
    end
end
fprintf('LS2: %d FEs\n', Global.evaluated - startFEs);
renderCurve(Global);
end