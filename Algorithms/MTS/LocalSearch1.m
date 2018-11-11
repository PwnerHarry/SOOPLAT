function [grade, xk, yk, SR, improve] = LocalSearch1(xk, yk, SR, improve, grade, BONUS1, BONUS2, Global)
startFEs = Global.evaluated;
if ~improve
    SR = SR / 2;
    if SR < 1e-15 %1e-15
        SR = 0.4 * (mean(Global.problem.upperbound) - mean(Global.problem.lowerbound));
    end
end
improve = false;
N = Global.problem.dimension;
cbf = Global.bestFitness;
for i = randperm(N)
    oldxk = xk;
    oldyk = yk;
    xk(i) = truncate2boundary(xk(i) - SR, i, Global);
    yk = Global.evaluate(xk);
    if yk < cbf
        grade = grade + BONUS1;
        cbf = yk;
    end
    if yk == oldyk
        xk = oldxk;
    else
        if yk > oldyk
            xk = oldxk;
            yk = oldyk;
            xk(i) = truncate2boundary(xk(i) + 0.5 * SR, i, Global);
            yk = Global.evaluate(xk);
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
renderCurve(Global);
end