function cratio = LS(Global)
global LS_RESTART_FLAG BICCA2_LS_STAGNANT_FLAG;
persistent FEs LS_RESTART_COUNTER SR Improved ofLocalSearch stagnant_threshold;
if LS_RESTART_FLAG || isempty(SR)
    SR = (0.3 + 0.1 * rand()) * mean(Global.problem.upperbound - Global.problem.lowerbound);
    Improved = true;
    stagnant_threshold = floor(log2(SR / 1e-15));
    BICCA2_LS_STAGNANT_FLAG = false;
    LS_RESTART_COUNTER = 0;
    FEs = 25 * Global.problem.dimension;
    LS_RESTART_FLAG = false;
end
fitness_before_LS = Global.bestFitness;
fprintf('LS: ');
improvement = zeros(1, Global.problem.dimension);
stagnant_counter = 0;
if BICCA2_LS_STAGNANT_FLAG == true
    fprintf('LS escaped from stagnant\n');
    cratio = 0;
    return;
end
startFEs = Global.evaluated;
while Global.evaluated - startFEs < FEs
    [new_improvement, ~, ~, SR, Improved] = LocalSearch1(Global.bestIndividual, Global.bestFitness, SR, Improved, improvement, Global);
    if ~Improved
        stagnant_counter = stagnant_counter + 1;
        if stagnant_counter >= stagnant_threshold
            fprintf('LS escaped from stagnant\n');
            BICCA2_LS_STAGNANT_FLAG = true;
            break;
        end
    else
        improvement = new_improvement;
        stagnant_counter = 0;
    end
end
fitness_after_LS = Global.bestFitness;
cratio = min(1, abs((fitness_before_LS - fitness_after_LS) / fitness_before_LS));
if cratio < 0.05
    LS_RESTART_COUNTER = LS_RESTART_COUNTER + 1;
    if LS_RESTART_COUNTER >= 2
        LS_RESTART_COUNTER = 0;
        fprintf('BICCA2: LS RESTART\n');
        LS_RESTART_FLAG = true;
    end
else
    LS_RESTART_COUNTER = 0;
    LS_RESTART_FLAG = false;
end
fprintf('CRATIO = %.3f%%\n', 100 * cratio);
hold on;
fill([0, 1, 1, 0], [fitness_after_LS, fitness_after_LS, fitness_before_LS, fitness_before_LS], [0.8, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', '0.2');
drawnow;
hold off;
end
function [improvement, xk, yk, SR, improve] = LocalSearch1(xk, yk, SR, improve, improvement, Global)
if ~improve
    SR = SR / 2;
    if SR < 1e-15
        SR = (0.3 + 0.1 * rand()) *  (mean(Global.problem.upperbound) - mean(Global.problem.lowerbound));
    end
end
improve = false;
cbf = Global.bestFitness;
[~, I] = sort(improvement, 'descend');
for i = I
    oldxk = xk;
    oldyk = yk;
    xk(i) = truncate2boundary(xk(i) - SR, i, Global);
    if isequal(xk, oldxk)
        yk = oldyk;
    else
        yk = Global.evaluate(xk);
    end
    if yk < cbf
        improvement(i) = yk - oldyk;
        cbf = yk;
    end
    if yk == oldyk
        xk = oldxk;
    else
        if yk > oldyk
            xk = oldxk;
            xk(i) = truncate2boundary(xk(i) + 0.5 * SR, i, Global);
            if isequal(xk, oldxk)
                yk = oldyk;
            else
                yk = Global.evaluate(xk);
            end
            if yk < cbf
                improvement(i) = yk - oldyk;
                cbf = yk;
            end
            if yk >= oldyk
                xk = oldxk;
                yk = oldyk;
            else
                improvement(i) = yk - oldyk;
                improve = true;
            end
        else
            improvement(i) = yk - oldyk;
            improve = true;
        end
    end
end
renderCurve(Global);
end
function x = truncate2boundary(x, i, Global)
span = Global.problem.upperbound(i) - Global.problem.lowerbound(i);
if x > Global.problem.upperbound(i)
    cycles = floor((x - Global.problem.lowerbound(i)) / span);
    x = x - cycles * span;
end
if x < Global.problem.lowerbound(i)
    cycles = floor((Global.problem.upperbound(i) - x) / span);
    x = x + cycles * span;
end
if x > Global.problem.upperbound(i) || x < Global.problem.lowerbound(i)
    x = span * rand() + Global.problem.lowerbound(i);
end
end