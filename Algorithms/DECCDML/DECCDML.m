% Author: Dr. Zhenyu Yang
% Modified by: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
% Description:
% This file is an implementation of cooperative co-evolution which
% uses SaNSDE algorithm as subcomponent optimizer.
% References:
% Omidvar, M.,  Li, X. and Yao, X. (2010), "Cooperative Co-evolution
% with Delta Grouping for Large Scale Non-separable Function Optimization",
% in Proceedings of Congress of Evolutionary Computation (CEC 2010), IEEE,
% p.1762 - 1769.
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Mohammad Nabi Omidvar
% e-mail: mn.omidvar AT gmail.com
% Copyright notice: (c) 2013 Mohammad Nabi Omidvar
function DECCDML(Global)
% the initial population
NP = 50;
pop = initPopulation([], NP, Global);
oldpop = pop;
Global.evaluate(pop);
%[bestval, ibest] = min(val);
%bestmem = pop(ibest, :);
%prev_best = bestmem;
prev_best_val = Global.bestFitness;
%oldpop = pop;
% the initial crossover rate for SaNSDE
ccm = 0.5;
dims = [5 10 25 50 100];
group_size = 100;
Cycle = 0;
while ~Global.terminated
    Cycle = Cycle + 1;
    fprintf('DECC-DML: CYCLE %d\n', Cycle);
    if Global.bestFitness == prev_best_val
        group_size = dims(randi(5));
    end
    prev_best_val = Global.bestFitness;
    delta = abs(oldpop-pop);
    oldpop = pop;
    group = delta_grouping(Global.problem.dimension, group_size, mean(delta));
    for i = 1: numel(group)
        dims = group{i};
        lb = Global.problem.lowerbound(dims);
        ub = Global.problem.upperbound(dims);
        OBJFUNC = @(X) Global.evaluate(combine(X, Global.bestIndividual, dims));
        [pop(:, dims), ccm] = SaNSDE2('-ub', ub, '-lb', lb, '-objfunc', OBJFUNC, '-ccm', ccm, '-NP', NP, '-initPop', pop(:, dims), '-xbest', Global.bestIndividual(dims), '-fbest', Global.bestFitness, '-maxGen', 1);
        renderCurve(Global);
    end
    Global.evaluate(pop);
end
end

function group = delta_grouping(dim, subdim, delta);
[S, dim_rand] = sort(delta);
group = {};
for i = 1:subdim:dim
    index = dim_rand(i:min(dim, i + subdim - 1));
    group = {group{1:end} index};
end
end
