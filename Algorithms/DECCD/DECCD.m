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
function DECCD(Global)
% the initial population
group_size = 100;
NP = 50;
pop = initPopulation([], NP, Global);
oldpop = pop;
Global.evaluate(pop);
% the initial crossover rate for SaNSDE
ccm = 0.5;
Cycle = 0;
while ~Global.terminated
    Cycle = Cycle + 1;
    fprintf('DECC-D: CYCLE %d\n', Cycle);
    delta = abs(oldpop-pop);
    oldpop = pop;
    group = delta_grouping(Global.problem.dimension, group_size, mean(delta));
    for i = 1: numel(group)
        dims = group{i};
        lb = Global.problem.lowerbound(dims);
        ub = Global.problem.upperbound(dims);
        OBJFUNC = @(X) Global.evaluate(combine(X, Global.bestIndividual, dims));
        [pop(:, dims), ccm] = SaNSDE2('-ub', ub, '-lb', lb, '-objfunc', OBJFUNC, '-ccm', ccm, '-NP', NP, '-initPop', pop(:, dims), '-xbest', Global.bestIndividual(dims), '-fbest', Global.bestFitness, '-maxGen', 1);
    end
    renderCurve(Global);
    [~, xbest_index] = ismember(Global.bestIndividual, pop, 'rows');
    evaluate_index = setdiff(1: NP, xbest_index);
    Global.evaluate(pop(evaluate_index, :));
end