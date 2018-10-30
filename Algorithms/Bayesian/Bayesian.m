function Bayesian(Global)
UseParallel = false;
AcquisitionFunctionName = 'expected-improvement-per-second-plus';
IsObjectiveDeterministic = true;
fhd = @Global.evaluate;
D = Global.problem.dimension;
varstr = '';
for d = 1: D
    eval(['v', num2str(d), ' = optimizableVariable(''x', num2str(d), ''', [Global.problem.lowerbound(' , num2str(d), '), Global.problem.upperbound(', num2str(d), ')]);']);
    if length(varstr) > 1
        varstr = [varstr, ', ', 'v',  num2str(d)];
    else
        varstr = ['v',  num2str(d)];
    end
end
eval(['vars = [', varstr, ']']);
bayesopt(fhd, vars, 'ExplorationRatio', 0.1, 'GPActiveSetSize', 2000, 'NumSeedPoints', 100 * Global.problem.dimension, 'MaxObjectiveEvaluations', Global.evaluation - Global.evaluated, 'IsObjectiveDeterministic', IsObjectiveDeterministic, 'AcquisitionFunctionName', AcquisitionFunctionName, 'UseParallel', UseParallel);
end