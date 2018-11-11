function directOptimize(varargin)
proStr = {'Global', 'optimizer'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
GNP = 50;
X = initPopulation(1: Global.problem.dimension, GNP, Global);
OBJFUNC = @(X) Global.evaluate(X);
fit = feval(OBJFUNC, X);
if strcmp(optimizer, 'SHADE')
	LSHADE('-monitor_flag', true, '-objfunc', OBJFUNC, '-lb', Global.problem.lowerbound, '-ub', Global.problem.upperbound,'-HM', [], '-fbest', Global.bestFitness, '-xbest', Global.bestIndividual, '-initPop', X, '-initFit', fit, '-maxFEs', Global.evaluation - Global.evaluated, '-NP', GNP);
else
    error('directOptimize unfinished');
end
end