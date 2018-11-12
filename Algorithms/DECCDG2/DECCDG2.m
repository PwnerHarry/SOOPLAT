function DECCDG2(varargin)
if nargin == 1
    Global = varargin{1};
    NP = 50;
    optimizer = 'SaNSDE';
else
    arginStrings = {'Global', 'NP', 'optimizer'};
    IsString = find(cellfun(@ischar, varargin));
    [~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], arginStrings, 'UniformOutput', false));
    for i = find(Loc ~= 0)
        eval([arginStrings{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
    end
    if ~exist('Global', 'var')
        error('no global');
    end
    if ~exist('NP', 'var')
        NP = 50;
    end
    if ~exist('optimizer', 'var')
        optimizer = 'SaNSDE';
    end
end

% the initial population
pop = initPopulation([], NP, Global);
Global.evaluate(pop);
% the initial crossover rate for SaNSDE
ccm = 0.5;
maxgen = 100;
Cycle = 0;
[group, seps] = DG2(Global);
if ~isempty(seps)
    group = [group, seps];
end
while ~Global.terminated
    Cycle = Cycle + 1;
    fprintf('DECC-DG2: CYCLE %d\n', Cycle);
	for i = 1: numel(group)
        dims = group{i};
        lb = Global.problem.lowerbound(dims);
        ub = Global.problem.upperbound(dims);
        OBJFUNC = @(X) Global.evaluate(combine(X, Global.bestIndividual, dims));
        [pop(:, dims), ccm] = SaNSDE2('-ub', ub, '-lb', lb, '-objfunc', OBJFUNC, '-ccm', ccm, '-NP', NP, '-initPop', pop(:, dims), '-xbest', Global.bestIndividual(dims), '-fbest', Global.bestFitness, '-maxGen', maxgen);
        renderCurve(Global);
    end
end
end