function directOptimize(varargin)
proStr = {'Global', 'optimizer'};
IsString = find(cellfun(@ischar, varargin(1: end - 1)));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], proStr, 'UniformOutput', false));
for i = 1: numel(IsString)
    eval([proStr{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if strcmp(optimizer, 'LSHADE')
    feval(optimizer, '-initNP', 1800, '-maxFEs', Global.evaluation, '-dims', 1: Global.problem.dimension,'-minNP', 4, '-Global', Global);
end
end

