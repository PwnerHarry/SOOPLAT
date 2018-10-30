function recordGroupings(varargin)
Global = [];
groups = [];
activeness = [];
separables = [];
arginStrings = {'Global', 'groups', 'activeness', 'separables'};
IsString = find(cellfun(@ischar, varargin));
[~, Loc]  = ismember(varargin(IsString), cellfun(@(S)['-', S], arginStrings, 'UniformOutput', false));
for i = find(Loc ~= 0)
    eval([arginStrings{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
end
if Global.evaluated == 0 %#ok<*BDSCI>
    Global.additionals.groupings.percentage = [];
    Global.additionals.groupings.groups = {};
    Global.additionals.groupings.separables = {};
    Global.additionals.groupings.accuracy = [];
    Global.additionals.activeness = [];
else
    if ~isempty(groups) || ~isempty(separables)
        Global.additionals.groupings.percentage = [Global.additionals.groupings.percentage; Global.evaluated / Global.evaluation];
        Global.additionals.groupings.groups = [Global.additionals.groupings.groups, {groups}];
        Global.additionals.groupings.separables = [Global.additionals.groupings.separables, separables];
        Global.additionals.groupings.accuracy = [Global.additionals.groupings.accuracy, groupingAccuracy(groups, Global.problem)];
        % fprintf('grouping accuracy: %.2f\n', 100 * Global.additionals.groupings.accuracy(end));
    end
    if ~isempty(activeness)
        Global.additionals.activeness = [Global.additionals.activeness; [Global.evaluated / Global.evaluation, mean(activeness)]];
    end
end
end

