function gatheredData = gatherData(varargin)
% gimme a lot of loaded dataframes, I will collect their information
algorithms = {};
functions = {};
keep_index = [];
for i = 1: nargin
    if ~isempty(varargin{i})
        keep_index = [keep_index, i];
    end
end
Varargin = {varargin{keep_index}};
Nargin = numel(Varargin);
for i = 1: Nargin
    dataFrames(i).obj = Varargin{i};
    if isempty(dataFrames(i).obj.algorithm)
        algoname = 'unknown';
    else
        algoname = dataFrames(i).obj.algorithm;
    end
    loc = find(ismember(algorithms, algoname));
    if numel(loc) > 0
        dataFrames(i).algoIndex = loc;
    else
        algorithms = [algorithms, algoname];
        dataFrames(i).algoIndex = numel(algorithms);
    end
    func = [func2str(dataFrames(i).obj.suite), '_F', num2str(dataFrames(i).obj.func)];
    loc = find(ismember(functions, func));
    if numel(loc) > 0
        dataFrames(i).funcIndex = loc;
    else
        functions = [functions, func];
        dataFrames(i).funcIndex = numel(functions);
    end
end
gatheredData.algorithm = [];
gatheredData.data = [];
gatheredData.runtime = [];
for i = 1: numel(algorithms)
    for j = 1: numel(functions)
        gatheredData(i, j).algorithm = algorithms{i};
        gatheredData(i, j).func = functions{j};
    end
end
for i = 1: Nargin
    if isempty(dataFrames(i).obj.algorithm)
        algoname = 'unknown';
    else
        algoname = dataFrames(i).obj.algorithm;
    end
    loc1 = find(ismember(algorithms, algoname));
    func = [func2str(dataFrames(i).obj.suite), '_F', num2str(dataFrames(i).obj.func)];
    loc2 = find(ismember(functions, func));
    gatheredData(loc1, loc2).data = [gatheredData(loc1, loc2).data, dataFrames(i).obj.bestFitness];
    gatheredData(loc1, loc2).runtime = [gatheredData(loc1, loc2).runtime, seconds(dataFrames(i).obj.runtime.total)];
end
end