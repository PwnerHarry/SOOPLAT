MAXPOINTS = 1001;
location = 'Results/LSO08/D10000/MTS';
dirOutput = dir(fullfile(location, '*'));
filenames = {dirOutput.name}';
reduce_index = [];
for i = 1: numel(filenames)
    expr = sprintf('.mat');
    if isempty(strfind(filenames{i}, expr))
        reduce_index = [reduce_index, i];
    end
end
filenames(reduce_index) = [];
for i = 1: numel(filenames)
    filename = filenames{i};
	load(fullfile(location, filename));
    try
        G.reduceTrace(MAXPOINTS);
    catch ME
        clear G;
        continue;
    end
    save(fullfile(location, filename), 'G');
    clear G;
end