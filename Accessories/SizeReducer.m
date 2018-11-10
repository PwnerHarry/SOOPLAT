MAXPOINTS = 1001;
location = 'Results/LSO08/D10000/BICCA2';
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
    G.reduceTrace(MAXPOINTS);
    save(fullfile(location, filename), 'G');
    clear G;
end