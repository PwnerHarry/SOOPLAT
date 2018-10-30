function [spercentage, sdata] = collectData(samples, fileFolder, algorithm, problem, dimension, evaluation)
dirOutput = dir(fullfile(fileFolder, '*'));
filenames = {dirOutput.name}';
reduce_index = [];
for i = 1: numel(filenames)
    expr = sprintf('%s_%s_D%d.mat', algorithm, problem, dimension);
    if isempty(strfind(filenames{i}, expr))
        reduce_index = [reduce_index, i];
    end
end
filenames(reduce_index) = [];
sdata = NaN(numel(filenames), samples);
spercentage = linspace(0, 1, samples);
for i = 1: numel(filenames)
	load(fullfile(fileFolder, filenames{i}));
    sdata(i, :) = curveSampler(spercentage, G.trace(:, 1), G.trace(:, 2));
    clear G;
end
end

