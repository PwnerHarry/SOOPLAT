function [DF, AT, RT] = makeStats(folder, algorithm, benchmark, milestones, samplepoints)
data_frame = [];
data_index = 1;
fileFolder = fullfile(folder);
dirOutput = dir(fullfile(fileFolder, '*'));
filenames = {dirOutput.name}';
filenames(1: 2) = [];
fnum = size(filenames, 1);
for i = 1: fnum
    fn = filenames{i};
    id = regexp(fn, '\[\S+\]', 'match');
    id = char(regexp(id{:}, '\d\S+\d', 'match'));
	func = char(regexp(fn, 'F\d+.mat', 'match'));
    func = str2num(func(2: end - 3));
    suite = char(regexp(fn, '\]\S+_', 'match'));
    suite = suite(2: end - 1);
    if ~strcmp(suite, benchmark)
       continue; 
    end
    load([folder, '\', filenames{i}]);
     if ~isempty(obj.algorithm) && ~strcmp(algorithm, obj.algorithm)
        continue;
     end
    data_frame(data_index).id = id;
    data_frame(data_index).suite = suite;
    data_frame(data_index).func = func;
    data_frame(data_index).algorithm = obj.algorithm;
    sp = linspace(0, 1, samplepoints);
    SP = sp * obj.evaluation;
    T = obj.mileStone(SP);
    data_frame(data_index).trace = T(:, 2);
    MS = obj.mileStone(milestones);
    data_frame(data_index).milestone = MS(:, 2);
    data_index = data_index + 1;
end
[~, index] = sort(cat(1, data_frame.func), 'ascend');
data_frame = data_frame(index);
[~, index] = sort(string({data_frame.suite}));
DF = data_frame(index);
funcs = unique(cat(1, data_frame.func));
for i = 1: numel(funcs)
    avg_data(i).func = funcs(i);
    avg_data(i).counter = 0;
    avg_data(i).trace = [];
    avg_data(i).avgtrace = [];
    avg_data(i).meanMS = [];
    avg_data(i).stdMS = [];
    avg_data(i).bestMS = [];
    avg_data(i).worstMS = [];
    avg_data(i).medianMS = [];
    avg_data(i).MS = [];
end
for i = 1: numel(data_frame)
    avg_data((funcs == data_frame(i).func)).trace = [avg_data((funcs == data_frame(i).func)).trace; data_frame(i).trace'];
    avg_data((funcs == data_frame(i).func)).MS = [avg_data((funcs == data_frame(i).func)).MS; data_frame(i).milestone'];
    avg_data((funcs == data_frame(i).func)).counter = avg_data((funcs == data_frame(i).func)).counter + 1;
end
for i = 1: numel(avg_data)
    avg_data(i).bestMS = min(avg_data(i).MS, [], 1);
    avg_data(i).worstMS = max(avg_data(i).MS, [], 1);
    avg_data(i).medianMS = median(avg_data(i).MS, 1);
    %avg_data(i).MS = sortrows(avg_data(i).MS, size(avg_data(i).MS, 2), 'ascend');
    avg_data(i).meanMS = mean(avg_data(i).MS, 1);
    avg_data(i).stdMS = std(avg_data(i).MS);
    %avg_data(i).trace = sortrows(avg_data(i).trace, size(avg_data(i).trace, 2), 'ascend');
    avg_data(i).avgtrace = mean(avg_data(i).trace, 1);
end
colName = cell(1, numel(avg_data));
for i = 1: numel(avg_data)
    colName{i} = sprintf('F%d', avg_data(i).func);
end
rowName = cell(1, numel(sp));
for i = 1: numel(sp)
    rowName{i} = ['P', num2str(100 * sp(i))];
end
AT = array2table(cat(1, avg_data.avgtrace)', 'VariableNames', colName, 'RowNames', rowName);
for i = 1: numel(funcs)
    compact_avg_data(i).func = avg_data(i).func;
    for j = 1: numel(milestones)
        compact_avg_data(i).milestone(j).best = avg_data(i).bestMS(j);
        compact_avg_data(i).milestone(j).median = avg_data(i).medianMS(j);
        compact_avg_data(i).milestone(j).worst = avg_data(i).worstMS(j);
        compact_avg_data(i).milestone(j).mean = avg_data(i).meanMS(j);
        compact_avg_data(i).milestone(j).std = avg_data(i).stdMS(j);
    end
end
result_table = zeros(5 * numel(milestones), numel(funcs));
for i = 1: size(result_table, 1)
    for j = 1: size(result_table, 2)
        ctg = mod(i, 5);
        msp = ceil(i / 5);
        if ctg == 1
            result_table(i, j) = compact_avg_data(j).milestone(msp).best;
        elseif ctg == 2
            result_table(i, j) = compact_avg_data(j).milestone(msp).median;
        elseif ctg == 3
            result_table(i, j) = compact_avg_data(j).milestone(msp).worst;
        elseif ctg == 4
            result_table(i, j) = compact_avg_data(j).milestone(msp).mean;
        elseif ctg == 0
            result_table(i, j) = compact_avg_data(j).milestone(msp).std;
        end
    end
end
colName = cell(1, numel(funcs));
for i = 1: numel(funcs)
    colName{i} = sprintf('F%d', avg_data(i).func);
end
rowName = cell(1, 5 * numel(milestones));
for i = 1: 5 * numel(milestones)
    msp = ceil(i / 5);
    ctg = mod(i, 5);
    if ctg == 1
        rowName{i} = sprintf('%.2e best', milestones(msp));
    elseif ctg == 2
        rowName{i} = sprintf('%.2e median', milestones(msp));
    elseif ctg == 3
        rowName{i} = sprintf('%.2e worst', milestones(msp));
    elseif ctg == 4
        rowName{i} = sprintf('%.2e mean', milestones(msp));
    elseif ctg == 0
        rowName{i} = sprintf('%.2e std', milestones(msp));
    end
end
RT = array2table(result_table, 'VariableNames', colName, 'RowNames', rowName);
writetable(AT, 'AT.xlsx', 'WriteRowNames', true);
writetable(RT, 'RT.xlsx', 'WriteRowNames', true);
end