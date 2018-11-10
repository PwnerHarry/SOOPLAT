function varargout = friedmanTest(gatheredData)
algorithms = {gatheredData(:, 1).algorithm};
functions = {gatheredData(1, :).func};
if numel(functions) == 1
    error('Friedman test should not be conducted on single test case');
end
fullData = nan(numel(functions), numel(algorithms));
trim_flag = false;
for i = 1: numel(algorithms)
    for j = 1: numel(functions)
        if isempty(gatheredData(i, j).data)
            trim_flag = true;
        else
            fullData(j, i) = mean(gatheredData(i, j).data);
        end
    end
end
participantTable = array2table(fullData, 'variableNames', algorithms, 'rowNames', functions);
fullTable = array2table(fullData, 'variableNames', algorithms, 'rowNames', functions);
while trim_flag
    warning('trimming initiated, please check participantTable and fullTable');
    trim_flag = false;
    remove_function_index = [];
    for i = 1: numel(functions)
       if sum(~isnan(participantTable{i, :})) <= 1 
           remove_function_index = [remove_function_index, i];
           trim_flag = true;
       end
    end
    participantTable(remove_function_index, :) = [];
    functions(remove_function_index) = [];
    remove_algorithm_index = [];
    for j = 1: numel(algorithms)
       if sum(~isnan(participantTable{:, j})) <= 1 
           remove_algorithm_index = [remove_algorithm_index, j];
           trim_flag = true;
       end
    end
    participantTable(:, remove_algorithm_index) = [];
    algorithms(remove_algorithm_index) = [];
end
[p, tbl, stats] = friedman(table2array(participantTable), 1, 'off');
Chisq = tbl{2, 5};
rankTable = array2table(stats.meanranks, 'variableNames', algorithms, 'rowNames', {'meanrank'});
if nargout == 2
    varargout = {rankTable, participantTable};
elseif nargout == 3
    varargout = {rankTable, participantTable, fullTable};
elseif nargout == 5
    varargout = {rankTable, participantTable, fullTable, Chisq, p};
end
end
