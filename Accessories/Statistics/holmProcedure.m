function varargout = holmProcedure(gatheredData)
    [rankTable, participantTable, ~] = friedmanTest(gatheredData);
    meanranks = rankTable{1, :};
    [n, k] = size(table2array(participantTable));
    z = NaN(k, k);
    p = NaN(k, k);
    for i = 1: k
        for j = i + 1: k
            z(i, j) = abs(meanranks(i) - meanranks(j)) * sqrt(6 * n / (k * (k + 1)));
            z(j, i) = z(i, j);
            p(i, j) = normpdf(z(i, j));
            p(j, i) = p(i, j);
        end
    end
    p = array2table(p, 'variableNames', rankTable.Properties.VariableNames, 'rowNames', rankTable.Properties.VariableNames);
    z = array2table(z, 'variableNames', rankTable.Properties.VariableNames, 'rowNames', rankTable.Properties.VariableNames);
    if nargout == 1
        varargout = {p};
    elseif nargout == 2
        varargout = {p, z};
    elseif nargout == 3
        varargout = {p, z, participantTable};
    end
end