
% set the search range, and call the decc algorithm

function bestval = runcompe(fname, func_num, D, NP, Max_Gen, runindex, fid);

    if(ismember(func_num, [1 4 7:9 12:14 17:20]))
      lb = -100 * ones(NP, D);
      ub = 100 * ones(NP, D);
   elseif (ismember(func_num, [2 5 10 15]))
      lb = -5 * ones(NP, D);
      ub = 5 * ones(NP, D);
   elseif (ismember(func_num,[3 6 11 16]))
      lb = -32 * ones(NP, D);
      ub = 32 * ones(NP, D);
   end

    % the main step, call decc(), see the decc.m
    bestval = decc(fname, func_num, D, lb, ub, NP, Max_Gen, runindex, fid);

end
