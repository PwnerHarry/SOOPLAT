% Author: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
%
% ------------
% Description:
% ------------
% This files is the entry point for running the differential gropuing algorithm.

clear all;

% Specify the functions that you want the differential grouping algorithm to 
% identify its underlying grouping structure.
func = [20:-1:1];
for i=func
    func_num = i

    t1 = [1 4 7 8 9 12 13 14 17 18 19 20];
    t2 = [2 5 10 15];
    t3 = [3 6 11 16];

    if (ismember(func_num, t1))
        lb = -100;
        ub = 100;
    elseif (ismember(func_num, t2))
        lb = -5;
        ub = 5;
    elseif (ismember(func_num, t3))
        lb = -32;
        ub = 32;
    end

    opts.lbound  = lb;
    opts.ubound  = ub;
    opts.dim     = 1000;
    opts.epsilon = 1e-3;

    addpath('cec2010');
    addpath('cec2010/datafiles');
    global initial_flag;
    initial_flag = 0;

    [seps, nonseps, FEs] = dg('benchmark_func', func_num, opts);

    filename = sprintf('./results/F%02d', func_num);
    save (filename, 'seps', 'nonseps', 'FEs', '-v7');
end

