% Author: Dr. Zhenyu Yang
% Modified by: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
%
% ------------
% Description:
% ------------
% This files is the entry point for the DECC-D algorithm which is
% a cooperative co-evolutionary DE based on the delta grouping.
% this algorithm also uses a simple scheme to dynamically change the 
% subcomponent sizes during the course of evolution.
% The details of this technique is described in the paper listed in
% the References section.
%
% -----------
% References:
% -----------
% Omidvar, M.,  Li, X. and Yao, X. (2010), "Cooperative Co-evolution
% with Delta Grouping for Large Scale Non-separable Function Optimization",
% in Proceedings of Congress of Evolutionary Computation (CEC 2010), IEEE,
% p.1762 - 1769.
%
% --------
% License:
% --------
% This program is to be used under the terms of the GNU General Public License 
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Mohammad Nabi Omidvar
% e-mail: mn.omidvar AT gmail.com
% Copyright notice: (c) 2013 Mohammad Nabi Omidvar


clear;
% set random seed
rand('state', sum(100*clock));
randn('state', sum(100*clock));

% problem dimension
D = 1000;

% population size
NP = 50;

% number of independent runs
runs = 25;

% number of fitness evaluations
FEs = 3e+6;

% number of generations
Max_Gen = FEs/NP;

% for the benchmark functions initialization
global initial_flag;

myfunc = [1:20];
addpath('benchmark');
addpath('benchmark/datafiles');

for fun=myfunc
    func_num = fun;
    bestval = [];
    for runindex = 1:runs
        % trace the fitness
        filename = sprintf('trace/tracef%02d_%02d.txt', func_num, runindex);
        fid = fopen(filename, 'w');
        
        initial_flag = 0;
        
        % the main step, call runcompe(), see the runcompe.m  
        val = runcompe('benchmark_func', func_num, D, NP, Max_Gen, runindex, fid);
        % bestval = [bestval; val];
        % fclose(fid);
    end
    
    % the best results of each independent run
    filename = sprintf('result/bestf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    % fprintf(fid, '%g\n', bestval);
    % fclose(fid);
    
    % mean
    filename = sprintf('result/meanf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    % fprintf(fid, '%g\n', mean(bestval));
    % fclose(fid);
    
    % std
    filename = sprintf('result/stdf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    % fprintf(fid, '%g\n', std(bestval));
    % fclose(fid);
end

