% Author: Dr. Zhenyu Yang
% Modified by: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
%
% ------------
% Description:
% ------------
% This file is an implementation of cooperative co-evolution which
% uses SaNSDE algorithm as subcomponent optimizer.
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


function bestval = decc(fname, func_num, dim, Lbound, Ubound, popsize, itermax, runindex, fid);

F_error = 0.000;

% for fitness trace
tracerst = [];

% the initial population
pop = Lbound + rand(popsize, dim) .* (Ubound-Lbound);

val = feval(fname, pop, func_num);
[bestval, ibest] = min(val);
bestmem = pop(ibest, :);
prev_best = bestmem;
prev_best_val = bestval;
oldpop = pop;

% the initial crossover rate for SaNSDE
group = {};
ccm = 0.5;
subdim = 100;
sansde_iter = 1;
Cycle = 0;
iter = 0;
delta = 0;
while (iter < itermax)
	Cycle = Cycle + 1;

    delta = abs(oldpop-pop);
    oldpop = pop;

    group = delta_grouping(dim, subdim, mean(delta));

	group_num = size(group, 2);

	for i = 1:group_num
		oneitermax = sansde_iter;
		if (iter + oneitermax >= itermax)
			oneitermax = itermax - iter;
		end
		if (oneitermax == 0)
			break;
		end

		dim_index = group{i};

        subpop = pop(:, dim_index); 
        subLbound = Lbound(:, dim_index);        
        subUbound = Ubound(:, dim_index);
        
        [subpopnew, bestmemnew, bestvalnew, tracerst, ccm] = sansde(fname, func_num, dim_index, subpop, bestmem, bestval, subLbound, subUbound, oneitermax, ccm);
	
        iter = iter + oneitermax;
        % fprintf(fid, '%g\n', tracerst);
        
        pop(:, dim_index) = subpopnew;
        bestmem = bestmemnew;
		bestval = bestvalnew;
    end
    
	val = feval(fname, pop, func_num);
	[best, ibest] = min(val);
	if (best < bestval)
		bestval = best;
		bestmem = pop(ibest, :);
	end
   
    fprintf(1, 'fun = %d, run = %d, Cycle = %d, Dim = %d, Subdim = %d\n', func_num, runindex, Cycle, dim, subdim);
    fprintf(1, 'popsize = %d, bestval = %g\n\n', popsize, bestval);        

end

