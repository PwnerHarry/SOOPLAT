% Author: Mohammad Nabi Omidvar
% email : mn.omidvar AT gmail.com
% ------------
% Description:
% ------------
% dg - This function runs the differential grouping
%      procedure to identify the non-separable groups
%      in cec'2010 benchmark problems.o
% -------
% Inputs:
% -------
%    fun        : the function suite for which the interaction structure
%                 is going to be identified in this case benchmark_func
%                 of cec'2010.
%
%    fun_number : the function number.
%
%    options    : this variable contains the options such as problem
%                 dimensionality, upper and lower bounds and the parameter
%                 epsilon used by differential grouping.
%    input3 - Description
% --------
% Outputs:
% --------
%    sep      : a vector of all separable variables.
%    allgroups: a cell array containing all non-separable groups.
%    FEs      : the total number of fitness evaluations used.
% --------
% License:
% --------
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Mohammad Nabi Omidvar
% e-mail: mn.omidvar AT gmail.com
% Copyright notice: (c) 2013 Mohammad Nabi Omidvar
function [seps, allgroups] = DG(Global)
% dim       = options.dim;
dim = Global.problem.dimension;
epsilon = 1e-3;
% r         = ub - lb;
dims = 1: dim;
seps  = [];
allgroups = {};
while length(dims) >= 1
    renderCurve(Global);
    n = length(dims);
    group = dims(1);
    group_ind = 1;
    p1 = Global.problem.lowerbound;
    p2 = p1;
    p2(dims(1)) = Global.problem.lowerbound(dims(1));
    % delta1 = feval(fun, p1, fun_number) - feval(fun, p2, fun_number);
    delta1 = Global.evaluate(p1) - Global.evaluate(p2);
    for i = 2: n
        p3 = p1;
        p4 = p2;
        temp = 0;
        p3(dims(i)) = temp;
        p4(dims(i)) = temp;
        % delta2 = feval(fun, p3, fun_number) - feval(fun, p4, fun_number);
        delta2 = Global.evaluate(p3) - Global.evaluate(p4);
        if abs(delta1 - delta2) > epsilon
            group = [group; dims(i)];
            group_ind = [group_ind; i];
        end
    end
    if length(group) == 1
        seps = [seps; group];
    else
        allgroups = {allgroups{1:end}, group};
    end
    if ~isempty(dims)
        dims(group_ind) = [];
    end
end
end