function [bestmemnew, bestval] = de_weight(bestmem, superval, Lbound, Ubound, popsize, itermax, group, Global)
%function [bestmemnew, bestval, flag] = de_weight(bestmem, superval, Lbound, Ubound, popsize, itermax, group, Global)
%flag = 0;
F = 0.5;
CR = 0.9;
bestmemnew = bestmem;
bestval = Global.evaluate(bestmem); % bestval = feval(fname, bestmem, func_num);
tmpmem = bestmem;
tmpmem(tmpmem == 0) = 1;

dim = size(group, 2);
subLbound = zeros(1, dim);
subUbound = zeros(1, dim);
for i = 1: dim
	tmpUbound = Ubound(group{i}) ./ tmpmem(group{i});
	tmpLbound = Lbound(group{i}) ./ tmpmem(group{i});
	index = find(tmpUbound < tmpLbound);
	tmp = tmpUbound(index);
	tmpUbound(index) = tmpLbound(index);
	tmpLbound(index) = tmp;
	subUbound(i) = min(tmpUbound);
	subLbound(i) = max(tmpLbound);
end
subUbound = ones(popsize, 1) * subUbound;
subLbound = ones(popsize, 1) * subLbound;
index = find(subLbound > subUbound, 1);
if ~isempty(index)
	%flag = 1;
    return;
end
pop = subLbound + rand(popsize, dim) .* (subUbound - subLbound);
randindex = randperm(popsize);
popold(randindex(1), :) = 1;
gpop = getgpop(pop, bestmem, group);
%val = feval(fname, gpop, func_num);
val = Global.evaluate(gpop);
[best, ibest] = min(val);
if (best < bestval)
	bestval = best;
	bestmemnew = gpop(ibest, :);
end
%pm1 = zeros(popsize, dim);
%pm2 = zeros(popsize, dim);
%pm3 = zeros(popsize, dim);
%ui = zeros(popsize, dim);
%mui = zeros(popsize, dim);
%mpo = zeros(popsize, dim);
%a1 = zeros(popsize);
%a2 = zeros(popsize);
%a3 = zeros(popsize);
rot = (0:1:popsize-1);
%ft = zeros(popsize);
%ind = zeros(4);
iter = 0;
while (iter < itermax) && ~Global.terminated
    popold = pop;
    ind = randperm(4);
    a1 = randperm(popsize);
    rt = rem(rot+ind(1), popsize);
    a2 = a1(rt+1);
    rt = rem(rot+ind(2), popsize);
    a3 = a2(rt+1);
    pm1 = popold(a1, :);
    pm2 = popold(a2, :);
    pm3 = popold(a3, :);
    mui = (rand(popsize, dim) < CR);
    mpo = (mui < 0.5);
    ui = pm3 + F*(pm1-pm2);
    ui = popold.*mpo + ui.*mui;
	index = find(ui > subUbound);
	ui(index) = subUbound(index) - mod((ui(index)-subUbound(index)), (subUbound(index)-subLbound(index)));
	index = find(ui < subLbound);
	ui(index) = subLbound(index) + mod((subLbound(index)-ui(index)), (subUbound(index)-subLbound(index)));
    gpopnew = getgpop(ui, bestmem, group);
    tempval = Global.evaluate(gpopnew); % tempval = feval(fname, gpopnew, func_num);
    index = find(tempval <= val);
    val(index) = tempval(index);
    pop(index, :) = ui(index, :);
	gpop(index, :) = gpopnew(index, :);
    [best, ibest] = min(val);
	if best < bestval
		bestval = best;
		bestmemnew = gpop(ibest, :);
    end
    if bestval < superval
        superval = bestval;
    end
    iter = iter + 1;
end
end
function gpop = getgpop(wpop, mem, group)
    dim = size(mem, 2);
    [popsize, wdim] = size(wpop);
    gpop = zeros(popsize, dim);
    for i = 1:wdim
        gpop(:, group{i}) = wpop(:, i) * ones(1, size(group{i}, 2));
    end
    gpop = gpop .* (ones(popsize, 1) * mem);
end