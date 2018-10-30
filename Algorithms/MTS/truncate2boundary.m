function x = truncate2boundary(x, i, Global)
span = Global.problem.upperbound(i) - Global.problem.lowerbound(i);
if x > Global.problem.upperbound(i)
    cycles = floor((x - Global.problem.lowerbound(i)) / span);
    x = x - cycles * span;
end
if x < Global.problem.lowerbound(i)
    cycles = floor((Global.problem.upperbound(i) - x) / span);
    x = x + cycles * span;
end
if x > Global.problem.upperbound(i) || x < Global.problem.lowerbound(i)
    x = span * rand() + Global.problem.lowerbound(i);
end
end