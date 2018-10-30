function group = grouping(method, Global)
if strcmp(method, 'DG')
    [seps, allgroups] = DG(Global);
    group = [allgroups, seps];
end
end

