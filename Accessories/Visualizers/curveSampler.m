function sdata = curveSampler(spercentage, percentage, curve)
sdata = NaN(size(spercentage));
for i = 1: numel(spercentage)
    [~, closest_index] = min(abs(percentage - spercentage(i) * ones(size(percentage))));
    sdata(i) = curve(closest_index);
end
end

