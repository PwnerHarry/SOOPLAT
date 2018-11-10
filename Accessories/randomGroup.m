function groups = randomGroup(dimension, GS)
d = 1: dimension;
d = d(randperm(length(d)));
groups = {};
for i = 1: ceil(dimension / GS)
    begin_index = (i - 1) * GS + 1;
    end_index = min(i * GS, dimension);
    groups{i} = d(begin_index: end_index);
end
end