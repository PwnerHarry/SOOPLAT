% ------------
% Description:
% ------------
% This fuction is used to implement delta grouping based on how the
% mean of the population moves in the decision space. The delta values are
% calculated in decc.m are passed to this fuction for forming the groups.
%
%--------
% Inputs:
%--------
%    dim : the dimensions of the problem.
%
%    subdim : the size of each subcomponent.
%
%    delta : population movement information that is calculated in decc.m

function group = delta_grouping(dim, subdim, delta);
   [~, dim_rand] = sort(delta);
   group = {};
   for i = 1: subdim: dim
      index = dim_rand(i: min(dim, i + subdim - 1));
      group = {group{1:end}, index};
   end
end
