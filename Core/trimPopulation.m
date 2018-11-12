function X = trimPopulation(varargin)
% X = trimPopulation(X, lb, ub)
% X = trimPopulation(X, lb, ub, mode)
% This function trims the population within a bounding box defined by lb and ub
% Note that if lb and ub are column vectors, the individuals in X should also be column vectors and vice versa.
% mode = 'bound': (default)
% pull the values back to the boundary
% mode = 'random':
% replace the violated values to somewhere random between lb and ub
if nargin == 3
    X = varargin{1};
    lb = varargin{2};
    ub = varargin{3};
    mode = 'bound';
elseif nargin == 4
    X = varargin{1};
    lb = varargin{2};
    ub = varargin{3};
    mode = varargin{4};
end
if iscolumn(lb) && iscolumn(ub)
    NP = size(X, 2);
    MIN = repmat(lb, 1, NP);
    MAX = repmat(ub, 1, NP);
elseif isrow(lb) && isrow(ub)
    NP = size(X, 1);
    MIN = repmat(lb, NP, 1);
    MAX = repmat(ub, NP, 1);
end
overflow_ub_index = X > MAX;
overflow_lb_index = X < MIN;
if strcmp(mode, 'bound')
    X(overflow_ub_index) = MAX(overflow_ub_index);
    X(overflow_lb_index) = MIN(overflow_lb_index);
elseif strcmp(mode, 'mod')
    X(overflow_lb_index) = MIN(overflow_lb_index) + mod((MIN(overflow_lb_index) - X(overflow_lb_index)), (MAX(overflow_lb_index) - MIN(overflow_lb_index)));
    X(overflow_ub_index) = MAX(overflow_ub_index) - mod((X(overflow_ub_index) - MAX(overflow_ub_index)), (MAX(overflow_ub_index) - MIN(overflow_ub_index)));
elseif strcmp(mode, 'random')
    if isrow(lb) && isrow(ub)
        R = ones(NP, 1) * lb + (ones(NP, 1) * (ub - lb)) .* rand(size(X));
    elseif iscolumn(lb) && iscolumn(ub)
        R = lb * ones(1, NP) + ((ub - lb) * ones(1, NP)) .* rand(size(X));
    end
    X(overflow_ub_index) = R(overflow_ub_index);
    X(overflow_lb_index) = R(overflow_lb_index);
end
end