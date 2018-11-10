classdef INDIVIDUAL < handle
    properties (SetAccess=?GLOBAL, GetAccess=public)
        solution % solution in the solution space
        fitness % fitness value
        evaluated
    end
    methods
        function obj = INDIVIDUAL(varargin)
            if nargin == 1
                obj.solution = varargin{1};
                obj.fitness = NaN;
                obj.evaluated = false;
            else
                error('wrong varargin');
            end
        end
    end
end
