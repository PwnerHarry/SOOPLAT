classdef PROBLEM < handle
    properties
        name
        dimension
        lowerbound
        upperbound
        functionhandle
        idealfitness
        idealsolution
        idealgroups
        idealseparables
    end
    methods
        function obj = PROBLEM()
            obj.name = 'UNDEFINED';
            obj.idealfitness = 0;
            obj.idealsolution = [];
            obj.dimension = NaN;
            obj.lowerbound = [];
            obj.upperbound = [];
            obj.functionhandle = [];
            obj.idealgroups = {};
            obj.idealseparables = [];
        end
        function shiftedX = shift(obj, X)
            if istable(X)
                X = table2array(X);
            end
            shiftedX = X - repmat(obj.idealsolution, size(X, 1), 1);
        end
        function E = geteps(obj)
            D = length(obj.idealsolution);
            E = 0;
            for d = 1: D
                x = obj.idealsolution;
                x(d) = x(d) + eps;
                E = max(E, obj.functionhandle(x));
                x = obj.idealsolution;
                x(d) = x(d) - eps;
                E = max(E, obj.functionhandle(x));
            end
            x = obj.idealsolution;
            x = x + eps;
            E = max(E, obj.functionhandle(x));
            x = obj.idealsolution;
            x = x - eps;
            E = max(E, obj.functionhandle(x));
            E = max(E, eps);
        end
    end
end