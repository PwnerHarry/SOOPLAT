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
    end
end