classdef PATTERN < handle
    properties
        core
        success
        fail
        scatters
        compactness
        additional
    end
    methods
        function obj = PATTERN(varargin)
            obj.additional = [];
            obj.success = 0;
            obj.fail = 0;
            obj.compactness = NaN;
            obj.scatters = [];
            if nargin == 1
                obj.core = varargin{1};
            elseif nargin == 3
                obj.core = varargin{1};
                obj.success = varargin{2};
                obj.fail = varargin{3};
            end
        end
        function a = activeness(obj)
            if obj.success > 0 || obj.fail > 0
                a = obj.success / (obj.success + obj.fail);
            else
                a = NaN;
            end
        end
    end
end

