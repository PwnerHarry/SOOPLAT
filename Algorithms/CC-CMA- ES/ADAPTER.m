classdef ADAPTER < handle
    properties
        actions
        pointer
        matrix
    end
    
    methods
        function obj = ADAPTER(actions, layers)
            obj.actions = actions;
            obj.matrix = ones(numel(actions), layers);
            obj.pointer = 1;
        end
        
        function choice = decide(obj)
            chance = softmax(sum(obj.matrix, 2));
            choicenum = rand();
            S = 0;
            k = 1;
            while k <= length(chance)
                S = S + chance(k);
                if choicenum <= S
                    choice = obj.actions{k};
                    return;
                end
                k = k + 1;
            end
        end
        
        function update(obj, choice, value)
            actnum = ismember(obj.actions, choice);
            obj.matrix(:, obj.pointer) = 1;
            obj.matrix(actnum, obj.pointer) = value;
            obj.pointer = mod(obj.pointer, size(obj.matrix, 2)) + 1;
        end
    end
end

