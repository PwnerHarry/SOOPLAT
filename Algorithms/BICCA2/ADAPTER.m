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
            for i = 1: numel(obj.actions)
                if sum(obj.matrix(i, :)) == size(obj.matrix, 2)
                    choice = obj.actions{i};
                    return;
                end
            end
            mask = find(obj.matrix ~= 1);
            temp_matrix = obj.matrix;
            temp_matrix(mask) = mapminmax(temp_matrix(mask)', 0.01, 0.99);
            chance = softmax(sum(temp_matrix, 2));
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
