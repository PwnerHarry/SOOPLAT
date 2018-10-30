classdef BOLTZADA < handle
    properties
        actions
        pointer
        matrix
    end
    methods
        function obj = BOLTZADA(actions, layers)
            obj.actions = actions;
            obj.matrix = zeros(numel(actions), layers);
            obj.pointer = 1;
        end
        function choice = decide(obj)
            for i = 1: numel(obj.actions)
                if sum(obj.matrix(i, :)) == 0 % size(obj.matrix, 2)
                    choice = obj.actions{i};
                    return;
                end
            end
            temp_matrix = obj.matrix;
            % temp_matrix(temp_matrix == -1) = 0;
            mask = find(obj.matrix ~= 0);
            temp_matrix(mask) = mapminmax(temp_matrix(mask)', 0.01, 0.99);
            chance = softmax(5 * mean(temp_matrix, 2));
            fprintf('Adapter Matrix:\t\t\t\t\t\tChance:\n')
            disp([temp_matrix, chance]);
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
            obj.matrix(:, obj.pointer) = 0;
            if value == 0
                value = NaN;
            end
            obj.matrix(actnum, obj.pointer) = value;
            obj.pointer = mod(obj.pointer, size(obj.matrix, 2)) + 1;
        end
    end
end
