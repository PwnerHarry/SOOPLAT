function parmain()
% This script is for parallel running
algorithms = {'CSO', 'FRG'};
function_index = 1;
runtimes = 1;
for R = 1: runtimes
    for i = function_index
        for j = 1: numel(algorithms)
            Global(i, j) = GLOBAL('-algorithm', algorithms{j}, '-evaluation', 5e6, '-dimension', 1000, '-problem', sprintf('LSO08F%d', i));
        end
    end
    title(sprintf("Run %d of %d", R, runtimes));
    drawnow;
    for i = function_index
        parfor j = 1: numel(algorithms)
            feval(algorithms{j}, Global(i, j));
            Global(i, j).draw();
        end
    end
end
end