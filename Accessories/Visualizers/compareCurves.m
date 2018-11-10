function compareCurves(varargin)
    for i = 1: nargin
        dataFrames(i) = varargin{i};
    end
    clf;
    for i = 1: numel(dataFrames)
        H(i) = semilogy(dataFrames(i).trace(:, 1), dataFrames(i).trace(:, 2), 'LineWidth', 2);
        if isempty(dataFrames(i).algorithm)
            algoname{i} = 'unknown';
        else
            algoname{i} = dataFrames(i).algorithm;
        end
        hold on;
    end
    L = legend(H, algoname);
    set(L, 'Fontname', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12); 
    hold off;
    drawnow;
end