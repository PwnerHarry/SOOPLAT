function [nonseps, seps] = DG2(Global)
[delta, lambda, evaluations] = ism(Global);
[nonseps, seps, theta, epsilon] = dsm(evaluations, lambda, Global);
end

function [delta, lambda, evaluations] = ism(Global)
    center = 0.5 * (Global.problem.upperbound + Global.problem.lowerbound);
    %ind0 = Global.problem.dimension + 1;
    f_archive = NaN(Global.problem.dimension, Global.problem.dimension);
    fhat_archive = NaN(Global.problem.dimension, 1);
    delta1 = NaN(Global.problem.dimension, Global.problem.dimension);
    delta2 = NaN(Global.problem.dimension, Global.problem.dimension);
    lambda = NaN(Global.problem.dimension, Global.problem.dimension);
    p1 = Global.problem.lowerbound;
    fp1 = Global.evaluate(p1);
    counter = 0;
    for i=1:Global.problem.dimension-1
        if(~isnan(fhat_archive(i)))
            fp2 = fhat_archive(i);
        else
            p2 = p1;
            p2(i) = center(i);
            fp2 = Global.evaluate(p2);

            fhat_archive(i) = fp2;
        end
        for j=i+1:Global.problem.dimension
            counter = counter + 1;
            if(~isnan(fhat_archive(j)))
                fp3 = fhat_archive(j);
            else
                p3 = p1;
                p3(j) = center(j);
                fp3 = Global.evaluate(p3);
                fhat_archive(j) = fp3;
            end
            p4 = p1;
            p4(i) = center(i);
            p4(j) = center(j);
            fp4 = Global.evaluate(p4);
            f_archive(i, j) = fp4;
            f_archive(j, i) = fp4;
            d1 = fp2 - fp1;
            d2 = fp4 - fp3;
            delta1(i, j) = d1;
            delta2(i, j) = d2;
            lambda(i, j) = abs(d1 - d2);
        end
        renderCurve(Global);
    end
    evaluations.base = fp1;
    evaluations.fhat = fhat_archive;
    evaluations.F = f_archive;
    delta.delta1 = delta1;
    delta.delta2 = delta2;
end

function [nonseps, seps, theta, epsilon] = dsm(evaluations, lambda, Global)
    fhat_archive = evaluations.fhat;
    f_archive = evaluations.F;
    fp1 = evaluations.base;
    F1 = ones(Global.problem.dimension, Global.problem.dimension) * fp1;
    F2 = repmat(fhat_archive', Global.problem.dimension, 1);
    F3 = repmat(fhat_archive, 1, Global.problem.dimension);
    F4 = f_archive;
    FS = cat(3, F1, F2, F3, F4);
    Fmax = max(FS, [], 3);
    %Fmin = min(FS, [], 3);
    FS = cat(3, F1 + F4, F2 + F3);
    Fmax_inf = max(FS, [], 3);
    theta = nan(Global.problem.dimension);
    %reliable_calcs = 0;
    muM = eps / 2;
    gamma = @(n)((n.*muM)./(1-n.*muM));
    errlb = gamma(2) * Fmax_inf;
    errub = gamma(Global.problem.dimension^0.5) * Fmax;
    I1 = lambda <= errlb;
    theta(I1) = 0;
    I2 = lambda >= errub;
    theta(I2) = 1;
    %si1 = sum(sum(I1));
    %si3 = sum(sum(I2));
    I0 = (lambda == 0);
    c0 = sum(sum(I0));
    count_seps = sum(sum(~I0 & I1));
    count_nonseps = sum(sum(I2));
    reliable_calcs = count_seps + count_nonseps;
    w1 = ((count_seps+c0) / (c0+reliable_calcs));
    w2 = (count_nonseps / (c0+reliable_calcs));
    epsilon = w1*errlb + w2*errub;
    %grayind = (lambda < errub) & (lambda > errlb);
    %grayindsum = sum(sum(grayind));
    AdjTemp = lambda > epsilon;
    idx = isnan(theta);
    theta(idx) = AdjTemp(idx);
    theta = theta | theta';
    theta(logical(eye(Global.problem.dimension))) = 1;
    components = findConnComp(theta);
    h = @(x)(length(x) == 1);
    sizeone = cellfun (h, components);
    seps = components(sizeone);
    seps = cell2mat(seps);
    components(sizeone) = [];
    nonseps = components;
end

function [components] = findConnComp(C)
% C - connection matrix
% labels =[1 1 1 2 2 3 3 ...]  lenght(labels)=L, label for each vertex
% labels(i) is order number of connected component, i is vertex number
% rts - roots, numbers of started vertex in each component

L=size(C,1); % number of vertex

% Breadth-first search:
labels=zeros(1,L); % all vertex unexplored at the begining
rts=[];
ccc=0; % connected components counter
while true
    ind=find(labels==0);
    if ~isempty(ind)
        fue=ind(1); % first unexplored vertex
        rts=[rts fue];
        list=[fue];
        ccc=ccc+1;
        labels(fue)=ccc;
        while true
            list_new=[];
            for lc=1:length(list)
                p=list(lc); % point
                cp=find(C(p,:)); % points connected to p
                cp1=cp(labels(cp)==0); % get only unexplored vertecies
                labels(cp1)=ccc;
                list_new=[list_new cp1];
            end
            list=list_new;
            if isempty(list)
                break;
            end
        end
    else
        break;
    end
end

group_num = max(labels);
allgroups = cell(1, group_num);
for i = 1:group_num
    allgroups{i} = find(labels == i); 
end
components = allgroups;
end