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
        fprintf('DG2: %.2f%%\n', 100 * i / (Global.problem.dimension - 1));
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
        Global.draw();
    end
    evaluations.base = fp1;
    evaluations.fhat = fhat_archive;
    evaluations.F = f_archive;
    delta.delta1 = delta1;
    delta.delta2 = delta2;
end