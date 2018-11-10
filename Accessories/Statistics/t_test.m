function string_versus = t_test(mu1_sequence, sigma1_sequence, mu2_sequence, sigma2_sequence)
    independent_runs = 20;
    better_counter = 0;
    worse_counter = 0;
    insignificant_counter = 0;
    S = numel(mu1_sequence);
    if S ~= numel(sigma1_sequence) || S~= numel(mu2_sequence) || S ~= numel(sigma2_sequence)
        error('wtf');
    end
    for i = 1: numel(mu1_sequence)
        mu1 = mu1_sequence(i);
        sigma1 = sigma1_sequence(i);
        mu2 = mu2_sequence(i);
        sigma2 = sigma2_sequence(i);
        data1 = normrnd(mu1, sigma1, [1, independent_runs]);
        data2 = normrnd(mu2, sigma2, [1, independent_runs]);
        [h, ~] = ttest2(data1, data2, 'Vartype', 'unequal');
        if h == 1
            if mu1 < mu2
                better_counter = better_counter + 1;
            else
                worse_counter = worse_counter + 1;
            end
        else
            insignificant_counter = insignificant_counter + 1;
        end
    end
    string_versus = sprintf('%d/%d/%d', better_counter, insignificant_counter, worse_counter);
end