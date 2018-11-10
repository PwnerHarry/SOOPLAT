function flag = stagCheck(P)
global espace_counter;
persistent M STD counter;
if isempty(P)
    M = [];
    STD = [];
    counter = 0;
end
currentM = mean(P, 1);
currentSTD = std(P, 0, 1);
if isequal(M, currentM) && isequal(STD, currentSTD)
    counter = counter + 1;
    if counter >= size(P, 2)
        flag = true;
        espace_counter = espace_counter + 1;
        fprintf('stagnant escaped for the %dth time\n', espace_counter);
        return;
    end
else
    M = currentM;
    STD = currentSTD;
end
flag = false;
end