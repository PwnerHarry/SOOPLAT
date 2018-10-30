function [nonseps, seps] = DG2(Global)
[delta, lambda, evaluations] = ism(Global);
[nonseps, seps, theta, epsilon] = dsm(evaluations, lambda, Global);
end