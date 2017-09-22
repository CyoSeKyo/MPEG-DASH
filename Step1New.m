function [fval, u, flag, w] = Step1New(r, TimeSlot, Cache, maxl, W, Segmentlen, Para)

f = zeros(1, maxl); A = zeros(1, maxl); b = zeros(1, 1); lb = zeros(maxl, 1);
for i = 1 : maxl  
    f(i) = -TimeSlot / r(i);
    A(1, i) = 1;
    lb(i) = max(0, r(i) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot));
end
b(1) = W;

Aeq = []; beq = [];

[w, fval, exitflag, output, lambda] = linprog(f, A, b, Aeq, beq, lb);
u = zeros(1, maxl + 1);

if (size(w, 1) ~= 0 && (abs(sum(w) - W) < 1e-2 || sum(w) < W))
    u = zeros(1, maxl + 1);
    flag = 1;
    for i = 1 : maxl
        u(1, i) = lambda.lower(i);
    end
    u(1, 1 + maxl) = lambda.ineqlin(1);
else
    flag = 0;
end
