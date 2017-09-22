function Ans = FindMinR(R, lr, x)
Ans = R(1);
for i = 1 : lr
    if (R(i) > x)
        Ans = R(i);
        break;
    end
end
