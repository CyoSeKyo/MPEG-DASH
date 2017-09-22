function r = Dp_back(Buffer, R, lenR, W, SerTag, Segmentlen, lack, maxl, alpha, LastR, Group, TimeSlot, Para)

INF = 1e9;
Dp = ones(maxl, max(lenR), W * 100 + maxl) * -INF;
Trace = zeros(maxl, max(lenR), W * 100 + maxl, 2);

Beta = zeros(maxl, max(lenR));
Get = zeros(1, maxl);
Lc = 0;
for i = 1 : maxl
    if (SerTag(i) ~= 1)
        continue;
    end
    Lc = Lc + 1;
    Get(Lc) = i;
    for j = 1 : lenR(Group(i))
        Beta(Lc, j) = max(0, R(Group(i), j) * (1 - (Buffer(i) - Segmentlen(i) * Para) / TimeSlot)) * lack(i);
        Beta(Lc, j) = floor(Beta(Lc, j) * 100) + 1;
    end
end

A = 200;
B = 150;

for i = 1 : lenR(Group(Get(1)))
    cr = R(Group(Get(1)), i);
    Dp(1, i, Beta(1, i)) = A * log(alpha + cr) - B * max(LastR(Get(1)) - cr, 0);
end

for i = 2 : Lc
    for l = 1 : lenR(Group(Get(i - 1)))
        for k = 1 : W * 100 + maxl
            if (Dp(i - 1, l, k) > -INF)
                for j = 1 : lenR(Group(Get(i)))
                    cr = R(Group(Get(i)), j);
                    if (k + Beta(i, j) < W * 100 + maxl)
                        if (Dp(i, j, k + Beta(i, j)) < Dp(i - 1, l, k) + A * log(alpha + cr) - B * max(LastR(Get(i)) - cr, 0))
                            Dp(i, j, k + Beta(i, j)) = Dp(i - 1, l, k) + A * log(alpha + cr) - B * max(LastR(Get(i)) - cr, 0);
                            Trace(i, j, k + Beta(i, j), 1) = l;
                            Trace(i, j, k + Beta(i, j), 2) = k;
                        end
                    end
                end
            end
        end
    end
end

r = zeros(1, maxl);

nowi = Lc; nowj = 1; nowk = 1;
for j = 1 : lenR(Group(Get(Lc)))
    for k = 1 : W * 100 + maxl
        if (Dp(Lc, j, k) > Dp(Lc, nowj, nowk))
            nowj = j; nowk = k;
        end
    end
end

while(1)
    if (nowi == 0)
        break;
    end
    r(Get(nowi)) = R(Group(Get(nowi)), nowj);
    tempj = Trace(nowi, nowj, nowk, 1);
    tempk = Trace(nowi, nowj, nowk, 2);
    nowj = tempj; nowk = tempk;
    nowi = nowi - 1;
end
