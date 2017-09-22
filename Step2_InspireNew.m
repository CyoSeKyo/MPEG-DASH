function [r, r0] = Step2_InspireNew(A, B, C, U, cu, LL, cl, Cache, maxl, R, lenR, W, Group, TimeSlot, alpha, StartDelay, Segmentlen, Para, LastR)

Upper = zeros(1, maxl);
for i = 1 : maxl
    Upper(i) = lenR(Group(i));
end

for q = 1 : cl
    while(1)
        Left = 0; Right = W * LL(q, maxl + 1);
        for i = 1 : maxl
            temp = -LL(q, i) + LL(q, maxl + 1);
            if (temp < 0)
                Left = Left + W * temp;
            end
            Right = Right - max(0,  R(Group(i), Upper(i)) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot)) * LL(q, i);
        end
        tl = Left; tr = Right;
        if (tl > tr)
            minx = -1e10; p = 0;
            for i = 1 : maxl
                if (Upper(i) > 1)
                    t1 = max(0, R(Group(i), Upper(i)) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot)) * LL(q, i);
                    t2 = max(0, R(Group(i), Upper(i) - 1) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot)) * LL(q, i);
                    if (t1 - t2 > minx)
                        minx = t1 - t2;
                        p = i;
                    end
                end
            end
            Upper(p) = Upper(p) - 1;
        else
            break;
        end
    end
end

stack = zeros(1, maxl);
cur = 1;

r0 = 1e10;
r = zeros(1, maxl);
Beta = zeros(1, maxl);

while(cur > 0)
    if (cur > maxl)
        Cr = zeros(1, maxl);
        avgr = 0;
        for i = 1 : maxl
            Cr(i) = R(Group(i), stack(i));
            Beta(i) = max(0, Cr(i) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot));
            avgr = avgr + Cr(i);
        end
        avgr = avgr / maxl;
        flag = 1;
        for q = 1 : cl
            Left = 0; Right = 0;
            for i = 1 : maxl
                temp = -LL(q, i) + LL(q, maxl + 1);
                if (temp < 0)
                    Left = Left + temp * W;
                end
                Right = Right - LL(q, i) * Beta(i);
            end
            if (Left > Right + LL(q, maxl + 1) * W)
                flag = 0;
                break;
            end
        end
        if (flag == 1)
            rr = -1e11;
            for p = 1 : cu
                sumr = -U(p, maxl + 1) * W;
                for i = 1 : maxl
                    temp = -TimeSlot / Cr(i) - U(p, i) + U(p, maxl + 1);
                    if (temp < 0)
                        sumr = sumr + temp * W;
                    end
                    %sumr = sumr - A * log(alpha + Cr(i)) + B * (avgr - Cr(i)) * (avgr - Cr(i)) - (Cache(i) + TimeSlot + U(p, i) * Beta(i));
                    sumr = sumr - A * log(alpha + Cr(i)) + B * (avgr - Cr(i)) * (avgr - Cr(i)) + C * max((LastR(i) - Cr(i)), 0) - (Cache(i) + TimeSlot + U(p, i) * Beta(i));
                    %sumr = sumr - A * log(alpha + Cr(i)) + B * max((LastR(i) - Cr(i)), 0) - Cache(i) + TimeSlot + U(p, i) * Beta(i);                    end
                end
                rr = max(rr, sumr);
            end
            if (rr < r0)
                r0 = rr;
                for i = 1 : maxl
                    r(i) = Cr(i);
                end
            end
        end
        cur = cur - 1;
    else
        while(cur > 0 && stack(cur) == Upper(cur))
            cur = cur - 1;
        end
        if (cur == 0)
            break;
        end
        stack(cur) = stack(cur) + 1;
        cur = cur + 1;
        while(cur <= maxl)
            stack(cur) = 1;
            cur = cur + 1;
        end
    end
end