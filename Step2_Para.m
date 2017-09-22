function [r, r0] = Step2_Para(U, cu, LL, cl, Cache, Tag, maxl, R, lenR, W, Group, TimeSlot, c, P, alpha, StartDelay, LastR, lack, GetR, Segmentlen, Para, B)

Beta = zeros(1, maxl);

bl = 0;
index = zeros(1, maxl);
for i = 1 : maxl
    if (Tag(i) == 1 && GetR(i) == -1)
        bl = bl + 1;
        index(bl) = i;
    end
end

%fuck = 2;
A = 200;
%B = 75;
C = 1;
stack = zeros(1, maxl);
cur = 1;

r0 = 1e10;
r = zeros(1, maxl);

while(cur > 0)
    if (cur > bl)
        Cr = zeros(1, maxl);
        for i = 1 : bl
            Cr(index(i)) = R(Group(index(i)), stack(i));
        end
        for i = 1 : maxl
            if (GetR(i) ~= -1)
                if (StartDelay(i) > 0)
                    %Cr(i) = R(Group(i), 1);
                %else
                    Cr(i) = GetR(i);
                end
            end
        end
        avgr = 0; cr = 0;
        for i = 1 : maxl
            if (Tag(i) == 1)
                Beta(i) = max(0, Cr(i) * (1 - (Cache(i) - Segmentlen(i) * Para) / TimeSlot)) * lack(i);
                avgr = avgr + Cr(i); cr = cr + 1;
            end
        end
        avgr = avgr / cr;
        flag = 1;
        for q = 1 : cl
            Left = 0; Right = 0;
            for i = 1 : maxl
                if (Tag(i) == 1)
                    temp = -LL(q, i) + LL(q, maxl + 1) + c(i) * LL(q, maxl + 2);
                    if (temp < 0)
                        Left = Left + temp * W;
                    end
                    Right = Right - LL(q, i) * Beta(i);
                end
            end
            if (Left > Right + LL(q, maxl + 1) * W + LL(q, maxl + 2) * P)
                flag = 0;
                break;
            end
        end
        if (flag == 1)
            rr = -1e11;
            for p = 1 : cu
                sumr = -U(p, maxl + 1) * W - U(p, maxl + 2) * P;
                for i = 1 : maxl
                    if (Tag(i) == 1)
                        temp = -TimeSlot / Cr(i) - U(p, i) + U(p, maxl + 1) + c(i) * U(p, maxl + 2);
                        if (temp < 0)
                            sumr = sumr + temp * W;
                        end
                        sumr = sumr - A * log(alpha + Cr(i)) + B * (avgr - Cr(i)) * (avgr - Cr(i)) - (Cache(i) + TimeSlot + U(p, i) * Beta(i)) * C;                    
                        %sumr = sumr - A * log(alpha + Cr(i)) + B * max((LastR(i) - Cr(i)), 0) - Cache(i) + TimeSlot + U(p, i) * Beta(i);
                    end
                end
                rr = max(rr, sumr);
            end
            if (rr < r0)
                r0 = rr;
                for i = 1 : maxl
                    if (Tag(i) == 1)
                        r(i) = Cr(i);
                    end
                end
            end
        end
        cur = cur - 1;
    else
        while(cur > 0 && stack(cur) == lenR(Group(index(cur))))
            cur = cur - 1;
        end
        if (cur == 0)
            break;
        end
        stack(cur) = stack(cur) + 1;
        cur = cur + 1;
        while(cur <= bl)
            stack(cur) = 1;
            cur = cur + 1;
        end
    end
end