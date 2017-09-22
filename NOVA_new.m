function [Rating, Cache, Bandwidth, Mean, lenD, Rmean, Rvar, DownRate, cP, StarveCnt] = NOVA_new(maxl, Data, ServeTime, Group, c, CacheSize, StartDelay, TimeSlot, W, P, Segmentlen);

%Bit Rate Set
R = zeros(6, 6); lenR = zeros(1, 5);
R(1, 1) = 0.1; R(1, 2) = 0.2; R(1, 3) = 0.4; R(1, 4) = 0.8; lenR(1) = 4;
R(2, 1) = 0.2; R(2, 2) = 0.5; R(2, 3) = 1; R(2, 4) = 2.5; lenR(2) = 4;
% R(2, 1) = 0.2; R(2, 2) = 0.4; R(2, 3) = 0.6; R(2, 4) = 0.8; R(2, 5) = 1.0; lenR(2) = 5;
R(2, 1) = 0.35; R(2, 2) = 0.6; R(2, 3) = 1.0; R(2, 4) = 2.0; R(2, 5) = 3.0; lenR(2) = 5;
R(3, 1) = 0.3; lenR(3) = 1;
R(4, 1) = 0.1; R(4, 2) = 0.3; R(4, 3) = 0.5; R(4, 4) = 0.7; lenR(4) = 4;
R(5, 1) = 0.1; R(5, 2) = 0.2; R(5, 3) = 0.3; lenR(5) = 3;

%Simulate Parameter
TimeSlot = 0.01;
Slotlen = 200000;
StopCr = 5;
Bandwidth = zeros(maxl, Slotlen);
Cache = zeros(maxl, Slotlen);
Rating = zeros(maxl, Slotlen);
Sl = zeros(1, maxl);
El = zeros(1, maxl);
cW = zeros(1, Slotlen);
cP = zeros(1, Slotlen);
Seeing = zeros(maxl, Slotlen);
accl = zeros(maxl, Slotlen);
preTl = zeros(1, maxl);
preRemain = ServeTime;
preState = zeros(1, maxl);
eps = 1e-3;
Download = zeros(maxl, max(ServeTime ./ Segmentlen));
DownRate = zeros(maxl, max(ServeTime ./ Segmentlen));
lenD = zeros(1, maxl);
Cflag = zeros(1, maxl);
Change = 0;

StarveCnt = 0;
nowSl = 0;
Tag = zeros(1, maxl); % 0 -- Not arriving; 1 -- Being served; -1 -- Ending
Remain = zeros(1, maxl); %Remaining Time need to be served later for each user

%Nova
xi = 0.05;
M_nova = zeros(maxl, Slotlen);
B_nova = zeros(maxl, Slotlen);
D_nova = zeros(maxl, Slotlen);
CalcTime = 2;
minB = 20;
Beta = 2;

while(1)
    nowSl = nowSl + 1;
    if (nowSl == 63336)
        -1;
    end
    CurrentTime = nowSl * TimeSlot;
    for i = 1 : maxl
        if (Tag(i) ~= 1)
            continue;
        end
        if (StartDelay(i) < eps && abs(preRemain(i) - Remain(i)) < eps)
            if (Cache(i, nowSl) - Segmentlen(i) < -eps)
                if (preState(i) == 1)
                    StarveCnt = StarveCnt + 1;
                end
                preState(i) = -1;
                %Remain(i) = Remain(i) + Segmentlen(i);
                if (Rating(i, nowSl - 1) ~= R(Group(i), 1))
                    accl(i, nowSl) = accl(i, nowSl) - Cache(i, nowSl);
                    Cache(i, nowSl) = 0;
                end
                StartDelay(i) = Segmentlen(i);
                preTl(i) = nowSl;
            end
        end
        if (Cache(i, nowSl) > Remain(i))
            Tag(i) = -1;
            El(i) = nowSl;
        end
    end
    
    %Is all Users have been served
    flag = 1;  %Tag of Ending the Simulation
    for i = 1 : maxl
        if (Tag(i) ~= -1)
            flag = 0;
        end
    end
    if (flag == 1) %End the Simulation
        break;
    end
    
    if (mod(CurrentTime, CalcTime) < eps)
        for i = 1 : maxl;
            if (Data(i) > CurrentTime && Data(i) < CurrentTime + CalcTime)
                Tag(i) = 1;
                Sl(i) = nowSl;
                Rating(i, nowSl) = R(Group(i), 1);
                Remain(i) = ServeTime(i);
                preRemain(i) = ServeTime(i);
                preTl(i) = nowSl - 1;
                M_nova(i, nowSl - 1) = round(lenR(Group(i)) / 2);
                B_nova(i, nowSl - 1) = 40 / 0.05;
                D_nova(i, nowSl - 1) = 1;
            end
        end
    end
    
    SerTag = zeros(1, maxl);
    %Arrange Rating
    for i = 1 : maxl
        if (Tag(i) ~= 1)
            continue;
        end
        M_nova(i, nowSl) = M_nova(i, nowSl - 1);
        B_nova(i, nowSl) = B_nova(i, nowSl - 1);
        D_nova(i, nowSl) = D_nova(i, nowSl - 1);
        if (mod(CurrentTime, CalcTime) < eps)
            B_nova(i, nowSl) = B_nova(i, nowSl) + xi * CalcTime / (1 + Beta);
        end
        if (StartDelay(i) > eps)
            SerTag(i) = 1;
            Rating(i, nowSl) = R(Group(i), 1);
            continue;
        end
        Rating(i, nowSl) = Rating(i, nowSl - 1);        
        f = zeros(1, maxl);
        if (Cache(i, nowSl) < StopCr * Segmentlen(i))
            temp = 0;
            SerTag(i) = 1;
            if (Cflag(i) == 0)
                continue;
            end
            Change = Change + 1;
            Cflag(i) = 0;
            Buffer = Cache(i, nowSl);
            while(Buffer - Segmentlen(i) > -eps)
                temp = temp + 1;
                Buffer = Buffer - Segmentlen(i);
            end
            accl(i, nowSl) = accl(i, nowSl) - Buffer;
            Cache(i, nowSl) = temp * Segmentlen(i);
            Buffer = Cache(i, nowSl) - (preRemain(i) - Remain(i));
            %new Rating
            qi = 1;
            temp = max(0, (B_nova(i, nowSl) - 20) / 0.05);
            f(i) = -0.005 * (B_nova(i, nowSl) / 0.05 + temp * temp);
            
            comq = qi - 0.05 * (qi - M_nova(i, nowSl - 1)) * (qi - M_nova(i, nowSl - 1));
            comq = comq - f(i) / (1 + Beta) * R(Group(i), qi) - 0.01 * 10 * D_nova(i, nowSl - 1) / 0.01 * R(Group(i), qi);
            for j = 2 : lenR(Group(i))
                comj = j - 0.05 * (j - M_nova(i, nowSl - 1)) * (j - M_nova(i, nowSl - 1));
                comj = comj - f(i) / (1 + Beta) * R(Group(i), j) - 0.01 * 10 * D_nova(i, nowSl - 1) / 0.01 * R(Group(i), j);
                if (comq < comj)
                    qi = j;
                    comq = comj;
                end
            end
            
            B_nova(i, nowSl) = max(B_nova(i, nowSl) - xi * Segmentlen(i), minB);
            M_nova(i, nowSl) = M_nova(i, nowSl) + xi * (qi - M_nova(i, nowSl - 1));
            D_nova(i, nowSl) = max(D_nova(i, nowSl) + xi * (Segmentlen(i) * R(Group(i), qi) - Segmentlen(i)), 1.5);%d0=1;
            Rating(i, nowSl) = R(Group(i), qi);
        else
            SerTag(i) = -1;
        end
    end
    
    %Deliver Width
    if (mod(CurrentTime, CalcTime) < eps)
        f = zeros(1, maxl); A = zeros(1, maxl); lb = zeros(1, maxl);
        for i = 1 : maxl
            if (SerTag(i) == 1)
                temp = max(0, (B_nova(i, nowSl - 1) - 20) / 0.05);
                f(i) = -0.005 * (B_nova(i, nowSl - 1) / 0.05 + temp * temp);
                A(i) = 1 / 0.01;
                lb(i) = 0.001;
            end
        end
        b = 1;
        Aeq = []; beq = [];
        [r, fval] = linprog(f, A, b, Aeq, beq, lb);
        
        for i = 1 : maxl
            if (SerTag(i) == 1)
                B_nova(i, nowSl) = B_nova(i, nowSl - 1) + xi * (TimeSlot / (1 + 0));
                Bandwidth(i, nowSl) = r(i) * 100 * W;
            end
        end
    else
        for i = 1 : maxl
            if (SerTag(i) == 1)
                Bandwidth(i, nowSl) = Bandwidth(i, nowSl - 1);
            end
        end
    end
    
    for i = 1 : maxl
        if (Tag(i) == 1 && StartDelay(i) < eps)
            if (StartDelay(i) < eps)
                Remain(i) = Remain(i) - TimeSlot;
            end
            Seeing(i, nowSl) = ServeTime(i) - Remain(i);
        end
    end
    
    for i = 1 : maxl
        if (SerTag(i) == 1)
            if (StartDelay(i) < eps && preRemain(i) - Remain(i) - Segmentlen(i) > -eps)
                preRemain(i) = Remain(i);
                Cache(i, nowSl) = Cache(i, nowSl) - Segmentlen(i);
            end
            Cache(i, nowSl + 1) = min(CacheSize(i), Cache(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl)));
            accl(i, nowSl + 1) = accl(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl));
            if (accl(i, nowSl + 1) - (lenD(i) + 1) * Segmentlen(i) > -eps)
                if (i == 1 && abs(accl(i, nowSl + 1) - 598) < eps)
                    accl(i, nowSl + 1)
                end
                lenD(i) = lenD(i) + 1;
                Download(i, lenD(i)) = (nowSl - preTl(i)) * TimeSlot;
                DownRate(i, lenD(i)) = Rating(i, nowSl);
                preTl(i) = nowSl;
                Cflag(i) = 1;
            end
            if (Cache(i, nowSl + 1) > StartDelay(i))
                StartDelay(i) = 0;
                preState(i) = 1;
            end
        end
        if (SerTag(i) == -1)
            preTl(i) = nowSl;
            if (StartDelay(i) < eps && preRemain(i) - Remain(i) - Segmentlen(i) > -eps)
                preRemain(i) = Remain(i);
                Cache(i, nowSl) = Cache(i, nowSl) - Segmentlen(i);
            end
            Cache(i, nowSl + 1) = Cache(i, nowSl);
            accl(i, nowSl + 1) = accl(i, nowSl);
            Rating(i, nowSl) = Rating(i, nowSl - 1);
            if (Cache(i, nowSl + 1) > StartDelay(i))
                StartDelay(i) = 0;
                preState(i) = 1;
            end
        end
    end
    
    for i = 1 : maxl
        if (Tag(i) == 1 && Cache(i, nowSl) < CacheSize(i))
            cW(nowSl) = cW(nowSl) + Bandwidth(i, nowSl);
            cP(nowSl) = cP(nowSl) + Bandwidth(i, nowSl) * c(i);
        end
    end
    cW(nowSl) = cW(nowSl) / W;
    cP(nowSl) = cP(nowSl) / P;
    
end

[Mean, Var, Rmean, Rvar] = Output(Cache, Rating, Sl, El, maxl);

Mean = zeros(1, maxl);
for i = 1 : maxl
    while(Cache(i, El(i)) > 0)
        View = min(Cache(i, El(i)), TimeSlot);
        Seeing(i, El(i) + 1) = Seeing(i, El(i)) + View;
        accl(i, El(i) + 1) = accl(i, El(i));
        Cache(i, El(i) + 1) = Cache(i, El(i)) - View;
        El(i) = El(i) + 1;
    end
    Mean(i) = mean(Cache(i, Sl(i) : El(i)));
end