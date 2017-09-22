function [Rating, Cache, Bandwidth, Mean, lenD, Rmean, Rvar, DownRate, cP, StarveCnt] = ELASTIC(maxl, Data, ServeTime, Group, c, CacheSize, StartDelay, TimeSlot, W, P, Segmentlen, Wrule)

%Fi(qi) = R(Group(i), qi;
R = zeros(6, 6); lenR = zeros(1, 5);
R(1, 1) = 0.1; R(1, 2) = 0.2; R(1, 3) = 0.4; R(1, 4) = 0.8; lenR(1) = 4;
R(2, 1) = 0.2; R(2, 2) = 0.5; R(2, 3) = 1; R(2, 4) = 2.5; lenR(2) = 4;
% R(2, 1) = 0.2; R(2, 2) = 0.4; R(2, 3) = 0.6; R(2, 4) = 0.8; R(2, 5) = 1.0; lenR(2) = 5;
R(2, 1) = 0.35; R(2, 2) = 0.6; R(2, 3) = 1.0; R(2, 4) = 2.0; R(2, 5) = 3.0; lenR(2) = 5;
R(3, 1) = 0.3; lenR(3) = 1;
R(4, 1) = 0.1; R(4, 2) = 0.3; R(4, 3) = 0.5; R(4, 4) = 0.7; lenR(4) = 4;
R(5, 1) = 0.1; R(5, 2) = 0.2; R(5, 3) = 0.3; lenR(5) = 3;

TimeSlot = 0.01;
StopCr = 5;
%Segmentlen = 5;
Slotlen = 200000;
Bandwidth = zeros(maxl, Slotlen); %Bandwith Record for each User
Cache = zeros(maxl, Slotlen); %Cache Record for each User
Rating = zeros(maxl, Slotlen); %Rating Record for each User
Sl = zeros(1, maxl); %Start time slot
El = zeros(1, maxl); %End time slot
cW = zeros(1, Slotlen);
cP = zeros(1, Slotlen);
preState = zeros(1, maxl);
preRemain = ServeTime;
preTl = zeros(1, maxl);
Seeing = zeros(maxl, Slotlen);
accl = zeros(maxl, Slotlen);
Download = zeros(maxl, max(ServeTime ./ Segmentlen));
DownRate = zeros(maxl, max(ServeTime ./ Segmentlen));
lenD = zeros(1, maxl);
QI = zeros(1, maxl);
eps = 1e-2;
Change = 0;
Cflag = zeros(1, maxl);

StarveCnt = 0;
CurrentTime = TimeSlot; %Time-Coordinate
nowSl = 0; %Sl - Coordinate
Tag = zeros(1, maxl); % 0 -- Not arriving; 1 -- Being served; -1 -- Ending
Remain = zeros(1, maxl); %Remaining Time need to be served later for each user
%Nova
Level = zeros(maxl, Slotlen);

while(1)
    nowSl = nowSl + 1;
    if (nowSl == 36)
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
    
    for i = 1 : maxl
        if (Data(i) > CurrentTime && Data(i) < CurrentTime + TimeSlot)
            Tag(i) = 1;
            Sl(i) = nowSl;
            Rating(i, nowSl) = R(Group(i), 1);
            Remain(i) = ServeTime(i);
            preTl(i) = nowSl;
        end
    end
    
    SerTag = zeros(1, maxl);
    for i = 1 : maxl
        if (Tag(i) ~= 1)
            continue;
        end
        if (StartDelay(i) > eps)
            SerTag(i) = 1;
            Rating(i, nowSl) = R(Group(i), 1);
            Level(i, nowSl) = 1;
            continue;
        end
        Rating(i, nowSl) = Rating(i, nowSl - 1);
        Level(i, nowSl) = Level(i, nowSl - 1);
        if (Cache(i, nowSl) < StopCr * Segmentlen(i))
            SerTag(i) = 1; 
            if (Cflag(i) == 0)
                continue;
            end          
            Cflag(i) = 0;
            Change = Change + 1;
            temp = 0;
            Buffer = Cache(i, nowSl);         
            while(Buffer - Segmentlen(i) > -eps)
                temp = temp + 1;
                Buffer = Buffer - Segmentlen(i);
            end
            Cache(i, nowSl) = temp * Segmentlen(i);
            accl(i, nowSl) = accl(i, nowSl) - Buffer;
            R_d = 1;
            R_q = Cache(i, nowSl);
            R_r = 0;
            for j = max(1, lenD(i) - 5) : (lenD(i) - 1)
                R_r = R_r + DownRate(i, j) * Segmentlen(i) / Download(i, j);
            end
            if (lenD(i) == 1) 
                R_r = DownRate(i, lenD(i)) * Segmentlen(i) / Download(i, lenD(i));
            else
                R_r = R_r / min(5, lenD(i) - 1) * 0.25 + DownRate(i, lenD(i)) * Segmentlen(i) / Download(i, lenD(i)) * 0.75;
            end
            QT = Segmentlen(i) * 5;
            QI(i) = QI(i) + Download(i, lenD(i)) * (R_q - QT);
            Level(i, nowSl) = max(1, floor(R_r / (R_d - R_q / 100 - QI(i) / 1000)));
            Level(i, nowSl) = min(Level(i, nowSl), lenR(Group(i)));
            Rating(i, nowSl) = R(Group(i), Level(i, nowSl));
           % if (Rating(i, nowSl) == Rating(i, nowSl - 1))
           %     Cache(i, nowSl) = Cache(i, nowSl) + Addback;
           %     accl(i, nowSl) = accl(i, nowSl) + Addback;
           % end
        else
            SerTag(i) = -1;
        end
    end
    
    sumR = 0; cnt = 0;
    for i = 1 : maxl
        if (SerTag(i) == 1)
            sumR = sumR + Rating(i, nowSl);
            cnt = cnt + 1;
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
            if (Wrule == 1)
                Bandwidth(i, nowSl) = W * Rating(i, nowSl) / sumR;
            else
                Bandwidth(i, nowSl) = W / cnt;
            end
            if (StartDelay(i) < eps && preRemain(i) - Remain(i) - Segmentlen(i) > -eps)
                preRemain(i) = Remain(i);
                Cache(i, nowSl) = Cache(i, nowSl) - Segmentlen(i);
            end
            Cache(i, nowSl + 1) = min(CacheSize(i), Cache(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl)));
            accl(i, nowSl + 1) = accl(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl));
            if (accl(i, nowSl + 1) - (lenD(i) + 1) * Segmentlen(i) > -eps)
                lenD(i) = lenD(i) + 1;
                Download(i, lenD(i)) = (nowSl - preTl(i)) * TimeSlot;
                DownRate(i, lenD(i)) = Rating(i, nowSl);
                preTl(i) = nowSl;
                Cflag(i) = 1;
            end
            %if (flag == 1)
            %    accl(i, nowSl + 1) = accl(i, nowSl + 1) - Segmentlen(i);
            %end
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
            accl(i, nowSl + 1) = accl(i, nowSl); %+ Cache(i, nowSl + 1) - Cache(i, nowSl);
            Rating(i, nowSl) = Rating(i, nowSl - 1);
            Bandwidth(i, nowSl) = Bandwidth(i, nowSl - 1);
            Level(i, nowSl) = Level(i, nowSl - 1);
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
