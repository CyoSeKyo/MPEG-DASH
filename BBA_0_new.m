function [Rating, Cache, Bandwidth, Mean, lenD, Rmean, Rvar, StarveCnt, DownRate, cP] = BBA_0_new(maxl, Data, ServeTime, Group, c, CacheSize, StartDelay, TimeSlot, W, P, Segmentlen, Wrule)
%Basic Parameter
R = zeros(6, 6); lenR = zeros(1, 5);
R(1, 1) = 0.1; R(1, 2) = 0.2; R(1, 3) = 0.4; R(1, 4) = 0.8; lenR(1) = 4;
R(2, 1) = 0.2; R(2, 2) = 0.5; R(2, 3) = 1; R(2, 4) = 2.5; lenR(2) = 4;
% R(2, 1) = 0.2; R(2, 2) = 0.4; R(2, 3) = 0.6; R(2, 4) = 0.8; R(2, 5) = 1.0; lenR(2) = 5;
R(2, 1) = 0.35; R(2, 2) = 0.6; R(2, 3) = 1.0; R(2, 4) = 2.0; R(2, 5) = 3.0; lenR(2) = 5;
R(3, 1) = 0.3; lenR(3) = 1;
R(4, 1) = 0.1; R(4, 2) = 0.3; R(4, 3) = 0.5; R(4, 4) = 0.7; lenR(4) = 4;
R(5, 1) = 0.1; R(5, 2) = 0.2; R(5, 3) = 0.3; lenR(5) = 3;

StopCr = 5;

TimeSlot = 0.01;
Slotlen = 200000;
Bandwidth = zeros(maxl, Slotlen); %Bandwith Record for each User
Cache = zeros(maxl, Slotlen); %Cache Record for each User
Rating = zeros(maxl, Slotlen); %Rating Record for each User
Sl = zeros(1, maxl); %Start time slot
El = zeros(1, maxl); %End time slot
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

%Serflag = zeros(1, maxl); % 1 -- User will be served in current TimeSlot; 0 -- Not
reservoir = 3; cushion = 5;
% reservoir = 5; cushion = 7;
StarveCnt = 0;
nowSl = 0; %Sl - Coordinate
Tag = zeros(1, maxl); % 0 -- Not arriving; 1 -- Being served; -1 -- Ending
Remain = zeros(1, maxl); %Remaining Time need to be served later for each user

while(1)
    nowSl = nowSl + 1;
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
            preTl(i) = nowSl - 1;
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
            continue;
        end
        Rating(i, nowSl) = Rating(i, nowSl - 1);
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
            R_l = 0; R_r = 0;
            if (nowSl > 1 && Rating(i, nowSl - 1) > 0)
                if (Rating(i, nowSl - 1) == R(Group(i), lenR(Group(i))))
                    R_r = R(Group(i), lenR(Group(i)));
                else
                    R_r = FindMinR(R(Group(i), :), lenR(Group(i)), Rating(i, nowSl - 1));
                end
                
                if (Rating(i, nowSl - 1) == R(Group(i), 1))
                    R_l = R(Group(i), 1);
                else
                    R_l = FindMaxR(R(Group(i), :), lenR(Group(i)), Rating(i, nowSl - 1));
                end
            end
            
            if (lenR(Group(i)) == 1)
                Fbr = R(Group(i), 1);
            else
                slant = (R(Group(i), lenR(Group(i))) - R(Group(i), 1)) / cushion;
                Fbr = (Buffer - reservoir) * slant + R(Group(i), 1);
            end
            if (Buffer <= reservoir)
                Rating(i, nowSl) = R(Group(i), 1);
            else if (Buffer >= reservoir + cushion)
                    Rating(i, nowSl) = R(Group(i), lenR(Group(i)));
                else if (Fbr >= R_r)
                        Rating(i, nowSl) = FindMaxR(R(Group(i), :), lenR(Group(i)), Fbr);
                    else if (Fbr <= R_l)
                            Rating(i, nowSl) = FindMinR(R(Group(i), :), lenR(Group(i)), Fbr);
                        else
                            Rating(i, nowSl) = Rating(i, nowSl - 1);
                        end
                    end
                end
            end
           % if (Rating(i, nowSl) == Rating(i, nowSl - 1))
           %     Cache(i, nowSl) = Cache(i, nowSl) + Addback;
           % end
        else
            SerTag(i) = -1;
        end
    end
    
    sumR = 0;
    cnt = 0;
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
                %   if (Cache(i, nowSl) < Segmentlen(i))
                %       StarveCnt = StarveCnt + 1;
                %       Remain(i) = Remain(i) + Segmentlen(i);
                %       if (Rating(i, nowSl - 1) ~= R(Group(i), 1))
                %           Cache(i, nowSl) = 0;
                %       end
                %       StartDelay(i) = Segmentlen(i);
                %       preTl(i) = nowSL;
                %   else
                preRemain(i) = Remain(i);
                Cache(i, nowSl) = Cache(i, nowSl) - Segmentlen(i);
                %  end
            end
            Cache(i, nowSl + 1) = min(CacheSize(i), Cache(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl)));            
            accl(i, nowSl + 1) = accl(i, nowSl) + TimeSlot * (Bandwidth(i, nowSl) / Rating(i, nowSl));
            if (accl(i, nowSl + 1) - (lenD(i) + 1) * Segmentlen(i) > -eps)
                if (i == 1)
                    accl(i, nowSl + 1)
                end
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
                % if (Cache(i, nowSl) < Segmentlen(i))
                %     StarveCnt = StarveCnt + 1;
                %     Remain(i) = Remain(i) + Segmentlen(i);
                %     if (Rating(i, nowSl - 1) ~= R(Group(i), 1))
                %         Cache(i, nowSl) = 0;
                %     end
                %     StartDelay(i) = Segmentlen(i);
                % else
                preRemain(i) = Remain(i);
                Cache(i, nowSl) = Cache(i, nowSl) - Segmentlen(i);
                % end
            end
            Cache(i, nowSl + 1) = Cache(i, nowSl);
            accl(i, nowSl + 1) = accl(i, nowSl); %+ Cache(i, nowSl + 1) - Cache(i, nowSl);
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

Color = 'ymcrgbwk';
figure(4)
hold on
for i = 1 : maxl
    who = i;
    %subplot(maxl, 1, who);
    %plot(Sl(who) : El(who), Seeing(who, Sl(who) : El(who)) /  ServeTime(who), 'r');
    %plot(Sl(who) : El(who), accl(who, Sl(who) : El(who)) /  ServeTime(who), 'b');
    plot((Sl(who) : El(who)) * TimeSlot, Rating(who, Sl(who) : El(who)), Color(i),'LineWidth',2);
    %xlim([0,400]);
end
hold off

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