function [Data, ServeTime, Group, c, CacheSize, StartDelay, Segmentlen] = GetInput(maxl)

Avg = 300;
%Users arrvial -- Possion distribution
Lemda = 0.1;          %Parameter of Poisson Process
Data = zeros(1, maxl); %Data
cur = 1;
Data(cur) = random('exponential', 1 / Lemda) + 10;
while(cur < maxl)
    Data(cur + 1) = Data(cur) + random('exponential', 1 / Lemda) + 10;
    cur = cur + 1;
end

cost = ones(1, 5) * 0.2;
%cost(1) = 1.0; cost(2) = 0.8; cost(3) = 0.6; cost(4) = 0.4; cost(5) = 0.2; %Problems

Group = zeros(1, maxl); %Group-ID for each User
ServeTime = zeros(1 ,maxl); %Service Time for each User
CacheSize = zeros(1, maxl); %Size of Cache for each User
StartDelay = zeros(1, maxl); %StartDelay for each User
c = zeros(1, maxl);
Segmentlen = zeros(1, maxl);

Segcnt = 10;

for i = 1 : maxl
    %Get Group-ID
    Segmentlen(i) = 4;
    CacheSize(i) = Segcnt * Segmentlen(i);
    StartDelay(i) = Segmentlen(i);
    temp = rand(1);
%     for j = 1 : 5
%         if (temp < j * 0.2)
%             Group(i) = j;
%             break;
%         end
%     end
    temp = rand(1);
    for j = 1 : 5
        if (temp < j * 0.2)
            c(i) = cost(j);
            break;
        end
    end
    c(i) = 0.1;
    Group(i) = 2;
    %Get Service Time
    if (Group(i) == 2 || Group(i) == 4)
        ServeTime(i) = Avg;
    else
        while(1)
            temp = rand(1);
            temp = -Avg * log(temp);
            if (temp > Avg / 1.5 && temp < Avg * 1.5)
                break;
            end
        end
        ServeTime(i) = round(temp);
    end
end