function [TFeature]=TFeatureExtractionNew(dataLead,QFeature,sampFreq,timeAxis)

%%% T Wave Detection

dataLeadTWave = dataLead;%dataLead;%-mean(dataLead);
dataLeadOrg = dataLead;

RInd = QFeature.RInd;
QRSOn = QFeature.QRSOn;
qrsOn = QFeature.qrsOn;

QRSOff = QFeature.QRSOff;
qrsOff = QFeature.qrsOff;

S = dataLead(QRSOff);


HTWave = filter([ 1 zeros(1,5) -1 0 -1 zeros(1,5) 1] , [1 -1],dataLeadTWave);
HTWave = HTWave(14:length(HTWave));
HTWave = [zeros(1,13) HTWave];

% HTWave = filter([1/9 2/9 3/9 2/9 1/9],1,HTWave);

% HDTWave = diff(HTWave);


RRInt = zeros(1,length(RInd) - 1);
for var1=1:1:length(RInd) - 1
    RRInt(var1) = RInd(var1 + 1) - RInd(var1);
end

% RRIntAvg = mean(RRInt);

    
samp700 = find((round(timeAxis*10000)/10000) == 700/sampFreq);
% samp500 = find((round(timeAxis*10000)/10000) == 500/sampFreq);
% samp120 = find((round(timeAxis*10000)/10000) == 120/sampFreq);
% samp100 = find((round(timeAxis*10000)/10000) == 100/sampFreq);

bw = zeros(1,length(RInd));
mw = zeros(1,length(RInd));
ew = zeros(1,length(RInd));

for var2=1:1:length(RInd)
    flagTPeak = 0;
    if (var2 >= 4)
        RRAvg = mean(RRInt(var2-3:var2-1));
    elseif (var2 == 3)
        RRAvg = mean(RRInt(var2-2:var2-1));
    else
        RRAvg = RRInt(1);
    end
    
    if (RRAvg > samp700)
        bw(var2) = 120;
        ew(var2) = 500;
    else
        bw(var2) = 100;
        ew(var2) = ceil(0.7*RRAvg);
    end
    
    mw(var2) = ceil(0.8*(ew(var2) - bw(var2)));
end

win = zeros(length(RInd),length(HTWave));

TShape = zeros(1,length(RInd));
TInd = zeros(1,length(RInd));
TPeak = zeros(1,length(RInd));
TOff = zeros(1,length(RInd));
TOffInd = zeros(1,length(RInd));

for var2=1:1:length(RInd)
    if ((RInd(var2) + ew(var2)) < length(HTWave))
        for var3=1:1:length(HTWave)
            if ((var3 >= RInd(var2) + bw(var2)) && (var3 < RInd(var2) + mw(var2)))
                win(var2,var3) = 1;
            elseif ((var3 >= RInd(var2) + mw(var2)) && (var3 <= RInd(var2) + ew(var2)))
                win(var2,var3) = 1 - ((var3 - RInd(var2) - mw(var2))/(ew(var2) - mw(var2)));
            else
                win(var2,var3) = 0;
            end
        end
    else
        win(var2,:) = zeros(1,length(HTWave));
    end
    
    GWin = win(var2,:).*HTWave;
    
    [GWMax locMax] = max(GWin);
    [GWMin locMin] = min(GWin);
    
    if (locMax < locMin)
        TShape(var2) = 1;  % Positive T Wave
        thTWave = GWMin/15;
        temp1 = find(GWin(locMax:locMin) < 0);
        TInd(var2) = locMax + temp1(1) -1;
        TPeak(var2) = dataLead(TInd(var2));
        temp2 = find(GWin(locMin:length(GWin)) > thTWave);
        TOffInd(var2) = locMin + temp2(1) - 1;
        TOff(var2)  = dataLead(TOffInd(var2));
    elseif(locMax > locMin)
        TShape(var2) = 3;           % Negative T Wave
        thTWave = GWMax/15;
        temp1 = find(GWin(locMin:locMax) > 0);
        TInd(var2) = locMin + temp1(1) -1;
        TPeak(var2) = dataLead(TInd(var2));
        temp2 = find(GWin(locMax:length(GWin)) < thTWave);
        TOffInd(var2) = locMax + temp2(1) - 1;
        TOff(var2)  = dataLead(TOffInd(var2));
    else
        TShape(var2) = 2;
        thTWave = 0;
        TInd(var2) = 0;
        TPeak(var2) = 0;
        TOff(var2) = 0;
    end
end


TInd = TInd(TInd~=0);
TPeak = TPeak(TPeak~=0);
TOffInd = TOffInd(TOffInd~=0);
TOff = TOff(TOff~=0);


% Ttime = timeAxis(TInd);
% TOtime = timeAxis(TOffInd);


if (length(QRSOn) == length(TOffInd))
    QTInt = timeAxis(TOffInd - QRSOn);
    QTPeak = timeAxis(TInd - QRSOn);
    
    QTIntInit = timeAxis(TOffInd - qrsOn);
    QTPeakInit = timeAxis(TInd - qrsOn);
else
    QTInt = timeAxis(TOffInd - QRSOn(1:length(TOffInd)));
    QTIntInit = timeAxis(TOffInd - qrsOn(1:length(TOffInd)));
    
    QTPeak = timeAxis(TInd - QRSOn(1:length(TInd)));
    QTPeakInit = timeAxis(TInd - qrsOn(1:length(TInd)));
end



heartRate = length(RInd)*60/(length(dataLead)/sampFreq);


STPointInd = zeros(1,length(RInd));
STPoint = zeros(1,length(RInd));
for var2=1:1:length(RInd)
    if (heartRate > 60 && (ceil(QRSOff(var2) + 60) < length(dataLeadOrg)))
        STPointInd(var2) = ceil(QRSOff(var2) + 60);
        STPoint(var2) = dataLeadOrg(STPointInd(var2));
    elseif(ceil(QRSOff(var2) + 80) < length(dataLeadOrg))
        STPointInd(var2) = ceil(QRSOff(var2) + 80);
        STPoint(var2) = dataLeadOrg(STPointInd(var2));
    end
end

% STPointTime = timeAxis(STPointInd);

STSlope = zeros(1,length(STPointInd));
for var2=1:1:length(STPointInd)
    STSlope(var2) = (STPoint(var2) - S(var2))/((STPointInd(var2) - QRSOff(var2))/sampFreq);
end

STPointInd = STPointInd(STPointInd~=0);
STPoint = STPoint(STPoint~=0);


STShape = zeros(1,length(STPointInd));

P1 = zeros(1,length(STPointInd));
P2 = zeros(1,length(STPointInd));
P3 = zeros(1,length(STPointInd));

for var2=1:1:length(STPointInd)
    P1(var2) = dataLeadOrg(ceil(QRSOff(var2) - 60));
    P2(var2) = dataLeadOrg(ceil(QRSOff(var2)));
    P3(var2) = dataLeadOrg(ceil(QRSOff(var2) + 60));

    if ((P1(var2) < P2(var2)) && (P1(var2) < P3(var2)) && (P2(var2) < P3(var2)))
        STShape(var2) = 1;
    elseif ((P1(var2) < P2(var2)) && (P3(var2) < P2(var2)))
        STShape(var2) = 2;
    elseif ((P1(var2) > P2(var2)) && (P3(var2) > P2(var2)))
        STShape(var2) = 4;
    elseif ((P1(var2) > P2(var2)) && (P1(var2) > P3(var2)) && (P2(var2) > P3(var2)))
        STShape(var2) = 5;
    else
        STShape(var2) = 3;
    end
end



TFeature.STPoint = STPoint;
TFeature.STPointInd = STPointInd;
TFeature.STSlope = STSlope;
TFeature.STShape = STShape;


TFeature.QTInt = QTInt;
TFeature.QTPeak = QTPeak;

TFeature.QTIntInit = QTIntInit;
TFeature.QTPeakInit = QTPeakInit;



TFeature.TShape = TShape;
TFeature.TInd = TInd;
TFeature.TPeak = TPeak;
TFeature.TOff = TOff;
TFeature.TOffInd = TOffInd;
