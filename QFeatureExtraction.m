function [QFeature]=QFeatureExtraction(dataLead,sampFreq,timeAxis)


%% 1 -40 Hz filter
dataLead40 = dataLead;%-mean(dataLead);
 
pb140 = 1/(sampFreq/2);
pb240 = 40/(sampFreq/2);
passBand40 = [pb140 pb240];
order = 4;
[b40,a40] = butter(order,passBand40);
bpf40 = filtfilt(b40,a40,dataLead40);
 
s40 = diff(bpf40)/2;

s40 = [s40 0];

qrsLen = find(timeAxis == 132/sampFreq);

zeroPadS40 = [zeros(1,floor(qrsLen/2)) s40 zeros(1,ceil(qrsLen/2))];

nonLinTrfm = zeros(1,length(s40));
for var1 = 1:1:length(s40)
    nonLinTrfm(var1) = sum(abs(zeroPadS40(var1:var1 + qrsLen)));
end


thresh = nonLinTrfm > (0.3*max(nonLinTrfm));
threshIndLeft = find(diff([0 thresh])==1);
threshIndRight = find(diff([thresh 0])==-1); 

shift100 = 100;
tallTWave = find(timeAxis == 160/sampFreq);

lMin = zeros(1,length(threshIndLeft));
lMax = zeros(1,length(threshIndLeft));
adaptThresh = zeros(1,length(threshIndLeft));
qrsOn = zeros(1,length(threshIndLeft));
qrsOff = zeros(1,length(threshIndLeft));

if(threshIndLeft(1) <= 100)
    lMin(1) = min(nonLinTrfm(threshIndLeft(1):threshIndRight(1)) + shift100);
    lMax(1) = max(nonLinTrfm(threshIndLeft(1):threshIndRight(1)) + shift100);
    adaptThresh(1) = lMin(1) + (lMax(1) - lMin(1))/2.4;
    temp1 = nonLinTrfm(threshIndLeft(1):threshIndRight(1)) > (adaptThresh(1));
    if(isempty(find(temp1, 1)) == 1)
        leftOne = 0;
        rightOne = 0;
    else
        leftOne = find(diff([0 temp1])==1);
        rightOne = find(diff([temp1 0])==-1); 
    end
    qrsOn(1) = threshIndLeft(1) + leftOne;
    qrsOff(1) = threshIndLeft(1) + rightOne;
else
    lMin(1) = min(nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) + shift100);
    lMax(1) = max(nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) + shift100);
    adaptThresh(1) = lMin(1) + (lMax(1) - lMin(1))/2.4;
    temp1 = nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) > (adaptThresh(1));
    qrsOn(1) = threshIndLeft(1) - shift100 + find(diff([0 temp1])==1);
    qrsOff(1) = threshIndLeft(1) - shift100 + find(diff([temp1 0])==-1); 
end

% lMin(1) = min(nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) + shift100);
% lMax(1) = max(nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) + shift100);
% 
% adaptThresh(1) = lMin(1) + (lMax(1) - lMin(1))/2.4;
% temp1 = nonLinTrfm(threshIndLeft(1) - shift100 :threshIndRight(1)) > (adaptThresh(1));
% qrsOn(1) = threshIndLeft(1) - shift100 + find(diff([0 temp1])==1);
% qrsOff(1) = threshIndLeft(1) - shift100 + find(diff([temp1 0])==-1); 

for noBeat=2:1:length(threshIndLeft)
    if (threshIndRight(noBeat) - threshIndRight(noBeat - 1) > tallTWave)
        lMin(noBeat) = min(nonLinTrfm(threshIndLeft(noBeat) - shift100 :threshIndRight(noBeat)) + shift100);
        lMax(noBeat) = max(nonLinTrfm(threshIndLeft(noBeat) - shift100 :threshIndRight(noBeat)) + shift100);
        adaptThresh(noBeat) = (adaptThresh(noBeat - 1) + lMin(noBeat) + (lMax(noBeat) - lMin(noBeat))/2.8)/2;
        temp1 = nonLinTrfm((threshIndLeft(noBeat) - shift100 ):threshIndRight(noBeat)) > (adaptThresh(noBeat));
        if(isempty(find(temp1, 1)) == 1)
            leftOne = 0;
            rightOne = 0;
        else
            leftOne = find(diff([0 temp1])==1);
            rightOne = find(diff([temp1 0])==-1); 
        end
        qrsOn(noBeat) = threshIndLeft(noBeat) - shift100 + min(leftOne) - 1;
        qrsOff(noBeat) = threshIndLeft(noBeat) - shift100 + max(rightOne)- 1;
    end
end




%% 1 -90 Hz filter

dataLead90 = dataLead;%-mean(dataLead);
 
pb190 = 1/(sampFreq/2);
pb290 = 90/(sampFreq/2);
passBand90 = [pb190 pb290];
order = 4;
[b90,a90] = butter(order,passBand90);
bpf90 = filtfilt(b90,a90,dataLead90);

Hdn = filter([1/2 0 -1/2], 1,bpf90);

qrsOn = qrsOn(qrsOn>0);
qrsOff = qrsOff(qrsOff>0);

qrsOntemp1 = qrsOn(qrsOn~=qrsOff);
qrsOff = qrsOff(qrsOn~=qrsOff);
qrsOn = qrsOntemp1;

dMax = zeros(1,length(qrsOn));
dMin = zeros(1,length(qrsOn));


for var1=1:1:length(qrsOn)
    dMax(var1) = max(abs(Hdn(qrsOn(var1):qrsOff(var1))));
    if (qrsOn(var1) <= 150)
        dMin(var1) = 0.5*max(abs(Hdn(qrsOn(var1):qrsOn(var1))));  %% instead of 30, we used 150.    
    else
        dMin(var1) = max(abs(Hdn(qrsOn(var1) - 150:qrsOn(var1))));  %% instead of 30, we used 150.
    end
end

dMax(1) = mean(dMax(2:length(dMax)));
dMin(1) = mean(dMin(2:length(dMax)));

dth = zeros(1,length(qrsOn));
for var1=1:1:length(qrsOn)
    dth(var1) = 1.05*max([dMax(var1)/20 dMin(var1)]);
end

pne = zeros(1,length(Hdn));

for var2 = 1:1:length(dth)
    for var1=qrsOn(var2):1:qrsOff(var2)
        if(Hdn(var1) >= dth(var2))
            pne(var1) = 1;
        elseif(Hdn(var1)<=-dth(var2))
            pne(var1) = -1;
        else
            pne(var1) = 0;
        end
    end
end

RPeak = zeros(1,length(qrsOn));
RPeakFilt = zeros(1,length(qrsOn));
RInd = zeros(1,length(qrsOn));
SPeak = zeros(1,length(qrsOn));
SPeakFilt = zeros(1,length(qrsOn));
SInd = zeros(1,length(qrsOn));
QPeak = zeros(1,length(qrsOn));
QPeakFilt = zeros(1,length(qrsOn));
QInd = zeros(1,length(qrsOn));  
for var1=1:1:length(qrsOn)
    flagR = 0;
    flagR1 = 0;
    for start=qrsOn(var1):1:qrsOff(var1)
        if (pne(start)>0 && flagR1==0)
            RLeft = start;
            flagR = 1;
            flagR1=1;
        end
        if ((pne(start)<0) && flagR==1)
            RRight = start;
            flagR=0;
            flagR1=1;
%         else
%             RRight = start;
        end
    end
    
    if ((flagR==0) && (flagR1==1))
        [RPeak(var1) RInd(var1)] = max(dataLead(RLeft:RRight));
        RInd(var1) = RInd(var1) + RLeft - 1;
        RPeakFilt(var1) = dataLead(RInd(var1));

        [SPeak(var1) SInd(var1)] = min(dataLead(RInd(var1):qrsOff(var1)));
        SInd(var1) = SInd(var1) + RInd(var1) - 1;
        SPeakFilt(var1) = dataLead(SInd(var1));

        [QPeak(var1) QInd(var1)] = min(dataLead(qrsOn(var1):RInd(var1)));
        QInd(var1) = QInd(var1) + qrsOn(var1) - 1;
        QPeakFilt(var1) = dataLead(QInd(var1));
    end
end

tempInd = find(QInd==0);

for var1=1:1:length(tempInd)
    qrsOn(tempInd(var1)) = 0;
    qrsOff(tempInd(var1)) = 0;
end

qrsOn = qrsOn(qrsOn~=0);
qrsOff = qrsOff(qrsOff~=0);

RPeak = RPeak(RPeak~=0);
RInd = RInd(RInd~=0);
RPeakFilt = RPeakFilt(RPeakFilt~=0);

SPeak = SPeak(SPeak~=0);
SInd = SInd(SInd~=0);
SPeakFilt = SPeakFilt(SPeakFilt~=0);

QPeak = QPeak(QPeak~=0);
QInd = QInd(QInd~=0);
QPeakFilt = QPeakFilt(QPeakFilt~=0);


QRSSign = zeros(1,length(qrsOn));
for var1=1:1:length(QPeakFilt)
    if (QPeakFilt(var1) >= 0)
        QRSSign(var1) = 1;
    else
        QRSSign(var1) = -1;
    end
        
end

QFeature.QRSSign = QRSSign;

PAndNLoc = zeros(1,length(pne));
PAndNNum = zeros(1,length(qrsOn));
for var1=1:1:length(qrsOn)
    flagPN1 = 0;
    flagPN2 = 0;
    for var2=QInd(var1):1:SInd(var1)
        if ((pne(var2) > 0) && (flagPN1 == 0))
            PAndNLoc(var2) = var2;
            PAndNNum(var1) = PAndNNum(var1) + 1;
            flagPN1 = 1;
            flagPN2 = 0;
        elseif ((pne(var2) < 0) && (flagPN2 == 0))
            PAndNLoc(var2) = var2;
            PAndNNum(var1) = PAndNNum(var1) + 1;
            flagPN2 = 1;
            flagPN1 = 0;
        end
    end
end

PAndNLoc = PAndNLoc(PAndNLoc~=0);

QFeature.PAndNLoc = PAndNLoc;
QFeature.PAndNNum = PAndNNum;

RRIntInd = zeros(1,length(RInd) - 1);
RRInt = zeros(1,length(RInd) - 1);
for var1=1:1:length(RInd) - 1
    RRIntInd(var1) = RInd(var1 + 1) - RInd(var1);
    RRInt(var1) = timeAxis(RInd(var1 + 1) - RInd(var1));
end

QFeature.RRInt = RRInt;
QFeature.RRIntInd = RRIntInd;


% Rtime = timeAxis(RInd);
% Stime = timeAxis(SInd);
% Qtime = timeAxis(QInd);

QFeature.QPeakFilt = QPeakFilt;
QFeature.RPeakFilt = RPeakFilt;
QFeature.SPeakFilt = SPeakFilt;


QFeature.QPeak = QPeak;
QFeature.RPeak = RPeak;
QFeature.SPeak = SPeak;


QFeature.QInd = QInd;
QFeature.RInd = RInd;
QFeature.SInd = SInd;




%% Smoothing Filter filter

dataLeadSF = dataLead;%-mean(dataLead);

pb1SF = 0.5/(sampFreq/2);
pb2SF = 40/(sampFreq/2);
passBandSF = [pb1SF pb2SF];
order = 4;
[bSF,aSF] = butter(order,passBandSF);
bpfSF = filtfilt(bSF,aSF,dataLeadSF);

% HSF1 = smooth(dataLeadSF,'lowess');

HSF = filter([1/9 2/9 3/9 2/9 1/9],1,bpfSF);

HFdn = diff(HSF);

HFdn = [HFdn 0];

thOn = zeros(1,length(qrsOn));
thOff = zeros(1,length(qrsOn));

fdOn = zeros(1,length(qrsOn));
fdOff = zeros(1,length(qrsOn));
for var1=1:1:length(qrsOn)
    [thOn(var1) fdOn(var1)] = max(abs(HFdn(qrsOn(var1):QInd(var1))));
    [thOff(var1) fdOff(var1)] = max(abs(HFdn(SInd(var1):qrsOff(var1))));
%     min1 = min(abs(HFdn(qrsOn(var1):QInd(var1))));
%     min2 = min(abs(HFdn(SInd(var1):qrsOff(var1))));
%     thOn(var1) = (thOn(var1) + min1)/2;
%     thOff(var1) = (thOff(var1) + min2)/2;
    
    thOn(var1) = thOn(var1)/100;
    thOff(var1) = thOff(var1)/100;
    fdOn(var1) = fdOn(var1) + qrsOn(var1);
    fdOff(var1) = fdOff(var1) + SInd(var1);
end


% qrsOn

QRSOn = zeros(1,length(qrsOn));
QRSOff = zeros(1,length(qrsOn));
for var1=1:1:length(qrsOn)
    flag1 = 0;
    flag2 = 0;
    for var2=fdOn(var1):-1:qrsOn(var1)
        if(HFdn(var2)>thOn(var1) && flag1 == 0)
            QRSOn(var1) = var2;
            flag1 = 1;
        elseif (flag1 == 0)
            QRSOn(var1) = qrsOn(var1);
        end
    end
    
    for var2=fdOff(var1):1:qrsOff(var1)
        if(HFdn(var2)<thOff(var1) && flag2 == 0)
            QRSOff(var1) = var2;
            flag2 = 1;
        elseif (flag2 == 0)
            QRSOff(var1) = qrsOff(var1);
        end
    end
end

% QRSOn
% QRSOff

% qrsOnTime = timeAxis(QRSOn);
% qrsOffTime = timeAxis(QRSOff);
% Q = bpf90(QRSOn);
% S = bpf90(QRSOff);dataLead
QRSOnAmp = dataLead(QRSOn);
QRSOffAmp = dataLead(QRSOff);


% qrsInitOnTime = timeAxis(qrsOn);
% qrsInitOffTime = timeAxis(qrsOff);

QRSOnInit = dataLead(qrsOn);
QRSOffInit = dataLead(qrsOff);


QRSDur = timeAxis(QRSOff - QRSOn);
QRSInitDur = timeAxis(qrsOff - qrsOn);

QRSPos = zeros(1,length(QRSOn));
QRSNeg = zeros(1,length(QRSOff));
for var2=1:1:length(QRSOn)
    QRSPos(var2) = max(dataLead(QRSOn(var2):QRSOff(var2)));
    QRSNeg(var2) = min(dataLead(QRSOn(var2):QRSOff(var2)));
end


QFeature.QRSPos = QRSPos;
QFeature.QRSNeg = QRSNeg;

QFeature.QRSDur = QRSDur;
QFeature.QRSInitDur = QRSInitDur;

QFeature.qrsOn = qrsOn;
QFeature.qrsOff = qrsOff;

QFeature.QInit = QRSOnInit;
QFeature.SInit = QRSOffInit;


QFeature.QRSOn = QRSOn;
QFeature.QRSOff = QRSOff;


QFeature.QRSOnAmp = QRSOnAmp;
QFeature.QRSOffAmp = QRSOffAmp;
