function [PFeature]=PFeatureExtractionNew(dataLead,QFeature,TFeature,sampFreq,timeAxis)


%% P Wave detection

dataLeadPWave = dataLead;%dataLead;%-mean(dataLead);

RRInt = QFeature.RRInt;
qrsOn = QFeature.qrsOn;
QRSOn = QFeature.QRSOn;

TOffInd = TFeature.TOffInd;

% dataLeadPWave = dataLeadPWave - mean(dataLeadPWave);

pb1P = 0.3/(sampFreq/2);
pb2P = 30/(sampFreq/2);
passBandP = [pb1P pb2P];
order = 4;
[bP,aP] = butter(order,passBandP);
bpfP = filtfilt(bP,aP,dataLeadPWave);

HPFdn = diff(bpfP);

HPFdn = [HPFdn 0];

winP = zeros(length(RRInt),length(HPFdn));
FPWin = zeros(length(RRInt),length(HPFdn));
OPWin = zeros(length(RRInt),length(HPFdn));

for var2=1:1:length(RRInt)
    for var3=1:1:length(HPFdn)
        if ((var3 >= TOffInd(var2)) && (var3 <= floor(TOffInd(var2) + 0.15*(qrsOn(var2+1) - TOffInd(var2)))))
            winP(var2,var3) = 1 - (floor(TOffInd(var2) + 0.15*(qrsOn(var2+1) - TOffInd(var2))) - var3)/(floor(TOffInd(var2) + 0.15*(qrsOn(var2+1) - TOffInd(var2))) - TOffInd(var2));
        elseif ((var3 > floor(TOffInd(var2) + 0.15*(qrsOn(var2+1) - TOffInd(var2)))) && (var3 <= qrsOn(var2+1)))
            winP(var2,var3) = 1;
        else
            winP(var2,var3) = 0;
        end
    end
    FPWin(var2,:) = winP(var2,:).*HPFdn;
    OPWin(var2,:) = winP(var2,:).*bpfP;
end

PShape = zeros(1,length(RRInt));

thOnP = zeros(1,length(RRInt));
thOffP = zeros(1,length(RRInt));

absMaxP = zeros(1,length(RRInt));
absMaxIndP = zeros(1,length(RRInt));

minRight = zeros(1,length(RRInt));
minRightInd = zeros(1,length(RRInt));

minLeft = zeros(1,length(RRInt));
minLeftInd = zeros(1,length(RRInt));

maxRight = zeros(1,length(RRInt));
maxRightInd = zeros(1,length(RRInt));

maxLeft = zeros(1,length(RRInt));
maxLeftInd = zeros(1,length(RRInt));

POnInd = zeros(1,length(RRInt));
POffInd = zeros(1,length(RRInt));

PAmp = zeros(1,length(RRInt));
PAmpInd = zeros(1,length(RRInt));

PShapeThreshold = 0.00028;

for var2=1:1:length(RRInt)
    [absMaxP(var2) absMaxIndP(var2)] = max(abs(FPWin(var2,TOffInd(var2):qrsOn(var2+1))));
    absMaxIndP(var2) = absMaxIndP(var2) + TOffInd(var2) -1;
    
    if (FPWin(var2,absMaxIndP(var2)) > 0)
        [minRight(var2) minRightInd(var2)] = min(FPWin(var2,absMaxIndP(var2):qrsOn(var2+1)));
        [minLeft(var2) minLeftInd(var2)] = min(FPWin(var2,TOffInd(var2):absMaxIndP(var2)));
        
        minRightInd(var2) = minRightInd(var2) + absMaxIndP(var2) - 1;
        minLeftInd(var2) = minLeftInd(var2) + TOffInd(var2) - 1;
        
        thOnP(var2) = minLeft(var2)/5;
        thOffP(var2) = minRight(var2)/5;
        
        flag1 = 0;
        flag2 = 0;
        for var4=minLeftInd(var2):-1:TOffInd(var2)
            if(FPWin(var2,var4)<thOnP(var2) && flag1 == 0)
                POnInd(var2) = var4;
                flag1 = 1;
            elseif(flag1 == 0)
                POnInd(var2) = TOffInd(var2) + 80;
            end
        end
%         if(flag1 == 0)
%           POnInd(var2) = TOffInd(var2) + 20;
%           flag1 = 1;
%         end
        for var4=minRightInd(var2):1:qrsOn(var2+1)
            if(FPWin(var2,var4)<thOffP(var2) && flag2 == 0)
                POffInd(var2) = var4;
                flag2 = 1;
            elseif(flag2 == 0)
                POffInd(var2) = qrsOn(var2+1) - 5;
            end
        end
%         if(flag2 == 0)
%           POffInd(var2) = qrsOn(var2+1) - 5;
%           flag2 = 1;
%         end
        
        
        if (OPWin(var2,minRightInd(var2)) < -PShapeThreshold*(OPWin(var2,absMaxIndP(var2))))
            PShape(var2) = 1;           % positive P wave
            [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if (OPWin(var2,minLeftInd(var2)) < -PShapeThreshold*(OPWin(var2,absMaxIndP(var2))))
            PShape(var2) = 2;               % negative P Shape
            [PAmp(var2) PAmpInd(var2)] = min(dataLead(POnInd(var2):POffInd(var2)));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if ((OPWin(var2,minRightInd(var2)) < -PShapeThreshold*(OPWin(var2,absMaxIndP(var2)))) && (OPWin(var2,minLeftInd(var2)) < -PShapeThreshold*(OPWin(var2,absMaxIndP(var2)))))
            PShape(var2) = 4;           % biphase II P wave
            [PAmp(var2) PAmpInd(var2)] = max(abs(dataLead(POnInd(var2):POffInd(var2))));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if (PAmp(var2) == 0)
            [~, minIndP] = min(FPWin(var2,POnInd(var2):POffInd(var2)));
            [~, maxIndP]= max(FPWin(var2,POnInd(var2):POffInd(var2)));
            
            minIndP = minIndP + POnInd(var2) -1;
            maxIndP = maxIndP + POnInd(var2) -1;
            if (maxIndP < minIndP)
                PShape(var2) = 1;
                [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            elseif (maxIndP > minIndP)
                PShape(var2) = 2;
                [PAmp(var2) PAmpInd(var2)] = min(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            else
                PShape(var2) = 0;
                [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            end
        end
        
    else
        [maxRight(var2) maxRightInd(var2)] = max(FPWin(var2,absMaxIndP(var2):qrsOn(var2+1)));
        [maxLeft(var2) maxLeftInd(var2)] = max(FPWin(var2,TOffInd(var2):absMaxIndP(var2)));
        
        maxRightInd(var2) = maxRightInd(var2) + absMaxIndP(var2) - 1;
        maxLeftInd(var2) = maxLeftInd(var2) + TOffInd(var2) - 1;
        
        thOnP(var2) = maxLeft(var2)/5;
        thOffP(var2) = maxRight(var2)/5;
        
        
        flag1 = 0;
        flag2 = 0;
        for var4=maxLeftInd(var2):-1:TOffInd(var2)
            if(FPWin(var2,var4)<thOnP(var2) && flag1 == 0)
                POnInd(var2) = var4;
                flag1 = 1;
            elseif(flag1 == 0)
                POnInd(var2) = TOffInd(var2) + 80;
            end
        end
%         if(flag1 == 0)
%           POnInd(var2) = TOffInd(var2) + 20;
%           flag1 = 1;
%         end
        for var4=maxRightInd(var2):1:qrsOn(var2+1)
            if(FPWin(var2,var4)<thOffP(var2) && flag2 == 0)
                POffInd(var2) = var4;
                flag2 = 1;
            elseif(flag2 == 0)
                POffInd(var2) = qrsOn(var2+1) - 5;
            end
        end
%         if(flag2 == 0)
%           POffInd(var2) = qrsOn(var2+1) - 5;
%           flag2 = 1;
%         end
        
        
        if (OPWin(var2,maxRightInd(var2)) > PShapeThreshold*(OPWin(var2,absMaxIndP(var2))))
            PShape(var2) = 2;           % negative P wave
            [PAmp(var2) PAmpInd(var2)] = min(dataLead(POnInd(var2):POffInd(var2)));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if (OPWin(var2,maxLeftInd(var2)) > PShapeThreshold*(OPWin(var2,absMaxIndP(var2))))
            PShape(var2) = 1;               % positive P Shape
            [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if ((OPWin(var2,maxRightInd(var2)) > PShapeThreshold*(OPWin(var2,absMaxIndP(var2)))) && (OPWin(var2,maxLeftInd(var2)) > PShapeThreshold*(OPWin(var2,absMaxIndP(var2)))))
            PShape(var2) = 3;           % biphase I P wave
            [PAmp(var2) PAmpInd(var2)] = max(abs(dataLead(POnInd(var2):POffInd(var2))));
            PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
        end
        
        if (PAmp(var2) == 0)
            [~, minIndP] = min(FPWin(var2,POnInd(var2):POffInd(var2)));
            [~, maxIndP]= max(FPWin(var2,POnInd(var2):POffInd(var2)));
            
            minIndP = minIndP + POnInd(var2) -1;
            maxIndP = maxIndP + POnInd(var2) -1;
            if (maxIndP < minIndP)
                PShape(var2) = 1;
                [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            elseif (maxIndP > minIndP)
                PShape(var2) = 2;
                [PAmp(var2) PAmpInd(var2)] = min(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            else
                PShape(var2) = 0;
                [PAmp(var2) PAmpInd(var2)] = max(dataLead(POnInd(var2):POffInd(var2)));
                PAmpInd(var2) = PAmpInd(var2) + POnInd(var2) -1;
            end
        end
    end
    
end


% POnTime = timeAxis(POnInd);
% POffTime = timeAxis(POffInd);

POn = dataLead(POnInd);
POff = dataLead(POffInd);

PDur = timeAxis(POffInd - POnInd);
PRInt = timeAxis(QRSOn(2:length(QRSOn)) - POnInd);
PRIntInit = timeAxis(qrsOn(2:length(QRSOn)) - POnInd);

PAmp = dataLead(PAmpInd)/2000; % converted to physical units by dividing with gain=2000

PFeature.POn = POn;
PFeature.POff = POff;

PFeature.PDur = PDur;
PFeature.PRInt = PRInt;
PFeature.PRIntInit = PRIntInit;

PFeature.POnInd = POnInd;
PFeature.POffInd = POffInd;

PFeature.PShape = PShape;

PFeature.PAmp = PAmp;
PFeature.PAmpInd = PAmpInd;



% POn = bpfP(POn);
% POff = bpfP(POff);
