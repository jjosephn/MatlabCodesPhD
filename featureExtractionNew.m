function [feature]=featureExtractionNew(dataLeadOrg,sampFreq)

len1 = length(dataLeadOrg);

pb1Filt = 1/(sampFreq/2);
pb2Filt = 100/(sampFreq/2);
passBandFilt = [pb1Filt pb2Filt];
order = 4;
[bFilt,aFilt] = butter(order,passBandFilt);
dataLeadFilt = filtfilt(bFilt,aFilt,dataLeadOrg);
dataLeadOrg = dataLeadFilt - mean(dataLeadFilt);

% dataLeadOrg = dataLeadOrg - mean(dataLeadOrg);

dataLead = dataLeadOrg;


% create a time axis

tt = 1/sampFreq:1/sampFreq:ceil(len1/sampFreq);
timeAxis = tt(1:len1);

feature.timeAxis = timeAxis;

QFeature = QFeatureExtractionNew(dataLead,sampFreq,timeAxis);


feature.RRInt = QFeature.RRInt;
feature.RRIntInd = QFeature.RRIntInd;

feature.PAndNLoc = QFeature.PAndNLoc;
feature.PAndNNum = QFeature.PAndNNum;

feature.QRSSign = QFeature.QRSSign;

% Rtime = timeAxis(RInd);
% Stime = timeAxis(SInd);
% Qtime = timeAxis(QInd);

feature.QPeakFilt = QFeature.QPeakFilt;
feature.RPeakFilt = QFeature.RPeakFilt;
feature.SPeakFilt = QFeature.SPeakFilt;


feature.QPeak = QFeature.QPeak;
feature.RPeak = QFeature.RPeak;
feature.SPeak = QFeature.SPeak;


feature.QInd = QFeature.QInd;
feature.RInd = QFeature.RInd;
feature.SInd = QFeature.SInd;

feature.QRSPos = QFeature.QRSPos;
feature.QRSNeg = QFeature.QRSNeg;

feature.QRSDur = QFeature.QRSDur;
feature.QRSInitDur = QFeature.QRSInitDur;

feature.qrsOn = QFeature.qrsOn;
feature.qrsOff = QFeature.qrsOff;

feature.QRSOnInit = QFeature.QInit;
feature.QRSOffInit = QFeature.SInit;


feature.QRSOn = QFeature.QRSOn;
feature.QRSOff = QFeature.QRSOff;


feature.QRSOnAmp = QFeature.QRSOnAmp;
feature.QRSOffAmp = QFeature.QRSOffAmp;

feature.DeltaWave = zeros(1,length(feature.RRInt));


TFeature = TFeatureExtractionNew(dataLead,QFeature,sampFreq,timeAxis);


feature.QTInt = TFeature.QTInt;
feature.QTPeak = TFeature.QTPeak;

feature.QTIntInit = TFeature.QTIntInit;
feature.QTPeakInit = TFeature.QTPeakInit;



feature.TShape = TFeature.TShape;
feature.TInd = TFeature.TInd;
feature.TPeak = TFeature.TPeak;
feature.TOff = TFeature.TOff;
feature.TOffInd = TFeature.TOffInd;


feature.STPoint = TFeature.STPoint;
feature.STPointInd = TFeature.STPointInd;
feature.STSlope = TFeature.STSlope;
feature.STShape = TFeature.STShape;


PFeature = PFeatureExtractionNew(dataLead,QFeature,TFeature,sampFreq,timeAxis);

feature.POn = PFeature.POn;
feature.POff = PFeature.POff;

feature.PDur = PFeature.PDur;
feature.PRInt = PFeature.PRInt;
feature.PRIntInit = PFeature.PRIntInit;

feature.POnInd = PFeature.POnInd;
feature.POffInd = PFeature.POffInd;

feature.PAmp = PFeature.PAmp;
feature.PAmpInd = PFeature.PAmpInd;

feature.PShape = PFeature.PShape;







