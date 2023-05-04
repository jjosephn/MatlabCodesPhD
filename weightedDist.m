function [WDD] = weightedDist(x,y,sampFreq)

xFeat = featureExtraction(x,sampFreq);

yFeat = featureExtraction(y,sampFreq);

deltaXY = zeros(1,18);

TWeight = [0 0.2 0.4;
           0.2 0 0.2;
           0.4 0.2 0];

PWeight = [0 0.2 0.2 0.2 0.3 0.2 0.2 0.2 0.4;
           0.2 0 0.1 0.1 0.3 0.4 0.4 0.4  0.2;
           0.2 0.1 0 0.2 0.3 0.4 0.4 0.4 0.2;
           0.2 0.1 0.1 0 0.3 0.4 0.4 0.4 0.2;
           0.3 0.3 0.3 0.3 0 0.3 0.3 0.3 0.3;
           0.2 0.4 0.4 0.4 0.3 0 0.1 0.1 0.2;
           0.2 0.4 0.4 0.4 0.3 0.1 0 0.1 0.2;
           0.2 0.4 0.4 0.4 0.3 0.1 0.1 0 0.2;
           0.4 0.2 0.2 0.2 0.3 0.2 0.2 0.2 0];
       
STWeight = [0 0.1 0.2 0.1 0.4;
            0.1 0 0.2 0.4 0.1;
            0.2 0.2 0 0.2 0.2;
            0.1 0.4 0.2 0 0.1;
            0.4 0.1 0.2 0.1 0];
       
lambda = [2 1 1 1 1 1 1 0.25 0.25 1 1 1 2 2 1 1 1 1];

Lambda = diag(lambda);

normWDD = zeros(1,length(xFeat.RRInt));
unNormWDD = zeros(1,length(xFeat.RRInt));

minLength = min([length(xFeat.RRInt) length(yFeat.RRInt)]);

xFeatVector = zeros(minLength,18);

yFeatVector = zeros(minLength,18);
       
for var1=1:1:minLength

    xFeatVector(var1,1) = xFeat.RRInt(var1);
    xFeatVector(var1,2) = xFeat.QRSDur(var1);
    xFeatVector(var1,3) = xFeat.QTInt(var1);
    xFeatVector(var1,4) = xFeat.QTPeak(var1);
    xFeatVector(var1,5) = xFeat.PDur(var1);
    xFeatVector(var1,6) = xFeat.PRInt(var1);
    xFeatVector(var1,7) = xFeat.PAndNNum(var1);
    xFeatVector(var1,8) = xFeat.QRSSign(var1);
    xFeatVector(var1,9) = xFeat.DeltaWave(var1);
    xFeatVector(var1,10) = xFeat.TShape(var1);
    xFeatVector(var1,11) = xFeat.PShape(var1);
    xFeatVector(var1,12) = xFeat.STShape(var1);
    xFeatVector(var1,13) = xFeat.QRSPos(var1);
    xFeatVector(var1,14) = xFeat.QRSNeg(var1);
    xFeatVector(var1,15) = xFeat.PAmp(var1);
    xFeatVector(var1,16) = xFeat.TPeak(var1);
    xFeatVector(var1,17) = xFeat.STPoint(var1);
    xFeatVector(var1,18) = xFeat.STSlope(var1);
    
    yFeatVector(var1,1) = yFeat.RRInt(var1);
    yFeatVector(var1,2) = yFeat.QRSDur(var1);
    yFeatVector(var1,3) = yFeat.QTInt(var1);
    yFeatVector(var1,4) = yFeat.QTPeak(var1);
    yFeatVector(var1,5) = yFeat.PDur(var1);
    yFeatVector(var1,6) = yFeat.PRInt(var1);
    yFeatVector(var1,7) = yFeat.PAndNNum(var1);
    yFeatVector(var1,8) = yFeat.QRSSign(var1);
    yFeatVector(var1,9) = yFeat.DeltaWave(var1);
    yFeatVector(var1,10) = yFeat.TShape(var1);
    yFeatVector(var1,11) = yFeat.PShape(var1);
    yFeatVector(var1,12) = yFeat.STShape(var1);
    yFeatVector(var1,13) = yFeat.QRSPos(var1);
    yFeatVector(var1,14) = yFeat.QRSNeg(var1);
    yFeatVector(var1,15) = yFeat.PAmp(var1);
    yFeatVector(var1,16) = yFeat.TPeak(var1);
    yFeatVector(var1,17) = yFeat.STPoint(var1);
    yFeatVector(var1,18) = yFeat.STSlope(var1);
    
    for var2=1:1:length(xFeatVector)
        if (var2 == 9)
            deltaXY(var2) = 0;
        elseif (var2 == 10)
            deltaXY(var2) = TWeight(xFeatVector(var1,var2),yFeatVector(var1,var2));
        elseif (var2 == 11)
            deltaXY(var2) = PWeight(xFeatVector(var1,var2),yFeatVector(var1,var2));
        elseif (var2 == 12)
            deltaXY(var2) = STWeight(xFeatVector(var1,var2),yFeatVector(var1,var2));
        else
            deltaXY(var2) =(abs(xFeatVector(var1,var2) - yFeatVector(var1,var2))/(max([abs(xFeatVector(var1,var2)) abs(yFeatVector(var1,var2))])));
        end
    end
    
    normWDD(var1) = deltaXY*(Lambda/Lambda')*deltaXY';
    unNormWDD(var1) = deltaXY*Lambda*deltaXY';
end

WDD.normalisedWDD = mean(normWDD);
WDD.unNormalisedWDD = mean(unNormWDD);
WDD.deltaXY = deltaXY;
WDD.xFeat = xFeatVector;
WDD.yFeat = yFeatVector;

WDD.xSignal = xFeat.dataLead;
WDD.ySignal = yFeat.dataLead;


