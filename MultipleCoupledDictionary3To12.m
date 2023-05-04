clc;
close all;
clearvars;

tic

% Update the path for data base.

addpath('../PTBDatabaseFull')
addpath('../ScSR')
addpath('../KSVD_Matlab_ToolBox')
addpath(genpath('ScSR/RegularizedSC'));

format('short');
leadNumber=15;
[num,txt,raw] = xlsread('dbListUpdated.xls');
[~,n]=size(txt);
samples=30000;

% sampFreq = 1000;
% datMax=77;

%% Memory allocation

patientID{549} = [];
dataName{290}=[];
recordName{7} = [];
leadName{leadNumber + 1}=[];
xl{leadNumber}=[];
leadCobination{60}=[];
xyProjection{60}=[];
xyError{60}=[];
thetaVal=zeros(1,60);
sampFreq = 1000;

%% Preprocessing

fid1 = fopen('ptb.txt');
var1 = 1;
tline = fgetl(fid1);
while ischar(tline)
    str = sprintf(tline);
    str = regexprep(str,'">.*','');
    loc = strfind(str,'<');
    str(loc(1)) = [];
    loc = strfind(str,'"');
    str(loc(1)) = [];
    str2 = regexprep(str, 'patient(\w*)/','');
    str = regexprep(str,'/.*','');
    str2 = strcat(str2,'m');

    patientID{1,var1} = str;
    patientID{2,var1} = str2;
    var1 = var1 + 1;

    tline = fgets(fid1);
end

fclose(fid1);
pNo = 1;
dNo = 0;
classNum = 1;
classesInit{15} = [];
patientClass{549} = [];
prevPatientName = 'patient000';
count = 1;
while (pNo <= length(patientID))
    patientName = patientID{1,pNo};

    if (strcmp(patientID(1,pNo),patientName)==1)
        s1 = strcat(patientID(2,pNo));
        indx1 = find(strcmp(txt(:,1),s1),1);
        s2 = strcat(txt(indx1,2));
        s3 = strcat(txt(indx1,3));
        patientClass{1,pNo} = s2{1};
        patientID{3,pNo} = s2{1};
        patientID{4,pNo} = s3{1};
        indx2 = find(strcmp(classesInit,s2));
        if isempty(indx2)
            classesInit{1,classNum} = s2{1};
            classNum = classNum + 1;
        end
        mi=load(cell2mat(s1));
        
        x=mi.val;
        xtemp=zeros(size(x(1,1:samples)));
        
        % preprocessing

        fs=1000; 

        for leads=1:1:leadNumber
            leadName{leads}=sprintf('x%d',leads);
            bline = LPFilter(x(leads,1:samples),.7/fs);
            xtemp=x(leads,1:samples)-bline;
            xtemp=xtemp-mean(xtemp);
            xl{leads}=(xtemp/max(abs(xtemp)));
        end
        
        if (strcmp(patientName,prevPatientName)==1)
            count = count + 1;
            recordName{count} = sprintf('record%d',count);
            for leads=1:1:leadNumber
                data.(dataName{dNo}).(recordName{count}).(leadName{leads})=xl{leads};
            end
        else
            count = 1;
            dNo = dNo + 1;
            dataName{dNo}=sprintf('%s',patientName);
            recordName{count} = sprintf('record%d',count);
            for leads=1:1:leadNumber
                data.(dataName{dNo}).(recordName{count}).(leadName{leads})=xl{leads};
            end
        end

        clear LPFilter;
    end
    prevPatientName = patientName;
    pNo = pNo + 1;
end

%% Initialization

limbLeads = {'x1','x2','x3','x4','x5','x6'};

precordialLeads = {'x7','x8','x9','x10','x11','x12'};

frankLeads = {'x13','x14','x15'};

featureName = {'RRInt', 'QRSDur', 'QTInt', 'QTPeak', 'PDur', 'PRInt', 'PAndNNum','QRSSign', 'DeltaWave', 'TShape', 'PShape', 'STShape', 'QRSPos', 'QRSNeg', 'PAmp', 'TPeak', 'STPoint', 'STSlope'};

predictorLeadNames = {'x1','x2','x9'};
responseLeadNames = {'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12'};

allLeads = [predictorLeadNames responseLeadNames];


leadNameOrg = {'LeadI', 'LeadII', 'LeadIII', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'X', 'Y', 'Z'};

initPredictorLeads = {'x1','x2'};
predictorLeadSet = {'x7','x8','x9','x10','x11','x12'};

wedd = zeros(1,length(patientID));
RSquare = zeros(1,length(patientID));
corCoeff = zeros(1,length(patientID));
RMSE = zeros(1,length(patientID));

weddPrec = zeros(1,length(patientID));
RSquarePrec = zeros(1,length(patientID));
corCoeffPrec = zeros(1,length(patientID));
RMSEPrec = zeros(1,length(patientID));

weddDictAll = zeros(length(responseLeadNames),length(patientID));
RSquareDictAll = zeros(length(responseLeadNames),length(patientID));
corCoeffDictAll = zeros(length(responseLeadNames),length(patientID));
RMSEDictAll = zeros(length(responseLeadNames),length(patientID));

level = 7;

maxDict = 25;
dictionarySize = 40;

selDataCount = 0;
selLeadIndex = zeros(1,length(inclDataRecords));

for dat = 1:1:length(fieldnames(data))

    count = 1;
    while (count <= length(fieldnames(data.(dataName{dat}))))

        selDataCount = selDataCount + 1;
        fprintf('%srecord%d \n',dataName{dat},count);

        %% Training Data

        trainSamplesDict = 20000;
        trainSamplesLeads = trainSamplesDict;
        skipSamples = 5000;
        trainStartSample = 1;
        trainStopSample = trainSamplesDict;

        orgRespLeadsTrain = zeros(length(responseLeadNames),trainSamplesDict);
        for leads = 1:1:length(responseLeadNames)
            orgRespLeadsTrain(leads,:) = data.(dataName{dat}).(recordName{count}).(responseLeadNames{leads})(trainStartSample:trainStopSample);
            orgRespLeadsTrain(leads,:) = (orgRespLeadsTrain(leads,:)/max(abs(orgRespLeadsTrain(leads,:))));
        end


        orgPredLeadsTrain = zeros(length(predictorLeadNames),trainSamplesDict);
        for leads = 1:1:length(predictorLeadNames)
            orgPredLeadsTrain(leads,:) = data.(dataName{dat}).(recordName{count}).(predictorLeadNames{leads})(trainStartSample:trainStopSample);
            orgPredLeadsTrain(leads,:) = (orgPredLeadsTrain(leads,:)/max(abs(orgPredLeadsTrain(leads,:))));
        end


        %% Data for Dictionary Learning

        trainData = [orgRespLeadsTrain ; orgPredLeadsTrain];

        windowLength = 1;
        eVar = 1;
        EnergyTrain = zeros(length(allLeads),length(trainData(1,:))/windowLength);
        for winStart = 1:windowLength:length(trainData(1,:))
            for leads = 1:1:length(allLeads)
                EnergyTrain(leads,eVar) = sum(abs(trainData(leads,winStart:winStart + windowLength - 1)).^2);
            end
            eVar = eVar + 1;                
        end

        sortEnergy = sort(EnergyTrain(1,:),'ascend');

        energyLength = length(sortEnergy);
        energyLevels = floor(energyLength/maxDict);

        energyValOn = zeros(1,maxDict);
        energyValOff = zeros(1,maxDict);

        for var1=1:1:maxDict
            if var1 == 1
                energyValOn(var1) = 0;
                energyValOff(var1) = sortEnergy(energyLevels);
            elseif var1 == maxDict
                energyValOn(var1) = energyValOff(var1 - 1);
                energyValOff(var1) = max(sortEnergy);
            else
                energyValOn(var1) = energyValOff(var1 - 1);
                energyValOff(var1) = sortEnergy(var1*energyLevels);
            end
        end

        for var1=1:1:maxDict
            eVar = 1;
            start = 1;
            stop = windowLength;
            for winStart = 1:windowLength:length(trainData(1,:))
                if (EnergyTrain(1,eVar) > energyValOn(var1) && EnergyTrain(1,eVar) <= energyValOff(var1))
                    trainDataDict.Dict{var1}(:,start:stop) = trainData(:,winStart:winStart + windowLength - 1);
                    start = start + windowLength;
                    stop = stop + windowLength;
                end
                eVar = eVar + 1;                
            end     
        end

        for var1=1:1:maxDict
            normtrainDataDict = sqrt(sum(trainDataDict.Dict{var1}.^2, 1)); 
            trainDataDict.Dict{var1} = trainDataDict.Dict{var1}./repmat(normtrainDataDict, size(trainDataDict.Dict{var1}, 1), 1);            
            trainDataDict.Dict{var1} = trainDataDict.Dict{var1}(:,randperm(length(trainDataDict.Dict{var1}(1,:))));
        end

        %% Sample Skipping

        skipStartSample = trainStopSample + 1;
        skipStopSample = trainStopSample + skipSamples;
        skipSamplesDict = skipStopSample - skipStartSample + 1;
        skipSampleLeads = skipSamplesDict;


        %% Testing Data

        testStartSample = skipStopSample + 1;
        testStopSample = samples;
        testSamples = samples - skipStopSample;
        testSamplesDict = samples - skipStopSample;
        testSampleLeads = testSamplesDict;


        orgRespLeadsTest = zeros(length(responseLeadNames),testSamplesDict);
        for leads = 1:1:length(responseLeadNames)
            orgRespLeadsTest(leads,:) = data.(dataName{dat}).(recordName{count}).(responseLeadNames{leads})(testStartSample:testStopSample);
            orgRespLeadsTest(leads,:) = (orgRespLeadsTest(leads,:)/max(abs(orgRespLeadsTest(leads,:)))); 
            ecgPlots.(dataName{dat}).(recordName{count}).(responseLeadNames{leads}).orgTestLead = 2000*orgRespLeadsTest(leads,:);
        end


        orgPredLeadsTest = zeros(length(predictorLeadNames),testSamplesDict);
        for leads = 1:1:length(predictorLeadNames)
            orgPredLeadsTest(leads,:) = data.(dataName{dat}).(recordName{count}).(predictorLeadNames{leads})(testStartSample:testStopSample);
            orgPredLeadsTest(leads,:) = (orgPredLeadsTest(leads,:)/max(abs(orgPredLeadsTest(leads,:))));
            ecgPlots.(dataName{dat}).(recordName{count}).(predictorLeadNames{leads}).orgTestLead = 2000*orgPredLeadsTest(leads,:);
            ecgPlots.(dataName{dat}).(recordName{count}).(predictorLeadNames{leads}).leadData = 2000*orgPredLeadsTest(leads,:);
        end

        %% Dictionary Learning CVX Optimization and and Reconstruction
        
        iterMax = 1;
        weddDict = zeros(iterMax,length(responseLeadNames));
        RSquareDict = zeros(iterMax,length(responseLeadNames));
        corCoeffDict = zeros(iterMax,length(responseLeadNames));
        RMSEDict = zeros(iterMax,length(responseLeadNames));

        recData = zeros(length(responseLeadNames),testSamplesDict);
        for iterNo = 1:1:iterMax

            % dictionary learning

            disp('Starting to  train the dictionary');

            param.K   = dictionarySize;          % dictionary size
            param.numIteration = 10;      % number of iterations
            param.L = 2;               % % number of elements in each linear combination.
            param.errorGoal = 0.0000001;
            param.errorFlag = 0;        % if 1, arbitrary number of coefficients is used for representation of each signal    
%             param.errorFlag = 1;        % if 1, arbitrary number of coefficients is used for representation of each signal    
            param.preserveDCAtom = 0;
            param.InitializationMethod = 'DataElements';
            param.displayProgress = 1;

            for var1=1:1:maxDict
                fprintf('Training Dictionary No.%d \n',var1);
                [multDict.Org{var1}.Dictionary,multDict.Org{var1}.output] = KSVD(trainDataDict.Dict{var1}, param);
                multDict.Org{var1}.dictResp = multDict.Org{var1}.Dictionary(1:length(responseLeadNames),:);
                multDict.Org{var1}.dictPred = multDict.Org{var1}.Dictionary(length(responseLeadNames)+1:end,:);
            end

            disp('Starting CVX');

            optParam.respLength = length(responseLeadNames);
            optParam.predLength = length(predictorLeadNames);
            optParam.dictionarySize = dictionarySize;
            optParam.K_target = param.L;

            % normalizing

            for var1=1:1:maxDict                    
                normDictResp = sqrt(sum(multDict.Org{var1}.dictResp.^2, 1)); 
                multDict.Org{var1}.dictResp = multDict.Org{var1}.dictResp./repmat(normDictResp, size(multDict.Org{var1}.dictResp, 1), 1);

                normDictPred = sqrt(sum(multDict.Org{var1}.dictPred.^2, 1)); 
                multDict.Org{var1}.dictPred = multDict.Org{var1}.dictPred./repmat(normDictPred, size(multDict.Org{var1}.dictPred, 1), 1);

                clear normDictResp normDictPred;   
            end

            % CVX Optimization

            for var1=1:1:maxDict
                fprintf('Optimization No.%d \n',var1);
                [multDict.OptOut{var1}] = optimizationCVXNew(multDict.Org{var1}.dictResp,multDict.Org{var1}.dictPred,trainDataDict.Dict{var1},optParam);
            end

            disp('Starting to  test the dictionary');

            % normalizing optimized data
            for var1=1:1:maxDict                    
                normDictResp = sqrt(sum(multDict.OptOut{var1}.dictResp.^2, 1)); 
                multDict.OptOut{var1}.dictResp = multDict.OptOut{var1}.dictResp./repmat(normDictResp, size(multDict.OptOut{var1}.dictResp, 1), 1);

                normDictPred = sqrt(sum(multDict.OptOut{var1}.dictPred.^2, 1)); 
                multDict.OptOut{var1}.dictPred = multDict.OptOut{var1}.dictPred./repmat(normDictPred, size(multDict.OptOut{var1}.dictPred, 1), 1);

                clear normDictResp normDictPred;   
            end

            eVar = 1;
            EnergyTest = zeros(length(predictorLeadNames),length(orgPredLeadsTest(1,:))/windowLength);
            for winStart = 1:windowLength:length(orgPredLeadsTest(1,:))
                for leads = 1:1:length(predictorLeadNames)
                    EnergyTest(leads,eVar) = sum(abs(orgPredLeadsTest(leads,winStart:winStart + windowLength - 1)).^2);
                end
                eVar = eVar + 1;                
            end

            energyValOff(maxDict) = max(EnergyTest(1,:));

            eVar = 1;
            start = 1;
            stop = windowLength;
            for winStart = 1:windowLength:length(orgPredLeadsTest(1,:))
                for var1=1:1:maxDict
                    if (EnergyTest(1,eVar) > energyValOn(var1) && EnergyTest(1,eVar) <= energyValOff(var1))
                        sparseRepTestDictPred = OMP(multDict.OptOut{var1}.dictPred,orgPredLeadsTest(:,winStart:winStart + windowLength - 1), optParam.K_target);                            
                        sparseRepFullTestDictPred(:,winStart:winStart + windowLength - 1) = full(sparseRepTestDictPred);

                        sparseRepTestDictResp = OMP(multDict.OptOut{var1}.dictResp,orgRespLeadsTest(:,winStart:winStart + windowLength - 1), optParam.K_target);
                        sparseRepFullTestDictResp(:,winStart:winStart + windowLength - 1) = full(sparseRepTestDictResp);

                        modSparseRepFullTestDictPred(:,winStart:winStart + windowLength - 1) = multDict.OptOut{var1}.mapMat*sparseRepFullTestDictPred(:,winStart:winStart + windowLength - 1);
                        recData(:,start:stop) = multDict.OptOut{var1}.dictResp*modSparseRepFullTestDictPred(:,winStart:winStart + windowLength - 1);
                    end
                end
                start = start + windowLength;
                stop = stop + windowLength;
                eVar = eVar + 1;                
            end


            for leads = 1:1:length(responseLeadNames)
                recData(leads,:) = smooth(recData(leads,:),5,'moving');
                recData(leads,:) = (recData(leads,:)/max(abs(recData(leads,:))));
            end

            level = 7;
            for leads = 1:1:length(responseLeadNames)
                ecgPlots.(dataName{dat}).(recordName{count}).(responseLeadNames{leads}).leadData = 2000*recData(leads,:);
                weddDict(iterNo,leads) = waveletDist(2000*orgRespLeadsTest(leads,:),2000*recData(leads,:),level)*100;
                RMSEDict(iterNo,leads) = sqrt(sum((2000*orgRespLeadsTest(leads,:) - 2000*recData(leads,:)).^2)/length(2000*orgRespLeadsTest(leads,:)));
                tempCorc = corc(2000*orgRespLeadsTest(leads,:),2000*recData(leads,:));
                corCoeffDict(iterNo,leads) = tempCorc(1,2);
                RSquareDict(iterNo,leads) = (1 - (sum((2000*recData(leads,:) - 2000*orgRespLeadsTest(leads,:)).^2)/sum(2000*orgRespLeadsTest(leads,:).^2)))*100;
            end

        end
        
        weddDictAll(:,selDataCount) = mean(weddDict,1)';
        RSquareDictAll(:,selDataCount) = mean(RSquareDict,1)';
        corCoeffDictAll(:,selDataCount) = mean(corCoeffDict,1)';
        RMSEDictAll(:,selDataCount) = mean(RMSEDict,1)';

        weddPrec(selDataCount) = (sum(weddDictAll(7:8,selDataCount)) + sum(weddDictAll(10:12,selDataCount)))/5;
        RSquarePrec(selDataCount) = (sum(RSquareDictAll(7:8,selDataCount)) + sum(RSquareDictAll(10:12,selDataCount)))/5;
        corCoeffPrec(selDataCount) = (sum(corCoeffDictAll(7:8,selDataCount)) + sum(corCoeffDictAll(10:12,selDataCount)))/5;
        RMSEPrec(selDataCount) = (sum(RMSEDictAll(7:8,selDataCount)) + sum(RMSEDictAll(10:12,selDataCount)))/5;

        wedd(selDataCount) = mean(weddDictAll(:,selDataCount));        
        RSquare(selDataCount) = mean(RSquareDictAll(:,selDataCount));
        corCoeff(selDataCount) = mean(corCoeffDictAll(:,selDataCount));
        RMSE(selDataCount) = mean(RMSEDictAll(:,selDataCount));
        count = count + 1;
    end
end

toc
        
savefile1 = sprintf('MultipleCoupledDictionary3To12');
save(savefile1,'ecgPlots','weddDictAll','RSquareDictAll','corCoeffDictAll','RMSEDictAll','wedd','RSquare','corCoeff','RMSE','weddPrec','RSquarePrec','corCoeffPrec','RMSEPrec','patientID');

