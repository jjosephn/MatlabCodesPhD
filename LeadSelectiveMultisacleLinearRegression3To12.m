clc;
close all;
clearvars;

tic

% Update the path for data base.

addpath('../PTBDatabaseFull')


format('short');
count=1;
leadNumber=15;
[num,txt,raw] = xlsread('dbListUpdated.xls');
[m,n]=size(txt);
samples=30000;
sampFreq = 1000;
% datMax=77;

%% Memory allocation

patientID{4,549} = [];
dataName{290}=[];
recordName{7} = [];
leadName{leadNumber + 1}=[];
xl{leadNumber}=[];
leadCobination{60}=[];
xyProjection{60}=[];
xyError{60}=[];
thetaVal=zeros(1,60);

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

predictionOrder{6-noInitPredLeads} = [];
predictionInd = 1;

defaultInitPredLeads = {'x1','x2'};
    
initSelPredLeads = {'x7','x8','x9','x10','x11','x12'};

leadNameOrg = {'LeadI', 'LeadII', 'LeadIII', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'X', 'Y', 'Z'};

independentLeads = [defaultInitPredLeads initSelPredLeads];

level = 7;

leads = 1;
while (leads <= length(independentLeads))
    corLeadWise.wav.(independentLeads{leads}) = 0;
    corLeadWise.lin.(independentLeads{leads}) = 0;
    rmseLeadWise.wav.(independentLeads{leads}) = 0;
    rmseLeadWise.lin.(independentLeads{leads}) = 0;
    weddLeadWise.wav.(independentLeads{leads}) = 0;
    weddLeadWise.lin.(independentLeads{leads}) = 0;
    RSquareLeadWise.wav.(independentLeads{leads}) = 0;
    RSquareLeadWise.lin.(independentLeads{leads}) = 0;
    leadCount.(independentLeads{leads}) = 0;
    leadCountWDD.(independentLeads{leads}) = 0;
    leads = leads + 1;
end

leadSelFreq = zeros(1,length(initSelPredLeads));
initRespLeadsAll{5} = [];
initRespLeadsAllOrg{5} = [];

% skipDataRecords = [9 16 17 39 49 78 89 122 127 138 162 173 182 190 208 223 230 262 283 284 306 313 315 318 324 326 329 330 338 356 359 362 364 367 376 396 404 405 ... 
%     407 418 440 444 448 452 468 499 513 534 537 538 543 545 546 547];

selDataCount = 0;
selLeadIndex = zeros(1,549 - length(skipDataRecords));
for dat = 1:1:length(fieldnames(data))
    
    count = 1;
    while (count <= length(fieldnames(data.(dataName{dat}))))

        selDataCount = selDataCount + 1;
        fprintf('patient-%d_record-%d \n',dat,count);

        trainSamples = 10000;
        trainSamplesLeads = trainSamples;
        skipSamples = 7000;
        trainStartSample = 1;
        trainStopSample = trainSamples;

        dataLead = data.(dataName{dat}).(recordName{count}).(leadName{1});
        len1 = length(dataLead);
        tt = 1/sampFreq:1/sampFreq:ceil(len1/sampFreq);
        timeAxis = tt(1:len1);

        initRespLeadsAll{5} = [];
        selLds = 1;
        while(selLds <= length(initSelPredLeads))
            selPredLead = initSelPredLeads{selLds};

            initRespLeadsAll{5} = [];
            var2=1;
            for var1=7:1:12
                if (strcmp(leadName{var1},selPredLead)==0)
                    initRespLeadsAll{var2} = leadName{var1};
                    var2 = var2 +1;
                end
            end

            lds = 1;
            initPredLeads = [defaultInitPredLeads selPredLead];
            while(lds <= length(initRespLeadsAll))
                initRespLeads = initRespLeadsAll(lds);
                newPredLeads = initPredLeads;
                totLeads = [newPredLeads initRespLeads];


                % Wavelet transform of predictor leads for training
                leads = 1;
                predXTr=zeros(trainSamplesLeads,length(newPredLeads));
                while (leads <= length(newPredLeads))
                    predXTrain.(newPredLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(trainStartSample:trainStopSample),level);
                    predXTr(:,leads) = data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(trainStartSample:trainStopSample)';
                    leads = leads + 1;
                end


                % Wavelet transform of response leads for training
                leads = 1;
                respXTr=zeros(trainSamplesLeads,length(initRespLeads));
                while (leads <= length(initRespLeads))
                    respXTrain.(initRespLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(trainStartSample:trainStopSample),level);
                    respXTr(:,leads) = data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(trainStartSample:trainStopSample)';
                    leads = leads + 1;
                end

                % linear regression 
                bLin.(initSelPredLeads{selLds}).(initRespLeads{1}) = regress(respXTr,[ones(size(predXTr(:,1))) predXTr]);

                % linear regression over wavelet coeff
                for l=1:1:level
                    pred=zeros(length(newPredLeads),length(predXTrain.(newPredLeads{1}).dc{l}));
                    for leads=1:1:length(newPredLeads)
                        pred(leads,:)=predXTrain.(newPredLeads{leads}).dc{l};
                    end
                    bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).dc{l} = regress(respXTrain.(initRespLeads{1}).dc{l}',[ones(size(pred(1,:)')) pred']);
                    clear pred;
                    if (l == level)
                        pred=zeros(length(newPredLeads),length(predXTrain.(newPredLeads{1}).ac{1}));
                        for leads=1:1:length(newPredLeads)
                            pred(leads,:)=predXTrain.(newPredLeads{leads}).ac{1};
                        end
                        bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).ac{1} = regress(respXTrain.(initRespLeads{1}).ac{1}',[ones(size(pred(1,:)')) pred']);
                        clear pred;
                    end
                end

                predictionOrder{predictionInd} = initRespLeads{1};
                predictionInd = predictionInd + 1;

                lds = lds + 1;        
            end
            clear initRespLeadsAll;
            selLds = selLds + 1;
        end


    %% Lead Selection Stage

        initRespLeadsAll{5} = [];
        weddWavAvg = zeros(1,length(initSelPredLeads));
        weddLinAvg = zeros(1,length(initSelPredLeads));

        RSquareWavAvg = zeros(1,length(initSelPredLeads));
        RSquareLinAvg = zeros(1,length(initSelPredLeads));

        corCoeffWavAvg = zeros(1,length(initSelPredLeads));
        corCoeffLinAvg = zeros(1,length(initSelPredLeads));

        RMSEWavAvg = zeros(1,length(initSelPredLeads));
        RMSELinAvg = zeros(1,length(initSelPredLeads));

        DSSWavAvg = zeros(1,length(initSelPredLeads));

        selLds = 1;
        while(selLds <= length(initSelPredLeads))
            selPredLead = initSelPredLeads{selLds};

            initRespLeadsAll{5} = [];
            var2=1;
            for var1=7:1:12
                if (strcmp(leadName{var1},selPredLead)==0)
                    initRespLeadsAll{var2} = leadName{var1};
                    var2 = var2 +1;
                end
            end

            skipStartSample = trainStopSample + 1;
            skipStopSample = trainStopSample + skipSamples;


            weddAvgWav = 0;
            weddAvgLin = 0;

            RSquareAvgWav = 0;
            RSquareAvgLin = 0;

            corCoeffAvgWav = 0;
            corCoeffAvgLin = 0;

            RMSEAvgWav = 0;
            RMSEAvgLin = 0;

            DSSAvgWav = 0;

            lds = 1;
            initPredLeads = [defaultInitPredLeads selPredLead];
            while(lds <= length(initRespLeadsAll))
                initRespLeads = initRespLeadsAll(lds);

                weddOrgWavLin = 0;
                weddOrgLin = 0;
                corOrgWavLin = zeros(2,2);
                corOrgLin = zeros(2,2);
                RMSEOrgWavLin = 0;
                RMSEOrgLin = 0;
%                     RSquareWavLin = 0;
%                     RSquareLin = 0;

                DSSInit = 0;

                newPredLeads = initPredLeads;
                totLeads = [newPredLeads initRespLeads];

                % Wavelet transform of predictor leads for testing
                leads = 1;
                intPredXTe=zeros(skipSamples,length(newPredLeads));
                while (leads <= length(newPredLeads))
                    intPredXTest.(newPredLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(skipStartSample:skipStopSample),level);
                    intPredXTe(:,leads) = data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(skipStartSample:skipStopSample)';
                    leads = leads + 1;
                end

                % Wavelet transform of response leads for testing
                leads = 1;
                intRespXTe=zeros(skipSamples,length(initRespLeads));
                while (leads <= length(initRespLeads))
                    intRespXTest.(initRespLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(skipStartSample:skipStopSample),level);
                    intRespXTe(:,leads) = data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(skipStartSample:skipStopSample)';
                    leads = leads + 1;
                end

                % linear regression reconstruction
                sum = bLin.(initSelPredLeads{selLds}).(initRespLeads{1})(1)*ones(skipSamples,1);
                for leads=1:1:length(newPredLeads)
                    sum = sum + bLin.(initSelPredLeads{selLds}).(initRespLeads{1})(leads+1)*intPredXTe(:,leads);
                end

                intRecXLin.(initRespLeads{1}) = sum';
                intRecXLin.(initRespLeads{1})=(intRecXLin.(initRespLeads{1})/max(abs(intRecXLin.(initRespLeads{1}))));


                intRecModData.Lin.(dataName{dat}).(recordName{count}).(initRespLeads{1}) = intRecXLin.(initRespLeads{1});

                clear sum;

                % linear regression Wavelet reconstruction

                for l=1:1:level
                    sum = bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).dc{l}(1);
                    for leads=1:1:length(newPredLeads)
                        sum = sum + bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).dc{l}(leads+1)*intPredXTest.(newPredLeads{leads}).dc{l}';
                    end
                    intPredXTest.(initRespLeads{1}).dc{l} = sum';
                    clear sum;
                    if (l == level)
                        sum = bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).ac{1}(1);
                        for leads=1:1:length(newPredLeads)
                            sum = sum + bWavLin.(initSelPredLeads{selLds}).(initRespLeads{1}).ac{1}(leads+1)*intPredXTest.(newPredLeads{leads}).ac{1}';
                        end
                        intPredXTest.(initRespLeads{1}).ac{1} = sum';
                        clear sum;
                    end
                end

                intRecXWavLin.(initRespLeads{1}) = invWltTfm(intPredXTest.(initRespLeads{1}),skipSamples,level);
                intRecXWavLin.(initRespLeads{1}) = (intRecXWavLin.(initRespLeads{1})/max(abs(intRecXWavLin.(initRespLeads{1}))));

                intRecModData.wavLin.(dataName{dat}).(recordName{count}).(initRespLeads{1}) = intRecXWavLin.(initRespLeads{1});

                % For Distortion Measure
                orgTestLead = 2000*data.(dataName{dat}).(recordName{count}).(initRespLeads{1})(skipStartSample:skipStopSample);
                recWavLin = 2000*intRecModData.wavLin.(dataName{dat}).(recordName{count}).(initRespLeads{1});
                recLin = 2000*intRecModData.Lin.(dataName{dat}).(recordName{count}).(initRespLeads{1});

                weddOrgWavLin = weddOrgWavLin + waveletDist(orgTestLead,recWavLin,level)*100;
                weddOrgLin = weddOrgLin + waveletDist(orgTestLead,recLin,level)*100;

                corOrgWavLin = corOrgWavLin + corc(orgTestLead,recWavLin);
                corOrgLin = corOrgLin + corc(orgTestLead,recLin);

                RMSEOrgWavLin = RMSEOrgWavLin + sqrt(sum((orgTestLead - recWavLin).^2)/length(orgTestLead));
                RMSEOrgLin = RMSEOrgLin + sqrt(sum((orgTestLead - recLin).^2)/length(orgTestLead));

                RSquareWavLin = (1 - (sum((recWavLin - orgTestLead).^2)/sum(orgTestLead.^2)))*100;
                RSquareLin = (1 - (sum((recLin - orgTestLead).^2)/sum(orgTestLead.^2)))*100;

                intWedd.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).wedd = weddOrgWavLin;
                intWedd.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).wedd = weddOrgLin;

                intRSquare.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).wedd = RSquareWavLin;
                intRSquare.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).wedd = RSquareLin;

                weddAvgWav = weddAvgWav +  weddOrgWavLin;
                weddAvgLin = weddAvgLin + weddOrgLin;

                RSquareAvgWav = RSquareAvgWav +  RSquareWavLin;
                RSquareAvgLin = RSquareAvgLin + RSquareLin;

                corCoeffAvgWav = corCoeffAvgWav + corOrgWavLin(1,2);
                corCoeffAvgLin = corCoeffAvgLin + corOrgLin(1,2);

                RMSEAvgWav = RMSEAvgWav + RMSEOrgWavLin;
                RMSEAvgLin = RMSEAvgLin + RMSEOrgLin;

                DSSAvgWav = DSSAvgWav + DSSInit;

                lds = lds + 1;    
            end

            weddWavAvg(selLds) = weddAvgWav/length(initRespLeadsAll);
            weddLinAvg(selLds) = weddAvgLin/length(initRespLeadsAll);

            RSquareWavAvg(selLds) = RSquareAvgWav/length(initRespLeadsAll);
            RSquareLinAvg(selLds) = RSquareAvgLin/length(initRespLeadsAll);

            corCoeffWavAvg(selLds) = corCoeffAvgWav/length(initRespLeadsAll);
            corCoeffLinAvg(selLds) = corCoeffAvgLin/length(initRespLeadsAll);

            RMSEWavAvg(selLds) = RMSEAvgWav/length(initRespLeadsAll);
            RMSELinAvg(selLds) = RMSEAvgLin/length(initRespLeadsAll);


            clear intRecXLin intRecXWavLin intRespXTe intRespXTest intPredXTe intPredXTest initRespLeadsAll;
            selLds = selLds + 1;
        end

        weightDSS = dec2bin(1:1:2^(3)-1)-'0';

        DSSVal = zeros(length(weightDSS(:,1)),length(initSelPredLeads));
        maxDSSVal = zeros(1,length(initSelPredLeads));
        maxDSSInd = zeros(1,length(initSelPredLeads));


        for var1 = 1:1:length(weightDSS(:,1))
            w1 = weightDSS(var1,1);
            w2 = weightDSS(var1,2);
            w3 = weightDSS(var1,3);

            DSSVal(var1,:) = w1*(100./weddWavAvg) + w2*corCoeffWavAvg + w3*RSquareWavAvg/100;
            [maxDSSVal(var1), maxDSSInd(var1)] = max(DSSVal(var1,:));

        end

        %% Testing

        clear intWedd intRecModData;

        DSSValue.(dataName{dat}).(recordName{count}).DSSVal = DSSVal;

        initRespLeadsAll{5} = [];

        selLeadLoc = mode(maxDSSInd);

        selLeadIndex(selDataCount) = selLeadLoc;

        leadSelFreq(selLeadLoc) = leadSelFreq(selLeadLoc) + 1;

        testStartSample = skipStopSample + 1;
        testStopSample = samples;
        testSamples = samples - skipStopSample;

        selPredLead = initSelPredLeads{selLeadLoc};

        initRespLeadsAll{5} = [];
        initRespLeadsAllOrg{5} = [];
        var2=1;
        for var1=7:1:12
            if (strcmp(leadName{var1},selPredLead)==0)
                initRespLeadsAll{var2} = leadName{var1};
                initRespLeadsAllOrg{var2} = leadNameOrg{var1};
                var2 = var2 +1;
            end
        end

        weddAvgWav = 0;
        weddAvgLin = 0;

        RSquareAvgWav = 0;
        RSquareAvgLin = 0;

        corCoeffAvgWav = 0;
        corCoeffAvgLin = 0;

        RMSEAvgWav = 0;
        RMSEAvgLin = 0;


        lds = 1;
        initPredLeads = [defaultInitPredLeads selPredLead];
        for leads = 1:1:length(initPredLeads)
            orgTestLead = 2000*data.(dataName{dat}).(recordName{count}).(initPredLeads{leads})(testStartSample:testStopSample);
            ecgPlots.(dataName{dat}).(recordName{count}).(initPredLeads{leads}).orgTestLead = orgTestLead;
            ecgPlots.(dataName{dat}).(recordName{count}).(initPredLeads{leads}).leadData = orgTestLead;

        end

        for leads = 3:1:6
            orgTestLead = 2000*data.(dataName{dat}).(recordName{count}).(leadName{leads})(testStartSample:testStopSample);
            ecgPlots.(dataName{dat}).(recordName{count}).(leadName{leads}).orgTestLead = orgTestLead;
            ecgPlots.(dataName{dat}).(recordName{count}).(leadName{leads}).leadData = orgTestLead;
        end

        while(lds <= length(initRespLeadsAll))
            initRespLeads = initRespLeadsAll(lds);

            initRespLeadsOrg = initRespLeadsAllOrg(lds);

            weddOrgWavLin = 0;
            weddOrgLin = 0;
            corOrgWavLin = zeros(2,2);
            corOrgLin = zeros(2,2);
            RMSEOrgWavLin = 0;
            RMSEOrgLin = 0;
%                 RSquareWavLin = 0;
%                 RSquareLin = 0;
            newPredLeads = initPredLeads;
            totLeads = [newPredLeads initRespLeads];

            % Wavelet transform of predictor leads for testing
            leads = 1;
            predXTe=zeros(testSamples,length(newPredLeads));
            while (leads <= length(newPredLeads))
                predXTest.(newPredLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(testStartSample:testStopSample),level);
                predXTe(:,leads) = data.(dataName{dat}).(recordName{count}).(newPredLeads{leads})(testStartSample:testStopSample)';
                leads = leads + 1;
            end

            % Wavelet transform of response leads for testing
            leads = 1;
            respXTe=zeros(testSamples,length(initRespLeads));
            while (leads <= length(initRespLeads))
                respXTest.(initRespLeads{leads}) = wltTfm(data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(testStartSample:testStopSample),level);
                respXTe(:,leads) = data.(dataName{dat}).(recordName{count}).(initRespLeads{leads})(testStartSample:testStopSample)';
                leads = leads + 1;
            end

            % linear regression reconstruction
            sum = bLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1})(1)*ones(testSamples,1);
            for leads=1:1:length(newPredLeads)
                sum = sum + bLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1})(leads+1)*predXTe(:,leads);
            end

            RecXLin.(initRespLeads{1}) = sum';
            RecXLin.(initRespLeads{1})=(RecXLin.(initRespLeads{1})/max(abs(RecXLin.(initRespLeads{1}))));


            recModData.Lin.(dataName{dat}).(recordName{count}).(initRespLeads{1}) = RecXLin.(initRespLeads{1});

            clear sum;

            % linear regression Wavelet reconstruction

            for l=1:1:level
                sum = bWavLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1}).dc{l}(1);
                for leads=1:1:length(newPredLeads)
                    sum = sum + bWavLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1}).dc{l}(leads+1)*predXTest.(newPredLeads{leads}).dc{l}';
                end
                RecXWavBand.(initRespLeads{1}).dc{l} = sum';
                clear sum;
                if (l == level)
                    sum = bWavLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1}).ac{1}(1);
                    for leads=1:1:length(newPredLeads)
                        sum = sum + bWavLin.(initSelPredLeads{selLeadLoc}).(initRespLeads{1}).ac{1}(leads+1)*predXTest.(newPredLeads{leads}).ac{1}';
                    end
                    RecXWavBand.(initRespLeads{1}).ac{1} = sum';
                    clear sum;
                end
            end

            if (isempty(find(skipDataRecords==selDataCount, 1)))
                leadCountWDD.(initRespLeads{1}) = leadCountWDD.(initRespLeads{1}) + 1;
            end

            leadCount.(initRespLeads{1}) = leadCount.(initRespLeads{1}) + 1;

            RecXWavLin.(initRespLeads{1}) = invWltTfm(RecXWavBand.(initRespLeads{1}),testSamples,level);
            RecXWavLin.(initRespLeads{1}) = (RecXWavLin.(initRespLeads{1})/max(abs(RecXWavLin.(initRespLeads{1}))));

            recModData.wavLin.(dataName{dat}).(recordName{count}).(initRespLeads{1}) = RecXWavLin.(initRespLeads{1});

            % For Distortion Measure
            orgTestLead = 2000*data.(dataName{dat}).(recordName{count}).(initRespLeads{1})(testStartSample:testStopSample);
            recWavLin = 2000*recModData.wavLin.(dataName{dat}).(recordName{count}).(initRespLeads{1});
            recLin = 2000*recModData.Lin.(dataName{dat}).(recordName{count}).(initRespLeads{1});

            ecgPlots.(dataName{dat}).(recordName{count}).(initRespLeads{1}).orgTestLead = orgTestLead;
            ecgPlots.(dataName{dat}).(recordName{count}).(initRespLeads{1}).leadData = recWavLin;

            weddOrgWavLin = weddOrgWavLin + waveletDist(orgTestLead,recWavLin,level)*100;
            weddOrgLin = weddOrgLin + waveletDist(orgTestLead,recLin,level)*100;

            corOrgWavLin = corOrgWavLin + corc(orgTestLead,recWavLin);
            corOrgLin = corOrgLin + corc(orgTestLead,recLin);

            RMSEOrgWavLin = RMSEOrgWavLin + sqrt(sum((orgTestLead - recWavLin).^2)/length(orgTestLead));
            RMSEOrgLin = RMSEOrgLin + sqrt(sum((orgTestLead - recLin).^2)/length(orgTestLead));

            RSquareWavLin = (1 - (sum((recWavLin - orgTestLead).^2)/sum(orgTestLead.^2)))*100;
            RSquareLin = (1 - (sum((recLin - orgTestLead).^2)/sum(orgTestLead.^2)))*100;

            corCoeff.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).cor = corOrgWavLin(1,2);
            corLeadWise.wav.(initRespLeads{1}) = corLeadWise.wav.(initRespLeads{1}) + corOrgWavLin(1,2);
            corCoeff.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).cor = corOrgLin(1,2);
            corLeadWise.lin.(initRespLeads{1}) = corLeadWise.lin.(initRespLeads{1}) + corOrgLin(1,2);

            RMSE.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).rmse = RMSEOrgWavLin;
            rmseLeadWise.wav.(initRespLeads{1}) = rmseLeadWise.wav.(initRespLeads{1}) + RMSEOrgWavLin;
            RMSE.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).rmse = RMSEOrgLin;
            rmseLeadWise.lin.(initRespLeads{1}) = rmseLeadWise.lin.(initRespLeads{1}) + RMSEOrgLin;

            wedd.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).wedd = weddOrgWavLin;
            weddLeadWise.wav.(initRespLeads{1}) = weddLeadWise.wav.(initRespLeads{1}) + weddOrgWavLin;
            wedd.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).wedd = weddOrgLin;
            weddLeadWise.lin.(initRespLeads{1}) = weddLeadWise.lin.(initRespLeads{1}) + weddOrgLin;

            RSquare.(dataName{dat}).(recordName{count}).wav.(initRespLeads{1}).rsquare = RSquareWavLin;
            RSquareLeadWise.wav.(initRespLeads{1}) = RSquareLeadWise.wav.(initRespLeads{1}) + RSquareWavLin;
            RSquare.(dataName{dat}).(recordName{count}).lin.(initRespLeads{1}).rsquare = RSquareLin;
            RSquareLeadWise.lin.(initRespLeads{1}) = RSquareLeadWise.lin.(initRespLeads{1}) + RSquareLin;

            weddAvgWav = weddAvgWav +  weddOrgWavLin;
            weddAvgLin = weddAvgLin + weddOrgLin;

            RSquareAvgWav = RSquareAvgWav +  RSquareWavLin;
            RSquareAvgLin = RSquareAvgLin + RSquareLin;

            corCoeffAvgWav = corCoeffAvgWav + corOrgWavLin(1,2);
            corCoeffAvgLin = corCoeffAvgLin + corOrgLin(1,2);

            RMSEAvgWav = RMSEAvgWav + RMSEOrgWavLin;
            RMSEAvgLin = RMSEAvgLin + RMSEOrgLin;
            dataNumber = dataNumber + 1;

            lds = lds + 1;
        end

        numRespLead = 5;

        wedd.(dataName{dat}).(recordName{count}).wav.weddAvg = weddAvgWav/numRespLead;
        wedd.(dataName{dat}).(recordName{count}).lin.weddAvg = weddAvgLin/numRespLead;

        RSquare.(dataName{dat}).(recordName{count}).wav.RSquareAvg = RSquareAvgWav/numRespLead;
        RSquare.(dataName{dat}).(recordName{count}).lin.RSquareAvg = RSquareAvgLin/numRespLead;

        corCoeff.(dataName{dat}).(recordName{count}).wav.corCoeffAvg = corCoeffAvgWav/numRespLead;
        corCoeff.(dataName{dat}).(recordName{count}).lin.corCoeffAvg = corCoeffAvgLin/numRespLead;

        RMSE.(dataName{dat}).(recordName{count}).wav.RMSEAvg = RMSEAvgWav/numRespLead;
        RMSE.(dataName{dat}).(recordName{count}).lin.RMSEAvg = RMSEAvgLin/numRespLead;
        
        count = count + 1;
        
    end             %%% while ends here
       
end    
% 
toc

weddAvgAll = zeros(1,selDataCount);
RSquareAvgAll = zeros(1,selDataCount);
corCoeffAvgAll = zeros(1,selDataCount);
RMSEAvgAll = zeros(1,selDataCount);
var1 = 1;
for dat = 1:1:length(dataName)
    for count = 1:1:length(fieldnames(data.(dataName{dat})))
        weddAvgAll(var1) = wedd.(dataName{dat}).(recordName{count}).wav.weddAvg;
        RSquareAvgAll(var1) = RSquare.(dataName{dat}).(recordName{count}).wav.RSquareAvg;
        corCoeffAvgAll(var1) = corCoeff.(dataName{dat}).(recordName{count}).wav.corCoeffAvg;
        RMSEAvgAll(var1) = RMSE.(dataName{dat}).(recordName{count}).wav.RMSEAvg;        
        var1 = var1 + 1;
    end
end



numRespLead = length(initSelPredLeads);
bandCorVal = zeros(numRespLead,numRespLead);
weddWavVal = zeros(1,numRespLead);
weddLinVal = zeros(1,numRespLead);
RSquareWavVal = zeros(1,numRespLead);
RSquareLinVal = zeros(1,numRespLead);
corWavVal = zeros(1,numRespLead);
corLinVal = zeros(1,numRespLead);
rmseWavVal = zeros(1,numRespLead);
rmseLinVal = zeros(1,numRespLead);

leads = 1;
while (leads <= length(initSelPredLeads))
    
    corLeadWise.wav.(initSelPredLeads{leads}) = corLeadWise.wav.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    corWavVal(leads) = corLeadWise.wav.(initSelPredLeads{leads});
    
    rmseLeadWise.wav.(initSelPredLeads{leads}) = rmseLeadWise.wav.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    rmseWavVal(leads) = rmseLeadWise.wav.(initSelPredLeads{leads});
    
    weddLeadWiseVal.wav.(initSelPredLeads{leads}) = weddLeadWise.wav.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    weddWavVal(leads) = weddLeadWiseVal.wav.(initSelPredLeads{leads});
    
    RSquareLeadWiseVal.wav.(initSelPredLeads{leads}) = RSquareLeadWise.wav.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    RSquareWavVal(leads) = RSquareLeadWiseVal.wav.(initSelPredLeads{leads});
    
    
    
    corLeadWise.lin.(initSelPredLeads{leads}) = corLeadWise.lin.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    corLinVal(leads) = corLeadWise.lin.(initSelPredLeads{leads});
    
    rmseLeadWise.lin.(initSelPredLeads{leads}) = rmseLeadWise.lin.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    rmseLinVal(leads) = rmseLeadWise.lin.(initSelPredLeads{leads});
    
    weddLeadWiseVal.lin.(initSelPredLeads{leads}) = weddLeadWise.lin.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    weddLinVal(leads) = weddLeadWiseVal.lin.(initSelPredLeads{leads});
    
    RSquareLeadWiseVal.lin.(initSelPredLeads{leads}) = RSquareLeadWise.lin.(initSelPredLeads{leads}) /leadCount.(initSelPredLeads{leads});
    RSquareLinVal(leads) = RSquareLeadWiseVal.lin.(initSelPredLeads{leads});
    
    leads = leads + 1;
end

savefile1 = sprintf('LeadSelectiveMultisacleLinearRegression3To12');
save(savefile1,'ecgPlots', 'wedd', 'weddAvgAll', 'RSquareAvgAll' ,'corCoeffAvgAll' ,'RMSEAvgAll','leadCount','leadCountWDD','selLeadIndex','DSSValue');

