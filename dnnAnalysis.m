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

predictorLeadNames = {'x1','x2','x9'};
responseLeadNames = {'x7','x8','x10','x11','x12'};
otherLeads = {'x1','x2','x3','x4','x5','x6','x9'};

allLeads = {'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12'};

totData = 549;

weddDNN = zeros(totData,length(responseLeadNames));
RMSEDNN = zeros(totData,length(responseLeadNames));
corCoeffDNN = zeros(totData,length(responseLeadNames));
RSquareDNN = zeros(totData,length(responseLeadNames));

% 

level = 7;
totData = 549;
dataCount = 1;
for dat = 1:1:length(fieldnames(data))
    count = 1;
    while (count <= length(fieldnames(data.(dataName{dat}))))

        filename = sprintf('/home/jjnalli/Documents/Python/DNN/dataOutLSTM16V3/dataOut%d.mat',dataCount-1);
%         filename = sprintf('/home/jjnalli/Documents/Python/DNN/dataOutGRU16V3/dataOut%d.mat',dataCount-1);
%         filename = sprintf('/home/jjnalli/Documents/Python/DNN/dataOutSimpleRNN16V3/dataOut%d.mat',dataCount-1);
        load(filename)

        startData = length(prediction) - 5000 + 1;
        stopData = length(prediction);
        testStartSample = 25001;
        testStopSample = 30000;
        
        for leads = 1:1:length(otherLeads)
            orgTestLeadTemp = 2000*data.(dataName{dat}).(recordName{count}).(otherLeads{leads})(testStartSample:testStopSample);
            ecgPlots.(dataName{dat}).(recordName{count}).(otherLeads{leads}).leadData = orgTestLeadTemp;
            ecgPlots.(dataName{dat}).(recordName{count}).(otherLeads{leads}).orgTestLead = orgTestLeadTemp;
        end

        for leads = 1:1:length(responseLeadNames)
            predictionTemp = double(2000*prediction(leads,startData:stopData));
            originalTemp = double(2000*original(leads,startData:stopData));
            ecgPlots.(dataName{dat}).(recordName{count}).(responseLeadNames{leads}).leadData = predictionTemp;
            ecgPlots.(dataName{dat}).(recordName{count}).(responseLeadNames{leads}).orgTestLead = originalTemp;
            weddDNN(dataCount,leads) = waveletDist(originalTemp, predictionTemp,level)*100;
            RMSEDNN(dataCount,leads) = sqrt(sum((originalTemp -  predictionTemp).^2)/length(originalTemp));
            tempCorc = corc(originalTemp,predictionTemp);
            corCoeffDNN(dataCount,leads) = tempCorc(1,2);
            RSquareDNN(dataCount,leads) = (1 - (sum(( predictionTemp - originalTemp).^2)/sum(originalTemp.^2)))*100;
        end
        
        dataCount = dataCount + 1;
        count = count + 1;
    end
end

weddDNN = weddDNN';
RMSEDNN = RMSEDNN';
corCoeffDNN = corCoeffDNN';
RSquareDNN = RSquareDNN';

weddDNNAvg = mean(weddDNN);
RMSEDNNAvg = mean(RMSEDNN);
corCoeffDNNAvg = mean(corCoeffDNN);
RSquareDNNAvg = mean(RSquareDNN);

% LSTM16V3
% GRU16V3
% SimpleRNN16V3

savefile1 = sprintf('dnnPythonTest3SimpleRNN16V3');
save(savefile1,'ecgPlots','weddDNN','RSquareDNN','corCoeffDNN','RMSEDNN','weddDNNAvg','RSquareDNNAvg','corCoeffDNNAvg','RMSEDNNAvg','patientID');



