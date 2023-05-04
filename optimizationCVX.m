function [OptOut] = optimizationCVX(dictResp,dictPred,trainData,optParam)

respLength = optParam.respLength;
predLength = optParam.predLength;
dictionarySize = optParam.dictionarySize;
K_target = optParam.K_target;


% normDictResp = sqrt(sum(dictResp.^2, 1)); 
% dictResp = dictResp./repmat(normDictResp, size(dictResp, 1), 1);
% 
% normDictPred = sqrt(sum(dictPred.^2, 1)); 
% dictPred = dictPred./repmat(normDictPred, size(dictPred, 1), 1);

sparseRepTrainDictResp = OMP(dictResp,trainData(1:respLength,:), K_target);
sparseRepFullTrainDictResp = full(sparseRepTrainDictResp);

sparseRepTrainDictPred = OMP(dictPred,trainData(respLength+1:end,:), K_target);
sparseRepFullTrainDictPred = full(sparseRepTrainDictPred);


lambdaMapMat = 0.002;
gammaMapMat = 0.9;
lambdaResp = 1;
lambdaPred = 1;

% cvx_begin
%     variable mapMat(dictionarySize,dictionarySize)
%     minimize( gammaMapMat*norm(mapMat*sparseRepFullTrainDictPred - sparseRepFullTrainDictResp,'fro') + lambdaMapMat*norm(mapMat,1) )
% cvx_end
% 
% OptOut.mapMat = full(mapMat);
% OptOut.dictResp = dictResp;
% OptOut.dictPred = dictPred;



cvx_begin
    variable mapMat(dictionarySize,dictionarySize);
    variable dictResp(respLength,dictionarySize);
    variable dictPred(predLength,dictionarySize);
    minimize( gammaMapMat*norm(mapMat*sparseRepFullTrainDictPred - sparseRepFullTrainDictResp,'fro') + lambdaMapMat*norm(mapMat,1) + ... 
        norm(trainData(1:respLength,:) - dictResp*sparseRepFullTrainDictResp,'fro') + lambdaResp*norm(sparseRepFullTrainDictResp,1) + ... 
        norm(trainData(respLength+1:end,:) - dictPred*sparseRepFullTrainDictPred,'fro') + lambdaPred*norm(sparseRepFullTrainDictPred,1) );
    subject to
        for var1 = 1:1:dictionarySize
            norm(dictResp(:,var1),2) <= 1;
            norm(dictPred(:,var1),2) <= 1;
%             norm(dictResp(:,var1),2) >= 1;
%             norm(dictPred(:,var1),2) >= 1;
%             norm(mapMat(:,var1),2) <= 1;
        end
cvx_end

OptOut.mapMat = full(mapMat);
OptOut.dictResp = full(dictResp);
OptOut.dictPred = full(dictPred);

