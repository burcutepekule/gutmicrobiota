clear all;close all;clc;

% PART 0 : READ FILES
trtTime      = 1000; %total time of treatment in days
trtCountVec  = 1:20; %number of treatment courses
resC0        = 0;
resC1        = 6;
hBin         = 1;
infCount     = 0;
trtInfoKeep  = [];
binSeqKeep   = [];
numDataPts   = [];
root_0       = './SIM_RESULTS/';
for k=1:length(trtCountVec)
    % for k=1
    trtCount   = trtCountVec(k); %once in every 30 days
    filname    = ['TL_' num2str(trtTime) '_N_' num2str(trtCount)];
    directory  = [root_0 filname '/'];
    trtInit =[];trtLen =[];samplePops=[];c0counter=[];com=[];
    for s=1:size(dir([directory '*.txt']),1)/4
        [k s]
        filname01 = [directory 'schedule_trtInit.txt'];
        filname02 = [directory 'schedule_trtLen.txt'];
        filname03 = [directory 'samplePops.txt'];
        filname04 = [directory 'extCounter.txt'];
        if((exist(filname01, 'file')+exist(filname02, 'file')+exist(filname03, 'file')+exist(filname04, 'file')==8))
            trtInit     = [trtInit;readTxtFile(['schedule_trtInit'], directory)];
            trtLen      = [trtLen; readTxtFile(['schedule_trtLen'], directory)];
            samplePops  = [samplePops; readTxtFile(['samplePops'], directory)];
            c0counter   = [c0counter;  readTxtFile(['extCounter'], directory)];
            maxSize     = min([size(trtInit,1),size(trtLen,1),size(samplePops,1),size(c0counter,1)]);
            trtInit     = trtInit(1:maxSize,:);
            trtLen      = trtLen(1:maxSize,:);
            samplePops  = samplePops(1:maxSize,:);
            c0counter   = c0counter(1:maxSize,:);
        end
        
    end
    numDataPts{k}  = size(samplePops,1);
    if(size(samplePops,1)>0)
        samplePopsTemp = samplePops;
        c0counterTemp  = c0counter;
        if(trtCount == 1)
            trtInitTemp = trtInit;
            trtLenTemp  = trtLen;
            trtFirstExp = trtInitTemp(:,1);
            trtLastExp  = trtInitTemp(:,end)+trtLenTemp(:,end);
            trtLastDur  = trtLenTemp(:,end);
            meanTrtLen  = trtLenTemp;
            stdTrtLen   = zeros(size(trtLenTemp));
                        [~,~,binSeqTemp] = trt2bin(trtInitTemp,trtLenTemp,trtTime,trtLastExp,samplePopsTemp);
            com=[];
            for cm=1:size(binSeqTemp,1)
                com(cm) = centerOfMass(binSeqTemp(cm,:));
            end
        else
            trtInitTemp = trtInit;
            trtLenTemp  = trtLen;
            trtFirstExp = trtInitTemp(:,1);
            trtLastExp  = trtInitTemp(:,end)+trtLenTemp(:,end);
            trtLastDur  = trtLenTemp(:,end);
            meanTrtLen  = mean(trtLenTemp,2);
            stdTrtLen   = std(trtLenTemp,[],2);
            [~,~,binSeqTemp] = trt2bin(trtInitTemp,trtLenTemp,trtTime,trtLastExp,samplePopsTemp);
            com=[];
            for cm=1:size(binSeqTemp,1)
                com(cm) = centerOfMass(binSeqTemp(cm,:));
            end
        end
        trtInfoKeep = [trtInfoKeep; repmat(trtCount,size(samplePopsTemp,1),1) sum(trtLenTemp,2) trtFirstExp trtLastExp trtLastDur meanTrtLen stdTrtLen com' samplePopsTemp];
        binSeqKeep  = [binSeqKeep; binSeqTemp samplePopsTemp];
    end
    size(trtInfoKeep)
end

obsTime   = trtTime+360+5;
popScale  = 1e13;
dt        = 0.01;
indsSC0   = 10:5:130; %%%%% CHANGE THIS IF NEW PREDICTORS ADDED %%%%%
indsRC0   = indsSC0+1;
indsSC1   = indsSC0+2;
indsRC1   = indsSC0+3;

binSeq    = binSeqKeep(:,1:trtTime);
prevs     = trtInfoKeep(:,indsRC0)./(trtInfoKeep(:,indsSC0)+trtInfoKeep(:,indsRC0));% prevelances
replic    = 25; %%TDF=0 added
maxtdf    = 15*(replic-1);
binSeqTdf = zeros(replic*size(binSeq,1),maxtdf+trtTime);
probsTdf  = zeros(replic*size(binSeq,1),1);
sz        = size(binSeq,1);
for tdf=1:25
    tdf
    binSeqTdf((tdf-1)*sz+1:tdf*sz,maxtdf-15*(tdf-1)+1:maxtdf-15*(tdf-1)+trtTime) = binSeq;
    binSeqTdf((tdf-1)*sz+1:tdf*sz,maxtdf-15*(tdf-1)+trtTime+1:end)       = -1*ones(sz,15*(tdf-1));
    probsTdf((tdf-1)*sz+1:tdf*sz)                                        = prevs(:,tdf)';
end
save('allData','-v7.3');

obsTime      = trtTime+360+5;
popScale     = 1e13;
dt           = 0.01;
indsSC0 = 10:5:130; %%%%% CHANGE THIS IF NEW PREDICTORS ADDED %%%%%
indsRC0 = indsSC0+1;
indsSC1 = indsSC0+2;
indsRC1 = indsSC0+3;
resInfLoadC0  = [trtInfoKeep(:,1:indsSC0(1)-2) trtInfoKeep(:,indsRC0)];
resInfLoadC1  = [trtInfoKeep(:,1:indsSC0(1)-2) trtInfoKeep(:,indsRC1)];
resInfProbsC0 = [trtInfoKeep(:,1:indsSC0(1)-2) trtInfoKeep(:,indsRC0)./(trtInfoKeep(:,indsSC0)+trtInfoKeep(:,indsRC0))];
resInfProbsC1 = [trtInfoKeep(:,1:indsSC0(1)-2) trtInfoKeep(:,indsRC1)./(trtInfoKeep(:,indsSC1)+trtInfoKeep(:,indsRC1))];
resInfProbsC0(isnan(resInfProbsC0))=0;

% PART 2 : SAVE PREDICTORS
tempProbsC0    = resInfProbsC0;
%%%PREDICTORS%%%
% trtCount sum(trtLenTemp,2) trtFirstExp trtLastExp trtLastDur samplePopsTemp
numTrts        = tempProbsC0(:,1); %predictor 1 -> num of treatments
inten          = tempProbsC0(:,2)./(tempProbsC0(:,4)-tempProbsC0(:,3)); %predictor 2 -> intensity
days2trt       = tempProbsC0(:,3); %predictor 3 -> t1
totalDays      = tempProbsC0(:,2); %predictor 4 -> sum d_{i}
durLastTrt     = tempProbsC0(:,5); %predictor 5 -> duration of last treatment
coeffOfVar     = tempProbsC0(:,7)./tempProbsC0(:,6);
com            = tempProbsC0(:,8);

TdfVec      = (0:15:360)';
probsMat    = tempProbsC0(:,indsSC0(1)-1:indsSC0(1)+length(TdfVec)-2)';
numTrts_2   = repmat(numTrts',length(TdfVec),1);
durLastTrt_2= repmat(durLastTrt',length(TdfVec),1);
inten_2     = repmat(inten',length(TdfVec),1);
totalDays_2 = repmat(totalDays',length(TdfVec),1);
days2trt_2  = repmat(days2trt',length(TdfVec),1);
coeffOfVar_2= repmat(coeffOfVar',length(TdfVec),1);
com_2       = repmat(com',length(TdfVec),1);
TdfVec_2    = repmat(TdfVec,size(numTrts,1),1);
predictorMat= [numTrts_2(:) totalDays_2(:) durLastTrt_2(:) days2trt_2(:) TdfVec_2(:) com_2(:) coeffOfVar_2(:) probsMat(:)];
save('predictorMat','predictorMat')

