clear all;close all;clc;
rng(3) % for  reproducibility
log10On = 0;
for pow = 3
    clearvars -except pow
    load('predictorMat')
    trainSubsetLen_REG_R = 5*10^pow; 
    
    folderName   = ['./R/PI_DATA_5E' num2str(pow)];
    mkdir(folderName);
    DATA             = predictorMat;
    DATA_POSITIVE    = predictorMat(predictorMat(:,end)>0,:);
    DATA_ZER0        = predictorMat(predictorMat(:,end)==0,:);
    
    if(size(DATA_POSITIVE,1)<trainSubsetLen_REG_R)
        trainInds_REG_POS        = randsample(1:size(DATA_POSITIVE,1),trainSubsetLen_REG_R,true);
    else
        trainInds_REG_POS        = randsample(1:size(DATA_POSITIVE,1),trainSubsetLen_REG_R);
    end
    if(size(DATA_ZER0,1)<trainSubsetLen_REG_R)
        trainInds_REG_ZER        = randsample(1:size(DATA_ZER0,1),trainSubsetLen_REG_R,true);
    else
        trainInds_REG_ZER        = randsample(1:size(DATA_ZER0,1),trainSubsetLen_REG_R);
    end
    DATA_REG_TRAIN_R_BLNCD   = [DATA_ZER0(trainInds_REG_ZER,:);DATA_POSITIVE(trainInds_REG_POS,:)];
    if(size(DATA_POSITIVE,1)<trainSubsetLen_REG_R)
        DATA_REG_TRAIN_R_POSONLY = DATA_POSITIVE(randsample(1:size(DATA_POSITIVE,1),trainSubsetLen_REG_R,true),:);
    else
        DATA_REG_TRAIN_R_POSONLY = DATA_POSITIVE(randsample(1:size(DATA_POSITIVE,1),trainSubsetLen_REG_R),:);
    end
    if(size(DATA,1)<trainSubsetLen_REG_R)
        DATA_REG_TRAIN_R         = DATA(randsample(1:size(DATA,1),trainSubsetLen_REG_R,true),:);
    else
        DATA_REG_TRAIN_R         = DATA(randsample(1:size(DATA,1),trainSubsetLen_REG_R),:);
    end
    [transdat,lambda]        = boxcox(DATA_REG_TRAIN_R_POSONLY(:,end));
    DATA_REG_TRAIN_R_BOXCOX  = [DATA_REG_TRAIN_R_POSONLY(:,1:end-1) transdat];
    DATA_REG_TRAIN_R_CLS     = DATA_REG_TRAIN_R;
    DATA_REG_TRAIN_R_CLS(:,end) = DATA_REG_TRAIN_R(:,end)>0;
    
    
    DATA_REG_TRAIN_R_CLS        = DATA_REG_TRAIN_R;
    DATA_REG_TRAIN_R_CLS(:,end) = DATA_REG_TRAIN_R(:,end)>0;
    
    DATA_CLS_R_min         = DATA_REG_TRAIN_R_CLS(DATA_REG_TRAIN_R_CLS(:,end)==0,:); % FIND MINORTY CLASS
    DATA_CLS_R_maj         = DATA_REG_TRAIN_R_CLS(DATA_REG_TRAIN_R_CLS(:,end)>0,:);% FIND MAJORITY CLASS
    if(size(DATA_CLS_R_maj,1)<size(DATA_CLS_R_min,1))
        subsampBalance         = randsample(1:size(DATA_CLS_R_maj,1),size(DATA_CLS_R_min,1),true); % UNDERSAMPLE MAJORITY CLASS
    else
        subsampBalance         = randsample(1:size(DATA_CLS_R_maj,1),size(DATA_CLS_R_min,1)); % UNDERSAMPLE MAJORITY CLASS
    end
    
    DATA_CLS_R_maj_sub     = DATA_CLS_R_maj(subsampBalance,:); % UNDERSAMPLED MAJORITY CLASS
    DATA_CLS_TRAIN_R_BLNCD = [DATA_CLS_R_min ;DATA_CLS_R_maj_sub]; % UNDERSAMPLED MAJORITY CLASS + MINORITY CLASS -> BALANCED CLASS TRAINING DATASET
    
    
    DATA_REG_TRAIN_R_POSONLY   = DATA_POSITIVE(trainInds_REG_POS,:); %# SAMPLES trainInds_REG_POS
    DATA_REG_TRAIN_R_ZERONLY   = DATA_ZER0(trainInds_REG_ZER,:); %# SAMPLES trainInds_REG_POS
    DATA_REG_TRAIN_R_BLNCD     = [DATA_REG_TRAIN_R_POSONLY(1:end/2,:);DATA_REG_TRAIN_R_ZERONLY(1:end/2,:)];
    DATA_REG_TRAIN_R_CLS_BLNCD = DATA_REG_TRAIN_R_BLNCD;
    
    
    DATA_REG_TRAIN_R_CLS_BLNCD(:,end) = DATA_REG_TRAIN_R_CLS_BLNCD(:,end)>0;
    
    
    DATA_REG_TRAIN_R_POSONLY       = [rand(size(DATA_REG_TRAIN_R_POSONLY,1),1) DATA_REG_TRAIN_R_POSONLY];
    DATA_REG_TRAIN_R_CLS_BLNCD     = [rand(size(DATA_REG_TRAIN_R_CLS_BLNCD,1),1) DATA_REG_TRAIN_R_CLS_BLNCD];
    
    COLUMNS={'rnd','N','sumdi','dN','t1','Tdf','COM','cv','prevelance'};
    
    DATA_REG_TRAIN_TABLE_POSONLY = array2table(DATA_REG_TRAIN_R_POSONLY,'VariableNames',COLUMNS);
    writetable(DATA_REG_TRAIN_TABLE_POSONLY,[folderName '/DATA_REG_TRAIN_POSONLY.csv'])
    
    DATA_REG_TRAIN_TABLE_CLS_BLNCD = array2table(DATA_REG_TRAIN_R_CLS_BLNCD,'VariableNames',COLUMNS);
    writetable(DATA_REG_TRAIN_TABLE_CLS_BLNCD,[ folderName '/DATA_REG_TRAIN_CLS_BLNCD.csv'])
    
end

