function [binSeq,probVec,allSeqs] = trt2bin(trtInit,trtLen,totalTime,trtLastExp,samplePopsTemp)
indsSC0 = 2:5:117;
indsRC0 = indsSC0+1;
binSeq = [];probVec=[];allSeqs = [];
days           = [15 30 45 60 90 120 240 360];
daysMax        = 360;
daysInds       = days./15;
for k=1:size(trtInit,1)
    tempSeq  = -1*ones(1,trtLastExp(k)-1);
    %     allProbs = 100*samplePopsTemp(k,indsRC0)./(samplePopsTemp(k,indsSC0)+samplePopsTemp(k,indsRC0));
    for i=1:size(trtInit,2)
        [k i]
        %         probsMat       = allProbs(:,daysInds)';
        tempSeq(trtInit(k,i):trtInit(k,i)+trtLen(k,i)-1)=1;
        
        %         binSeq  = [binSeq;   repmat(binSequences,length(days),1) days'];
        %         probVec = [probVec;  probsMat(:)];
        %
        %         for l=1:(size(samplePopsTemp,2)/5)
        %             l
        %             binSeq(ix,:) = [zeros(1,totalTime-trtLastExp(k)) tempSeq 15*l];
        %             probVec(ix)  = allProbs(l);
        %             ix = ix +1;
        %         end
    end
    binSequences   = [zeros(1,totalTime-trtLastExp(k)+1) tempSeq];
    allSeqs        = [allSeqs;  binSequences];
end
end

