function [com,offset] = centerOfMass(binSeqTemp)
inds=find(binSeqTemp~=0);
tempSequence = binSeqTemp(inds);
tempSequence = double((tempSequence+1)>0);%map -1s to 0, 1s stay as 1s
stats        = regionprops(tempSequence);
com          = stats.Centroid(1);
offset       = inds(1);
end

