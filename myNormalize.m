function [norm_data] = myNormalize(bla)
inds0=find(bla==0);
norm_data = (bla - min(bla)) / ( max(bla) - min(bla) );
norm_data(inds0)=0;
end

