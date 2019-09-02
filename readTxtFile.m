function [ variable_out] = readTxtFile( varname, directory )
longname  = [directory varname '.txt'];
if isfile(longname)
    [variable_out,~]=importdata(longname);
else
    variable_out=[];
end
end

