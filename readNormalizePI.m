function [ variable_out] = readNormalizePI( varname, directory, len )
longname  = [directory varname '.txt'];
if isfile(longname)
    [variable_out,~]=importdata(longname);
    variable_out = variable_out.data;
    if(max(variable_out)>0)
         variable_out = myNormalize(variable_out);
        variable_out = variable_out./max(variable_out);
    else
        variable_out = variable_out;
    end
else
    variable_out=NaN(1,8);
end
variable_out = num2cell(variable_out);
end
