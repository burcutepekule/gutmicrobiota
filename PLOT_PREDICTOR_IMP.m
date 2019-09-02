clear all;close all;clc;

setTrees  = 50:50:200;
pow       = 3;
treeSizePick =200; %pick the tree size for plotting
predNames = {'num_trees','rnd','N','sum_di','d_N','t_1', 'T_df','COM','c_v'};
varTypes  = repmat({'double'},1,length(predNames));

directory = './R/PI_RESULTS/';
figname   = ['./R/PI_RESULTS/PI_T' num2str(treeSizePick) '_S_5E' num2str(pow)]; %filename for the figure

tblNormal_REG = table('Size',[length(setTrees) length(predNames)],'VariableTypes',varTypes,'VariableNames',predNames);
tblNormal_CLS = table('Size',[length(setTrees) length(predNames)],'VariableTypes',varTypes,'VariableNames',predNames);
tblCond_REG   = table('Size',[length(setTrees) length(predNames)],'VariableTypes',varTypes,'VariableNames',predNames);
tblCond_CLS   = table('Size',[length(setTrees) length(predNames)],'VariableTypes',varTypes,'VariableNames',predNames);

for s = 1:length(setTrees)
    tblNormal_REG(s,1)     = num2cell(setTrees(s));
    tblNormal_CLS(s,1)     = num2cell(setTrees(s));
    tblCond_REG(s,1)       = num2cell(setTrees(s));
    tblCond_CLS(s,1)       = num2cell(setTrees(s));
    tblNormal_REG(s,2:end) = readNormalizePI(['tableNormal_T' num2str(setTrees(s)) '_REG_S_5E' num2str(pow)], directory);
    tblNormal_CLS(s,2:end) = readNormalizePI(['tableNormal_T' num2str(setTrees(s)) '_CLS_S_5E' num2str(pow)], directory);
    tblCond_REG(s,2:end)   = readNormalizePI(['tableCond_T' num2str(setTrees(s)) '_REG_S_5E' num2str(pow)], directory);
    tblCond_CLS(s,2:end)   = readNormalizePI(['tableCond_T' num2str(setTrees(s)) '_CLS_S_5E' num2str(pow)], directory);
end
tableCond_REG = table2array(tblCond_REG(tblCond_REG.num_trees==treeSizePick,2:end));
tableCond_CLS = table2array(tblCond_CLS(tblCond_CLS.num_trees==treeSizePick,2:end));
plotPI(tableCond_CLS,tableCond_REG,figname,0)
