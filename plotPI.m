function [] = plotPI_noHM(tableCond_CLS,tableCond_REG,figname,notdf)
%%
if(notdf==1)
    predictorNames = {'rnd','$N$','$\sum_{i}d_{i}$','$d_{N}$','$(t_1-T_{I})$','$COM$','$c_{v}$'};
else 
    predictorNames = {'rnd','$N$','$\sum_{i}d_{i}$','$d_{N}$','$(t_1-T_{I})$','$T_{df}$','$COM$','$c_{v}$'};
end
dr = [0.8431,0,0.1490];
bb = [0.3020,0.3922,0.5529]; %blueberry
impsCond = [tableCond_CLS;tableCond_REG]';

% sz   = [800 350];
% screensize = get(0,'ScreenSize');
% set(0,'defaultAxesFontSize',22)
% xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the
% ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the
% h=figure('Position', [xpos , ypos, sz(1), sz(2)]);
% [ha, pos]=tight_subplot(1,1,[.01 .01],[.01 .01],[.01 .01]);
% set(ha(1),'position',[0.09 0.18 0.9 0.8])
% axes(ha(1))

figure;
axes1 = axes('Position',[0.12 0.198 0.85 0.65]);

h=bar(impsCond);
h(1).FaceColor = bb;
h(2).FaceColor = dr;
grid on;
xticklabels(predictorNames)
% leg={'CLS','REG'};
leg={'Classification','Regression'};
ylabel('Relative Importance','interpreter','latex')
xlabel('Predictors','interpreter','latex')
legend(leg,'interpreter','latex','location','northwest','Fontsize',22);
set(gca,'TickLabelInterpreter','latex')
axis([0.5 length(predictorNames)+0.5 -0.1 1.05])

print(figname,'-dpng')
print(figname,'-depsc')

end

