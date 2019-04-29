clearvars

root=pwd;
filename = [root,'\clinical_behavioral_091517.txt'];
tb = readtable(filename);

%% CAPS histogram
caps = tb.caps_total_pm(tb.isExcluded_behavior==0 & tb.isGain==1);
subjnum = length(caps(~isnan(caps)))

fig = figure
hist= histogram(caps,10)

%axis property
ax = gca;
ax.XLabel.String = {'CAPS Total'};
ax.XLabel.FontSize = 35;

% ax.XTick = [1,2,3,4];
ax.Box = 'off';
ax.FontSize = 25;
ax.LineWidth =3;
ax.YLabel.String = 'subjects count'; 
ax.YLabel.FontSize = 35;
% ax.YLim = [-1.2,0.6];

title('')


reexp = tb.R_fi_pm(tb.isExcluded_behavior==0 & tb.isGain==1);
histogram(reexp,10);


avoid = tb.A_fi_pm(tb.isExcluded_behavior==0 & tb.isGain==1);
histogram(avoid,10);
