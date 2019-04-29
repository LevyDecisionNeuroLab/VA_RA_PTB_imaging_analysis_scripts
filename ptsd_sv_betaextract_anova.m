clearvars

root=pwd;
file = 'BetaExtracts_CV_withExclustionGLM_ROI_Bartra_vStr.xlsx';
filename = fullfile(root,file);
tbold = readtable(filename);
%exclude females
exclude = [120 1216 1325 1338 50];
tb = tbold(~ismember(tbold.Subject,exclude),:);
tb.id = tb.Subject;

covfile =fullfile('D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral', filesep, 'covariates_101118.txt');
covold = readtable(covfile);
%exclude subject
cov = covold(~ismember(covold.id,exclude), :);

% combine table 
tbbeta = join(tb,cov,'Keys','id');


% creat data table that is not cross-tabed, for anova
analytbNames = {'id', 'group', 'isGain', 'isRisk', 'svbeta'};
id = repmat(tbbeta.Subject,4,1);
group = repmat(tbbeta.group,4,1);
isGain = [ones(height(tbbeta)*2,1);zeros(height(tbbeta)*2,1)];
isRisk = [zeros(height(tbbeta),1); ones(height(tbbeta),1); zeros(height(tbbeta),1); ones(height(tbbeta),1)];
svBeta = [tbbeta.Amb_gains_DisplayXP1; tbbeta.Risk_gains_DisplayXP1; tbbeta.Amb_loss_DisplayXP1; tbbeta.Risk_loss_DisplayXP1];

analytb = table(id, group, isGain, isRisk, svBeta);

%% Avona
% gain/loss, risk/ambig, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R')),...
    {analytb.isGain(~strcmp(analytb.group, 'R'));...
    analytb.isRisk(~strcmp(analytb.group, 'R'));...
    analytb.group(~strcmp(analytb.group, 'R'))},'full');

% gain, risk/ambig, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==1),...
    {analytb.isRisk(~strcmp(analytb.group, 'R')& analytb.isGain ==1);...
    analytb.group(~strcmp(analytb.group, 'R')& analytb.isGain ==1)},'full');

% gain, ambig, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==1 & analytb.isRisk == 0),...
    {analytb.group(~strcmp(analytb.group, 'R')& analytb.isGain ==1 & analytb.isRisk == 0)},'full');

% gain, ambig, CC/PTSD/Remitted
[p,tbl,stats,terms]=anovan(analytb.svBeta(analytb.isGain ==1 & analytb.isRisk == 0),...
    {analytb.group(analytb.isGain ==1 & analytb.isRisk == 0)},'full');

% gain, risk, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==1 & analytb.isRisk == 1),...
    {analytb.group(~strcmp(analytb.group, 'R')& analytb.isGain ==1 & analytb.isRisk == 1)},'full');

% gain, risk, CC/PTSD/Remitted
[p,tbl,stats,terms]=anovan(analytb.svBeta(analytb.isGain ==1 & analytb.isRisk == 1),...
    {analytb.group(analytb.isGain ==1 & analytb.isRisk == 1)},'full');

% loss, risk/ambig, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==0),...
    {analytb.isRisk(~strcmp(analytb.group, 'R')& analytb.isGain ==0);...
    analytb.group(~strcmp(analytb.group, 'R')& analytb.isGain ==0)},'full');

% loss, ambig, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==0 & analytb.isRisk == 0),...
    {analytb.group(~strcmp(analytb.group, 'R') & analytb.isGain ==0 & analytb.isRisk == 0)},'full');

% loss, risk, CC/PTSD
[p,tbl,stats,terms]=anovan(analytb.svBeta(~strcmp(analytb.group, 'R') & analytb.isGain ==0 & analytb.isRisk == 1),...
    {analytb.group(~strcmp(analytb.group, 'R')& analytb.isGain ==0 & analytb.isRisk == 1)},'full');



%% Bargraph
% each row represents: AG ,RG,AL,RL(CC),  AG ,RG,AL,RL(PTSD),
% each column represents: CC, PTSD
plotmean = [nanmean(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
    nanmean(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));,...
    nanmean(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    nanmean(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')));...
    nanmean(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
    nanmean(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));...    
    nanmean(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    nanmean(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')))...
    ];
plotstd = [nanstd(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
    nanstd(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));...  
    nanstd(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    nanstd(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')));...   
    nanstd(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
    nanstd(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));,...    
    nanstd(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    nanstd(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')))...
    ];
plotcount = [length(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
    length(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));...    
    length(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    length(analytb.svBeta(analytb.isGain == 1 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')));...    
    length(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'CC'))),...
     length(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==0 & strcmp(analytb.group, 'PTSD')));...   
    length(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'CC'))),...
    length(analytb.svBeta(analytb.isGain == 0 & analytb.isRisk ==1 & strcmp(analytb.group, 'PTSD')))...
    ];
plotsem = plotstd ./ sqrt(plotcount);

% figure, beta bar plots for four conditions

fig = figure
set(fig, 'Position', [90 200 1120 700])
bplot = bar(plotmean);
hold on
errorbar([1,2,3,4]-0.145, plotmean(:,1), plotsem(:,2), '.','Color',[0,0,0],'LineWidth',2);
hold on
errorbar([1,2,3,4]+0.145,plotmean(:,2),plotsem(:,2),'.','Color',[0,0,0],'LineWidth',2);


%bar color
bplot(1).FaceColor = [104,160,66]/255;
bplot(1).EdgeColor = [104,160,66]/255;
bplot(2).FaceColor = [237,125,49]/255;
bplot(2).EdgeColor = [237,125,49]/255;
bplot(1).BarWidth = 0.9;

%axis property
ax = gca;
% ax.XTickLabel = {'AG','RG','AL','RL'};
ax.XTickLabel = '';
ax.Box = 'off';
ax.FontSize = 35;
ax.LineWidth =4;
% ax.YLabel.String = 'beta';
% ax.YLabel.FontSize = 18;
ax.YLim = [-0.05, 0.08];

title(file, 'FontSize',16)

leg = legend('CC','PTSD');
leg.FontSize = 25;
