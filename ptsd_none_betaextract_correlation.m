clearvars

%%
root=pwd;

file = 'BetaExtracts_none_femaleIn_noRemit_Day2_GLM_ROI_none_femaleIn_noRimt_Day2_AllCond_capsotalCovar_p0.005Alphasim_vmPFC.xlsx';
filename = fullfile(root,file);
tbold = readtable(filename);

%exclude females
% exclude = [120 1216 1325 1338 50];
% exclude remitted
exclude = [3, 8, 108, 111, 122, 1005, 1203, 1223, 1246, 1254, 1278, 1290, 1328, 1354];

tb = tbold(~ismember(tbold.Subject,exclude),:);
tb.id = tb.Subject;

covfile =fullfile('D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral', filesep, 'covariates_012419.txt');
covold = readtable(covfile);
%exclude subject 1300
cov = covold(~ismember(covold.id,exclude), :);

% combine table 
tbbeta = join(tb,cov,'Keys','id');

% % tb and cov should have a same subject list in the same order
tbbeta.All = (tbbeta.Amb_gains_Display + tbbeta.Risk_gains_Display + tbbeta.Amb_loss_Display + tbbeta.Risk_loss_Display) ./ 4;
tbbeta.Gain = (tbbeta.Amb_gains_Display + tbbeta.Risk_gains_Display) ./ 2;
tbbeta.Loss = (tbbeta.Amb_loss_Display + tbbeta.Risk_loss_Display) ./ 2;
tbbeta.Risk = (tbbeta.Risk_loss_Display + tbbeta.Risk_gains_Display) ./ 2;
tbbeta.Ambig = (tbbeta.Amb_loss_Display + tbbeta.Amb_gains_Display) ./ 2;

%% Exclude remitted PTSD
tbbeta = tbbeta(~strcmp(tbbeta.group,'R'),:);

%% correlation 
betaNames = {'Amb_gains_Display', 'Risk_gains_Display', 'Amb_loss_Display', 'Risk_loss_Display','All', 'Gain', 'Loss'};
clusterNames = {'capsR', 'capsA', 'capsN', 'capsDA', 'capsAA', 'caps'};
attNames = {'alpha_gain', 'beta_gain', 'alpha_loss', 'beta_loss', 'risk_gain','ambig_gain','risk_loss','ambig_loss'};
fitTypes = {'ordinary', 'robust'};

% individual plot and fitting
betaName = 'Amb_gains_Display';
covName = 'comp1_noRemit';
fitType = 'ordinary';

%% Scatter plot and regression fitting
% plot regression with symptom or behavioral attitude
regrName = covName;

y = tbbeta.(betaName);
x = tbbeta.(regrName);
% x = log(x);

figure    
scatter(x,y, 'filled', 'LineWidth',5);
hold on

ax=gca;
ax.FontSize = 20;
ax.LineWidth = 3;
% ax.YLim = [-0.3,0.3];
% ax.YTick = [ -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3];
% ax.XLim = [0,120]; 
                
if strcmp(fitType,'robust')
    % robust regression
    [b,robustmdl1]= robustfit(x,y); 
    linex = linspace(min(x),max(x));
    liney = b(2)*linex+b(1);
    plot(linex, liney, 'color','k', 'LineWidth', 2);
    [corrmat, pmat] = corrcoef([x,y],'rows','complete');
    
    
    % print text of coeff and p value
    txt1 = ['regression coeff = ', num2str(b(2))];
    txt2 = ['p = ', num2str(round(robustmdl1.p(2),4,'significant'))];

    txt3 = ['correlation coeff =',num2str(corrmat(1,2))];
    txt4 = ['p = ', num2str(round(pmat(1,2),4,'significant'))];
    xlab = xlim;
    ylab = ylim;
    txt = {txt1;txt2;[];txt3;txt4};
    text(xlab(2)-(xlab(2)-xlab(1))/4, ylab(2)-(ylab(2)-ylab(1))/8, txt, 'FontSize',8)

    title([regrName ' with ' betaName ' Robust'])

% elseif strcmp(fitType, 'pearson')
%     [corrmat, pmat] = corrcoef([x,y],'rows','complete')
    
elseif strcmp(fitType, 'ordinary')
    % ordinary linear regression
    mdl1 = LinearModel.fit(x,y); % creates a linear model of the responses y to a tb matrix x
    coeff = table2array(mdl1.Coefficients);
%     linex = linspace(min(x)-(max(x)-min(x))/30,max(x)+(max(x)-min(x))/20);
    linex = linspace(min(x),max(x));
    liney = coeff(2,1)*linex+coeff(1,1);
    plot(linex, liney, 'color','k', 'LineWidth', 2);
    
    % print text of coeff, r2 and p value
    txt1 = ['R^{2} = ',num2str(mdl1.Rsquared.Ordinary)];
    txt2 = ['p = ', num2str(round(coeff(2,4),4,'significant'))];
    txt3 = ['coeff = ', num2str(coeff(2,1))];
    xlab = xlim;
    ylab = ylim;
    txt = {txt3;txt1;txt2};
    text(xlab(2)-(xlab(2)-xlab(1))/1.3, ylab(2), txt, 'FontSize',8)

    title([regrName ' with ' betaName ' OLS'])
    
    % residual plot
    figure
    scatter(x,mdl1.Residuals.Raw)
    hold on
    title([regrName ' with ' betaName ' OLS ' 'Residuals'])

end

%% Multilinear model
% cluster2fit = {'R_fi_pm', 'A_fi_pm', 'N_fi_pm', 'DA_fi_pm', 'AA_fi_pm', 'caps_total_pm'};
% cluster2fit = {'caps_total_pm'};
% cluster2fit = {'N_fi_pm', 'caps_total_pm'};
% cluster2fit = {'N_fi_pm'};
% cluster2fit = {'DA_fi_pm', 'caps_total_pm'};
% 
% xMulti = zeros(length(cov.caps_total_pm), length(cluster2fit)+1);
% xMulti(:,1) = ones(length(cov.caps_total_pm),1);
% for i = 1:length(cluster2fit)
%     xMulti(:,i+1) = cov.(cluster2fit{i});    
% end
% 
% [b,bint,r,rint,stats] = regress(y,xMulti);


lm = fitlm(tbbeta, 'Loss~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'All~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'Gain~capsR+capsA+capsN+capsDA+capsAA')


lm = fitlm(tbbeta, 'Risk~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'Ambig~capsR+capsA+capsN+capsDA+capsAA')

lm = fitlm(tbbeta, 'Amb_gains_Display~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'Risk_gains_Display~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'Amb_loss_Display~capsR+capsA+capsN+capsDA+capsAA')
lm = fitlm(tbbeta, 'Risk_loss_Display~capsR+capsA+capsN+capsDA+capsAA')


lm = fitlm(tbbeta, 'All ~ bdi + stai_x1 + stai_x2')

%%
lm = fitlm(tbbeta, 'Gain~R_fi_pm+A_fi_pm+N_fi_pm+DA_fi_pm+AA_fi_pm')
% plot mean and se for the predictors
plotmean = lm.Coefficients.Estimate(2:6);
plotsem = lm.Coefficients.SE(2:6);
fig = figure
bplot = bar(plotmean);
hold on
errorbar(plotmean, plotsem, '.','Color',[0,0,0],'LineWidth',2);

ax = gca;
ax.XTickLabel = '';
ax.Box = 'off';
ax.FontSize = 25;
ax.LineWidth =3;

%%
lm = fitlm(tbbeta, 'Loss~R_fi_pm+A_fi_pm+N_fi_pm+DA_fi_pm+AA_fi_pm')
plotmean = lm.Coefficients.Estimate(2:6);
plotsem = lm.Coefficients.SE(2:6);
fig = figure
bplot = bar(plotmean);
hold on
errorbar(plotmean, plotsem, '.','Color',[0,0,0],'LineWidth',2);

ax = gca;
ax.XTickLabel = '';
ax.Box = 'off';
ax.FontSize = 25;
ax.LineWidth =3;

lm = fitlm(tbbeta, 'All~caps_total_pm')
lm = fitlm(tbbeta, 'All~N_fi_pm')
lm = fitlm(tbbeta, 'All~caps_total_pm+N_fi_pm')

lm = fitlm(tbbeta, 'Amb_gains_Display~caps_total_pm')
lm = fitlm(tbbeta, 'Amb_gains_Display~N_fi_pm')
lm = fitlm(tbbeta, 'Amb_gains_Display~caps_total_pm+N_fi_pm')

lm = fitlm(tbbeta, 'Risk_gains_Display~caps_total_pm')
lm = fitlm(tbbeta, 'Risk_gains_Display~N_fi_pm')
lm = fitlm(tbbeta, 'Risk_gains_Display~caps_total_pm+N_fi_pm')

lm = fitlm(tbbeta, 'Amb_loss_Display~caps_total_pm')
lm = fitlm(tbbeta, 'Amb_loss_Display~N_fi_pm')
lm = fitlm(tbbeta, 'Amb_loss_Display~caps_total_pm+N_fi_pm')

lm = fitlm(tbbeta, 'Risk_loss_Display~caps_total_pm')
lm = fitlm(tbbeta, 'Risk_loss_Display~N_fi_pm')
lm = fitlm(tbbeta, 'Risk_loss_Display~caps_total_pm+N_fi_pm')


%% symptom, attitudes, and activation
lm = fitlm(tbbeta, 'All ~ caps + alpha_gain')
lm = fitlm(tbbeta, 'All ~ caps + beta_gain')
lm = fitlm(tbbeta, 'All ~ caps + alpha_loss')
lm = fitlm(tbbeta, 'All ~ caps + beta_loss')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ caps + alpha_gain')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ caps + beta_gain')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ caps + alpha_loss')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ caps + beta_loss')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ comp1_femaleIn_noRemit + alpha_gain')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ comp1_femaleIn_noRemit + beta_gain')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ comp1_femaleIn_noRemit + alpha_loss')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ comp1_femaleIn_noRemit + beta_loss')

    
    
    
    
    
    
    
    
    
    
    
    
