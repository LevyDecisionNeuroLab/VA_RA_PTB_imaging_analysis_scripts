clearvars
root=pwd;
%% Read data from none GLM beta extract

file_none = 'BetaExtracts_none_femaleIn_noRemit_Day2_GLM_ROI_none_femaleIn_noRimt_Day2_AllCond_capsotalCovar_p0.005Alphasim_vmPFC.xlsx';
filename_none = fullfile(root,file_none);
tbold_none = readtable(filename_none);

%exclude females
% exclude = [120 1216 1325 1338 50];
% exclude remitted
exclude = [3, 8, 108, 111, 122, 1005, 1203, 1223, 1246, 1254, 1278, 1290, 1328, 1354];

tb_none = tbold_none(~ismember(tbold_none.Subject,exclude),:);
tb_none.id = tb_none.Subject;

covfile =fullfile('D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral', filesep, 'covariates_012419.txt');
covold = readtable(covfile);
%exclude subject 1300
cov = covold(~ismember(covold.id,exclude), :);

% combine table 
tbbeta_none = join(tb_none,cov,'Keys','id');

% % tb and cov should have a same subject list in the same order
tbbeta_none.All = (tbbeta_none.Amb_gains_Display + tbbeta_none.Risk_gains_Display + tbbeta_none.Amb_loss_Display + tbbeta_none.Risk_loss_Display) ./ 4;
tbbeta_none.Gain = (tbbeta_none.Amb_gains_Display + tbbeta_none.Risk_gains_Display) ./ 2;
tbbeta_none.Loss = (tbbeta_none.Amb_loss_Display + tbbeta_none.Risk_loss_Display) ./ 2;
tbbeta_none.Risk = (tbbeta_none.Risk_loss_Display + tbbeta_none.Risk_gains_Display) ./ 2;
tbbeta_none.Ambig = (tbbeta_none.Amb_loss_Display + tbbeta_none.Amb_gains_Display) ./ 2;

%% Read data from sv GLM beta extract

file_sv = 'BetaExtracts_SV_femaleIn_noRemit_GLM_ROI_SV_femaleIn_noRimt_RiskLoss_Contrast_p0.001Alphasim_lLentiform564.xlsx';
filename_sv = fullfile(root,file_sv);
tbold_sv = readtable(filename_sv);

tb_sv = tbold_sv(~ismember(tbold_sv.Subject,exclude),:);
tb_sv.id = tb_sv.Subject;


tb_sv.Amb_gainsMinusloss = tb_sv.Amb_gains_DisplayXP1 - tb_sv.Amb_loss_DisplayXP1;
tb_sv.Risk_gainsMinusloss = tb_sv.Risk_gains_DisplayXP1 - tb_sv.Risk_loss_DisplayXP1;

tb_sv_beta = tb_sv(:,[3,5,7,9,12,13:14]); % columns of SV betas and id

% change column name in tb_sv
tb_sv_beta.Properties.VariableNames = {'Amb_gains_sv', 'Risk_gains_sv', 'Amb_loss_sv', 'Risk_loss_sv', 'id', 'Amb_gainsMinusloss_sv', 'Risk_gainsMinusloss_sv'};

%% final table for analysis
tbbeta = join(tbbeta_none,tb_sv_beta,'Keys','id');

%% symptom, attitudes, and activation
lm = fitlm(tbbeta, 'All ~ caps')
lm = fitlm(tbbeta, 'Risk_gains_Display ~ caps')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ caps')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ caps')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ caps')

lm = fitlm(tbbeta, 'All ~ caps + alpha_gain')
lm = fitlm(tbbeta, 'All ~ caps + beta_gain')
lm = fitlm(tbbeta, 'All ~ caps + alpha_loss')
lm = fitlm(tbbeta, 'All ~ caps + beta_loss')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ caps + alpha_gain')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ caps + beta_gain')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ caps + alpha_loss')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ caps + beta_loss')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ caps')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ caps')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ caps')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ caps')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ caps + Risk_gains_sv')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ caps + Amb_gains_sv')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ caps + Risk_loss_sv')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ caps + Amb_loss_sv')

lm = fitlm(tbbeta, 'Amb_loss_sv ~ caps')
lm = fitlm(tbbeta, 'Risk_loss_sv ~ caps')


lm = fitlm(tbbeta, 'Risk_gains_Display ~ comp1_femaleIn_noRemit + alpha_gain')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ comp1_femaleIn_noRemit + beta_gain')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ comp1_femaleIn_noRemit + alpha_loss')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ comp1_femaleIn_noRemit + beta_loss')

lm = fitlm(tbbeta, 'Risk_gains_Display ~ comp1_femaleIn_noRemit + Risk_gains_sv')
lm = fitlm(tbbeta, 'Amb_gains_Display ~ comp1_femaleIn_noRemit + Amb_gains_sv')
lm = fitlm(tbbeta, 'Risk_loss_Display ~ comp1_femaleIn_noRemit + Risk_loss_sv')
lm = fitlm(tbbeta, 'Amb_loss_Display ~ comp1_femaleIn_noRemit + Amb_loss_sv')


fig = figure
scatter3(tbbeta.Amb_loss_sv, tbbeta.caps, tbbeta.Amb_loss_Display, 'filled')
ax = gca;
ax.XLabel.String = 'SV neural representation';
ax.YLabel.String = 'Clinical measurement';
ax.ZLabel.String = 'General activation';
