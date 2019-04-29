% This scripts create COV files of CAPS and AlphaBeta from empty COV files
% created from BV on the first hand

% IMPORTANT NOTICE:
% .cov does not contain subject ID, so need to make sure the list of
% subject matches between excel and BV, because BV use a 'String' fashion
% to sort ID

clearvars

%% import the .txt file with all subjects' covariates in the right sequence as in BV
covFile = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\covariates_012419.txt';

covTable = readtable(covFile);

%% covariates for the 3 factors out of factor analysis
cov = xff('COV_67subj_empty_3.cov');

covVal = cov.Covariates;

covVal = table2array(covTable(:,35:37));
cov.Covariates(1:67,1:3) = covVal(1:67,1:3);

cov.SaveAs('COV_67subj_3factors.cov');
% cov = xff('COV_67subj_3factors.cov');

%% CAPS measurement
% including 5 factor and total
caps = xff('COV_67subj_empty_6.cov');

% nrOfSubj = caps.NrOfSubjectRows;
% nrOfCov = caps.NrOfCovariateColumns;
covCaps = caps.Covariates;
% runTimeVars = caps.RunTimeVars;


covCaps = table2array(covTable(:,15:20));
caps.Covariates(1:67,1:6) = covCaps(1:67,1:6);

caps.SaveAs('COV_67subj_caps.cov');
% caps = xff('COV_67subj_caps.cov');

%% model based and model free attitudes

ab = xff('COV_67subj_empty_8.cov');

covAb = ab.Covariates;
covAb = table2array(covTable(:,3:10));

ab.Covariates(1:67,1:8) = covAb(1:67,1:8);

ab.SaveAs('COV_67subj_att.cov');
ab = xff('COV_67subj_att.cov');


%% create covariate file in the format of .mat file
covTable2Include = covTable(:,[2:13,15:53]);
names = covTable2Include.Properties.VariableNames';
names(1) = {'subject_num'};
variates = table2array(covTable2Include);
save('COV_67subj_all.mat','names','variates');

% for day1, some subject were not included
day1exclude = [1232, 1278];
variates1 = variates(~ismember(variates(:,1),day1exclude),:);
save('COV_65subj_day1_all.mat','names','variates1');



