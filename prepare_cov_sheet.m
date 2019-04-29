%% prepare the .txt file to create .cov for imaging analysis
clearvars

%% read data
covFile = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\covariates_101118.txt';
compFileNoRemit = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\pca score gain_femaleIn_noRemitted_01172019.csv';
compFileRemit = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\pca score gain_femaleIn_Remitted_01172019.csv';
root='D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral';
clinicalfilename = fullfile(root, 'all clinical_allSubj.xlsx');

covTable = readtable(covFile);
compNoRemitTable = readtable(compFileNoRemit);
compRemitTable = readtable(compFileRemit);
clinicaltable = readtable(clinicalfilename);

%% PTSD and control

subj_mask_comp_noRemit = ismember(compNoRemitTable.id, covTable.id);
[sorted_id_noRemit, sort_idx_noRemit] = sort(compNoRemitTable.id, 'ascend');
comp1NoRemit = compNoRemitTable.comp1(sort_idx_noRemit & subj_mask_comp_noRemit);
comp2NoRemit = compNoRemitTable.comp2(sort_idx_noRemit & subj_mask_comp_noRemit);
comp3NoRemit = compNoRemitTable.comp3(sort_idx_noRemit & subj_mask_comp_noRemit);

% find which rows in the covTable needs to be fill in
subj_mask_noRemit = ismember(covTable.id, sorted_id_noRemit);
[sorted_id, sort_idx] = sort(covTable.id,'ascend');
subj_mask_noRemit_sorted = subj_mask_noRemit(sort_idx);
sort_idx_mask = sort_idx(subj_mask_noRemit_sorted);

covTable.comp1_femaleIn_noRemit(sort_idx_mask) = comp1NoRemit;
covTable.comp2_femaleIn_noRemit(sort_idx_mask) = comp2NoRemit;
covTable.comp3_femaleIn_noRemit(sort_idx_mask) = comp3NoRemit;


%% Remitted

subj_mask_comp_Remit = ismember(compRemitTable.id, covTable.id);
[sorted_id_Remit, sort_idx_Remit] = sort(compRemitTable.id, 'ascend');
comp1Remit = compRemitTable.comp1(sort_idx_Remit & subj_mask_comp_Remit);
comp2Remit = compRemitTable.comp2(sort_idx_Remit & subj_mask_comp_Remit);
comp3Remit = compRemitTable.comp3(sort_idx_Remit & subj_mask_comp_Remit);

% find which rows in the covTable needs to be fill in
subj_mask_Remit = ismember(covTable.id, sorted_id_Remit);
[sorted_id, sort_idx] = sort(covTable.id,'ascend');
subj_mask_Remit_sorted = subj_mask_Remit(sort_idx);
sort_idx_mask = sort_idx(subj_mask_Remit_sorted);

covTable.comp1_femaleIn_noRemit(sort_idx_mask) = comp1Remit;
covTable.comp2_femaleIn_noRemit(sort_idx_mask) = comp2Remit;
covTable.comp3_femaleIn_noRemit(sort_idx_mask) = comp3Remit;

%%  Clean and save covariate file
covTable.comp1_femaleIn_noRemit(covTable.comp1_femaleIn_noRemit == 0) = NaN;
covTable.comp2_femaleIn_noRemit(covTable.comp2_femaleIn_noRemit == 0) = NaN;
covTable.comp3_femaleIn_noRemit(covTable.comp3_femaleIn_noRemit == 0) = NaN;

writetable(covTable, "D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\covariates_012419.txt", 'Delimiter', '\t')

