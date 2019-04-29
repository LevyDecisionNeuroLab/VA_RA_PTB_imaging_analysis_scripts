%% prepare the .txt file to create .cov for imaging analysis
clearvars

%% read data
day1_file = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Analysis Ruonan\Fitpar files\Behavior data fitpar_020519\day1_par_nonpar.txt';
day2_file = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Analysis Ruonan\Fitpar files\Behavior data fitpar_020519\day2_par_nonpar.txt';
compFileNoRemit = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\pca score gain_femaleIn_noRemitted_01172019.csv';
compFileRemit = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\pca score gain_femaleIn_Remitted_01172019.csv';
clinical_file = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Clinical and behavioral\question scores EFA_09152018.xlsx';


day1 = readtable(day1_file);
day2 = readtable(day2_file);
compNoRemitTable = readtable(compFileNoRemit);
compRemitTable = readtable(compFileRemit);
clinical = readtable(clinical_file);

%% combine table
% clean clinical table, only need the clinical measurement
clinical_imaging = clinical(clinical.isGain == 1 & clinical.isExcluded_imaging == 0, [1, 2, 23:24, 33:53]);
clinical_name = clinical_imaging.Properties.VariableNames';

% cleance pca table, include only imaging subjects
% pcaNoRemit_imaging = compNoRemitTable(compNoRemitTable.isExcluded_imaging == 0, {'id',...
%     'bdiii_total','R_F_I_PastMonth_','A_F_I_PastMonth_','N_F_I_PastMonth_','DA_F_I_PastMonth_','AA_F_I_PastMonth_','total_pm',...
%     'ctq_physicalAbuse','ctq_emotionalAbuse','ctq_sexualAbuse','ctq_emotionalNeglect','ctq_physicalNeglect','ctq_total',...
%     'stai_x1_total','stai_x2_total','ces_total','des_taxon_sum','des_depersonalization_sum','des_amnsetic_sum','des_absorption_sum','des_total_sum',...
%     'comp1','comp2','comp3'});
% 
% pcaRemit_imaging = compRemitTable(compRemitTable.isExcluded_imaging == 0, {'id',...
%     'bdiii_total','R_F_I_PastMonth_','A_F_I_PastMonth_','N_F_I_PastMonth_','DA_F_I_PastMonth_','AA_F_I_PastMonth_','total_pm',...
%     'ctq_physicalAbuse','ctq_emotionalAbuse','ctq_sexualAbuse','ctq_emotionalNeglect','ctq_physicalNeglect','ctq_total',...
%     'stai_x1_total','stai_x2_total','ces_total','des_taxon_sum','des_depersonalization_sum','des_amnsetic_sum','des_absorption_sum','des_total_sum',...
%     'comp1','comp2','comp3'});

pcaNoRemit_imaging = compNoRemitTable(compNoRemitTable.isExcluded_imaging == 0, {'id',...
    'comp1','comp2','comp3'});

pcaRemit_imaging = compRemitTable(compRemitTable.isExcluded_imaging == 0, {'id',...
    'comp1','comp2','comp3'});

% stack remitted and non-remitted data
pca_imaging = vertcat(pcaNoRemit_imaging, pcaRemit_imaging);
pca_names = pca_imaging.Properties.VariableNames';

% include only imaging subjects
day1_imaging = day1(day1.isExcluded_imaging == 0, :);
day2_imaging = day2(day2.isExcluded_imaging == 0, :);

% subjects who are included in imaginag analysis but did not have PCA component score because of incomplete clinical data
subj_noPCA = day1_imaging.id(~ismember(day1_imaging.id, pca_imaging.id));

% add these subjects into PCA table, with empty entries
pca_imaging.id(length(pca_imaging.id)+1:length(pca_imaging.id)+3) = subj_noPCA;
pca_imaging(length(pca_imaging.id)-2:length(pca_imaging.id), 2:width(pca_imaging)) = {NaN};

% join pca and clinical table
pca_clinical_imaging = join(clinical_imaging, pca_imaging, 'Key', 'id');

% join table
day1_all = join(day1_imaging, pca_clinical_imaging, 'Key', 'id');
day2_all = join(day2_imaging, pca_clinical_imaging, 'Key', 'id');

day1_all.Properties.VariableNames'

% select variable and subjects for covariate file
day1_cov = day1_all(:, {'id','bdiii_total',...
    'R_F_I_PastMonth_','A_F_I_PastMonth_','N_F_I_PastMonth_','DA_F_I_PastMonth_','AA_F_I_PastMonth_','total_pm',...
    'ctq_physicalAbuse','ctq_emotionalAbuse','ctq_sexualAbuse','ctq_emotionalNeglect','ctq_physicalNeglect','ctq_total',...
    'stai_x1_total','stai_x2_total','ces_total',...
    'des_taxon_sum','des_depersonalization_sum','des_amnsetic_sum','des_absorption_sum','des_total_sum',...
    'comp1','comp2','comp3',...
    'G_risk25','G_risk50','G_risk75','G_amb24','G_amb50','G_amb74',...
    'L_risk25','L_risk50','L_risk75','L_amb24','L_amb50','L_amb74',...
    'isExcluded_behavior','isExcluded_imaging',...
    'age','isMale','kbit',...
    'G_risk','G_amb','G_amb_risk50','L_risk','L_amb','L_amb_risk50',...
    'G_alpha','G_beta','L_alpha','L_beta'});

day2_cov = day2_all(:, {'id','bdiii_total',...
    'R_F_I_PastMonth_','A_F_I_PastMonth_','N_F_I_PastMonth_','DA_F_I_PastMonth_','AA_F_I_PastMonth_','total_pm',...
    'ctq_physicalAbuse','ctq_emotionalAbuse','ctq_sexualAbuse','ctq_emotionalNeglect','ctq_physicalNeglect','ctq_total',...
    'stai_x1_total','stai_x2_total','ces_total',...
    'des_taxon_sum','des_depersonalization_sum','des_amnsetic_sum','des_absorption_sum','des_total_sum',...
    'comp1','comp2','comp3',...
    'G_risk25','G_risk50','G_risk75','G_amb24','G_amb50','G_amb74',...
    'L_risk25','L_risk50','L_risk75','L_amb24','L_amb50','L_amb74',...
    'isExcluded_behavior','isExcluded_imaging',...
    'age','isMale','kbit',...
    'G_risk','G_amb','G_amb_risk50','L_risk','L_amb','L_amb_risk50',...
    'G_alpha','G_beta','L_alpha','L_beta'});


%% create covariate file in the format of .mat file, with id sorted based onthe GLM order

% load id order in GLM
id_order_file = 'Order of subject for imaging covariate.xlsx';
order_id_day1 = xlsread(id_order_file, '65 subjects for day1');
order_id_day2 = xlsread(id_order_file, '67 subjects for day2 and both');


% day2
names_day2 = day2_cov.Properties.VariableNames';
names_day2(1) = {'subject_num'};
var_day2 = table2array(day2_cov);

[sorted_reference_id, reference_sort_idx] = sort(order_id_day2);
[sorted_id, sort_idx] = sort(var_day2(:,1));
reorder_idx(reference_sort_idx) = sort_idx;

var_day2_sorted = var_day2(reorder_idx, :);
clear reorder_idx reference_sort_idx sorted_reference_id sorted_id sort_idx

save('COV_67subj_day2_separate_fit.mat','names_day2','var_day2_sorted');


% day1, some subject were not included
day1exclude = [1232, 1278];
names_day1 = day1_cov.Properties.VariableNames';
names_day1(1) = {'subject_num'};
var_day1 = table2array(day1_cov);
var_day1 = var_day1(~ismember(var_day1(:,1),day1exclude),:);

[sorted_reference_id, reference_sort_idx] = sort(order_id_day1);
[sorted_id , sort_idx] = sort(var_day1(:,1));
reorder_idx(reference_sort_idx) = sort_idx;

var_day1_sorted = var_day1(reorder_idx, :);
clear reorder_idx reference_sort_idx sorted_reference_id sorted_id sort_idx

save('COV_65subj_day1_separate_fit.mat','names_day1','var_day1_sorted');

