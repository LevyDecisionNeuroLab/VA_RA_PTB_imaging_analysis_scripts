clearvars

%%
glm = xff('RA_PTSD_ppi_SV_gains_ROI_Bartra_original_vStr_GLM_SV.glm');
vmp = xff('dummy_map.vmp');
cov = load('COV_67subj_all.mat');

caps = cov.variates(:, strcmp(cov.names, 'caps'));


glm.Help('RFX_conmaps')

glm.SubjectPredictors

subj2include = [1206, 1063];

mapopts = struct('subsel', subj2include);
con = [0,0,0,0,0,0,0,0,0,0,1,0]; % ppi term
con = [0,1,0,1,0,0,0,0,0,0,0,0]; % subjective value, gain

c = glm.RFX_conmaps(con,mapopts);

c_size = size(c); %[87, 60, 69, 67]

caps_cov = reshape(caps, [1,1,1,c_size(4)]);
caps_cov = repmat(caps_cov, [c_size(3), c_size(1), c_size(2)]);

glm.Help('VOIBetas')