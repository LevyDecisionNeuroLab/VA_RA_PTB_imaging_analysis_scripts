clearvars
%% 

% First, we need to load the mdm for our model - VTC + PRT preferable
mdm = xff('RA_PTSD_none_prt_forppi.mdm');
mdm = xff('RA_PTSD_SV_prt_forppi.mdm');
mdm = xff('RA_PTSD_SV_day1_prt_forppi.mdm');
mdm = xff('RA_PTSD_SV_day2_prt_forppi.mdm');
mdm = xff('RA_PTSD_RewardValue_prt.mdm');
mdm = xff('RA_PTSD_RiskLevel_prt.mdm');

% VOI
seedvoi = xff('Bartra13_SV_ROI_original_fromZhihao_vStr.voi');
seedvoi = xff('Bartra13_SV_ROI_original_fromZhihao_vmPFC.voi');

% psycho term
psy_term =  {{'Amb_gains_Display x p1 + Risk_gains_Display x p1'}};
psy_term = {{'Amb_gains_Display + Risk_gains_Display > Amb_loss_Display + Risk_loss_Display'}};
psy_term = {{'Risk_gains_Display x p1 + Risk_loss_Display x p1'}};

psy_term = {{'Risk_gains_Display x p1 + Risk_loss_Display x p1'}};
psy_term = {{'Amb_gains_Display x p1 + Amb_loss_Display x p1 '}};

psy_term =  {{'Amb_gains_Display x p1'}}; % ambig gain
psy_term =  {{'Risk_gains_Display x p1'}}; % risk gain
psy_term = {{'Amb_loss_Display x p1'}}; % ambig loss
psy_term = {{'Risk_loss_Display x p1'}}; % risk loss



% output name
glm_name1 = ['RA_PTSD_ppi_SV_gains_' 'ROI_Bartra_original_vStr_' 'GLM_SV.glm'];
glm_name2 = ['RA_PTSD_ppi_none_gainsloss_contrast_' 'ROI_Bartra_original_vStr_' 'GLM_none.glm'];
glm_name3 = ['RA_PTSD_ppi_SV_loss_' 'ROI_Bartra_original_vStr_' 'GLM_SV.glm'];
glm_name4 = ['RA_PTSD_ppi_SV_gains_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV.glm'];
glm_name5 = ['RA_PTSD_ppi_SV_loss_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV.glm'];
glm_name6 = ['RA_PTSD_ppi_SV_day1_gain_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day1.glm'];
glm_name7 = ['RA_PTSD_ppi_SV_day1_gain_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day1.glm'];
glm_name8 = ['RA_PTSD_ppi_SV_day2_gain_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day2.glm'];
glm_name9 = ['RA_PTSD_ppi_SV_day1_gain_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV_day1.glm'];
glm_name10 = ['RA_PTSD_ppi_SV_day2_gain_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV_day2.glm'];
glm_name11 = ['RA_PTSD_ppi_RV_gain_' 'ROI_Bartra_original_vStr_' 'GLM_RV.glm'];
glm_name12 = ['RA_PTSD_ppi_RV_gain_' 'ROI_Bartra_original_vmPFC_' 'GLM_RV.glm'];
glm_name13 = ['RA_PTSD_ppi_RV_loss_' 'ROI_Bartra_original_vStr_' 'GLM_RV.glm'];
glm_name14 = ['RA_PTSD_ppi_RV_loss_' 'ROI_Bartra_original_vmPFC_' 'GLM_RV.glm'];
glm_name15 = ['RA_PTSD_ppi_SV_day1_loss_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day1.glm'];
glm_name16 = ['RA_PTSD_ppi_SV_day1_loss_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV_day1.glm'];
glm_name17 = ['RA_PTSD_ppi_SV_day2_loss_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day2.glm'];
glm_name18 = ['RA_PTSD_ppi_SV_day2_loss_' 'ROI_Bartra_original_vmPFC_' 'GLM_SV_day2.glm'];

glm_name19 = ['RA_PTSD_ppi_Risk_gain_' 'ROI_Bartra_original_vStr_' 'GLM_RiskLevel.glm'];
glm_name20 = ['RA_PTSD_ppi_Risk_loss_' 'ROI_Bartra_original_vStr_' 'GLM_RiskLevel.glm'];
glm_name21 = ['RA_PTSD_ppi_Risk_gain_' 'ROI_Bartra_original_vmPFC_' 'GLM_RiskLevel.glm'];
glm_name22 = ['RA_PTSD_ppi_Risk_loss_' 'ROI_Bartra_original_vmPFC_' 'GLM_RiskLevel.glm'];

glm_name23 = ['RA_PTSD_ppi_SV_day1_ambiggain_' 'ROI_Bartra_original_vStr_' 'GLM_SV_day1.glm'];

%% for paralell original bartra ROI

mdm_day1 = 'RA_PTSD_SV_day1_prt_forppi.mdm';
mdm_day2 = 'RA_PTSD_SV_day2_prt_forppi.mdm';

seedvoi_vstr = 'Bartra13_SV_ROI_original_fromZhihao_vStr.voi';
seedvoi_vmpfc = 'Bartra13_SV_ROI_original_fromZhihao_vmPFC.voi';

psy_term_ag =  {{'Amb_gains_Display x p1'}}; % ambig gain
psy_term_rg =  {{'Risk_gains_Display x p1'}}; % risk gain
psy_term_al = {{'Amb_loss_Display x p1'}}; % ambig loss
psy_term_rl = {{'Risk_loss_Display x p1'}}; % risk loss


mdm_all = {mdm_day1, mdm_day2};
seed_all = {seedvoi_vstr, seedvoi_vmpfc};
psy_term_all = {psy_term_ag, psy_term_rg, psy_term_al, psy_term_rl};

glm_name_all = {...
    ['RA_PTSD_ppi_SV_day1_AG_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_RG_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_AL_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_RL_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_AG_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_RG_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_AL_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_RL_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_AG_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_RG_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_AL_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_RL_' 'ROI_Bartra_original_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_AG_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_RG_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_AL_' 'ROI_Bartra_original_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_RL_' 'ROI_Bartra_original_vmPFC.glm'],... 
    };

pooljob = parpool('local', 8); 

for mdm_idx = 1:length(mdm_all)
    for seed_idx = 1:length(seed_all)
        parfor psy_term_idx = 1:length(psy_term_all)
            
            mdm = xff(mdm_all{mdm_idx});
            seedvoi = xff(seed_all{seed_idx});
            psy_term = psy_term_all{psy_term_idx};
            
            glm_name_idx = ...
                (mdm_idx-1)*(length(seed_all)*length(psy_term_all)) + ...
                (seed_idx-1)*length(psy_term_all) + ...
                psy_term_idx
            
            opt = struct( ...
                'ppivoi',    seedvoi, ...
                'ppicond',   psy_term, ...
                'tfilter',   'Inf', ...
                'tfilttype', 'dct',...    
                'outfile', glm_name_all{glm_name_idx});  
            
            
            glm = mdm.ComputeGLM(opt);
              
              fprintf('GLM number%d \n', glm_name_idx)
              fprintf([mdm_all{mdm_idx},' ', seed_all{seed_idx}, ' ', psy_term_all{psy_term_idx}{1}{1} '\n'])
              fprintf('GLM computed: %s \n\n', glm_name_all{glm_name_idx})
 
        end
    end
end

%% for paralell, peak sphere bartra ROI

mdm_day1 = 'RA_PTSD_SV_day1_prt_forppi.mdm';
mdm_day2 = 'RA_PTSD_SV_day2_prt_forppi.mdm';

seedvoi_vstr = 'Bartra13_SV_ROI_peakcoord_5mm_fromZhihao_vStr.voi';
seedvoi_vmpfc = 'Bartra13_SV_ROI_peakcoord_5mm_fromZhihao_vmPFC.voi';

psy_term_ag =  {{'Amb_gains_Display x p1'}}; % ambig gain
psy_term_rg =  {{'Risk_gains_Display x p1'}}; % risk gain
psy_term_al = {{'Amb_loss_Display x p1'}}; % ambig loss
psy_term_rl = {{'Risk_loss_Display x p1'}}; % risk loss


mdm_all = {mdm_day1, mdm_day2};
seed_all = {seedvoi_vstr, seedvoi_vmpfc};
psy_term_all = {psy_term_ag, psy_term_rg, psy_term_al, psy_term_rl};

glm_name_all = {...
    ['RA_PTSD_ppi_SV_day1_AG_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_RG_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_AL_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_RL_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day1_AG_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_RG_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_AL_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day1_RL_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_AG_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_RG_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_AL_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_RL_' 'ROI_Bartra_peakcoord_vStr.glm'],...
    ['RA_PTSD_ppi_SV_day2_AG_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_RG_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_AL_' 'ROI_Bartra_peakcoord_vmPFC.glm'],...
    ['RA_PTSD_ppi_SV_day2_RL_' 'ROI_Bartra_peakcoord_vmPFC.glm'],... 
    };

pooljob = parpool('local', 8); 

for mdm_idx = 1:length(mdm_all)
    for seed_idx = 1:length(seed_all)
        parfor psy_term_idx = 1:length(psy_term_all)
            
            mdm = xff(mdm_all{mdm_idx});
            seedvoi = xff(seed_all{seed_idx});
            psy_term = psy_term_all{psy_term_idx};
            
            glm_name_idx = ...
                (mdm_idx-1)*(length(seed_all)*length(psy_term_all)) + ...
                (seed_idx-1)*length(psy_term_all) + ...
                psy_term_idx
            
            opt = struct( ...
                'ppivoi',    seedvoi, ...
                'ppicond',   psy_term, ...
                'tfilter',   'Inf', ...
                'tfilttype', 'dct',...    
                'outfile', glm_name_all{glm_name_idx});  
            
            
            glm = mdm.ComputeGLM(opt);
              
              fprintf('GLM number%d \n', glm_name_idx)
              fprintf([mdm_all{mdm_idx},' ', seed_all{seed_idx}, ' ', psy_term_all{psy_term_idx}{1}{1} '\n'])
              fprintf('GLM computed: %s \n\n', glm_name_all{glm_name_idx})
 
        end
    end
end


%% compute single GLM

% computation option setup
opt = struct( ...
    'ppivoi',    seedvoi, ...
    'ppicond',   psy_term, ...
    'tfilter',   'Inf', ...
    'tfilttype', 'dct',...    
    'outfile', glm_name23);

% compute GLM
glm = mdm.ComputeGLM(opt);

