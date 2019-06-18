%%  VTC::VOITimeCourse
%   VTC::VOITimeCourse  - extract VOI time course data
%  
%   FORMAT:       voitc [, uvec] = vtc.VOITimeCourse(voi [, weight, fliplr])
%  
%   Input fields:
%  
%         voi         VOI file or coordinates (e.g. from VOI::BVCoords)
%         weight      (cell array of) Nx1 vector(s) with voxel weights
%                     give scalar 0 for unique, scalar 1 for no weighting,
%                     a scalar -1 for SVD after z-transform, or
%                     a scalar [Inf] to get a cell array of TxV arrays
%         fliplr      flip left/right (Z axes) for radiological convention
%  
%   Output fields:
%  
%         voitc       TxV time course of voi(s)
%         uvec        unique VTC voxel indices within VOI, so that
%                     voi.VOI(i).Voxels(uvec{i}) leads to those coordinates
%   
%   E.g. % load vtc and voi interactively
%          vtc = BVQXfile('*.vtc');
%          voi = BVQXfile('*.voi');

         % get each VOI's voxels time course and their unique indices
%          [voitc, voiuvec] = vtc.VOITimeCourse(voi, Inf);

%% find all VTC files of a single subject
clear all

root = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Imaging analysis\Imaging_anallysis_082018';
cd(root)
% load vois
voi1 = BVQXfile('PTSD_none_resp_corrected_femaleOut_noRemit_Day2_AllCond_capstotalCovar_p0.001Alphasim.voi'); 
voi2 = BVQXfile('Gilaie-DotanEtAl_2014_rPPC_7mmSphere.voi');

% subjs = [1063;1069;1072;115;1206;1208;1244;1266;1273;1284;1291;1304;1305;1309;1340;1344;1345;1346;30;38;53;56;58;60;75;83;95;96;99;105;1074;110;119;1205;1232;1237;1245;125;1280;1285;1350;45;82;85;87;88;93;98];
%subjs = [1063;1069;1072;115;1206;1208;1244;1266;1273;1284;1291;1304]
% subjs = [1305;1309;1340;1344;1345;1346;30;38;53;56]
% subjs = [58;60;75;83;95;96;99;105;1074;110]
subjs = [119;1205;1232;1237;1245;125;1280;1285;1350;45];
%subjs = [82;85;87;88;93;98];
corr = zeros(length(subjs),2);

% pooljob = parpool('local', 8);

for i=1:length(subjs)
    subj = subjs(i);
    cd(fullfile(root, 'BV_VTCs\'));
    vtcfiles = dir([num2str(subj),'_*.vtc']); % all vtc files for one subject
    vtc1 = zeros(490,length(vtcfiles));
    vtc2 = zeros(490,length(vtcfiles));
    
    corr_subj = zeros(1,2);
    
    for j = 1:length(vtcfiles)
       vtc = BVQXfile(vtcfiles(j).name);
       % extract VOI vtc
       [voitc1,voiuvec1] = vtc.VOITimeCourse(voi1,inf); % using [inf] will give you a V X T matrix of vtc, and you can average across all voxels
       [voitc2,voiuvec2] = vtc.VOITimeCourse(voi2,inf);

       vtc1(:,j) = mean(voitc1{1},2);
       vtc2(:,j) = mean(voitc2{1},2);
    end
    
    % link all sessions' VTCs together
    vtc1all = vtc1(:,1);
    vtc2all = vtc2(:,1);

    for n = 2:length(vtcfiles)
        vtc1all = [vtc1all;vtc1(:,n)];
        vtc2all = [vtc2all;vtc2(:,n)];
    end
    
    %correlation coefficient for each subject
    r = corrcoef(vtc1all, vtc2all);
    corr_subj(1,1) = subjs(i);
    corr_subj(1,2) = r(1,2);
    
    corr(i,:) = corr_subj
    
    % clear vtcfiles vtc1 vtc2 vtc voitc1 voitc2 vtc1all vtc2all r voiuvec1 voiuvec2
end

% delete(pooljob)