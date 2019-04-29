clearvars

%% Elman et al 2009 bilateral striatum
voi = xff('new:voi');
 
% add ROI
voi.AddSphericalVOI([24, 0, 2], 6);
voi.AddSphericalVOI([-2, 8, 11], 6);

voi.SaveAs('elman2009_bilateral_str.voi');

%% Sailer et al 2008
voi = xff('new:voi');

talCoord = [-23, 31, 54;...
            -53, -52, 28;...
            -4, 7, -2;...
            5, 9, 1;...
            -8, 56, -13;...
            -7, -11, 5];
        
for i=1:size(talCoord, 1)
    voi.AddSphericalVOI(talCoord(i,:), 5);
end

voi.SaveAs('sailer2008_patients-controls.voi');

