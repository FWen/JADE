function [azm, ele] = find_max_peak(spec, azm_range, ele_range)
    
[m1, im1] = max(spec);
[m2, im2] = max(m1);
azm = azm_range(im2);
ele = ele_range(im1(im2));

% [AZ,EL] = meshgrid(azm_range, ele_range);
% figure(1); mesh(AZ,EL,spec); xlabel('Azimuth (degree)'); ylabel('Elevation (degree)');
% ylim([azm_range(1), azm_range(end)]); zlim([0 1.1]);%view(90,0);

end