function [azm, ele] = az_el_estimation(b, azm_range, ele_range, array_pos)
    
d2pi = pi/180;
    
for ke=1:length(ele_range)
    rho = [cos(azm_range*d2pi).*cos(ele_range(ke)*d2pi); sin(azm_range*d2pi).*cos(ele_range(ke)*d2pi)];
    azel(ke,:) = abs(b'*exp(array_pos*rho));
end

[m1, im1] = max(azel);
[m2, im2] = max(m1);
azm = azm_range(im2);
ele = ele_range(im1(im2));

% [AZ,EL] = meshgrid(azm_range, ele_range);
% figure(1); mesh(AZ,EL,azel/max(azel(:))); xlabel('Azimuth (degree)'); ylabel('Elevation (degree)');
% ylim([azm_range(1), azm_range(end)]); zlim([0 1.1]);%view(90,0);

end