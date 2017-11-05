function a = arst(K_pos, azm, ele)

% array steering matrix

d2pi = pi/180;

rho = [cos(azm*d2pi).*cos(ele*d2pi); sin(azm*d2pi).*cos(ele*d2pi)];
a = exp(K_pos*rho);

end