% This program is used to calculate the polarization diffusivity x, 
% which characterizes how quickly the direction of polarization changes 
% with increasing depth of penetration for linearly polarized light inputs

function [r_pol ] = r_cal(s11, s33, N, mus, g)

mu = cos(pi/(N-1)*((1:N)-1));

d1 = (s11-s33)/sum(s11.*abs(gradient(mu)))/2;

sigma_d = sum(d1.*abs(gradient(mu)))*mus;

r_pol = sqrt((sigma_d*2+mus*(1-g))*mus*(1-g));

end