% This program is used to calculate various parameters of the scattering medium

function [g, mus, musp, s11, s12, s33, s34, s1, s2] = mie_cal(dia, lambda, n_water, N, fv, m_c)

mu = cos(pi/(N-1)*((1:N)-1));
xx = pi*dia/lambda*n_water; 
s1=zeros(1,N);
s2=zeros(1,N);
for i=1:N
    result_s=Mie_S12(m_c,xx,mu(i));
    s1(i)=result_s(1);
    s2(i)=result_s(2);
end

% Scattering matrix elements
s11=(abs(s2).^2+abs(s1).^2)/2;
s12=(abs(s2).^2-abs(s1).^2)/2;
s33=(conj(s1).*s2+conj(s2).*s1)/2;
s34=1i*(conj(s1).*s2-conj(s2).*s1)/2;

result_mie=Mie(m_c,xx);
qsca = result_mie(5);           % Scattering efficiency
qext = result_mie(4);           % Transport efficiency

Vsphere = 4/3*pi*(dia/2)^3;     % Particle volume
rho     = fv/Vsphere;           % Particle concentration, #/um^3

g = result_mie(8);              % Forward scattering coefficient
A       = pi*dia^2/4;           % Geometrical cross-sectional area, um^2
sigma_s = qsca*A;               % Scattering cross-section, um^2
mus     = sigma_s*rho*1e4;      % Scattering coefficient, cm^-1
musp    = mus*(1-g);            % Reduced scattering coefficient, cm^-1

end