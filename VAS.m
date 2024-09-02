%% Main program of the VAS model
clc
clear
%% Predifing Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

%% Basic parameters
% scattering medium parameters
n_water = 1.33;                    % Medium refractive index
dia = 0.5;                         % Particle diameter, μm
mr = 1.59;                         % Real part of the particle refractive index
mi = 0;                            % Imag part of the particle refractive index
m_c0 = mr+1i*mi;                   % Particle refractive index
m_c = m_c0/n_water;                % Relative refractive index
fv = 0.05/16.5;                    % Dilution, %

% input light field parameters
lambda = 0.532;                    % Wavelength, μm
N_obj = [1000,1000];               % Light field size
k = 2*pi/lambda;                   % Wave vector

% calculation parameters
n0 = m_c0*fv+n_water*(1-fv);       % Mean refractive index
dx_pixel = lambda/4;               % Pixel size
N_diffuser = 40;                   % Layers of the scattering medium
d = 50;                            % Distance between layers


%% Mie calculation
N = 10001;                         % Number of angular subdivisions

% Basic parameters of the scattering medium
% g-forward scattering coefficient( g=<cosθ> )
% mus-scattering coefficient，musp-reduced scattering coefficient
[g, mus, musp, s11, s12, s33, s34, s1, s2] = mie_cal(dia, lambda, n_water, N, fv, m_c);

ls = 1/mus*1e4;                    % Scattering mean free path，μm

% Calculation of single layer polarization transformation ratio
r_pol = r_cal(s11, s33, N, mus, g);


%% Scattering phase screen calculation
% Root mean square of the scattering phase distribution
sigma_p = sqrt(d/ls);

% Calculation of single scattering angular spectral distribution
Xsize = N_obj(2)*dx_pixel;
du = 2*pi/(Xsize);
umax = pi/(dx_pixel);
u = -umax:2*umax/N_obj(2):umax-2*umax/N_obj(2);
[U,V] = meshgrid(u,u);
k2 = (U.^2+V.^2);
eva = double(k2<(2*pi/lambda)^2);
theta = asin(sqrt(k2/k^2).*eva);
ntheta = fix(theta/pi*(N-1))+1;
[fai,rr] = cart2pol(U,V);
Fx = (s11(ntheta)+s12(ntheta).*cos(2*fai))/4/pi.*((exp(-d/ls)-exp(-d/ls./cos(theta)))./((1-cos(theta))+(1-eva)*1e-4)).*eva;
Fx(N_obj(1)/2+1,N_obj(2)/2+1) = s11(1)/4/pi*(d/ls)*exp(-d/ls);
Fx = Fx/Fx(N_obj(1)/2+1,N_obj(2)/2+1);
Fx = sqrt(Fx);
Fy = Fx';

% Calculation of energy loss for single-layer backscattering
mu = cos(pi/(N-1)*((1:N)-1));
mu1 = cos(pi/(N-1)*((N/2+1:N)-1));
fai1 = 0:2*pi/1000:2*pi-2*pi/1000;
[Mu,Fai] = meshgrid(mu,fai1);
[Mu1,Fai1] = meshgrid(mu1,fai1);
[dM,dN] = gradient(Mu);
[dM1,dN1] = gradient(Mu1);
ntheta1 = repmat(1:N,[1000,1]);
ntheta2 = repmat(floor(N/2+1:N),[1000,1]);
PP0 = s11(ntheta1)+s12(ntheta1).*cos(2*Fai);
PP1 = sum(PP0);
PP2 = PP1/sum(PP1.*abs(gradient(mu)));
PP = PP2(floor(N/2+2):end);

Los = zeros(1,floor(N/2));
for ii = 1:floor(N/2)
    mu0 = mu(ii);
    Los(ii) = sum(PP/2.*(1-exp(-d/ls/mu0+d/ls./mu1))./(mu0-mu1)*mu0.*mu1.*gradient(mu1));
end
LO = Los(ntheta);
Res = (1-LO).*eva;


% Calculate the phase screens for all layers
o_slice_x = zeros(N_obj(1),N_obj(2),N_diffuser);
o_slice_y = zeros(N_obj(1),N_obj(2),N_diffuser);

for m = 1:N_diffuser


    ph_seed_x = normrnd(0,sigma_p,N_obj);
    ph_seed_y = normrnd(0,sigma_p,N_obj);

    ang_x=angle(Ft(F(exp(1i*ph_seed_x)).*Fx));
    ang_y=angle(Ft(F(exp(1i*ph_seed_y)).*Fy));

    ph_mask_x=(ang_x-mean(ang_x(:)))/std2(ang_x)*sigma_p;
    ph_mask_y=(ang_y-mean(ang_y(:)))/std2(ang_y)*sigma_p;

    o_slice_x(:,:,m)=exp(1i*ph_mask_x);
    o_slice_y(:,:,m)=exp(1i*ph_mask_y);

end


%% Defining the input light field
i0 = zeros(N_obj(1),N_obj(2));
i0(:,:) = 1;


%% Transmission simulation
% Initialize the computational matrices
% The polarization direction is x here
phix1 = i0;
phiy1 = zeros(N_obj(1),N_obj(2));

phix2 = i0;
phiy2 = zeros(N_obj(1),N_obj(2));

psix = i0;
psiy = zeros(N_obj(1),N_obj(2));

psix1 = i0;
psiy1 = zeros(N_obj(1),N_obj(2));

psix2 = i0;
psiy2 = zeros(N_obj(1),N_obj(2));

ptix = zeros(N_obj(1),N_obj(2),N_diffuser+1);
ptiy = zeros(N_obj(1),N_obj(2),N_diffuser+1);
ptix(:,:,1) = i0;

bal = i0;

% Calculation of the angular spectral transport operator H
x = -N_obj(1)/2:N_obj(1)/2-1;
y = -N_obj(2)/2:N_obj(2)/2-1;
LX = N_obj(1)*dx_pixel;
LY = N_obj(2)*dx_pixel;
u = lambda*x/LX;
v = lambda*y/LY;
[uu,vv] = meshgrid(u,v);

k2 = (uu.^2+vv.^2);
eva = double(k2/lambda^2<(1/lambda)^2);

H = exp(1i*k*d*n0*sqrt((1-uu.^2-vv.^2).*eva));

% Layered transmission, repeating the steps for each layer:
% Compose phase mask → free space transmission → polarization transformation
for m = 2:(N_diffuser+1)

    % Calculate the polarization transformation factor of this layer
    c_pol1 = sqrt(0.5*(1+exp(-r_pol*d*1e-4*(m-1))));
    c_pol2 = sqrt(0.5*(1+exp(-r_pol*d*1e-4)));

    % Separation of ballistic and scattered light components
    phix1 = bal;
    phix2 = ptix(:,:,m-1)-phix1;
    phiy2 = ptiy(:,:,m-1)-phiy1;

    % Compose phase mask and ballistic light attenuation
    psix1 = phix1.*o_slice_x(:,:,m-1);
    psix2 = phix2.*o_slice_x(:,:,m-1);
    psiy1 = phiy1.*o_slice_y(:,:,m-1);
    psiy2 = phiy2.*o_slice_y(:,:,m-1);

    bal = bal*sqrt(exp(-m*d/ls)/exp(-(m-1)*d/ls));

    % Polarization transformation
    psix1_x=(psix1-bal)*c_pol1+bal;
    psix1_y=(psix1-bal)*(1-c_pol1);
    psiy1_y=(psiy1-mean2(psiy1))*c_pol1+mean2(psiy1);
    psiy1_x=(psiy1-mean2(psiy1))*(1-c_pol1);

    psix2_x=psix2*c_pol2;
    psix2_y=psix2*(1-c_pol2);
    psiy2_y=psiy2*c_pol2;
    psiy2_x=psiy2*(1-c_pol2);

    psix = psix1_x + psiy1_x + psix2_x + psiy2_x;
    psiy = psix1_y + psiy1_y + psix2_y + psiy2_y;

    % Free space transmission
    ptix(:,:,m) = Ft((F(psix)).*H.*Res);
    ptiy(:,:,m) = Ft((F(psiy)).*H.*Res);
    bal = Ft((F(bal)).*H);

end

% Output field
zz=5000;  % location of the imaging plane
H1=exp(1i*k*zz*n0*sqrt((1-uu.^2-vv.^2).*eva));
outputx=Ft((F(ptix(:,:,end))).*H1);
outputy=Ft((F(ptiy(:,:,end))).*H1);
figure
imagesc(x,y,abs(outputx).^2);title('amplitude of x output field');
figure
imagesc(x,y,abs(outputy).^2);title('amplitude of y output field');


%% Calculate simulated g and ls
fg=abs(F(ptix(:,:,round(ls/d+1))))/sum(sum(abs(i0).^2));
fg1=fg.^2;
g_cal=sum(sum(fg1.*cos(theta).*eva))
ls_error=fg1(N_obj(1)/2+1,N_obj(1)/2+1)/exp(-1)-1

