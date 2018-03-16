function [im,psf,y]=GenerateData3D(noise,level)
%----------------------------------------------
% function [im,psh,y]=GenerateData3D(noise,level)
%
% Generate data for 3D deconvolution example
%
% Input:   noise -> type of noise
%                    - 'Gaussian'  level is PSNR
%                    - 'Poisson'   level is average number of photons per pixel
%
% Outputs: im    -> ground truth (star like object)
%          psf   -> psf 
%          y     -> blurred and noisy data
%
%  Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%                     F. Soulez ferreol.soulez@univ-lyon1.fr
%----------------------------------------------

%% Ground truth
N=128;sz=[N N N];                       % Image size
im=StarLikeSample(3,N,6,20,5,0.7);     % Star-like object (help StarLikeSample to see the details of parameters)

%% PSF 
lamb=561;                % Illumination wavelength
res=60;                  % Resolution (nm)
Na=1.4;                  % Objective numerica aperture
fc=2*Na/lamb*res;        % cut-off frequency
ll=linspace(-0.5,0,sz(1)/2+1);
lr=linspace(0,0.5,sz(1)/2);
[X,Y,Z]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)],[ll,lr(2:end)]);
[th,ph,rho]=cart2sph(X,Y,Z);
OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
psf=real(ifftn(OTF));

%% Data
H=LinOpConv(OTF);
y_noNoise=H*im;
if strcmp(noise,'Gaussian')
    y = y_noNoise + max(y_noNoise(:)) .* 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise)));
elseif strcmp(noise,'Poisson')
    factor = level./mean(y_noNoise(:)) ;
    y_noNoise = y_noNoise.* factor;
    im = im.*factor;
    y = random('Poisson',y_noNoise);
else
    error('Wrong type of noise');
end

end