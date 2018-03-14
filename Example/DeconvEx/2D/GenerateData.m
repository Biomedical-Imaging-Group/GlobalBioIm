function [im,psf,y]=GenerateData(noise)
%----------------------------------------------
% function [im,psh,y]=GenerateData(noise)
%
% Generate data for 2D deconvolution examples
%
% Input:   noise -> type of noise
%                    - 'Gaussian'
%                    - 'Poisson'
%
% Outputs: im    -> ground truth (star like object)
%          psf   -> psf 
%          y     -> blurred and noisy data
%
%  Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%----------------------------------------------

%% Ground truth
N=256;sz=[N N];                       % Image size
im=StarLikeSample(2,N,25,20,3,0.7);   % Star-like object (help StarLikeSample to see the details of parameters)

%% PSF 
lamb=561;                % Illumination wavelength
res=30;                  % Resolution (nm)
Na=1.4;                  % Objective numerica aperture
fc=2*Na/lamb*res;        % cut-off frequency
ll=linspace(-0.5,0,sz(1)/2+1);
lr=linspace(0,0.5,sz(1)/2);
[X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
[th,rho]=cart2pol(X,Y);
OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
psf=real(ifft2(OTF));

%% Data
H=LinOpConv(OTF);
y_noNoise=H*im;
if strcmp(noise,'Gaussian')
    sigN=3e-2;
    y=y_noNoise+sigN*randn(sz);
elseif strcmp(noise,'Poisson')
    photBud=110;
    y=poissrnd(round(y_noNoise*photBud))/photBud;
else
    error('Wrong type of noise');
end
disp(['SNR data : ',num2str(20*log10(norm(y_noNoise(:))/norm(y_noNoise(:)-y(:))))]);

end