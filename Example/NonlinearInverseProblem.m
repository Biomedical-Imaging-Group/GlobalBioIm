
%     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

clear; %close all; 

disp('DECONVOLUTION');
% Data 
saturn = fitsread('saturn.fits');

% PSF
psf = fitsread('saturn_psf.fits');
psf = padarray(psf,[230,230]);
psf = psf(1:end-1, 1:end-1);
psf = psf/ sum(sum(psf));


% DECONVOLUTION with missing pixels
% x =argmin_x || H.x - y||_W^2 + mu  || D.x||_2^2  
% normal equation  (H'WH + mu D'D) x = H'.y
%                               A x = b
MissingFraction = 0.99;
Missing = (saturn>0)&(rand(size(saturn))>MissingFraction); 
data = saturn .* Missing;
W = Diagonal(double(Missing));

% Non stationnary gaussian noise accounting for photon noise
% w = zeros(size(saturn));
% w(Missing) = 100./(saturn(Missing) + 100);
% w(Missing) = w(Missing)./ mean(mean(w(Missing)));
% W = Diagonal(double(w));


% convolution operator
H = Convolution(fftshift(psf));
% Finite difference operator
D = Grad(size(saturn));

B = Identity(size(data));
zProx = JointL1(3); % JointL1( D*x)  == Total variation 
%zProx = L2();
tProx = NonNegativity();
rho1 = 1e-3;
rho2 =1e-3;

mu =.1; %hyperparameter
x0 = zeros(size(data));
cgmaxiter = 5;
maxiter =100;
x=ADMM_Restore(H,D,B,W,data, zProx, tProx,mu, rho1, rho2,x0,maxiter,cgmaxiter);


figure;
subplot(1, 2, 1);
imagesc(data,[0,max(max(data))]);
title(['Data (',num2str(100*MissingFraction),'% Missing pixels)']);
subplot(1, 2, 2);
imagesc(x,[0,max(max(x))]);
title('Deconvolution');
