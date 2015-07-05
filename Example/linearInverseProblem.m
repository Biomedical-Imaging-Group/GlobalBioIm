
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

clear; close all; clc;

disp('DECONVOLUTION');
% Data 
saturn = fitsread('saturn.fits');

% PSF
psf = fitsread('saturn_psf.fits');
psf = padarray(psf,[230,230]);
psf = psf(1:end-1, 1:end-1);
psf = psf/ sum(sum(psf));

% convolution operator
H = Convolution(fftshift(psf));
% Finite difference operator
D = Grad(size(saturn));

mu =1e-2; %hyperparameter

% DECONVOLUTION
% x =argmin_x || H.x - y||_2^2 + mu  || D.x||_2^2  
% normal equation  (H'H + mu D'D) x = H'.y
%                               A x = b

A = H'*H + mu * D'*D;

A = OneToMany({H,D},[1, mu]);

b = H'* saturn;
x = ConjGrad(A,b,  zeros(size(saturn)),100,{1,1});

figure;
subplot(1, 2, 1);
imagesc(saturn,[0,max(max(saturn))]);
title('Data');
subplot(1, 2, 2);
imagesc(x,[0,max(max(x))]);
title('Deconvolution');

input('DECONVOLUTION with missing pixels');
% DECONVOLUTION with missing pixels
% x =argmin_x || H.x - y||_W^2 + mu  || D.x||_2^2  
% normal equation  (H'WH + mu D'D) x = H'.y
%                               A x = b
MissingFraction = 0.99;
Missing = (saturn>0)&(rand(size(saturn))>MissingFraction); 
data = saturn .* Missing;
W = Diagonal(double(Missing));

mu =1e-1; %hyperparameter

A = H'*W*H + mu * D'*D;
b = H'* W*data;
x = ConjGrad(A,b,  zeros(size(saturn)),100);

figure;
subplot(1, 2, 1);
imagesc(data,[0,max(max(data))]);
title(['Data (',num2str(100*MissingFraction),'% Missing pixels)']);
subplot(1, 2, 2);
imagesc(x,[0,max(max(x))]);
title('Deconvolution');
