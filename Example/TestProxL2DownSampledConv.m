%------------------------------------------------------------------------
% This script test the implementation of the CostL2Composition prox when
% the comnposed Operator is of the form SH where S is a LinOpDownsample and
% H a LinOpConv. The direct computation of the prox [1] is compared with a
% resolution using the Conjugate Gradient algorithm.
%
% [1] Zhao Ningning et al. "Fast Single Image Super-Resolution Using a New 
%     Analytical Solution for l2-l2 Problems".
%     IEEE Transactions on Image Processing, 25(8), 3683-3697 (2016).
%------------------------------------------------------------------------
clear; close all; clc;
help TestProxL2DownSampledConv
%--------------------------------------------------------------
%  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------

%% Load Data
load('DeconvEx/2D/StarLikeSample');
load('DeconvEx/2D/psf');idx=129:384;
psf=fftshift(psf); psf=psf(idx,idx);

%% Operators
downFact=[2 2];                          % Downsampling factor along each direction          
H=LinOpConv(fft2(fftshift(psf)));        % Convolution operator
S=LinOpDownsample(H.sizein,downFact);    % Downsmpling operator
SH=S*H;

%% Generation data
y=H*im;
ys=S*y;
figure; 
subplot(1,3,1); imagesc(im); axis image;  colormap gray;title('Original Image');
subplot(1,3,2); imagesc(y); axis image;  colormap gray;title('Convolved Image');
subplot(1,3,3); imagesc(ys); axis image; colormap gray;title('Downsampled Convolved Image');

%% Cost L2 prox computation
% -- Cost L2
L2=CostL2([],ys);
F=L2*SH;
F.doPrecomputation=true;
% -- Parameters to call the prox
x=ones(H.sizein);
alpha=1e4;
% -- Prox cost
Fprox=alpha*F+CostL2([],x);

% -- Test the direct computation of the prox
tic; xp1=F.applyProx(x,alpha); t1=toc;
disp(['Computation time for direct method ',num2str(t1),' s']);

% -- Conjugate Gradient to solve the prox
A=SH'*SH + 1/alpha*LinOpIdentity(H.sizein);
b=SH'*ys + 1/alpha*x;
CG=OptiConjGrad(A,b,[],OutputOpti(0,[],0));
CG.maxiter=150;
CG.ItUpOut=2;
CG.run(zeros(size(x)));
xp2=CG.xopt;
t2=CG.time;
disp(['Computation time for CG method ',num2str(t2),' s']);
disp(['Speedup x',num2str(round(t2/t1))]);
for ii=1:length(CG.OutOp.evolxopt)
    evolFprox(ii)=Fprox*CG.OutOp.evolxopt{ii};
end

% -- Display
figure;
semilogy(CG.OutOp.iternum,evolFprox,'LineWidth',2); hold all; set(gca,'FontSize',14);grid;
semilogy(CG.OutOp.iternum,ones(size(evolFprox))*(Fprox*xp1),'LineWidth',2);
xlabel('Iterations');ylabel('Cost function');
legend('CG','Direct');
figure;
subplot(1,3,1);imagesc(xp1); axis image; title(['Direct : ',num2str(t1),' s'])
subplot(1,3,2);imagesc(xp2); axis image; title(['CG : ',num2str(t2),' s'])
subplot(1,3,3);imagesc(abs(xp1-xp2)); axis image; colorbar; title('Error');