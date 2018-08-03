%------------------------------------------------------------------------
% This script outlines the performances of a direct computation of the
% proximal operator of f(x)= 0.5||SHx - y||_2^2 where x is a 3D volume, S
% a downsampling operator (LinOpDownsample) the slices and H a convolution
% operator (LinOpConv). The direct computation of the prox [1,2] is compared with a
% resolution using the Conjugate Gradient algorithm.
%
% **References**
% [1] Zhao Ningning et al. "Fast Single Image Super-Resolution Using a New 
%     Analytical Solution for l2-l2 Problems".
%     IEEE Transactions on Image Processing, 25(8), 3683-3697 (2016).
% [2] Emmanuel Soubies and Michael Unser. "Computational Super-Sectioning for Single-Slice
%     Structured-Illumination Microscopy" (2018)
%------------------------------------------------------------------------
close all;
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

%% Parameters
downFact=[2 2];    % Downsampling factor along each direction        
regTikh0=1e-5;     % Regularization parameter Tikhonov order 0 (i.e. alpha in prox equal 1/regTikh0)

%% Load Data
[im,psf,y]=GenerateData('Gaussian',0);

%% Operators  
sz=size(psf);
H=LinOpConv(fft2(psf));                  % Convolution operator
S=LinOpDownsample(H.sizein,downFact);    % Downsmpling operator
fwd=S*H;                                 % Forward model
In=LinOpIdentity(S.sizein);              % Identity operator of sizein
Im=LinOpIdentity(S.sizeout);             % Identity operator of sizeout

%% Data generation 
y=H*im;
data=S*y;
figure; 
subplot(1,3,1); imagesc(im); axis image;  colormap gray;title('Original Image');
subplot(1,3,2); imagesc(y); axis image;  colormap gray;title('Convolved Image');
subplot(1,3,3); imagesc(data); axis image; colormap gray;title('Downsampled Convolved Image');

%% Cost L2 prox computation
disp('=============================================================');
disp('===== Test Prox (Direct inversion with Tikhonov order 0) ====');
disp('=============================================================');
% - Cost function 
%   0.5 ||SHx - data||_2^2 + regTikh0/2 ||x - x0||_2^2
% which corresponds to the prox of 0.5 ||SHx - y||_2^2 for alpha=1/regTikh0.
x0=zeros_(size(im));
b=fwd'*data/regTikh0 + x0;
cost=CostL2([],data)*fwd+regTikh0*CostL2([],x0);

% - Direct inversion (using Woodbury Formula)
tic;A=(In - fwd'*(regTikh0*Im + fwd*fwd')^(-1)*fwd);
xTikh0=A*b;t=toc;
disp(['Direct computation using Woodbury: ',num2str(t),' s']);

% -  Direct inversion with Less FFTs (Corollary III.2 in [2])
Lamb=LinOpDiag([],H.mtf);
P=LinOpSumPatches(S.sizein,S.sizein./S.df);
D=LinOpDiag([],P*abs(H.mtf).^2+regTikh0*prod(S.df));
tic;xTikh0_bis=real(ifftn((In - Lamb'*P'*D^(-1)*P*Lamb)*fftn(b)));t=toc;
disp(['Direct computation using Proposition 1: ',num2str(t),' s']);
disp(['Error between the two direct computations: ',num2str(norm(xTikh0(:)-xTikh0_bis(:)))]);

% - Comparison with the prox of CostL2Composition (implemented according to Corollary III.2 in [2])
cc=CostL2([],data)*fwd;cc.doPrecomputation=1;
tic;xTikh0_prox=cc.applyProx(x0,1/regTikh0);t=toc;
disp(['Use of the implemented prox: ',num2str(t),' s']);
tic;xTikh0_prox=cc.applyProx(x0,1/regTikh0);t=toc;
disp(['Second prox call (effect of precomputation): ',num2str(t),' s']);
disp(['Error between prox and direct computation : ',num2str(norm(xTikh0(:)-xTikh0_prox(:)))]);

% - Comparison with CG
AA=fwd'*fwd+regTikh0*In;
CG=OptiConjGrad(AA,b*regTikh0);
CG.CvOp=TestCvgStepRelative(1e-6);
CG.verbose=true;CG.OutOp.iterVerb=20;CG.OutOp.saveXopt=1;
CG.maxiter=100;
CG.ItUpOut=1;
CG.run(x0);
t2=CG.time;
disp(['Iterative computation with CG: ',num2str(CG.time),' s']);
disp(['Speedup x',num2str(round(t2/t))]);
disp('=============================================================');

% -- Display
figure;
semilogy(CG.OutOp.iternum,cellfun(@(x) cost*x,CG.OutOp.evolxopt),'LineWidth',2); hold all; set(gca,'FontSize',14);grid;
semilogy(CG.OutOp.iternum,ones(1,length(CG.OutOp.iternum))*(cost*xTikh0_prox),'LineWidth',2);
xlabel('Iterations');ylabel('Cost function');
legend('CG','Direct');
figure;
subplot(1,3,1);imagesc(xTikh0_prox); axis image; title(['Direct : ',num2str(t),' s'])
subplot(1,3,2);imagesc(CG.xopt); axis image; title(['CG : ',num2str(t2),' s'])
subplot(1,3,3);imagesc(abs(xTikh0_prox-CG.xopt)); axis image; colorbar; title('Diff');
colormap gray;

%% For Unitary Tests
global generateDataUnitTests stateTest
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valxTikh0=xTikh0;save('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0','valxTikh0');
        valxTikh0_bis=xTikh0_bis;save('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0_bis','valxTikh0_bis');
        valxTikh0_prox=xTikh0_prox;save('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0_prox','valxTikh0_prox');
        xopt=CG.xopt;save('Util/UnitTest/Data/TestProxL2DownSampledConv_xopt','xopt');
    else
        load('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0');
        load('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0_bis');
        load('Util/UnitTest/Data/TestProxL2DownSampledConv_xTikh0_prox');
        load('Util/UnitTest/Data/TestProxL2DownSampledConv_xopt');
        stateTest=isequal(xopt,CG.xopt) &&  isequal(valxTikh0,xTikh0) &&  isequal(valxTikh0_bis,xTikh0_bis) &&  isequal(valxTikh0_prox,xTikh0_prox);
    end
end