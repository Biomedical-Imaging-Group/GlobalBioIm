%--------------------------------------------------------------------------
% This script outlines the performances of a direct computation of the
% proximal operator of f(x)= 0.5||SHx - y||_2^2 where x is a 3D volume, S
% an operator that sums (LinOpSum) the slices and H a convolution operator
% (LinOpConv)  that performs a 2D convolution for each slice.
%
% There is two parts in the script:
%    - Test of the direct prox computation and comparison with iterative
%    Conjugate Gradient.
%    - Use this prox in an iterative algorithm to perform 2D deconvolution
%    with the consideration of "outr-of-focus slices".
%
% **Reference**
% [1] Emmanuel Soubies and Michael Unser. "Computational Super-Sectioning for Single-Slice
%     Structured-Illumination Microscopy" (2018)
%--------------------------------------------------------------------------
close all;
help TestProxL2SumConv
%--------------------------------------------------------------
% Copyright (C) 2018, E. Soubies emmanuel.soubies@epfl.ch
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
nbSl=1;        % Number of slices used above and below the focal plane
regTikh0=1e-5;

%% Load Data
[im,psf,y]=GenerateData('Gaussian',0);
psf=repmat(psf,[1 1 2*nbSl+1]);            % Form a 3D psf by replicating the 2D one (but works also for different psf at each plane)
sz=size(psf);
tmp=zeros_(sz);
tmp(:,:,nbSl+1)=im;
tmp(:,:,1)=im(:,[41:sz(2),1:40]);
tmp(:,:,end)=im(:,[sz(2)-40:sz(2),1:sz(2)-40-1]);
im=tmp;

%% 0perators
mtf=Sfft(psf,3);   % 2D ffts of each PSF slices
H=LinOpConv(mtf,1,[1 2]);                  % Convolution operator slice-by-slice 2D
S=LinOpSum(sz,3);
fwd=S*H;                                   % Forward model
In=LinOpIdentity(sz);                      % Identity operator of sizein
Im=LinOpIdentity(sz(1:2));                 % Identity operator of sizeout

%% Generation data
data=fwd*im;
figure;
for ii=1:sz(3)
    subplot(2,sz(3),ii);
    imagesc(im(:,:,ii)); axis image;  colormap gray;title(['GT Slice #',num2str(ii)]);
end
subplot(2,sz(3),sz(3)+nbSl+1); imagesc(data); axis image;  colormap gray;title('Data');

%% Test Prox (Direct inversion with Tikhonov order 0)
disp('=============================================================');
disp('===== Test Prox (Direct inversion with Tikhonov order 0) ====');
disp('=============================================================');
% - Cost function 
%   0.5 ||SHx - data||_2^2 + regTikh0/2 ||x - x0||_2^2
% which corresponds to the prox of 0.5 ||SHx - y||_2^2 for alpha=1/regTikh0.
x0=zeros_(H.sizein);x0(:,:,nbSl+1)=im(:,:,nbSl+1);
b=fwd'*data/regTikh0 + x0;
cost=CostL2([],data)*fwd + regTikh0*CostL2([],x0);

% - Direct inversion (using Woodbury Formula)
tic;A=(In - fwd'*(regTikh0*Im + fwd*fwd')^(-1)*fwd);
xTikh0=A*b;t=toc;
disp(['Direct method (manual derivation of Woodbury Formula) : ',num2str(t),' s']);

% -  Direct inversion with Less FFTs (Corollary III.4 in [1])
Lamb=LinOpDiag([],mtf);
D=LinOpDiag([],sum(abs(mtf).^2,3)+regTikh0);
tic;xTikh0_bis=real(iSfft((In - Lamb'*S'*D^(-1)*S*Lamb)*(Sfft(b,3)),3));t=toc;
disp(['Direct method bis (Corollary III.4 in [1]) : ',num2str(t),' s']);
disp(['Error between the two manual computation : ',num2str(norm(xTikh0(:)-xTikh0_bis(:)))]);

% - Comparison with the prox of CostL2Composition
cc=CostL2([],data)*fwd;cc.doPrecomputation=1;
tic;xTikh0_prox=cc.applyProx(x0,1/regTikh0);t=toc;
disp(['Direct method (with prox of CostL2Composition) : ',num2str(t),' s']);
tic;xTikh0_prox=cc.applyProx(x0,1/regTikh0);t=toc;
disp(['Second prox call (effect of precomputation): ',num2str(t),' s']);
disp(['Error between prox and manual computation : ',num2str(norm(xTikh0(:)-xTikh0_prox(:)))]);

% - Comparison with CG
AA=fwd'*fwd+regTikh0*In;
CG=OptiConjGrad(AA,b*regTikh0);
CG.maxiter=100;
CG.CvOp=TestCvgStepRelative(1e-6);
CG.verbose=true;CG.OutOp.iterVerb=20;CG.OutOp.saveXopt=1;
CG.ItUpOut=1;
CG.run(x0);
disp(['With CG : ',num2str(CG.time),' s']);
disp(['Speedup x',num2str(round(CG.time/t))]);
disp('=============================================================');

% -- Display
figure;set(gca,'FontSize',14); grid; hold all;set(gca,'yscale','log');
semilogy(CG.OutOp.iternum,cellfun(@(x) cost*x,CG.OutOp.evolxopt),'LineWidth',2);
semilogy(CG.OutOp.iternum,ones(1,length(CG.OutOp.iternum))*(cost*xTikh0_prox),'LineWidth',2);
legend('CG','Direct');
figure;
for ii=1:sz(3)
    subplot(3,sz(3),ii);
    imagesc(xTikh0_prox(:,:,ii)); axis image;  colormap gray;title(['Direct Slice #',num2str(ii)]);
    subplot(3,sz(3),sz(3)+ii);
    imagesc(CG.xopt(:,:,ii)); axis image;  colormap gray;title(['ConjGrad Slice #',num2str(ii)]);
    subplot(3,sz(3),2*sz(3)+ii);
    imagesc(abs(CG.xopt(:,:,ii)-xTikh0_prox(:,:,ii))); axis image;  colorbar; colormap gray;title(['Diff Slice #',num2str(ii)]);
end

%% For Unitary Tests
global generateDataUnitTests stateTest
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valxTikh0=xTikh0;save('Util/UnitTest/Data/TestProxL2SumConv_xTikh0','valxTikh0');
        valxTikh0_bis=xTikh0_bis;save('Util/UnitTest/Data/TestProxL2SumConv_xTikh0_bis','valxTikh0_bis');
        valxTikh0_prox=xTikh0_prox;save('Util/UnitTest/Data/TestProxL2SumConv_xTikh0_prox','valxTikh0_prox');
        xopt=CG.xopt;save('Util/UnitTest/Data/TestProxL2SumConv_xopt','xopt');
    else
        load('Util/UnitTest/Data/TestProxL2SumConv_xTikh0');
        load('Util/UnitTest/Data/TestProxL2SumConv_xTikh0_bis');
        load('Util/UnitTest/Data/TestProxL2SumConv_xTikh0_prox');
        load('Util/UnitTest/Data/TestProxL2SumConv_xopt');
        stateTest=isequal(xopt,CG.xopt) &&  isequal(valxTikh0,xTikh0) &&  isequal(valxTikh0_bis,xTikh0_bis) &&  isequal(valxTikh0_prox,xTikh0_prox);
    end
end