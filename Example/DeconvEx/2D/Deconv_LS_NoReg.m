%-----------------------------------------------------------
% Deconv_LS_NoReg script: Deconvolution by minimizing the
% Least-Squares function without any regularizer:
%     0.5 ||Hx - y||^2
% using Gradient Descent
%
% See LinOp, LinOpConv, Cost, CostL2, CostL2Composition, Opti,
% OptiGradDsct, OutpuOpti
%------------------------------------------------------------
clear; close all;
help Deconv_LS_NoReg
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

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
load('GT');                % Load ground truth (variable im)
load('psf');               % Load psf (variable psf)
imdisp(im,'Input Image (GT)',1);


% -- Convolution Operator definition
H=LinOpConv(fft2(psf));
% If true, applyHtH method will save its result. Hence for two consecutive HtH*x with the
% same x, no computation is done for the second call
H.memoizeOpts.applyHtH=1;  

% -- Generate data
load('data');    % load data (variable y)
imdisp(y,'Convolved and noisy data',1);

% -- Function definition
LS=CostL2([],y);  % Least-Sqaures data term
F=LS*H;
% For the CostL2, the precomputation save Hty=H'*y and then the gradient is
% computed using H.HtH(x) - Hty and the apply is computed using the HtH
% method
F.doPrecomputation=1;

% Note: For this example, the activation of the memoize for the applyHtH
% method together with the activation of the doPrecomputation flag leads to
% the faster execution time. This is because, both the gradient and the
% apply make use of a fast HtH.

% -- Gradient Descent LS
GD=OptiGradDsct(F);
GD.OutOp=OutputOpti(1,im,40);
GD.ItUpOut=1;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=200;         % max number of iterations
GD.run(zeros(size(y))); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

% -- Display
imdisp(GD.OutOp.evolxopt{end},'LS (GD)',1);
figure;plot(GD.OutOp.iternum,GD.OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');title('Cost evolution');
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(GD.OutOp.iternum,GD.OutOp.evolsnr,'LineWidth',1.5); 
legend('LS (GD)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[GD.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS (GD)'});set(gca,'XTickLabelRotation',45)
