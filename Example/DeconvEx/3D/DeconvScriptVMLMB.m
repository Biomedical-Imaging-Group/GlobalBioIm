%--------------------------------------------------------------
% This script performs 3D deconvolution using VMLMB with
%    - Data-term: Least-Squares or Kullback-Leibler
%    - regul: TV 
%--------------------------------------------------------------
close all;
help DeconvScriptVMLMB
%--------------------------------------------------------------
%  Copyright (C) 2024 E. Soubies emmanuel.soubies@irit.fr
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

% -- To run on GPU (0: CPU / 1: Matlab Parrallel Computing Toolbox / 2: CudaMat) 
%    Note: when using Matlab Parrallel Computing Toolbox (useGPU(1)), resolution with HessianSchatten Norm 
%          may nor be optimal because svd computations are done with a mex function (see Cost/CostUtils/HessianSchatten) 
%          on the cpu (so lost of time in transfert CPU/GPU).
useGPU(0)

%% Parameters
lamb=1e-4;        % Hyperparameter for initial deconvolution
maxIt=30;         % Max iterations
DataTerm=1;       % 1 for LS, 2 for KL

%% fix the random seed (for reproductibility)
rng(1);

%% Input image and psf
if DataTerm==1
    [im,psf,y]=GenerateData3D('Gaussian',25);
elseif DataTerm==2
    [im,psf,y]=GenerateData3D('Poisson',5);
end
Orthoviews(im,[],'Input Image (GT)');
Orthoviews(y,[],'Convolved and noisy data');
sz=size(y);

%% Conversion CPU/GPU is necessary
psf=gpuCpuConverter(psf);
im=gpuCpuConverter(im);
y=gpuCpuConverter(y);

%% Common operators and costs
H=LinOpConv(fftn(psf));                          % Convolution Operator
L2=CostL2(H.sizeout,y);                          % L2 cost function
KL=CostKullLeib(sz,y,1e-6);                      % Kullback-Leibler divergence data term

%% Deconvolution
% -- Functions definition
Freg=CostHyperBolic([sz,3],1e-7,4);    % Smooth TV regularizer: Mixed Norm 2-1
Opreg=LinOpGrad(sz);                   % Smooth TV regularizer: Operator Gradient

% -- VMLMB
if DataTerm==1
    C = L2*H + lamb*Freg*Opreg;
else
    C = KL*H + lamb*Freg*Opreg;
end
VMLMB=OptiVMLMB(C,0.,[]);  
VMLMB.OutOp=OutputOptiSNR(1,im,round(maxIt/10));
VMLMB.CvOp=TestCvgCombine('CostRelative',1e-4, 'StepRelative',1e-4); % identical to VMLMB.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4),TestCvgStepRelative(1e-5)); 
VMLMB.ItUpOut=round(maxIt/10); 
VMLMB.maxiter=maxIt;                             % max number of iterations
VMLMB.m=3;                                     % number of memorized step in hessian approximation
VMLMB.run(y); 

%% Displays
Orthoviews(VMLMB.xopt,[],'Deconvolved image');

% -- Plot Evolution SNR, cost  and Running Time for TV-Reg-Pos methods
figure;subplot(1,3,1); grid; hold all;
plot(VMLMB.OutOp.iternum,VMLMB.OutOp.evolcost,'LineWidth',1.5);  
set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('VMLMB');title('Cost evolution');
subplot(1,3,2); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(VMLMB.OutOp.iternum,VMLMB.OutOp.evolsnr,'LineWidth',1.5);
legend('VMLMB','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,3,3);hold on; grid; title('Runing Time');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[VMLMB.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3 4]);ylabel('Time (s)');
set(gca,'xticklabels',{'VMLMB'});set(gca,'XTickLabelRotation',45);
