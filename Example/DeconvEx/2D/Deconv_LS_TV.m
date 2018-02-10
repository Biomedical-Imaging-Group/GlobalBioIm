%-----------------------------------------------------------
% Deconv_LS_TV script: Deconvolution by minimizing the 
% Least-Squares function plus the TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% using 
%      - Chambolle-Pock
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostL2,   
% CostMixNorm12, Opti, OptiChambPock, OptiADMM, OutpuOpti
%------------------------------------------------------------
clear; close all;
help Deconv_LS_TV
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

% -- To run on GPU (0: CPU / 1: Matlab Parrallel Computing Toolbox / 2: CudaMat) 
useGPU(1)

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
load('StarLikeSample');    % Load image (variable im)
load('psf');               % Load psf (variable psf)
psf=gpuCpuConverter(psf);
im=gpuCpuConverter(im);
imdisp(im,'Input Image',1);

% -- Image padding
impad=zeros_(512); idx=129:384;
impad(idx,idx)=im;

% -- Convolution Operator definition
H=LinOpConv(fft2(psf));

% -- Generate data
load('data');    % load data (variable y)
y=gpuCpuConverter(y);
imdisp(y(idx,idx),'Convolved and noisy data',1);
sz=size(y);

% -- Functions definition
LS=CostL2([],y);            % Least-Squares data term
F=LS*H;
F.doPrecomputation=1;
R_N12=CostMixNorm21([sz,2],3);  % Mixed Norm 2-1
G=LinOpGrad(size(y));       % Operator Gradient
lamb=1e-3;                  % Hyperparameter

% -- Chambolle-Pock  LS + TV
OutCP=OutputOpti(1,impad,40);
CP=OptiChambPock(lamb*R_N12,G,F,OutCP);
CP.tau=15;                            % algorithm parameters
CP.sig=1/(CP.tau*G.norm^2)*0.99;      %
CP.ItUpOut=10;                        % call OutputOpti update every ItUpOut iterations
CP.maxiter=200;                       % max number of iterations
CP.run(y);                            % run the algorithm 

% -- ADMM LS + TV
Fn={lamb*R_N12};
Hn={G};rho_n=[1e-1];
% Here no solver needed in ADMM since the operator H'*H + alpha*G'*G is invertible
OutADMM=OutputOpti(1,impad,40);
ADMM=OptiADMM(F,Fn,Hn,rho_n,[],OutADMM);
ADMM.ItUpOut=10;       % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;      % max number of iterations
ADMM.run(y);           % run the algorithm 

% -- Display
imdisp(CP.xopt(idx,idx),'LS + TV (CP)',1);
imdisp(ADMM.xopt(idx,idx),'LS + TV (ADMM)',1);
figure;plot(OutCP.iternum,OutCP.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all;plot(OutADMM.iternum,OutADMM.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutCP.iternum,OutCP.evolsnr,'LineWidth',1.5); 
semilogy(OutADMM.iternum,OutADMM.evolsnr,'LineWidth',1.5);
legend('LS+TV (CP)','LS+TV (ADMM)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[CP.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[ADMM.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+TV (CP)','LS+TV (ADMM)'});set(gca,'XTickLabelRotation',45)  
