%-----------------------------------------------------------
% Deconv_KL_TV_NonNeg script: Deconvolution by minimizing
% the Kullback-Leibler divergence plus the NonNegativity constraint 
% with TV  regularizer:
%     \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + Tv(x)
% using:
%    - Primal-Dual Condat
%    - ADMM 
%    - RichardsonLucy-TV
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostKullLeib, CostNonNeg, 
% CostMixNorm12, Opti, OptiADMM, OptiRichLucy, OutpuOpti
% OptiPrimalDualCondat.
%------------------------------------------------------------
clear all; close all; clc;
help Deconv_KL_TV_NonNeg
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
load('StarLikeSample');    % Load image (variable im)
load('psf');               % Load psf (variable psf)
imdisp(im,'Input Image',1);

% -- Image padding
impad=zeros(512); idx=129:384;
impad(idx,idx)=im;

% -- Convolution Operator definition
H=LinOpConv(fft2(psf));

% -- Generate data
load('data');    % load data (variable y)
imdisp(y(idx,idx),'Convolved and noisy data',1);
sz=size(y);

% -- Functions definition
F=CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
R_POS=CostNonNeg(sz);          % Non-Negativity
R_N12=CostMixNorm21([sz,2],3); % Mixed Norm 2-1
G=LinOpGrad(sz);               % Operator Gradient
lamb=1e-2;                     % Hyperparameter  

% -- ADMM KL + TV + NonNeg
Fn={CostKullLeib([],y,1e-6),lamb*R_N12,R_POS};
Hn={H,G,LinOpDiag(sz)};
rho_n=[1e-2,1e-2,1e-2];
ADMM=OptiADMM([],Fn,Hn,rho_n);
ADMM.OutOp=MyOutputOpti(1,impad,40);
ADMM.ItUpOut=10;                                  % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;                                 % max number of iterations
ADMM.run(y);                                      % run the algorithm

% -- PrimalDual Condat KL + TV + NonNeg
Fn={lamb*R_N12};
Hn={G};
PDC=OptiPrimalDualCondat(F*H,R_POS,Fn,Hn);
PDC.OutOp=MyOutputOpti(1,impad,40);
PDC.tau=1e-2;          % set algorithm parameters
PDC.sig=1;             %
PDC.rho=1.95;          %
PDC.ItUpOut=10;        % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;       % max number of iterations
PDC.run(y);            % run the algorithm 

% -- Richardson-Lucy-TV  KL + TV + NonNeg (implicit)
RLTV=OptiRichLucy(F*H,1,lamb);
RLTV.OutOp=MyOutputOpti(1,impad,40);
RLTV.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
RLTV.maxiter=200;  % max number of iterations
RLTV.run(y);       % run the algorithm 

% -- Display
imdisp(ADMM.OutOp.evolxopt{end}(idx,idx),'KL+TV+POS (ADMM)',1);
imdisp(PDC.OutOp.evolxopt{end}(idx,idx),'KL+TV+POS (Condat)',1);
imdisp(RLTV.OutOp.evolxopt{end}(idx,idx),'KL+TV+POS (RL-TV)',1);

figure;plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5); grid;  
hold all;plot(PDC.OutOp.iternum,PDC.OutOp.evolcost,'LineWidth',1.5); 
plot(RLTV.OutOp.iternum,RLTV.OutOp.evolcost,'LineWidth',1.5); 
set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat','RL-TV');title('Cost evolution');

% -- Plot Evolution SNR and Running Time for Hess-Reg-Pos methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5);
semilogy(PDC.OutOp.iternum,PDC.OutOp.evolsnr,'LineWidth',1.5);
semilogy(RLTV.OutOp.iternum,RLTV.OutOp.evolsnr,'LineWidth',1.5);
legend('KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[RLTV.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)'});set(gca,'XTickLabelRotation',45)

