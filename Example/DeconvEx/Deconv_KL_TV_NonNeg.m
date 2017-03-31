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
clear all; close all; clc;%warning('off');
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
H=LinOpConv(psf);

% -- Generate data
load('data');    % load data (variable y)
imdisp(y(idx,idx),'Convolved and noisy data',1);

% -- Functions definition
F_KL=CostKullLeib(H,y,1e-6);% Kullback-Leibler divergence data term
R_POS=CostNonNeg();         % Non-Negativity
R_N12=CostMixNorm12([3]);   % Mixed Norm 2-1
G=LinOpGrad(size(y));       % Operator Gradient
lamb=1e-2;                  % Hyperparameter  

% -- ADMM KL + TV + NonNeg
Fn={CostKullLeib([],y,1e-6),MultScalarCost(R_N12,lamb),R_POS};
Hn={H,G,LinOpIdentity(size(impad))};
rho_n=[1e-2,1e-2,1e-2];
lap=zeros(size(impad)); lap(1,1)=4; lap(1,2)=-1;lap(2,1)=-1; lap(1,end)=-1;lap(end,1)=-1; Flap=fft2(lap);
solver = @(z,rho,x) real(ifft2((rho(1)*conj(H.mtf).*fft2(z{1}) + fft2(rho(2)*G'*z{2} + rho(3)*z{3}) )./(rho(1)*abs(H.mtf).^2 + rho(2)*Flap + rho(3))));  % solver to solve the x update
OutADMM=OutputOpti(1,impad,40);
ADMM=OptiADMM([],[],Fn,Hn,rho_n,solver,OutADMM);
ADMM.ItUpOut=10;                                  % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;                                 % max number of iterations
ADMM.run(y);                                      % run the algorithm

% -- PrimalDual Condat KL + TV + NonNeg
Fn={MultScalarCost(R_N12,lamb)};
Hn={G};
OutPDC=OutputOpti(1,impad,40);
PDC=OptiPrimalDualCondat(F_KL,R_POS,Fn,Hn,OutPDC);
PDC.tau=1e-2;          % set algorithm parameters
PDC.sig=1;             %
PDC.rho=1.95;          %
PDC.ItUpOut=10;        % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;       % max number of iterations
PDC.run(y);            % run the algorithm 

% -- Richardson-Lucy-TV  KL + TV + NonNeg (implicit)
OutRLTV=OutputOpti(1,impad,40);
RLTV=OptiRichLucy(F_KL,1,lamb,OutRLTV);
RLTV.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
RLTV.maxiter=200;  % max number of iterations
RLTV.run(y);       % run the algorithm 

% -- Display
imdisp(OutADMM.evolxopt{end}(idx,idx),'KL+TV+POS (ADMM)',1);
imdisp(OutPDC.evolxopt{end}(idx,idx),'KL+TV+POS (Condat)',1);
imdisp(OutRLTV.evolxopt{end}(idx,idx),'KL+TV+POS (RL-TV)',1);

figure;plot(OutADMM.iternum,OutADMM.evolcost,'LineWidth',1.5); grid;  
hold all;plot(OutPDC.iternum,OutPDC.evolcost,'LineWidth',1.5); 
plot(OutRLTV.iternum,OutRLTV.evolcost,'LineWidth',1.5); 
set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat','RL-TV');title('Cost evolution');

% -- Plot Evolution SNR and Running Time for Hess-Reg-Pos methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutADMM.iternum,OutADMM.evolsnr,'LineWidth',1.5);
semilogy(OutPDC.iternum,OutPDC.evolsnr,'LineWidth',1.5);
semilogy(OutRLTV.iternum,OutRLTV.evolsnr,'LineWidth',1.5);
legend('KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[RLTV.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)'});set(gca,'XTickLabelRotation',45)

