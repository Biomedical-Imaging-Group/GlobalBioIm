%-----------------------------------------------------------
% Deconv_LS_TV_NonNeg script: Deconvolution by minimizing the 
% Least-Squares function plus the NonNegativity constraint 
% with TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% using 
%      - Primal-Dual Condat
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostL2, CostNonNeg,  
% CostMixNorm12, Opti, OptiPrimalDualCondat, OptiADMM, OutpuOpti
%------------------------------------------------------------
clear all; close all; clc;
help Deconv_LS_TV_NonNeg
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
H.memoizeOpts.applyHtH=true;

% -- Generate data
load('data');    % load data (variable y)
imdisp(y(idx,idx),'Convolved and noisy data',1);
sz=size(y);

% -- Functions definition
LS=CostL2([],y);                 % Least-Squares data term
F=LS*H;
F.doPrecomputation=1;
R_N12=CostMixNorm21([sz,2],3);   % Mixed Norm 2-1
G=LinOpGrad(sz);                 % Operator Gradient
R_POS=CostNonNeg(sz);            % Non-Negativity
lamb=7e-4;                       % Hyperparameter

% -- ADMM LS + TV + NonNeg
Fn={lamb*R_N12,R_POS};
Hn={G,LinOpIdentity(sz)};
rho_n=[1e-1,1e-1];
OutADMM=MyOutputOpti(1,impad,40);
ADMM=OptiADMM(F,Fn,Hn,rho_n,[],OutADMM);
ADMM.ItUpOut=10;        % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;       % max number of iterations
ADMM.run(y);            % run the algorithm 

% -- PrimalDual Condat LS + TV + NonNeg
Fn={lamb*R_N12};
Hn={G};
OutPDC=MyOutputOpti(1,impad,40);
PDC=OptiPrimalDualCondat(F,R_POS,Fn,Hn,OutPDC);
PDC.tau=1;                                   % set algorithm parameters
PDC.sig=(1/PDC.tau-F.lip/2)/G.norm^2*0.9; %
PDC.rho=1.95;                                %
PDC.ItUpOut=10;                              % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;                             % max number of iterations
PDC.run(y);                                  % run the algorithm 


%% -- VMLMB LS + hyperbolicTV + NonNeg
hyperB = CostHyperBolic(G.sizeout,   1e-7,  3)*G;
hyperB.doPrecomputation=1;
hyperB.memoizeOpts.apply=true;
C = F+ lamb*hyperB;
OutVMLMB=MyOutputOpti(1,impad,40);
VMLMB=OptiVMLMB(C,0., [],OutVMLMB);                            %
VMLMB.ItUpOut=10;   
VMLMB.maxiter=200;                             % max number of iterations
VMLMB.run(y);                                  % run the algorithm 


% -- Display
imdisp(OutADMM.evolxopt{end}(idx,idx),'LS+TV+POS (ADMM)',1);
imdisp(OutPDC.evolxopt{end}(idx,idx),'LS+TV+POS (Condat)',1);
imdisp(OutVMLMB.evolxopt{end}(idx,idx),'LS+TV+POS (VMLMB)',1);
figure; plot(OutADMM.iternum,OutADMM.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all;plot(OutPDC.iternum,OutPDC.evolcost,'LineWidth',1.5);set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
plot(OutVMLMB.iternum,OutVMLMB.evolcost,'LineWidth',1.5);set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat','VMLMB');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutADMM.iternum,OutADMM.evolsnr,'LineWidth',1.5); 
semilogy(OutPDC.iternum,OutPDC.evolsnr,'LineWidth',1.5);
semilogy(OutVMLMB.iternum,OutVMLMB.evolsnr,'LineWidth',1.5);
legend('LS+TV+POS (ADMM)','LS+TV+POS (Condat)','LS+TV+POS (VMLMB)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[VMLMB.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+TV+POS (ADMM)','LS+TV+POS (Condat)','LS+TV+POS (VMLMB)'});set(gca,'XTickLabelRotation',45)

