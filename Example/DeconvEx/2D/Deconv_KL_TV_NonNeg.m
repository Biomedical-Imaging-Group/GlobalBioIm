%-----------------------------------------------------------
% Deconv_KL_TV_NonNeg script: Deconvolution by minimizing
% the Kullback-Leibler divergence plus the NonNegativity constraint 
% with TV  regularizer:
%     \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + Tv(x)
% using:
%    - Primal-Dual Condat
%    - ADMM 
%    - RichardsonLucy-TV
%    - VMLMB
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostKullLeib, CostNonNeg, 
% CostMixNorm12, Opti, OptiADMM, OptiRichLucy, OutpuOpti
% OptiPrimalDualCondat, OptiVMLMB.
%------------------------------------------------------------
close all; 
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
[im,psf,y]=GenerateData('Poisson',100);
imdisp(im,'Input Image (GT)',1);
imdisp(y,'Convolved and noisy data',1);
sz=size(y);

% -- Convolution Operator definition
H=LinOpConv(fft2(psf));
H.memoizeOpts.applyHtH = true;

% -- Functions definition
F=CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
R_POS=CostNonNeg(sz);          % Non-Negativity
R_N12=CostMixNorm21([sz,2],3); % Mixed Norm 2-1
G=LinOpGrad(sz);               % Operator Gradient
lamb=5e-3;                     % Hyperparameter  

% -- ADMM KL + TV + NonNeg
Fn={CostKullLeib([],y,1e-6),lamb*R_N12,R_POS};
Hn={H,G,LinOpDiag(sz)};
rho_n=[1e-3,1e-3,1e-3];
ADMM=OptiADMM([],Fn,Hn,rho_n);
ADMM.OutOp=OutputOptiSNR(1,im,40,[1 2]);
ADMM.ItUpOut=2;                                  % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;                                 % max number of iterations
ADMM.run(y);                                      % run the algorithm

% -- PrimalDual Condat KL + TV + NonNeg
Fn={lamb*R_N12,F};
Hn={G,H};
PDC=OptiPrimalDualCondat([],R_POS,Fn,Hn);
PDC.OutOp=OutputOptiSNR(1,im,40,[2 3]);
PDC.tau=100;          % set algorithm parameters
PDC.sig=1e-2;         %
PDC.rho=1.2;          %
PDC.ItUpOut=2;        % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;      % max number of iterations
PDC.run(y);           % run the algorithm 

% -- Richardson-Lucy-TV  KL + TV + NonNeg (implicit)
RLTV=OptiRichLucy(F*H,1,lamb);
RLTV.OutOp=OutputOptiSNR(1,im,40);
RLTV.ItUpOut=2;   % call OutputOpti update every ItUpOut iterations
RLTV.maxiter=200;  % max number of iterations
RLTV.run(y);       % run the algorithm 

%% -- VMLMB KL + hyperbolicTV + NonNeg
hyperB = CostHyperBolic(G.sizeout,   1e-7,  3)*G;
C = F*H+ lamb*hyperB; 
C.memoizeOpts.apply=true;
VMLMB=OptiVMLMB(C,0.,[]);  
VMLMB.OutOp=OutputOptiSNR(1,im,40);
VMLMB.ItUpOut=2; 
VMLMB.maxiter=200;                             % max number of iterations
VMLMB.m=3;                                     % number of memorized step in hessian approximation
VMLMB.run(y);                                  % run the algorithm 



% -- Display
imdisp(ADMM.xopt,'KL+TV+POS (ADMM)',1);
imdisp(PDC.xopt,'KL+TV+POS (Condat)',1);
imdisp(RLTV.xopt,'KL+TV+POS (RL-TV)',1);
imdisp(VMLMB.xopt,'VMLMB+TV+POS (RL-TV)',1);

figure;plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5); grid;  
hold all;plot(PDC.OutOp.iternum,PDC.OutOp.evolcost,'LineWidth',1.5); 
plot(RLTV.OutOp.iternum,RLTV.OutOp.evolcost,'LineWidth',1.5);  
plot(VMLMB.OutOp.iternum,VMLMB.OutOp.evolcost,'LineWidth',1.5); 
set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat','RL-TV','VMLMB');title('Cost evolution');

% -- Plot Evolution SNR and Running Time for TV-Reg-Pos methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5);
semilogy(PDC.OutOp.iternum,PDC.OutOp.evolsnr,'LineWidth',1.5);
semilogy(RLTV.OutOp.iternum,RLTV.OutOp.evolsnr,'LineWidth',1.5);
semilogy(VMLMB.OutOp.iternum,VMLMB.OutOp.evolsnr,'LineWidth',1.5);
legend('KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)','KL+TV+POS (VMLMB)','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[RLTV.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
bar(4,[VMLMB.time],'FaceColor',orderCol(4,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3 4]);ylabel('Time (s)');
set(gca,'xticklabels',{'KL+TV+POS (ADMM)','KL+TV+POS (Condat)','KL+TV+POS (RL-TV)','KL+TV+POS (VMLMB)'});set(gca,'XTickLabelRotation',45)

%% For Unitary Tests
global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valADMM=ADMM.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_ADMM','valADMM');
        valPDC=PDC.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_PDC','valPDC');
        valVMLMB=VMLMB.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_VMLMB','valVMLMB');
        valRLTV=RLTV.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_RLTV','valRLTV');
    else
        load('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_VMLMB');
        load('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_PDC');
        load('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_ADMM');
        load('Util/UnitTest/Data/Deconv_KL_TV_NonNeg_RLTV');
        stateTest=isequal(valADMM,ADMM.OutOp.evolcost) &&  isequal(valPDC,PDC.OutOp.evolcost) &&  isequal(valVMLMB,VMLMB.OutOp.evolcost) &&  isequal(valRLTV,RLTV.OutOp.evolcost);
        message=['Max error between costs evolution : ',num2str(max([norm(valGD-GD.OutOp.evolcost),...
            norm(valADMM-ADMM.OutOp.evolcost),norm(valVMLMB-VMLMB.OutOp.evolcost),norm(valRLTV-RLTV.OutOp.evolcost)]))];
    end
end