%-----------------------------------------------------------
% Deconv_LS_TV_NonNeg script: Deconvolution by minimizing the 
% Least-Squares function plus the NonNegativity constraint 
% with TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% using 
%      - Primal-Dual Condat
%      - ADMM 
%      - VMLMB
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostL2, CostNonNeg,  
% CostMixNorm12, Opti, OptiPrimalDualCondat, OptiADMM, OutpuOpti
% CostHyperBolic, OptiVMLMB
%------------------------------------------------------------
close all;
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

% -- To run on GPU (0: CPU / 1: Matlab Parrallel Computing Toolbox / 2: CudaMat) 
useGPU(0)

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
[im,psf,y]=GenerateData('Gaussian',20);
imdisp(im,'Input Image (GT)',1);
imdisp(y,'Convolved and noisy data',1);
sz=size(y);

% -- Convolution Operator definition
H=LinOpConv(fft2(psf));
H.memoizeOpts.applyHtH=true;

% -- Functions definition
LS=CostL2([],y);                 % Least-Squares data term
F=LS*H;
F.doPrecomputation=1;
F.memoizeOpts.apply=true;
R_N12=CostMixNorm21([sz,2],3);   % Mixed Norm 2-1
G=LinOpGrad(sz);                 % Operator Gradient
R_POS=CostNonNeg(sz);            % Non-Negativity
lamb=1e-3;                       % Hyperparameter

%% -- ADMM LS + TV + NonNeg
Fn={lamb*R_N12,R_POS};
Hn={G,LinOpIdentity(sz)};
rho_n=[1e-1,1e-1];
ADMM=OptiADMM(F,Fn,Hn,rho_n);
ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
% STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4,TestCvgADMM());  
ADMM.ItUpOut=2;             % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;           % max number of iterations
ADMM.run(y);   % run the algorithm 

%% -- PrimalDual Condat LS + TV + NonNeg
Fn={lamb*R_N12};
Hn={G};
PDC=OptiPrimalDualCondat(F,R_POS,Fn,Hn);
PDC.OutOp=OutputOptiSNR(1,im,40,[1 3]);
% STOP when the sum successives C = F*x + Fn*Hn*x is lower than 1e-4 or when the distance between two successive step is lower than 1e-5
PDC.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 3])  , 'StepRelative',1e-4); 
PDC.tau=1;                                   % set algorithm parameters
PDC.sig=(1/PDC.tau-F.lip/2)/G.norm^2*0.9;    %
PDC.rho=1.95;                                %
PDC.ItUpOut=2;                               % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;                             % max number of iterations
PDC.run(y);                     % run the algorithm 


%% -- VMLMB LS + hyperbolicTV + NonNeg
hyperB = CostHyperBolic(G.sizeout,   1e-7,  3)*G;
C = F+ lamb*hyperB; 
C.memoizeOpts.apply=true;
VMLMB=OptiVMLMB(C,0.,[]);  
VMLMB.OutOp=OutputOptiSNR(1,im,10);
VMLMB.CvOp=TestCvgCombine('CostRelative',1e-4, 'StepRelative',1e-4); % identical to VMLMB.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4),TestCvgStepRelative(1e-5)); 
VMLMB.ItUpOut=2; 
VMLMB.maxiter=200;                             % max number of iterations
VMLMB.m=3;                                     % number of memorized step in hessian approximation
VMLMB.run(y);                                  % run the algorithm 


%% -- Display
imdisp(ADMM.xopt,'LS+TV+POS (ADMM)',1);
imdisp(PDC.xopt,'LS+TV+POS (Condat)',1);
imdisp(VMLMB.xopt,'LS+TV+POS (VMLMB)',1);
figure; plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);
hold all;plot(PDC.OutOp.iternum,PDC.OutOp.evolcost,'LineWidth',1.5);set(gca,'FontSize',12);
plot(VMLMB.OutOp.iternum,VMLMB.OutOp.evolcost,'LineWidth',1.5);set(gca,'FontSize',12);
xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat','VMLMB');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(PDC.OutOp.iternum,PDC.OutOp.evolsnr,'LineWidth',1.5);
semilogy(VMLMB.OutOp.iternum,VMLMB.OutOp.evolsnr,'LineWidth',1.5);
legend('LS+TV+POS (ADMM)','LS+TV+POS (Condat)','LS+TV+POS (VMLMB)','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time ');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
%text( [1],[ADMM.time],[num2str(ADMM.OutOp.iternum(end)), ' iter'],'vert','bottom','horiz','center','FontSize',12);
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[VMLMB.time],'FaceColor',orderCol(3,:),'EdgeColor','k');

set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+TV+POS (ADMM)','LS+TV+POS (Condat)','LS+TV+POS (VMLMB)'});set(gca,'XTickLabelRotation',45)

%% For Unitary Tests
global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valADMM=ADMM.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_ADMM','valADMM');
        valPDC=PDC.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_PDC','valPDC');
        valVMLMB=VMLMB.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_VMLMB','valVMLMB');
    else
        load('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_VMLMB');
        load('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_PDC');
        load('Util/UnitTest/Data/Deconv_LS_TV_NonNeg_ADMM');
        stateTest=isequal(valADMM,ADMM.OutOp.evolcost) &&  isequal(valPDC,PDC.OutOp.evolcost) &&  isequal(valVMLMB,VMLMB.OutOp.evolcost);
        message=['Max error between costs evolution : ',num2str(max([norm(valPDC-PDC.OutOp.evolcost),...
            norm(valADMM-ADMM.OutOp.evolcost),norm(valVMLMB-VMLMB.OutOp.evolcost)]))];
    end
end