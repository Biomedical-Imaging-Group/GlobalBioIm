%-----------------------------------------------------------
% Deconv_LS_HessSchatt_NonNeg script: Deconvolution by minimizing the 
% Least-Squares function plus the NonNegativity constraint 
% with Hessian-Schatten regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*||Hess*x||_{S_p,1}
% using 
%      - Primal-Dual Condat
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpHess, Cost, CostL2, CostNonNeg,  
% CostMixNorm1Schatt, Opti, OptiPrimalDualCondat, OptiADMM, OutpuOpti
%------------------------------------------------------------
close all;
help Deconv_LS_HessSchatt_NonNeg
%--------------------------------------------------------------
%  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%                     F. Soulez ferreol.soulez@univ-lyon1.fr
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
[im,psf,y]=GenerateData('Gaussian',20);
imdisp(im,'Input Image (GT)',1);
imdisp(y,'Convolved and noisy data',1);
sz=size(y);

% -- Convolution Operator definition
H=LinOpConv(fft2(psf));

% -- Functions definition
LS=CostL2([],y);                 % Least-Sqaures data term
F=LS*H;
F.doPrecomputation=1;
Hess=LinOpHess(sz);                  % Hessian Operator
R_1sch=CostMixNormSchatt1([sz,3],1); % Mixed Norm 1-Schatten (p=1)
R_POS=CostNonNeg(sz);                % Non-Negativity
lamb=5e-3;                           % Hyperparameter

% -- ADMM LS + ShattenHess + NonNeg
Fn={lamb*R_1sch,R_POS};
Hn={Hess,LinOpIdentity(size(im))};
rho_n=[1e-1,1e-1];
ADMM=OptiADMM(F,Fn,Hn,rho_n);
ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 2]), 'StepRelative',1e-4); 
ADMM.OutOp=OutputOptiSNR(1,im,10,[1 2]);
ADMM.ItUpOut=1;           % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;          % max number of iterations
ADMM.run(zeros(size(y)));  % run the algorithm 

% -- PrimalDual Condat LS + ShattenHess + NonNeg
Fn={lamb*R_1sch};
Hn={Hess};
PDC=OptiPrimalDualCondat(F,R_POS,Fn,Hn);
PDC.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 3]), 'StepRelative',1e-4); 
PDC.OutOp=OutputOptiSNR(1,im,40,[1 3]);
PDC.tau=1;                % set algorithm parameters
PDC.sig=1e-2;             %
PDC.rho=1.7;             %
PDC.ItUpOut=1;           % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;          % max number of iterations
PDC.run(zeros(size(y)));  % run the algorithm 


% -- Display
imdisp(ADMM.xopt,'LS+HESS+POS (ADMM)',1);
imdisp(PDC.xopt,'LS+HESS+POS (Condat)',1);
figure; plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5); grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all;plot(PDC.OutOp.iternum,PDC.OutOp.evolcost,'LineWidth',1.5); grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(PDC.OutOp.iternum,PDC.OutOp.evolsnr,'LineWidth',1.5);
legend('LS+HESS+POS (ADMM)','LS+HESS+POS (Condat)','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+HESS+POS (ADMM)','LS+HESS+POS (Condat)'});set(gca,'XTickLabelRotation',45)

%% For Unitary Tests
global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valADMM=ADMM.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_HessSchatt_NonNeg_ADMM','valADMM');
        valPDC=PDC.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_HessSchatt_NonNeg_PDC','valPDC');
    else
        load('Util/UnitTest/Data/Deconv_LS_HessSchatt_NonNeg_PDC');
        load('Util/UnitTest/Data/Deconv_LS_HessSchatt_NonNeg_ADMM');
        stateTest=isequal(valADMM,ADMM.OutOp.evolcost) &&  isequal(valPDC,PDC.OutOp.evolcost);
        message=['Max error between costs evolution : ',num2str(norm(valADMM-ADMM.OutOp.evolcost),...
            norm(valPDC-PDC.OutOp.evolcost))];
    end
end