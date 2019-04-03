%-----------------------------------------------------------
% Deconv_KL_HessSchatt_NonNeg script: Deconvolution by minimizing the 
% Kullback-Leibler divergence plus the NonNegativity constraint 
% with Hessian-Schatten regularizer:
%    \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + lamb*||Hess*x||_{1,S_p}
% using 
%      - Primal-Dual Condat
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpHess, Cost, CostKullLeib, CostNonNeg,  
% CostMixNorm1Schatt, Opti, OptiPrimalDualCondat, OptiADMM, OutpuOpti
%------------------------------------------------------------
close all; 
help Deconv_KL_HessSchatt_NonNeg
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

% -- Functions definition
F=CostKullLeib([],y,1e-6);           % Kullback-Leibler divergence data term
Hess=LinOpHess(sz);                  % Hessian Operator
R_1sch=CostMixNormSchatt1([sz,3],1); % Mixed Norm 1-Schatten (p=1)
R_POS=CostNonNeg(sz);                % Non-Negativity
lamb=5e-3;                           % Hyperparameter

% -- ADMM KL + ShattenHess + NonNeg
Fn={CostKullLeib([],y,1e-6),lamb*R_1sch,R_POS};
Hn={H,Hess,LinOpDiag(sz)};
rho_n=[1e-2,1e-2,1e-2];
ADMM=OptiADMM([],Fn,Hn,rho_n,[]);
ADMM.OutOp=OutputOptiSNR(1,im,40,[1 2]);
ADMM.ItUpOut=5;        % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;       % max number of iterations
ADMM.run(y);            % run the algorithm 


% -- PrimalDual Condat KL + ShattenHess + NonNeg
Fn={lamb*R_1sch,F};
Hn={Hess,H};
PDC=OptiPrimalDualCondat([],R_POS,Fn,Hn);
PDC.OutOp=OutputOptiSNR(1,im,40,[2 3]);
PDC.tau=100;          % set algorithm parameters
PDC.sig=1e-2;            %
PDC.rho=1.2;          %
PDC.ItUpOut=5;        % call OutputOpti update every ItUpOut iterations
PDC.maxiter=200;       % max number of iterations
PDC.run(y);            % run the algorithm 

% -- Display
imdisp(ADMM.xopt,'KL+HESS+POS (ADMM)',1);
imdisp(PDC.xopt,'KL+HESS+POS (Condat)',1);
figure; plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5);
hold all; plot(PDC.OutOp.iternum,PDC.OutOp.evolcost,'LineWidth',1.5);
 grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('ADMM','Condat');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(PDC.OutOp.iternum,PDC.OutOp.evolsnr,'LineWidth',1.5);
legend('KL+HESS+POS (ADMM)','KL+HESS+POS (Condat)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[ADMM.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[PDC.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'KL+HESS+POS (ADMM)','KL+HESS+POS (Condat)'});set(gca,'XTickLabelRotation',45)

%% For Unitary Tests
global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valADMM=ADMM.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_HessSchatt_NonNeg_ADMM','valADMM');
        valPDC=PDC.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_KL_HessSchatt_NonNeg_PDC','valPDC');
    else
        load('Util/UnitTest/Data/Deconv_KL_HessSchatt_NonNeg_PDC');
        load('Util/UnitTest/Data/Deconv_KL_HessSchatt_NonNeg_ADMM');
        stateTest=isequal(valADMM,ADMM.OutOp.evolcost) &&  isequal(valPDC,PDC.OutOp.evolcost);
        message=['Max error between costs evolution : ',num2str(max(norm(valADMM-ADMM.OutOp.evolcost),...
            norm(valPDC-PDC.OutOp.evolcost)))];
    end
end