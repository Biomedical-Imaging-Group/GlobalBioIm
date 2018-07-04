%-----------------------------------------------------------
% Deconv_Ls_NonNeg_NoReg script: Deconvolution by minimizing 
% the Least-Squares function plus the NonNegativity constraint 
% without any regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x)
% using FISTA
%
% See LinOp, LinOpConv, Cost, CostL2, CostNonNeg, Opti, 
% OptiFBS, OutpuOpti
%------------------------------------------------------------
clear; close all; 
help Deconv_Ls_NonNeg_NoReg
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
LS=CostL2([],y);         % Least-Sqaures data term
R_POS=CostNonNeg(sz);    % Non-Negativity
F=LS*H;
F.doPrecomputation=1;

% -- FISTA LS + NonNeg
FBS=OptiFBS(F,R_POS);
FBS.OutOp=OutputOpti(1,im,20);
%FBS.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);  
FBS.ItUpOut=1;          % call OutputOpti update every ItUpOut iterations
FBS.fista=true;         % activate fista
FBS.maxiter=200;        % max number of iterations
FBS.run(zeros(size(y)));% run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 


%% -- VMLMB LS +  NonNeg
 
H.memoizeOpts.applyHtH=true;
VMLMB=OptiVMLMB(F,0.,[]);  
VMLMB.OutOp=OutputOpti(1,im,10);
%VMLMB.CvOp=TestCvgCombine('CostRelative',1e-4, 'StepRelative',1e-4); % identical to VMLMB.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4),TestCvgStepRelative(1e-5)); 
VMLMB.ItUpOut=2; 
VMLMB.maxiter=200;                             % max number of iterations
VMLMB.m=3;                                     % number of memorized step in hessian approximation
VMLMB.run(y);                                  % run the algorithm 

%% -- Display
imdisp(FBS.OutOp.evolxopt{end},'LS + NonNeg (FISTA)',1);
imdisp(VMLMB.OutOp.evolxopt{end},'LS+POS (VMLMB)',1);
figure;plot(FBS.OutOp.iternum,FBS.OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);
hold all;loglog(VMLMB.OutOp.iternum,VMLMB.OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);
xlabel('Iterations');ylabel('Cost');legend('LS+POS (FISTA)','LS+POS (VMLMB)');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(FBS.OutOp.iternum,FBS.OutOp.evolsnr,'LineWidth',1.5);
semilogy(VMLMB.OutOp.iternum,VMLMB.OutOp.evolsnr,'LineWidth',1.5);
xlabel('Iterations');ylabel('SNR (dB)');legend('LS+POS (FISTA)','LS+POS (VMLMB)','Location','southeast');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[FBS.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[VMLMB.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+POS (FISTA)','LS+POS (VMLMB)'});set(gca,'XTickLabelRotation',45)
