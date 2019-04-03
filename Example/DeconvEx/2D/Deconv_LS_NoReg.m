%-----------------------------------------------------------
% Deconv_LS_NoReg script: Deconvolution by minimizing the
% Least-Squares function without any regularizer:
%     0.5 ||Hx - y||^2
% using Gradient Descent, VMLMB and Conjugate Gradient
%
% See LinOp, LinOpConv, Cost, CostL2, CostL2Composition, Opti,
% OptiGradDsct, OptiConjGrad, OptiVMLMB, OutpuOpti
%------------------------------------------------------------
close all;
help Deconv_LS_NoReg
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
% If true, applyHtH method will save its result. Hence for two consecutive HtH*x with the
% same x, no computation is done for the second call
H.memoizeOpts.applyHtH=1;  
% -- Function definition
LS=CostL2([],y);  % Least-Sqaures data term
F=LS*H;
% For the CostL2, the precomputation save Hty=H'*y and then the gradient is
% computed using H.HtH(x) - Hty and the apply is computed using the HtH
% method
F.doPrecomputation=1;

% Note: For this example, the activation of the memoize for the applyHtH
% method together with the activation of the doPrecomputation flag leads to
% the faster execution time. This is because, both the gradient and the
% apply make use of a fast HtH.

% -- Gradient Descent LS
GD=OptiGradDsct(F);
GD.OutOp=OutputOptiSNR(1,im,40);
GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=200;         % max number of iterations
GD.run(y); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

%% -- VMLMB LS 
VMLMB=OptiVMLMB(F,[],[]);  
VMLMB.OutOp=OutputOptiSNR(1,im,40);
VMLMB.ItUpOut=2; 
VMLMB.maxiter=200;                             % max number of iterations
VMLMB.m=2;                                     % number of memorized step in hessian approximation (one step is enough for quadratic function)
VMLMB.run(y);                                  % run the algorithm 


%% -- ConjGrad LS 
A = H.makeHtH();
b = H'*y;
CG=OptiConjGrad(A,b);  
CG.OutOp=OutputOptiConjGrad(1,dot(y(:),y(:)),im,40);  
CG.ItUpOut=2; 
CG.maxiter=200;                             % max number of iterations
CG.run(y);                                  % run the algorithm 


%% -- Display
imdisp(GD.xopt,'LS (GD)',1);
imdisp(VMLMB.xopt,'LS (VMLMB)',1);
imdisp(CG.xopt,'LS (CG)',1);
figure;plot(GD.OutOp.iternum,GD.OutOp.evolcost,'LineWidth',1.5);
hold all; plot(VMLMB.OutOp.iternum,VMLMB.OutOp.evolcost,'LineWidth',1.5);
plot(CG.OutOp.iternum,CG.OutOp.evolcost,'LineWidth',1.5);
set(gca,'FontSize',12);
legend('GD','VMLMB');
grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');title('Cost evolution');
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(GD.OutOp.iternum,GD.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(VMLMB.OutOp.iternum,VMLMB.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(CG.OutOp.iternum,CG.OutOp.evolsnr,'LineWidth',1.5); 
legend('GD','VMLM','CG','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[GD.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[VMLMB.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(3,[CG.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'GD','VMLM','CG'});set(gca,'XTickLabelRotation',45)

%% For Unitary Tests
global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    if generateDataUnitTests
        valGD=GD.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_NoReg_GD','valGD');
        valVMLMB=VMLMB.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_NoReg_VMLMB','valVMLMB');
        valCG=CG.OutOp.evolcost;save('Util/UnitTest/Data/Deconv_LS_NoReg_CG','valCG');
    else
        if (exist('Util/UnitTest/Data/Deconv_LS_NoReg_GD.mat','file')==2) && ... 
                (exist('Util/UnitTest/Data/Deconv_LS_NoReg_VMLMB.mat','file')==2) && ...
                (exist('Util/UnitTest/Data/Deconv_LS_NoReg_CG.mat','file')==2)

        load('Util/UnitTest/Data/Deconv_LS_NoReg_GD');
        load('Util/UnitTest/Data/Deconv_LS_NoReg_VMLMB');
        load('Util/UnitTest/Data/Deconv_LS_NoReg_CG');
        stateTest=isequal(valGD,GD.OutOp.evolcost) &&  isequal(valVMLMB,VMLMB.OutOp.evolcost) &&  isequal(valCG,CG.OutOp.evolcost);
        message=['Max error between costs evolution : ',num2str(max([norm(valGD-GD.OutOp.evolcost),...
            norm(valVMLMB-VMLMB.OutOp.evolcost),norm(valCG-CG.OutOp.evolcost)]))];
        else
            error('you must generate data before setting : set generateDataUnitTests=1 in UnitScript.m');
        end
    end
end