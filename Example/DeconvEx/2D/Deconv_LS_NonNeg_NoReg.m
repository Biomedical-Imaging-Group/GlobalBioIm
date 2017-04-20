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
clear all; close all; clc;
help Deconv_Ls_NonNeg_NoReg
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
F_LS=CostL2(H,y);      % Least-Sqaures data term
R_POS=CostNonNeg();    % Non-Negativity

% -- FISTA LS + NonNeg
OutOp=OutputOpti(1,impad,40);
FBS=OptiFBS(F_LS,R_POS,OutOp);
FBS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
FBS.fista=true;   % activate fista
FBS.maxiter=200;  % max number of iterations
FBS.run(y);       % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

% -- Display
[v,n]=max(OutOp.evolsnr(:));
imdisp(OutOp.evolxopt{n}(idx,idx),'LS + NonNeg (FISTA)',1);
figure;plot(OutOp.iternum,OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');legend('LS+POS (FISTA)');title('Cost evolution');
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutOp.iternum,OutOp.evolsnr,'LineWidth',1.5);
xlabel('Iterations');ylabel('SNR (dB)');legend('LS+POS (FISTA)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[FBS.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+POS (FISTA)'});set(gca,'XTickLabelRotation',45)
