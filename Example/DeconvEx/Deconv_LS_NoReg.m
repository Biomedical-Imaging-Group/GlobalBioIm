%-----------------------------------------------------------
% Deconv_LS_NoReg script: Deconvolution by minimizing the
% Least-Squares function without any regularizer:
%     0.5 ||Hx - y||^2
% using Gradient Descent
%
% See LinOp, LinOpConv, Cost, CostL2, Opti, OptiGradDsct
% OutpuOpti
%------------------------------------------------------------
clear all; close all; clc;warning('off');
help Deconv_LS_NoReg
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

% -- Function definition
F_LS=CostL2(H,y);  % Least-Sqaures data term

% -- Gradient Descent LS
OutOp=OutputOpti(1,impad,40);
GD=OptiGradDsct(F_LS,OutOp);
GD.ItUpOut=10;     % call OutputOpti update every ItUpOut iterations
GD.maxiter=200;    % max number of iterations
GD.run(y);         % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

% -- Display
[v,n]=max(OutOp.evolsnr(:));
imdisp(OutOp.evolxopt{n}(idx,idx),'LS (GD)',1);
figure;plot(OutOp.iternum,OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');title('Cost evolution');
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutOp.iternum,OutOp.evolsnr,'LineWidth',1.5); 
legend('LS (GD)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[GD.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
set(gca,'xtick',[1]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS (GD)'});set(gca,'XTickLabelRotation',45)
