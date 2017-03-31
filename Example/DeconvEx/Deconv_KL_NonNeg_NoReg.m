%-----------------------------------------------------------
% Deconv_KL_NonNeg_NoReg script: Deconvolution by minimizing
% the Kullback-Leibler divergence  plus the NonNegativity constraint 
% without any regularizer:
%     \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x)
% using:
%    - FISTA
%    - RichardsonLucy
%
% See LinOp, LinOpConv, Cost, CostKullLeib, CostNonNeg, Opti, 
% OptiFBS, OptiRichLucy, OutpuOpti
%------------------------------------------------------------
clear all; close all; clc;warning('off');
help Deconv_KL_NonNeg_NoReg
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
F_KL=CostKullLeib(H,y,1e-6);     % Kullback-Leibler divergence data term
R_POS=CostNonNeg();              % Non-Negativity

% -- FISTA KL + NonNeg
OutFista=OutputOpti(1,impad,40);
FBS=OptiFBS(F_KL,R_POS,OutFista);
FBS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
FBS.fista=true;   % activate fista
FBS.maxiter=200;  % max number of iterations
FBS.gam=1e-2;     % set gamma parameter
FBS.run(y);       % run the algorithm 

% -- Richardson-Lucy KL + NonNeg (implicit)
OutRL=OutputOpti(1,impad,40);
RL=OptiRichLucy(F_KL,[],[],OutRL);
RL.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
RL.maxiter=200;  % max number of iterations
RL.run(y);       % run the algorithm 

% -- Display
[v,n]=max(OutFista.evolsnr(:));
imdisp(OutFista.evolxopt{n}(idx,idx),'KL + NonNeg (FISTA)',1);
[v,n]=max(OutRL.evolsnr(:));
imdisp(OutRL.evolxopt{n}(idx,idx),'KL + NonNeg (RL)',1);

figure;plot(OutFista.iternum,OutFista.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all; plot(OutRL.iternum,OutRL.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');title('Cost evolution');
legend('FISTA','RL');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutFista.iternum,OutFista.evolsnr,'LineWidth',1.5); 
semilogy(OutRL.iternum,OutRL.evolsnr,'LineWidth',1.5); 
legend('KL+POS (FBS)','KL+POS (RL)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[FBS.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[RL.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'KL+POS (FBS)','KL+POS (RL)'});set(gca,'XTickLabelRotation',45)
