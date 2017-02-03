%-----------------------------------------------------------
% Deconv_LS_HessSchatt script: Deconvolution by minimizing the
% Least-Squares function  plus the Hessian-Schatten regularizer:
%     0.5 ||Hx - y||^2  + lamb*||Hess*x||_{1,S_p}
% using 
%      - Chambolle-Pock
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpHess, Func, FuncLeastSquares,   
% FuncMixNorm1Schatt, Opti, OptiChambPock, OptiADMM, OutpuOpti
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;warning('off');
help Deconv_LS_HessSchatt

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
F_LS=FuncLeastSquares(y,H);      % Least-Sqaures data term
Hess=LinOpHess(size(impad));     % Hessian Operator
R_1sch=FuncMixNorm1Schatt([],1); % Mixed Norm 1-Schatten (p=1)
lamb=2e-3;                       % Hyperparameter

% -- Chambolle-Pock  LS + ShattenHess
OutCP=OutputOpti(1,impad,40);
CP=OptiChambPock(FuncMultScalar(R_1sch,lamb),Hess,F_LS,OutCP);
CP.tau=1;        % algorithm parameters
CP.sig=0.02;     %
CP.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
CP.maxiter=1000;  % max number of iterations
CP.run(y);       % run the algorithm 

% -- ADMM LS + ShattenHess
Fn={FuncLeastSquares(y),FuncMultScalar(R_1sch,lamb)};
Hn={H,Hess};rho_n=[1e-1,1e-1];
OutADMM=OutputOpti(1,impad,40);
ADMM=OptiADMM([],[],Fn,Hn,rho_n,[],OutADMM);
ADMM.maxiterCG=2;  % 2 CG iterations are sufficient for this example
ADMM.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;  % max number of iterations
ADMM.run(y);       % run the algorithm 

% -- Display
imdisp(OutCP.evolxopt{end}(idx,idx),'LS + Hess (CP)',1);
imdisp(OutADMM.evolxopt{end}(idx,idx),'LS + Hess (ADMM)',1);
figure; plot(OutCP.iternum,OutCP.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all;plot(OutADMM.iternum,OutADMM.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutCP.iternum,OutCP.evolsnr,'LineWidth',1.5); 
semilogy(OutADMM.iternum,OutADMM.evolsnr,'LineWidth',1.5);
legend('LS+TV (CP)','LS+TV (ADMM)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[CP.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[ADMM.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+HESS (CP)','LS+HESS (ADMM)'});set(gca,'XTickLabelRotation',45)


