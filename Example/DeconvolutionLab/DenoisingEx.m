%-----------------------------------------------------------
% -- Description
% Script which performs some denoising experiments done to compare 
% the Chambolle-Pock algorithm and the ADMM.
%
% In the second part of the script the Hessian Shatten Norm
% regularization is applied to the denoising problem also using
% the Chambolle-Pock and ADMM algorithms
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;warning('off');

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
im=load('StarLikeSample'); im=im.im;   % Load image
sz=size(im); im=im/max(im(:));        
imdisp(im,'Input Image',1);

% -- Generate data
sig_n=1e-1;  % gaussian noise level
y=im + sig_n*randn(sz);
imdisp(y,'Noisy data',1);

% -- Least Squares Functional definition
F=FuncLeastSquares(y);  % Func LS
% -- Regul term
G=LinOpGrad(size(y));      % Operator Gradient
R_N21=FuncMixNorm21([3]);  % Mixed Norm 2-1
lamb=8;                    % Hyperparameter (on the Least Squares functional)

% -- Get precise solution to plot convergence rates
disp('Compute precise solution (5000 iterations of Chambolle-Pock 1/N^2)...');
Out=OutputOpti();CP=OptiChambPock(R_N21,G,FuncMultScalar(F,lamb),Out);
CP.tau=1/G.norm;CP.sig=1/(CP.tau*G.norm^2);CP.gam=0.7*lamb;CP.maxiter=5000;CP.xtol=1e-10; 
CP.run(y);xstar=CP.xopt;
disp('... done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Denoising with TV regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% -- Chambolle-Pock  LS + TV
OutCP=OutputOpti(1,im,5);
CP=OptiChambPock(R_N21,G,FuncMultScalar(F,lamb),OutCP);
CP.tau=1e-2;CP.xtol=1e-10;
CP.sig=1/(CP.tau*G.norm^2);
CP.ItUpOut=100;   % call OutputOpti update every ItUpOut iterations
CP.maxiter=1000;  % max number of iterations
CP.run(y);        % run the algorithm 
subplot(2,2,1);imdisp(CP.xopt,'Chambolle-Pock',0);
err=[];for n=1:length(OutCP.evolxopt), err(n)=norm(xstar(:)-OutCP.evolxopt{n}(:)); end;
subplot(2,2,4);loglog(OutCP.iternum,err/err(1),'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Error');

% -- Chambolle-Pock  LS + TV
OutCP_acc=OutputOpti(1,im,5);
CPacc=OptiChambPock(R_N21,G,FuncMultScalar(F,lamb),OutCP_acc);
CPacc.tau=1/G.norm;CPacc.xtol=1e-10;
CPacc.sig=1/(CPacc.tau*G.norm^2);
CPacc.gam=0.7*lamb;
CPacc.ItUpOut=100;   % call OutputOpti update every ItUpOut iterations
CPacc.maxiter=1000;  % max number of iterations
CPacc.run(y);        % run the algorithm 
subplot(2,2,2);imdisp(CPacc.xopt,'Chambolle-Pock (1/N^2 version)',0);
err=[];for n=1:length(OutCP_acc.evolxopt), err(n)=norm(xstar(:)-OutCP_acc.evolxopt{n}(:)); end;
subplot(2,2,4);hold all;loglog(OutCP_acc.iternum,err/err(1),'LineWidth',1.5); 

% -- ADMM LS + TV
Fn={FuncMultScalar(F,lamb),R_N21};
Hn={LinOpIdentity(sz),G};rho_n=[10,10];
OutADMM=OutputOpti(1,im,5);
ADMM=OptiADMM([],[],Fn,Hn,rho_n,[],OutADMM);
ADMM.
ADMM.xtol=1e-10;
ADMM.ItUpOut=100;   % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=1000;  % max number of iterations
ADMM.run(y);        % run the algorithm 
subplot(2,2,3);imdisp(ADMM.xopt,'ADMM',0);
err=[];for n=1:length(OutADMM.evolxopt), err(n)=norm(xstar(:)-OutADMM.evolxopt{n}(:)); end;
subplot(2,2,4);loglog(OutADMM.iternum,err/err(1),'LineWidth',1.5);
loglog(1:200,1./[1:200],'--','LineWidth',1.5);loglog(1:200,1./[1:200].^2,'--','LineWidth',1.5);
legend('CP (1/N)','CP (1/N^2)','ADMM','1/N','1/N^2');

% -- Plot Evolution SNR and Running Time for TV-Reg methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(OutCP.iternum,OutCP.evolsnr,'LineWidth',1.5); 
semilogy(OutCP_acc.iternum,OutCP_acc.evolsnr,'LineWidth',1.5);
semilogy(OutADMM.iternum,OutADMM.evolsnr,'LineWidth',1.5);
legend('CP (1/N)','CP (1/N^2)','ADMM');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[CP.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[CPacc.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[ADMM.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3]);ylabel('Time (s)');
set(gca,'xticklabels',{'CP (1/N)','CP (1/N^2)','ADMM'});set(gca,'XTickLabelRotation',45)


