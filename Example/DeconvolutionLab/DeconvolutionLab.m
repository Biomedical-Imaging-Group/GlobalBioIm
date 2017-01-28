%-----------------------------------------------------------
%
%
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
im=load('StarLikeSample'); im=im.im;   % Load image
sz=size(im); im=im/max(im(:));                    
fc=0.15;                                % Build psf
f = linspace(-0.5,0.5,2*sz(1));
[Fx,Fy] = meshgrid(f,f);[THETA,RHO]=cart2pol(Fx,Fy); 
OTF = 1/pi*(2*(acos(RHO./fc)) - sin(2*acos(RHO/fc)));
OTF = abs(OTF); OTF(RHO > fc) = 0;OTF=fftshift(OTF);
psf=real(ifft2(OTF));
psf=psf/sum(psf(:));
imdisp(im,'Input Image',1);
imdisp(fftshift(psf),'PSF',1);colormap parula;

% -- Image padding
impad=zeros(512);
idx=129:384;
impad(idx,idx)=im;

% -- Convolution Operator definition
H=LinOpConv(psf);

% -- Generate data
sig_n=3e-2;  % gaussian noise level
phobud=1e5;  % photon budget
tmp=H*impad*phobud; 
y=zeros(512); y(idx,idx)=poissrnd(tmp(idx,idx))/phobud + sig_n*randn(256);
y=max(y,0);
imdisp(y(idx,idx),'Convolved and noisy data',1);

% -- Least Squares Functional definition
F_LS=FuncLeastSquares(y,H);  % Func LS
F_KL=FuncKullLeib(y,H);      % Func KL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non Regularized optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% -- Gradient Descent LS
VU_GDls=VerbUpdate(1,impad);
GD_LS=OptiGradDsct(F_LS,VU_GDls);
GD_LS.verb=20;        % verbose upate every verb iterations
GD_LS.maxiter=200;    % max number of iterations
GD_LS.run(y);         % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VU_GDls.evolerr(:));
subplot(1,2,1);imdisp(VU_GDls.evolxopt{n}(idx,idx),'LS (GD)',0);
subplot(1,2,2);loglog(VU_GDls.iternum,VU_GDls.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NonNegativity Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Func NonNeg definition
R_POS=FuncNonNeg();   

figure;
% -- FISTA LS + NonNeg
VU_FSBlsPos=VerbUpdate(1,impad);
FBS_LSPOS=OptiFBS(F_LS,R_POS,VU_FSBlsPos);
FBS_LSPOS.verb=20;      % verbose upate every verb iterations
FBS_LSPOS.fista=true;   % activate fista
FBS_LSPOS.maxiter=200;  % max number of iterations
FBS_LSPOS.run(y);       % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VU_FSBlsPos.evolerr(:));
subplot(2,3,1);imdisp(VU_FSBlsPos.evolxopt{n}(idx,idx),'LS + NonNeg (FISTA)',0);
subplot(2,3,4);loglog(VU_FSBlsPos.iternum,VU_FSBlsPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Error');

% -- FISTA KL + NonNeg
VU_FSBklPos=VerbUpdate(1,impad);
FBS_KLPOS=OptiFBS(F_KL,R_POS,VU_FSBklPos);
FBS_KLPOS.verb=20;      % verbose upate every verb iterations
FBS_KLPOS.fista=true;   % activate fista
FBS_KLPOS.maxiter=200;  % max number of iterations
FBS_KLPOS.gam=1e-2;
FBS_KLPOS.run(y);      % run the algorithm 
[v,n]=min(VU_FSBklPos.evolerr(:));
subplot(2,3,2);imdisp(VU_FSBklPos.evolxopt{n}(idx,idx),'KL + NonNeg (FBS)',0);
subplot(2,3,5);loglog(VU_FSBklPos.iternum,VU_FSBklPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');


% -- Richardson-Lucy KL + NonNeg (implicit)
VU_RLklPos=VerbUpdate(1,impad);
RL_KLPOS=OptiRichLucy(F_KL,VU_RLklPos);
RL_KLPOS.verb=20;      % verbose upate every verb iterations
RL_KLPOS.maxiter=200;  % max number of iterations
RL_KLPOS.run(y);      % run the algorithm 
[v,n]=min(VU_RLklPos.evolerr(:));
subplot(2,3,3);imdisp(VU_RLklPos.evolxopt{n}(idx,idx),'KL + NonNeg (RL)',0);
subplot(2,3,5);hold all; loglog(VU_RLklPos.iternum,VU_RLklPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('FISTA','RL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Regul term
G=LinOpGrad(size(y));      % Operator Gradient
R_N21=FuncMixNorm21([3]);  % Mixed Norm 2-1
lamb=5e-3;                 % Hyperparameter

figure;
% -- Chambolle-Pock  LS + TV
VU_CPlstv=VerbUpdate(1,impad);
CP_LSTV=OptiChambPock(FuncMultScalar(R_N21,lamb),G,F_LS,VU_CPlstv);
CP_LSTV.tau=10;
CP_LSTV.sig=1/(CP_LSTV.tau*G.norm^2)-eps;
CP_LSTV.verb=20;      % verbose upate every verb iterations
CP_LSTV.maxiter=200;  % max number of iterations
CP_LSTV.run(y);       % run the algorithm 
subplot(1,3,1);imdisp(VU_CPlstv.evolxopt{end}(idx,idx),'LS + TV (CP)',0);
subplot(1,3,3);loglog(VU_CPlstv.iternum,VU_CPlstv.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- ADMM LS + TV
Fn={FuncLeastSquares(y),FuncMultScalar(R_N21,lamb)};
Hn={H,G};
rho_n=[1e-1,1e-1];
VU_ADMMlstv=VerbUpdate(1,impad);
ADMM_LSTV=OptiADMM([],[],Fn,Hn,rho_n,[],VU_ADMMlstv);
ADMM_LSTV.verb=20;
ADMM_LSTV.maxiter=200;
ADMM_LSTV.run(y);
subplot(1,3,2);imdisp(VU_ADMMlstv.evolxopt{end}(idx,idx),'LS + TV (ADMM)',0);
subplot(1,3,3);hold all;loglog(VU_ADMMlstv.iternum,VU_ADMMlstv.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hessian Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Chambolle-Pock  LS + ShattenHess

% -- ADMM LS + ShattenHess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wavelet regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- FISTA LS + Wavelet

% -- ADMM LS + Wavelet

% -- Chambolle-Pock LS + Wavelet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Positivity + TV/Hessian/Wavelet Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% -- ADMM LS + TV + NonNeg
lamb=2e-3; 
Fn={FuncLeastSquares(y),FuncMultScalar(R_N21,lamb),R_POS};
Hn={H,G,LinOpIdentity(size(impad))};
rho_n=[1e-1,1e-1,1e-1];
VU_ADMMlstvPos=VerbUpdate(1,impad);
ADMM_LSTVPOS=OptiADMM([],[],Fn,Hn,rho_n,[],VU_ADMMlstvPos);
ADMM_LSTVPOS.verb=20;
ADMM_LSTVPOS.maxiter=200;
ADMM_LSTVPOS.run(y);
subplot(2,3,1);imdisp(VU_ADMMlstvPos.evolxopt{end}(idx,idx),'LS + TV + POS (ADMM)',0);
subplot(2,3,5);hold all;loglog(VU_ADMMlstvPos.iternum,VU_ADMMlstvPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- PrimalDual Condat LS + TV + NonNeg

% -- ADMM KL + TV + NonNeg
lamb=3e-2;                 % Hyperparameter
Fn={FuncKullLeib(y,[],0),FuncMultScalar(R_N21,lamb),R_POS};
Hn={H,G,LinOpIdentity(size(impad))};
rho_n=[1e-1,1e-1,1e-1];
VU_ADMMkltvPos=VerbUpdate(1,impad);
ADMM_KLTVPOS=OptiADMM([],[],Fn,Hn,rho_n,[],VU_ADMMkltvPos);
ADMM_KLTVPOS.verb=20;
ADMM_KLTVPOS.maxiter=200;
ADMM_KLTVPOS.run(y);
subplot(2,3,3);imdisp(VU_ADMMkltvPos.evolxopt{end}(idx,idx),'KL + TV + POS (ADMM)',0);
subplot(2,3,6);hold all;loglog(VU_ADMMkltvPos.iternum,VU_ADMMkltvPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- PrimalDual Condat KL + TV + NonNeg

% -- Richardson-Lucy-TV  KL + TV + NonNeg (implicit)


% -- ADMM LS + ShattenHess + NonNeg

% -- PrimalDual Condat LS + ShattenHess + NonNeg

% -- ADMM KL + ShattenHess + NonNeg

% -- PrimalDual Condat KL + ShattenHess + NonNeg



% -- ADMM LS + Wavelet + NonNeg

% -- PrimalDual Condat LS + Wavelet + NonNeg

% -- ADMM KL + Wavelet + NonNeg

% -- PrimalDual Condat KL + Wavelet + NonNeg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Comparison Error Ground Thruth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; grid; hold all; title('Error Ground Truth');
semilogy(VU_GDls.iternum,VU_GDls.evolerr,'LineWidth',1.5); set(gca,'FontSize',12);
semilogy(VU_FSBlsPos.iternum,VU_FSBlsPos.evolerr,'LineWidth',1.5); set(gca,'FontSize',12);
semilogy(VU_FSBklPos.iternum,VU_FSBklPos.evolerr,'LineWidth',1.5); set(gca,'FontSize',12);
semilogy(VU_RLklPos.iternum,VU_RLklPos.evolerr,'LineWidth',1.5); set(gca,'FontSize',12);
semilogy(VU_CPlstv.iternum,VU_CPlstv.evolerr,'LineWidth',1.5,'LineStyle','--'); set(gca,'FontSize',12);
semilogy(VU_ADMMlstv.iternum,VU_ADMMlstv.evolerr,'LineWidth',1.5,'LineStyle','--'); set(gca,'FontSize',12);
legend('LS (GD)','LS+POS (FBS)','KL+POS (FBS)','KL+POS (RL)','LS+TV (CP)','LS+TV (ADMM)');
