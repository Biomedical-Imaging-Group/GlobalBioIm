%-----------------------------------------------------------
%
%
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;

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
phobud=1e4;  % photon budget
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
VerbUpGDls=VerbUpdate(1,impad);
GD_LS=OptiGradDsct(F_LS,VerbUpGDls);
GD_LS.verb=20;        % verbose upate every verb iterations
GD_LS.maxiter=200;    % max number of iterations
GD_LS.run(y);       % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VerbUpGDls.evolerr(:));
subplot(2,1,1);imdisp(VerbUpGDls.evolxopt{n}(idx,idx),'LS (GD)',0);
subplot(2,1,2);loglog(VerbUpGDls.iternum,VerbUpGDls.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NonNegativity Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Func NonNeg definition
R_POS=FuncNonNeg();    

figure;
% -- FISTA LS + NonNeg
VerbUpFSBlsPos=VerbUpdate(1,impad);
FBS_LSPOS=OptiFBS(F_LS,R_POS,1,VerbUpFSBlsPos);
FBS_LSPOS.verb=20;      % verbose upate every verb iterations
FBS_LSPOS.fista=true;   % activate fista
FBS_LSPOS.maxiter=200;  % max number of iterations
FBS_LSPOS.run(y);      % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VerbUpFSBlsPos.evolerr(:));
subplot(2,3,1);imdisp(VerbUpFSBlsPos.evolxopt{n}(idx,idx),'LS + NonNeg (FBS)',0);
subplot(2,3,4);loglog(VerbUpFSBlsPos.iternum,VerbUpFSBlsPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- FISTA KL + NonNeg
VerbUpFSBklPos=VerbUpdate(1,impad);
FBS_KLPOS=OptiFBS(F_KL,R_POS,1,VerbUpFSBklPos);
FBS_KLPOS.verb=20;      % verbose upate every verb iterations
FBS_KLPOS.fista=true;   % activate fista
FBS_KLPOS.maxiter=200;  % max number of iterations
FBS_KLPOS.gam=1e-2;
FBS_KLPOS.run(y);      % run the algorithm 
[v,n]=min(VerbUpFSBklPos.evolerr(:));
subplot(2,3,2);imdisp(VerbUpFSBklPos.evolxopt{n}(idx,idx),'KL + NonNeg (FBS)',0);
subplot(2,3,5);loglog(VerbUpFSBklPos.iternum,VerbUpFSBklPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- Richardson-Lucy KL + NonNeg (implicit)
VerbUpRLklPos=VerbUpdate(1,impad);
RL_KLPOS=OptiRichLucy(F_KL,VerbUpRLklPos);
RL_KLPOS.verb=20;      % verbose upate every verb iterations
RL_KLPOS.maxiter=200;  % max number of iterations
RL_KLPOS.run(y);      % run the algorithm 
[v,n]=min(VerbUpRLklPos.evolerr(:));
subplot(2,3,3);imdisp(VerbUpRLklPos.evolxopt{n}(idx,idx),'KL + NonNeg (RL)',0);
subplot(2,3,5);hold all; loglog(VerbUpRLklPos.iternum,VerbUpRLklPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('FBS','RL');

%-- Comparison Error Ground Thruth between non regularized methods
figure; grid; hold all; title('Error Ground Thruth No-Regul Methods');
semilogy(VerbUpGDls.iternum,VerbUpGDls.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
semilogy(VerbUpFSBlsPos.iternum,VerbUpFSBlsPos.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
semilogy(VerbUpFSBklPos.iternum,VerbUpFSBklPos.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
semilogy(VerbUpRLklPos.iternum,VerbUpRLklPos.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
legend('LS (GD)','LS+POS (FBS)','KL+POS (FBS)','KL+POS (RL)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Chambolle-Pock  LS + TV

% -- ADMM LS + TV

% -- ADMM KL + TV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hessian Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Chambolle-Pock  LS + ShattenHess

% -- ADMM LS + ShattenHess

% -- ADMM KL + ShattenHess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wavelet regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- FISTA LS + Wavelet

% -- ADMM LS + Wavelet

% -- Chambolle-Pock LS + Wavelet

% -- FISTA KL + Wavelet

% -- ADMM KL + Wavelet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Positivity + TV/Hessian/Wavelet Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- ADMM LS + TV + NonNeg

% -- PrimalDual Condat LS + TV + NonNeg

% -- ADMM KL + TV + NonNeg

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


