%-----------------------------------------------------------
%
%
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;

% -- Input image and psf
im=double(imread('cameraman.tif'));              % Load image
im=im/max(im(:));
ll=linspace(-1,1,512);                           % Build psf
[X,Y]=meshgrid(ll,ll);
sig_psf=5e-3;
psf=fftshift(exp(-(X.^2+Y.^2)/(2*sig_psf^2)));
psf=psf/sum(psf(:));
imdisp(im,'Input Image',1);

% -- Image padding
impad=zeros(512);
idx=129:384;
impad(idx,idx)=im;

% -- Convolution Operator definition
H=LinOpConv(psf);

% -- Generate data
sig_n=1e-2;  % gaussian noise level
phobud=1e5;  % photon budget
tmp=H*impad*phobud; 
y=zeros(512); y(idx,idx)=poissrnd(tmp(idx,idx))/phobud + sig_n*randn(256);
imdisp(y(idx,idx),'Convolved and noisy data',1);

% -- Least Squares Functional definition
F_LS=FuncLeastSquares(y,H);  % Func LS

% -- Initial point
x0=zeros(size(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non Regularized optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% -- Gradient Descent LS
VerbUpGDls=VerbUpdate(1,impad);
GD_LS=OptiGradDsct(F_LS,VerbUpGDls);
GD_LS.verb=2;        % verbose upate every verb iterations
GD_LS.maxiter=30;    % max number of iterations
GD_LS.run(x0);       % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VerbUpGDls.evolerr(:));
subplot(2,2,1);imdisp(VerbUpGDls.evolxopt{n}(idx,idx),'LS (GD)',0);
subplot(2,2,3);loglog(VerbUpGDls.iternum,VerbUpGDls.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- Gradient Descent KL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NonNegativity Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Func NonNeg definition
R_POS=FuncNonNeg();    

figure;
% -- FISTA LS + NonNeg
VerbUpFSBlsPos=VerbUpdate(1,impad);
FBS_LSPOS=OptiFBS(F_LS,R_POS,1,VerbUpFSBlsPos);
FBS_LSPOS.verb=2;       % verbose upate every verb iterations
FBS_LSPOS.fista=true;   % activate fista
FBS_LSPOS.maxiter=30;   % max number of iterations
FBS_LSPOS.run(x0);      % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=min(VerbUpFSBlsPos.evolerr(:));
subplot(2,3,1);imdisp(VerbUpFSBlsPos.evolxopt{n}(idx,idx),'LS + NonNeg (FBS)',0);
subplot(2,3,4);loglog(VerbUpFSBlsPos.iternum,VerbUpFSBlsPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- FISTA KL + NonNeg

% -- Richardson-Lucy KL + NonNeg (implicit)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison Error Ground Thruth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; grid; hold all; title('Error Ground Thruth');
semilogy(VerbUpGDls.iternum,VerbUpGDls.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
semilogy(VerbUpFSBlsPos.iternum,VerbUpFSBlsPos.evolerr,'LineWidth',1.5);grid; set(gca,'FontSize',12);
legend('LS (GD)','LS+POS (FBS)');
