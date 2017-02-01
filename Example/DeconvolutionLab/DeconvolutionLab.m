%-----------------------------------------------------------
%
%
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;warning('off');

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
fftHty=conj(fft2(psf)).*fft2(y);
imdisp(y(idx,idx),'Convolved and noisy data',1);

% -- Least Squares Functional definition
F_LS=FuncLeastSquares(y,H);  % Func LS
F_KL=FuncKullLeib(y,H);      % Func KL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non Regularized optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% -- Gradient Descent LS
Out_GDls=OutputOpti(1,impad,40);
GD_LS=OptiGradDsct(F_LS,Out_GDls);
GD_LS.ItUpOut=10;     % call OutputOpti update every ItUpOut iterations
GD_LS.maxiter=200;    % max number of iterations
GD_LS.run(y);         % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=max(Out_GDls.evolsnr(:));
subplot(1,2,1);imdisp(Out_GDls.evolxopt{n}(idx,idx),'LS (GD)',0);
subplot(1,2,2);plot(Out_GDls.iternum,Out_GDls.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NonNegativity Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Func NonNeg definition
R_POS=FuncNonNeg();   

figure;
% -- FISTA LS + NonNeg
Out_FSBlsPos=OutputOpti(1,impad,40);
FBS_LSPOS=OptiFBS(F_LS,R_POS,Out_FSBlsPos);
FBS_LSPOS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
FBS_LSPOS.fista=true;   % activate fista
FBS_LSPOS.maxiter=200;  % max number of iterations
FBS_LSPOS.run(y);       % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
[v,n]=max(Out_FSBlsPos.evolsnr(:));
subplot(2,3,1);imdisp(Out_FSBlsPos.evolxopt{n}(idx,idx),'LS + NonNeg (FISTA)',0);
subplot(2,3,4);plot(Out_FSBlsPos.iternum,Out_FSBlsPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- FISTA KL + NonNeg
Out_FSBklPos=OutputOpti(1,impad,40);
FBS_KLPOS=OptiFBS(F_KL,R_POS,Out_FSBklPos);
FBS_KLPOS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
FBS_KLPOS.fista=true;   % activate fista
FBS_KLPOS.maxiter=200;  % max number of iterations
FBS_KLPOS.gam=1e-2;     % set gamma parameter
FBS_KLPOS.run(y);       % run the algorithm 
[v,n]=max(Out_FSBklPos.evolsnr(:));
subplot(2,3,2);imdisp(Out_FSBklPos.evolxopt{n}(idx,idx),'KL + NonNeg (FISTA)',0);
subplot(2,3,5);plot(Out_FSBklPos.iternum,Out_FSBklPos.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');


% -- Richardson-Lucy KL + NonNeg (implicit)
Out_RLklPos=OutputOpti(1,impad,40);
RL_KLPOS=OptiRichLucy(F_KL,Out_RLklPos);
RL_KLPOS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
RL_KLPOS.maxiter=200;  % max number of iterations
RL_KLPOS.run(y);       % run the algorithm 
[v,n]=max(Out_RLklPos.evolsnr(:));
subplot(2,3,3);imdisp(Out_RLklPos.evolxopt{n}(idx,idx),'KL + NonNeg (RL)',0);
subplot(2,3,5);hold all; plot(Out_RLklPos.iternum,Out_RLklPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('FISTA','RL');

% -- Plot Evolution SNR and Running Time for No-Reg methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(Out_GDls.iternum,Out_GDls.evolsnr,'LineWidth',1.5); 
semilogy(Out_FSBlsPos.iternum,Out_FSBlsPos.evolsnr,'LineWidth',1.5);
semilogy(Out_FSBklPos.iternum,Out_FSBklPos.evolsnr,'LineWidth',1.5); 
semilogy(Out_RLklPos.iternum,Out_RLklPos.evolsnr,'LineWidth',1.5); 
legend('LS (GD)','LS+POS (FBS)','KL+POS (FBS)','KL+POS (RL)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[GD_LS.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[FBS_LSPOS.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
bar(3,[FBS_KLPOS.time],'FaceColor',orderCol(3,:),'EdgeColor','k');
bar(4,[RL_KLPOS.time],'FaceColor',orderCol(4,:),'EdgeColor','k');
set(gca,'xtick',[1 2 3 4]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS (GD)','LS+POS (FBS)','KL+POS (FBS)','KL+POS (RL)'});set(gca,'XTickLabelRotation',45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Regul term
G=LinOpGrad(size(y));      % Operator Gradient
R_N12=FuncMixNorm12([3]);  % Mixed Norm 2-1
lamb=5e-4;                 % Hyperparameter

figure;
% -- Chambolle-Pock  LS + TV
Out_CPlstv=OutputOpti(1,impad,40);
CP_LSTV=OptiChambPock(FuncMultScalar(R_N12,lamb),G,F_LS,Out_CPlstv);
CP_LSTV.tau=15;
CP_LSTV.sig=1/(CP_LSTV.tau*G.norm^2)*0.99;
CP_LSTV.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
CP_LSTV.maxiter=200;  % max number of iterations
CP_LSTV.run(y);       % run the algorithm 
subplot(1,3,1);imdisp(Out_CPlstv.evolxopt{end}(idx,idx),'LS + TV (CP)',0);
subplot(1,3,3);plot(Out_CPlstv.iternum,Out_CPlstv.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- ADMM LS + TV
Fn={FuncMultScalar(R_N12,lamb)};
Hn={G};rho_n=[1e-1];
lap=zeros(size(impad)); lap(1,1)=4; lap(1,2)=-1;lap(2,1)=-1; lap(1,end)=-1;lap(end,1)=-1; Flap=fft2(lap);
solver = @(z,rho) real(ifft2((fftHty + rho(1)*fft2(G'*z{1}))./(abs(H.mtf).^2 + rho(1)*Flap)));
Out_ADMMlstv=OutputOpti(1,impad,40);
ADMM_LSTV=OptiADMM(FuncLeastSquares(y),H,Fn,Hn,rho_n,solver,Out_ADMMlstv);
ADMM_LSTV.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
ADMM_LSTV.maxiter=200;  % max number of iterations
ADMM_LSTV.run(y);       % run the algorithm 
subplot(1,3,2);imdisp(Out_ADMMlstv.evolxopt{end}(idx,idx),'LS + TV (ADMM)',0);
subplot(1,3,3);hold all;plot(Out_ADMMlstv.iternum,Out_ADMMlstv.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');

% -- Plot Evolution SNR and Running Time for TV-Reg methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(Out_CPlstv.iternum,Out_CPlstv.evolsnr,'LineWidth',1.5); 
semilogy(Out_ADMMlstv.iternum,Out_ADMMlstv.evolsnr,'LineWidth',1.5);
legend('LS+TV (CP)','LS+TV (ADMM)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[CP_LSTV.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[ADMM_LSTV.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+TV (CP)','LS+TV (ADMM)'});set(gca,'XTickLabelRotation',45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hessian Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Hessian Schatten Regul
Hess=LinOpHess(size(impad));     % Hessian Operator
R_1sch=FuncMixNorm1Schatt([],1); % Mixed Norm 1-Schatten (p=1)
lamb=2e-3;                       % Hyperparameter

figure
% -- Chambolle-Pock  LS + ShattenHess
Out_CPlshess=OutputOpti(1,impad,40);
CP_LSHESS=OptiChambPock(FuncMultScalar(R_1sch,lamb),Hess,F_LS,Out_CPlshess);
CP_LSHESS.tau=15;
CP_LSHESS.sig=0.001;
CP_LSHESS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
CP_LSHESS.maxiter=200;  % max number of iterations
CP_LSHESS.run(y);       % run the algorithm 
subplot(1,3,1);imdisp(Out_CPlshess.evolxopt{end}(idx,idx),'LS + Hess (CP)',0);
subplot(1,3,3);plot(Out_CPlshess.iternum,Out_CPlshess.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- ADMM LS + ShattenHess
Fn={FuncLeastSquares(y),FuncMultScalar(R_1sch,lamb)};
Hn={H,Hess};rho_n=[1e-1,1e-1];
Out_ADMMlshess=OutputOpti(1,impad,40);
ADMM_LSHESS=OptiADMM([],[],Fn,Hn,rho_n,[],Out_ADMMlshess);
ADMM_LSHESS.maxiterCG=2;  % 2 CG iterations are sufficient for this example
ADMM_LSHESS.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
ADMM_LSHESS.maxiter=200;  % max number of iterations
ADMM_LSHESS.run(y);       % run the algorithm 
subplot(1,3,2);imdisp(Out_ADMMlshess.evolxopt{end}(idx,idx),'LS + Hess (ADMM)',0);
subplot(1,3,3);hold all;plot(Out_ADMMlshess.iternum,Out_ADMMlshess.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');

% -- Plot Evolution SNR and Running Time for Hess-Reg methods
figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(Out_CPlshess.iternum,Out_CPlshess.evolsnr,'LineWidth',1.5); 
semilogy(Out_ADMMlshess.iternum,Out_ADMMlshess.evolsnr,'LineWidth',1.5);
legend('LS+TV (CP)','LS+TV (ADMM)');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time (200 iterations)');set(gca,'FontSize',12);
orderCol=get(gca,'ColorOrder');
bar(1,[CP_LSHESS.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[ADMM_LSHESS.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+HESS (CP)','LS+HESS (ADMM)'});set(gca,'XTickLabelRotation',45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wavelet regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- FISTA LS + Wavelet

% -- ADMM LS + Wavelet

% -- Chambolle-Pock LS + Wavelet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Positivity + TV/Hessian/Wavelet Regul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure;
% -- ADMM LS + TV + NonNeg
lamb=2e-3; 
Fn={FuncLeastSquares(y),FuncMultScalar(R_N21,lamb),R_POS};
Hn={H,G,LinOpIdentity(size(impad))};
rho_n=[1e-1,1e-1,1e-1];
VU_ADMMlstvPos=VerbUpdate(1,impad);
ADMM_LSTVPOS=OptiADMM([],[],Fn,Hn,rho_n,[],VU_ADMMlstvPos);
ADMM_LSTVPOS.verb=10;
ADMM_LSTVPOS.maxiter=200;
ADMM_LSTVPOS.run(y);
subplot(2,3,1);imdisp(VU_ADMMlstvPos.evolxopt{end}(idx,idx),'LS + TV + POS (ADMM)',0);
subplot(2,3,5);hold all;plot(VU_ADMMlstvPos.iternum,VU_ADMMlstvPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');

% -- PrimalDual Condat LS + TV + NonNeg

% -- ADMM KL + TV + NonNeg
lamb=3e-2;                 % Hyperparameter
Fn={FuncKullLeib(y,[],0),FuncMultScalar(R_N21,lamb),R_POS};
Hn={H,G,LinOpIdentity(size(impad))};
rho_n=[1e-1,1e-1,1e-1];
VU_ADMMkltvPos=VerbUpdate(1,impad);
ADMM_KLTVPOS=OptiADMM([],[],Fn,Hn,rho_n,[],VU_ADMMkltvPos);
ADMM_KLTVPOS.verb=10;
ADMM_KLTVPOS.maxiter=200;
ADMM_KLTVPOS.run(y);
subplot(2,3,3);imdisp(VU_ADMMkltvPos.evolxopt{end}(idx,idx),'KL + TV + POS (ADMM)',0);
subplot(2,3,6);hold all;plot(VU_ADMMkltvPos.iternum,VU_ADMMkltvPos.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
%}
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

