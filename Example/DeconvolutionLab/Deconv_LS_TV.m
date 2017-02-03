%-----------------------------------------------------------
% Deconv_LS_TV script: Deconvolution by minimizing the 
% Least-Squares function plus the TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% using 
%      - Chambolle-Pock
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpGrad, Func, FuncLeastSquares,   
% FuncMixNorm12, Opti, OptiChambPock, OptiADMM, OutpuOpti
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------------------
clear all; close all; clc;warning('off');
help Deconv_LS_TV

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
fftHty=conj(H.mtf).*fft2(y);

% -- Functions definition
F_LS=FuncLeastSquares(y,H);  % Least-Sqaures data term
R_N12=FuncMixNorm12([3]);    % Mixed Norm 2-1
G=LinOpGrad(size(y));        % Operator Gradient
lamb=1e-3;                   % Hyperparameter

% -- Chambolle-Pock  LS + TV
OutCP=OutputOpti(1,impad,40);
CP=OptiChambPock(FuncMultScalar(R_N12,lamb),G,F_LS,OutCP);
CP.tau=15;                            % algorithm parameters
CP.sig=1/(CP.tau*G.norm^2)*0.99;      %
CP.ItUpOut=10;                        % call OutputOpti update every ItUpOut iterations
CP.maxiter=200;                       % max number of iterations
CP.run(y);                            % run the algorithm 

% -- ADMM LS + TV
Fn={FuncMultScalar(R_N12,lamb)};
Hn={G};rho_n=[1e-1];
lap=zeros(size(impad)); lap(1,1)=4; lap(1,2)=-1;lap(2,1)=-1; lap(1,end)=-1;lap(end,1)=-1; Flap=fft2(lap);  
solver = @(z,rho) real(ifft2((fftHty + rho(1)*fft2(G'*z{1}))./(abs(H.mtf).^2 + rho(1)*Flap)));            % solver to solve the x update
OutADMM=OutputOpti(1,impad,40);
ADMM=OptiADMM(FuncLeastSquares(y),H,Fn,Hn,rho_n,solver,OutADMM);
ADMM.ItUpOut=10;   % call OutputOpti update every ItUpOut iterations
ADMM.maxiter=200;  % max number of iterations
ADMM.run(y);       % run the algorithm 


% -- Display
imdisp(OutCP.evolxopt{end}(idx,idx),'LS + TV (CP)',1);
imdisp(OutADMM.evolxopt{end}(idx,idx),'LS + TV (ADMM)',1);
figure;plot(OutCP.iternum,OutCP.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
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
set(gca,'xticklabels',{'LS+TV (CP)','LS+TV (ADMM)'});set(gca,'XTickLabelRotation',45)  
