%-----------------------------------------------------
% 3D-Deconvolution with Inverse Problem Library
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%----------------------------------------------------
clear all; close all; clc;

%% Load Data
load('acq');  % Acquisitions
load('psf');  % PSF
figure; imagesc(max(acq,[],3)); axis image; axis off; colormap gray;title('Max Intensity Data');

%% Some precomputation
fftacq=fftn(acq);

%% Building the Cost-Function
F_LS=CostL2([],acq);                  % Least-Sqaures functional
H=LinOpConv(fftshift(psf));           % Convolution Operator
R_N12=CostMixNorm12([4]);             % TV regularizer: Mixed Norm 2-1
G=LinOpGrad(size(acq));               % TV regularizer: Operator Gradient
R_POS=CostNonNeg();                   % Non-Negativity: Indicator function
Id=LinOpIdentity(size(psf));          % Non-Negativity: Identity Operator 
R_1sch=CostMixNorm1Schatt([],1);      % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
Hess=LinOpHess(size(acq));            % Hessian-Shatten: Hessian Operator  

%% Parameters
lamb=[1e-4];
maxiter=50;

%% Fourier of the filter G'G (Laplacian)
fGtG=fftn(G.fHtH);
%% Fourier of the filter Hess'Hess
fHesstHess=fftn(Hess.fHtH);

for ll=1:length(lamb)   
    %% TV Regul (ADMM minimization)
    Fn={lamb(ll)*R_N12,R_POS};   % Functionals F_n
    Hn={G,Id};                                   % Associated operators H_n
    rho_n=[1e-3,1e-3];                           % Multipliers rho_n
    solver = @(z,rho,x) real(ifftn((conj(H.mtf).*fftacq + fftn(rho(1)*G'*z{1} + rho(2)*z{2}) )./(abs(H.mtf).^2 + rho(1)*fGtG + rho(2))));
    OutADMM=OutputOpti(1,[],10);
    ADMM=OptiADMM(F_LS,H,Fn,Hn,rho_n,solver,OutADMM);
    ADMM.ItUpOut=10;            % call OutputOpti update every ItUpOut iterations
    ADMM.maxiter=maxiter;       % max number of iterations
    ADMM.run(acq);              % run the algorithm
    figure; imagesc(max(ADMM.xopt,[],3)); axis image; axis off; colormap gray;title('Max Intensity TV');
    
    %% Hessian-Shatten Regul (ADMM minimization)
    Fn={lamb(ll)*R_1sch,R_POS};  % Functionals F_n
    Hn={Hess,Id};                                % Associated operators H_n
    rho_n=[5e-3,5e-3];                           % Multipliers rho_n
    solver = @(z,rho,x) real(ifftn((conj(H.mtf).*fftacq + fftn(rho(1)*Hess'*z{1} + rho(2)*z{2}) )./(abs(H.mtf).^2 + rho(1)*fHesstHess + rho(2))));
    OutADMM=OutputOpti(1,[],10);
    ADMM=OptiADMM(F_LS,H,Fn,Hn,rho_n,solver,OutADMM);
    ADMM.ItUpOut=10;                % call OutputOpti update every ItUpOut iterations
    ADMM.maxiter=maxiter;           % max number of iterations
    ADMM.run(acq);                  % run the algorithm
    figure; imagesc(max(ADMM.xopt,[],3)); axis image; axis off; colormap gray;title('Max Intensity Hessian');
end