%----------------------------------------------------
% This script performs 3D deconvolution using ADMM with
%    - Data-term: Least-Squares or Kullback-Leibler
%    - regul: TV or Hessian-Schatten norm
%
% Before to run it set the following papareters:
% -- General 
%  outFolder=      -> Folder to save results
%  dataname=       -> file name data image
%  psfname=        -> file name psf
% -- Deconvolution
%  lamb=       -> Hyperparameter for final deconvolution (can be an array to loop)
%  maxIt=      -> Max iterations
%  Reg=        -> Choice regul: 1 for TV, 2 for Hessian-Schatten
%  DataTerm=   -> Choice data-term: 1 for LS, 2 for KL
%  usemask=    -> if true, use mask -> ||SHx - y||^2 where S is a LinOpSelectorPatch
%  padsz=      -> padding size the size of reconstructed x will be sizedata
%                 + 2*padsz (when usemask=true)
%                 WARNING: The new size will be given to the function
%                 fft_best_dim in order to get the closest size which is
%                 "optimal" in terms of fft computation. Hence the final
%                 used reconstruction size may differ from sizedata + 2*padsz
%
% WARNING : This script uses the functions loadtiff and saveastiff to read
% and save 3D tiff stack. These functions can be downloaded at:
% https://ch.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack?focused=7519470&tab=function
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%----------------------------------------------------

%% Reading data
psf=double(loadtiff(psfname));psf=psf/sum(psf(:));
y=double(loadtiff(dataname));maxy=max(y(:));y=y/maxy;
sz=size(y);
if usemask
    if DataTerm==2, error('Kullback-Leibler with mask not implemented'); end;
    sznew=fft_best_dim(sz+2*padsz);
    fprintf('Reconstruction size: [%d %d %d]. If Ok press enter. \n \n',sznew(1),sznew(2),sznew(3));
    pause;
    halfPad=(sznew-sz)/2;
    psf=padarray(psf,halfPad,0,'both');
end

%% Common operators and costs
H=LinOpConv(fftn(fftshift(psf)));                          % Convolution Operator
if usemask
    S=LinOpSelectorPatch(sznew,halfPad+1,sznew-halfPad);   % Selector Operatr if usemak activated
    L2=CostL2(S.sizeout,y);                                % L2 cost function
    LS=L2*S;                                               % Least-Sqaures data term
    pos=CostNonNeg(sznew);                                 % Non-Negativity: Indicator function
    Id=LinOpIdentity(sznew);                               % Identity Operator
else
    L2=CostL2(H.sizeout,y);                                % L2 cost function
    LS=L2*H;                                               % Least-Sqaures data term
    pos=CostNonNeg(sz);                                    % Non-Negativity: Indicator function
    Id=LinOpIdentity(sz);                                  % Identity Operator
    KL=CostKullLeib(sz,y,1e-3);                            % Kullback-Leibler divergence data term
end

%% Deconvolution
% -- Functions definition
if Reg==1
    if usemask
        Freg=CostMixNorm21([sznew,3],4);       % TV regularizer: Mixed Norm 2-1
        Opreg=LinOpGrad(sznew);                % TV regularizer: Operator Gradient
    else
        Freg=CostMixNorm21([sz,3],4);          % TV regularizer: Mixed Norm 2-1
        Opreg=LinOpGrad(sz);                   % TV regularizer: Operator Gradient
    end
elseif Reg==2
    if usemask
        Freg=CostMixNormSchatt1([sznew,6],1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
        Opreg=LinOpHess(sznew);                % Hessian-Shatten: Hessian Operator
    else
        Freg=CostMixNormSchatt1([sz,6],1);     % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
        Opreg=LinOpHess(sz);                   % Hessian-Shatten: Hessian Operator
    end
end

% -- Run the algorithm
for ii=1:length(lamb)
    costOld=Inf;
    if usemask
        xopt=S'*max(0,y);
    else
        xopt=max(0,y);
    end
    % -- ADMM
    if DataTerm==1
        if ~usemask
            Fn={lamb(ii)*Freg,pos};           % Functionals F_n
            Hn={Opreg,Id};                    % Associated operators H_n
            rho_n=[1e-3,1e-3];                % Multipliers rho_n
            ADMM=OptiADMM(LS,Fn,Hn,rho_n,[],MyOutputOpti(1,[],round(maxIt/10)));
        else
            Fn={CostL2(S.sizeout,y)*S,lamb(ii)*Freg,pos}; % Functionals F_n
            Hn={H,Opreg,Id};                    % Associated operators H_n
            rho_n=[1e-1,1e-3,1e-3];             % Multipliers rho_n
            ADMM=OptiADMM([],Fn,Hn,rho_n,[],MyOutputOpti(1,[],round(maxIt/10)));
        end
    else
        Fn={KL,lamb(ii)*Freg,pos};          % Functionals F_n
        Hn={H,Opreg,Id};                    % Associated operators H_n
        rho_n=[1e-1,1e-3,1e-3];             % Multipliers rho_n
        ADMM=OptiADMM([],Fn,Hn,rho_n,[],MyOutputOpti(1,[],round(maxIt/10)));
    end
    ADMM.ItUpOut=round(maxIt/10);
    ADMM.maxiter=maxIt;
    ADMM.run(xopt);
    xopt=ADMM.xopt;
    if usemask
        xopt=S*xopt;
    end
    % -- save
    saveastiff(xopt,[outFolder,'/lamb',num2str(lamb(ii)),'.tif']);
end