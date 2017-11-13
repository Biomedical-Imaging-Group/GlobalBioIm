clear;close all; clc;

%% Parameters
% -- General 
dataname='~/Desktop/InvPbLib/Example/DeconvEx/3D/Y.tif';   % file name data image
psfname='~/Desktop/InvPbLib/Example/DeconvEx/3D/H.tif';    % file name psf
outFolder='~/Desktop/InvPbLib/Example/DeconvEx/3D';
% -- Deconvolution
lamb=1e-5;        % Hyperparameter for initial deconvolution
maxIt=20;         % Max iterations
Reg=1;            % 1 for TV, 2 for Hessian-Schatten 
DataTerm=1;       % 1 for LS, 2 for KL
usemask=false;  
padsz=[0,0,64];

%% Run script
run ../../../setMatlabPath         
run DeconvScript