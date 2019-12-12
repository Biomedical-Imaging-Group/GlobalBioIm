%----------------------------------------------
% Deconvolution
%
% Description:
%    Deconvolution Pipeline
%
% Generated on 12-Dec-2019 using the GUI of the GlobalBioIm library
% Link to GlobalBioIm: <https://biomedical-imaging-group.github.io/GlobalBioIm/>
%----------------------------------------------
useGPU(0);

%% Load/Read Parameters
% - Operator-1 : SelectorPatch
H1_InputSize = [512 512];
H1_idxmin = [1 1];
H1_idxmax = [256 256];
% - Operator-2 : Conv
H2_InputSize = [512 512];
tmp=load("/home/esoubies/Bureau/GitHub/GlobalBioIm/GUI/psf.mat"); fd=fields(tmp);
H2_PSF = tmp.(fd{1});
% - Data-fidelity : L2
DF_y = double(imread("/home/esoubies/Bureau/GitHub/GlobalBioIm/GUI/data.png"));
% - Algorithm Fwd-Bkwd Splitting
Opt_gam = 1e-1;

%% GPU/CPU converter
H2_PSF = gpuCpuConverter(H2_PSF);
DF_y = gpuCpuConverter(DF_y);

%% Instanciate the Forward Model
% - Operator-1 : SelectorPatch
H1 = LinOpSelectorPatch(H1_InputSize,H1_idxmin,H1_idxmax,1);
% - Operator-2 : Conv
H2 = LinOpConv('PSF',H2_PSF,1,[],'Pad',H2_InputSize,0);

%% Instanciate the Cost function
% - Data-Fidelity : L2
DF = CostL2(H1.sizeout,DF_y);

%% Instanciate and Run the Optimization method
% - Algorithm Fwd-Bkwd Splitting
F = DF*(H1*H2);
P = CostNonNeg(F.sizein);
Opt = OptiFBS(F,P);

