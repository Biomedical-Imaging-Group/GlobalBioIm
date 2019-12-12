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
H2_index = [];
% - Data-fidelity : L2
DF_y = double(imread("/home/esoubies/Bureau/GitHub/GlobalBioIm/GUI/data.png"));
DF_Weight = [];
% - Regularization-1 : Total-Variation
R1_lambda = 0.2;
R1_index = [];
R1_BC = 'circular';
R1_res = [];

%% Instanciate the Forward Model
% - Operator-1 : SelectorPatch
H1 = LinOpSelectorPatch(H1_InputSize,H1_idxmin,H1_idxmax,1);
% - Operator-2 : Conv
H2 = LinOpConv('PSF',H2_PSF,1,H2_index,'Pad',H2_InputSize,0);

%% Instanciate the Cost function
% - Data-Fidelity : L2
DF = CostL2(H1.sizeout,DF_y,DF_Weight);
% - Regularization-1 : Total-Variation
OpReg1 = LinOpGrad(H2.sizein,R1_index,R1_BC,R1_res);
CostReg1 = R1_lambda * CostMixNorm21(OpReg1.sizeout,length(OpReg1.sizeout));

%% Instanciate and Run the Optimization method

