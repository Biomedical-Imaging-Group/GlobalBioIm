%----------------------------------------------
% Deconvolution
%
% Description:
%    Deconvolution Pipeline
%
% Generated on 17-Dec-2019 using the GUI of the GlobalBioIm library
% Link to GlobalBioIm: <https://biomedical-imaging-group.github.io/GlobalBioIm/>
%----------------------------------------------
useGPU(0);

%% Load/Read Parameters
% - Operator-1 : SelectorPatch
H1_InputSize = [300 300];
H1_idxmin = [1 1];
H1_idxmax = [256 256];
% - Operator-2 : Conv
H2_InputSize = [300 300];
tmp=load("/home/esoubies/Bureau/GitHub/GlobalBioIm/GUI/psf.mat"); fd=fields(tmp);
H2_PSF = tmp.(fd{1});
% - Data-fidelity : L2
DF_y = double(imread("/home/esoubies/Bureau/GccitHub/GlobalBioIm/GUI/data.png"));
% - Regularization-1 : Smooth-Total-Variation
R1_lambda = 0.3;
% - Algorithm Fwd-Bkwd Splitting
Opt_gam = 1;
Opt_fista = 1;
Opt_TolCost = 1e-4;
Opt_TolStep = 1e-4;
% - Path to save results
resultPath = '/home/esoubies/Bureau/GitHub/GlobalBioIm/GUI/deconv';

%% GPU/CPU converter
H2_PSF = gpuCpuConverter(H2_PSF);
DF_y = gpuCpuConverter(DF_y);

%% Instanciate the Forward Model
% - Operator-1 : SelectorPatch
H1 = LinOpSelectorPatch(H1_InputSize,H1_idxmin,H1_idxmax,1);
% - Operator-2 : Conv
H2 = LinOpConv('PSF',H2_PSF,1,[],'Centered','Pad',H2_InputSize,0);

%% Instanciate the Cost function
% - Data-Fidelity : L2
DF = CostL2(H1.sizeout,DF_y);
% - Regularization-1 : Smooth-Total-Variation
OpReg1 = LinOpGrad(H2.sizein);
CostReg1 = R1_lambda * CostHyperBolic(OpReg1.sizeout,[],length(OpReg1.sizeout));

%% Instanciate and Run the Optimization method
% - Algorithm Fwd-Bkwd Splitting
F = DF*(H1*H2) + CostReg1*(OpReg1);
P = CostNonNeg(F.sizein);
Opt = OptiFBS(F,P);
Opt.gam = Opt_gam;
Opt.fista = Opt_fista;
Opt.OutOp = OutputOpti(1,round(Opt.maxiter/10));
Opt.ItUpOut = round(Opt.maxiter/10);
Opt.CvOp = TestCvgCombine('CostRelative',Opt_TolCost, 'StepRelative',Opt_TolStep);
Opt.run(zeros(Opt.cost.sizein));

%% Display and Save Results
imdisp(Opt.xopt,['Deconvolution-Result'],1);
save([resultPath,'_OptiCell'],'Opt');
imwrite(uint8(Opt.xopt/max(Opt.xopt(:))*255),[resultPath,'.png']);

