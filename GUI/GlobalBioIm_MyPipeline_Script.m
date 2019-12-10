%----------------------------------------------
% MyPipeline
%
% Description:
%    Describe your pipeline
%
% Generated on 10-Dec-2019 using the GUI of the GlobalBioIm library
% Link to GlobalBioIm: <https://biomedical-imaging-group.github.io/GlobalBioIm/>
%----------------------------------------------

%% Load/Read Parameters
% - H1 : Grad
H1_InputSize = [];
H1_Index = [];

H1_res = [];
% - H2 : Conv
H2_InputSize = [];
H2_PSF = [];
H2_index = [];

%% Instanciate the Forward Model
% - H1 : Grad
H1 = LinOpGrad(H1_InputSize,index,H1_BC,H1_res);
% - H2 : Conv
H2 = LinOpConv('PSF',H2_PSF,1,H2_index,'Pad',H2_InputSize,0);
% - Complete Forward Model
H = H1*H2;

%% Instanciate the Cost function

%% Instanciate and Run the Optimization method

