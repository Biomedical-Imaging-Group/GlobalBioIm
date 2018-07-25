function useGPU(v)
%--------------------------------------------------------------
% function useGPU(v)
% function that allows to enable or disable GPU computation
%   v=0 : no use of gpu 
%   v=1 : uses the parrallel computing toolbox of Matlab (gpuArray) [1]
%   v=2 : uses the CudaMat Library [2]
%
% [1] https://ch.mathworks.com/help/distcomp/gpu-computing-in-matlab.html
% [2] https://nanoimaging.de/research/software/cudamat/#Related_Software
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------

global isInitCuda isGPU cuda_enabled gpuDev; isGPU=v;
if  ~isempty(cuda_enabled) && isGPU~=2
    disableCuda();
    cuda_enabled=[];
end
if ~isempty(gpuDev) && isGPU~=1
    reset(gpuDev);
    gpuDev=[];
end
if isGPU==1 && isempty(gpuDev)
    gpuDev=gpuDevice();
elseif isGPU==2 && isempty(cuda_enabled)
    assert(exist('initCuda')>0,'useGPU(2) requires CudaMat which seems to be not installed or not on Matlab path. See https://nanoimaging.de/research/software/cudamat/');
    if isempty(isInitCuda)
        isInitCuda=1;
        initCuda();
    end
    enableCuda();
end
end