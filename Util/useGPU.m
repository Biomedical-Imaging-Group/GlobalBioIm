function useGPU(v)
%--------------------------------------------------------------
% function useGPU(v)
% function that allows to enable or disable GPU computation
%   v=0 : no use of gpu 
%   v=1 : uses the parrallel computing toolbox of Matlab (gpuArray) [1]
%   v=2 : uses the CudaMat Library [2]
%
% [1]
% [2]
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------

global isGPU cuda_enabled; isGPU=v;
if isGPU==2
    assert(exist('initCuda')>0,'useGPU(2) requires CudaMat which seems to be not installed or not on Matlab path. See https://nanoimaging.de/research/software/cudamat/');
    initCuda();enableCuda();
else
    if ~isempty(cuda_enabled)
        disableCuda();
    end
end
end