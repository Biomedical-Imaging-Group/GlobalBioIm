function out=gpuCpuConverter(in)
%--------------------------------------------------------------
% function res=gpuCpuConverter(in)
%
% Convert the input variable to the right type according to the 
% use (or not) of GPU.
%
% Seel also: useGPU, isGPU
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------

global isGPU

if isempty(isGPU) || isGPU==0
    % No GPU
    out=double(in);
elseif isGPU==1
    % gpuArray of Matlab
    out=gpuArray(double(in));
elseif isGPU==2
    % CudaMat
    out=cuda(double(in));
end
end