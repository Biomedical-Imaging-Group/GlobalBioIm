function res=ones_(varargin)
%--------------------------------------------------------------
% function res=ones_(varargin)
%
% Generate a matrix of ones on the CPU or directly on the GPU if
% GPU computation is activated .
%
% See also: useGPU, zeros_
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------
global isGPU;

if isempty(isGPU) || isGPU==0
    % No GPU
    res=ones(varargin{:});
elseif isGPU==1
    % gpuArray of Matlab
    res = ones(varargin{:},'double','gpuArray');
elseif isGPU==2
    % CudaMat
    res=ones_cuda(varargin{:});
end
end