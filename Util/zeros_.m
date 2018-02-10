function res=zeros_(varargin)
%--------------------------------------------------------------
% function zeros_(varargin)
%
% Generate a matrix of zeros on the CPU or directly on the GPU if
% GPU computation is activated .
%
% See also: useGPU, ones_
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------
global isGPU;

if isempty(isGPU) || isGPU==0
    % No GPU
    res=zeros(varargin{:});
elseif isGPU==1
    % gpuArray of Matlab
    res = zeros(varargin{:},'double','gpuArray');
elseif isGPU==2
    % CudaMat
    res=zeros_cuda(varargin{:});
end
end