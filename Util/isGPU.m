function isGPU()
%--------------------------------------------------------------
% function isGPU()
%
% Returns the current CPU/GPU activation
%
% Seel also: useGPU, isGPU
%
% Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------
global isGPU

if  isempty(isGPU) || isGPU==0
    disp('Execution is made on CPU');
elseif isGPU==1
    disp('Execution is made on GPU with Matlab parrallel Computing Toolbox');
elseif isGPU==2
    disp('Execution is made on GPU with CudaMat');
end
end