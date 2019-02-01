function buildHessianSchatten(varargin)
%% buildHessianSchatten function
%   build the HessianSchatten norm mexgl files for CostMixNormSchatt1
%
%   You can give as a parameter of this function the path to your GCC
%   compiler. Ex: buildHessianSchatten('/usr/bin/gcc-6')

%     Copyright (C) 2018 F. Soulez ferreol.soulez@epfl.ch
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin==1
    pathToGcc=[' GCC=''',varargin{1},''' '];
else
    pathToGcc='';
end

[mpath,~,~] = fileparts(which('buildHessianSchatten'));
disp('Installing HessianSchatten');
disp('WARNING: depending on your system and compiler,  OPENMP may not be activated leading to slow computation.');
pth = cd;
cd(mpath);
MexOpt= ['-DUSE_BLAS_LIB ' '-DNEW_MATLAB_BLAS ' '-DINT_64BITS '  '-largeArrayDims ' 'COMPFLAGS=''$COMPFLAGS -Wall -mtune=native  -fomit-frame-pointer -O2 -fopenmp '''];
eval(['mex ',pathToGcc,' svd2D_recomp.cpp ',MexOpt]);
eval(['mex ',pathToGcc,' svd2D_decomp.cpp ',MexOpt]);
eval(['mex ',pathToGcc,' svd3D_recomp.cpp ',MexOpt]);
eval(['mex ',pathToGcc,' svd3D_decomp.cpp ',MexOpt]);
cd(pth);
end