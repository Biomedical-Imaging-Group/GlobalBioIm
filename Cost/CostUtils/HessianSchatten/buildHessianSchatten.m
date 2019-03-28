function buildHessianSchatten(options)
%% buildHessianSchatten function
%   build the HessianSchatten norm mexgl files for CostMixNormSchatt1
%
%   You can give as a parameter of this function the path to your GCC
%   compiler. Ex: buildHessianSchatten('GCC=/usr/bin/gcc-6')

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
if nargin==0
    options=[];
end

disp('Installing HessianSchatten');
get_architecture;
if linux
   options = [ options, ' CXXFLAGS='' -fopenmp ''',' LDFLAGS=''$LDFLAGS -fopenmp '''];
else
    disp('On your system and compiler,  OPENMP is desactivated leading to slow computation. This can be tuned using the options parameter:');
    disp('Example: options =  CXXFLAGS=  -fopenmp ');
end

[mpath,~,~] = fileparts(which('buildHessianSchatten'));
pth = cd;
cd(mpath);
MexOpt= ['-DUSE_BLAS_LIB ' '-DNEW_MATLAB_BLAS ' '-DINT_64BITS '  '-largeArrayDims ' ,options,  ' CXXFLAGS=''$CXXFLAGS -fPIC -Wall -mtune=native  -fomit-frame-pointer -O2  '''  ' LDFLAGS=''$LDFLAGS '''];
eval(['mex ',' svd2D_recomp.cpp ',MexOpt]);
eval(['mex ',' svd2D_decomp.cpp ',MexOpt]);
eval(['mex ',' svd3D_recomp.cpp ',MexOpt]);
eval(['mex ',' svd3D_decomp.cpp ',MexOpt]);
cd(pth);
end