function installOptimPack(options)
%% installOptimPack function
%   Install OptimPackLegacy and its matlab API from https://github.com/emmt/OptimPackLegacy
%
%   You can give as a parameter of this function the path to your GCC
%   compiler. Ex: installOptimPack('/usr/bin/gcc-6')

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


[mpath,~,~] = fileparts(which('installOptimPack'));
disp('Installing OptimPackLegacy');
pth = cd;
cd(mpath);
websave('OptimPackLegacy.zip','https://github.com/emmt/OptimPackLegacy/archive/master.zip');
unzip('OptimPackLegacy');
matDir = 'OptimPackLegacy-master/matlab/';
srcDir = 'OptimPackLegacy-master/src/';
MexOpt= [ '-largeArrayDims ' '-DUSE_BLAS_LIB ' '-DNEW_MATLAB_BLAS ' '-DINT_64BITS '  '-largeArrayDims ' ,options, ' CFLAGS=''-fPIC  -Wall -mtune=native  -fomit-frame-pointer -O2 ''' ];
CFiles =  [srcDir,'opl_vmlmb.c ',srcDir,'opl_algebra.c ',srcDir,'opl_lnsrch.c ',srcDir,'opl_utils.c ','-IOptimPackLegacy-master/src/ '];
eval(['mex ',matDir,'m_opl_vmlmb_get_reason.c ', CFiles,MexOpt]);
eval(['mex ',matDir,'m_opl_vmlmb_create.c ', CFiles,MexOpt]);
eval(['mex ',matDir,'m_opl_vmlmb_iterate.c ', CFiles,MexOpt]);
eval(['mex  ',matDir,'m_opl_vmlmb_restore.c ', CFiles,MexOpt]);
movefile('OptimPackLegacy-master/matlab/*.m');
delete('makeOptimPack.m');
delete('OptimPackLegacy.zip');
rmdir('OptimPackLegacy-master','s');
cd(pth);

end