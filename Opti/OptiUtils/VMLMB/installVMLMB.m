function installVMLMB()
%% installOptimPack function
%   Install VMLMB for Matlab from https://github.com/emmt/VMLMB
%

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

[mpath,~,~] = fileparts(which('installVMLMB'));
disp('Installing VMLMB for Matlab');
pth = cd;
cd(mpath);
websave('VMLMB.zip','https://github.com/emmt/VMLMB/archive/refs/heads/main.zip');
unzip('VMLMB.zip');
movefile('VMLMB-main/*');
rmdir('VMLMB-main','s');
addpath(genpath('matlab/'));
delete('VMLMB.zip');
% Remove VMLMB package in the other languages
rmdir('python','s');
rmdir('yorick','s');

cd('matlab/src/')
files = dir('.');
for ii = 3:length(files)
    if isunix || ismac
        system(['sed -i ''s/zeros/zeros_/g'' ',files(ii).name]);
        system(['sed -i ''s/ones/ones_/g'' ',files(ii).name]);
    elseif ispc
        system(['get-content ',files(ii).name,' | %{$_ -replace "zeros","zeros_"}']);
        system(['get-content ',files(ii).name,' | %{$_ -replace "ones","ones_"}']);
    end
end
cd(pth);
end