function mexCheckNCompile(filename, varargin)
%% mexCheckNCompile function
% Determine whether the mex function FILENAME exist. If not it compiles it
% with the option varargin

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

if exist(filename)~=3
    disp([filename,'is not compiled']);
    [mpath,name,~] = fileparts(which(filename));
    pth = cd;
    cd(mpath);
    if nargin==1
        varargin={''};
    end
    mex(varargin{:},[name,'.c'])
    cd(pth);
end
end