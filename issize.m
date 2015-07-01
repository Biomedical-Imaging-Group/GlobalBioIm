function out = issize(in)
%% ISSIZE function
% Determine whether the input describe a size (nameley a vector of strictly
% positive integers).
%

%     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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

if isnumeric(in)
    if isvector(in) && all(mod(in,1)==0) && all(in>0)
        out=true;
    else
        out=false;
    end
else
    out = false;
end
end

