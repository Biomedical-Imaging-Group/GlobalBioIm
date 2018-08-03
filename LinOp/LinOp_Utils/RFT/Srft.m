function y = Srft(x, Notindex)
%% Srft function
% Recursive function for sliced FFT. Computed the FFT along all dimension
% of x but those indexed by Notindex;
%
% Example:
% y = Srft(x,[3,4])
% will compute 2D FFTs along dims 1 and 2 only such
%  y(:,:,m,n) = fftn(x(:,:,m,n)) for all (m,n)
%
% See also iSfft fftn



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


if numel(Notindex)~=0
    tdirs = ones(1,ndims(x));
    tdirs(Notindex) = 0;  % do NOT transform these directions
    y = rft(x,tdirs);    
else
    y = rft(x);
end
