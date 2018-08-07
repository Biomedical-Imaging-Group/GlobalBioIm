function y = iSfft(x, Notindex, pad)
%% iSfft function
% Recursive function for sliced inverse FFT. Computed the inverse FFT along all dimension
% of x but those indexed by Notindex;
% 
% Example:
% y = iSfft(x,[3,4]) 
% will compute 2D iFFTs along dims 1 and 2 only such 
%  y(:,:,m,n) = ifftn(x(:,:,m,n} for all (m,n)
%
% See also Sfft fftn

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


sz=size(x);
if nargin < 2, Notindex=[]; end
if nargin < 3 || isempty(pad), pad=sz; end
ndms=length(sz);
index=setdiff(1:ndms,Notindex);
t=pad>=sz;

elems = repmat({':'}, 1,ndms);
if t(index(1))
    y=ifft(x,pad(index(1)),index(1));
else
    y=ifft(x,[],index(1));
    elems{index(1)}=1:pad(index(1));
end
for n=2:length(index)
    if t(index(n))
        y=ifft(y,pad(index(n)),index(n));
    else
        y=ifft(y,[],index(n));
        elems{index(n)}=1:pad(index(n));
    end
end
y=y(elems{:});