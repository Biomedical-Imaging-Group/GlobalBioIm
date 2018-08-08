function y = Sfft(x, Notindex,pad)
%% Sfft function
% Recursive function for sliced FFT. Computed the FFT along all dimension
% of x but those indexed by Notindex;
%
% Example:
% y = Sfft(x,[3,4])
% will compute 2D FFTs along dims 1 and 2 only such
%  y(:,:,m,n) = fftn(x(:,:,m,n)) for all (m,n)
%
% See also iSfft fftn



%     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch, E. Soubies emmanuel.soubies@epfl.ch
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
    y=fft(x,pad(index(1)),index(1));
else
    y=fft(x,[],index(1));
    elems{index(1)}=1:pad(index(1));
end
for n=2:length(index)
    if t(index(n))
        y=fft(y,pad(index(n)),index(n));
    else
        y=fft(y,[],index(n));
        elems{index(n)}=1:pad(index(n));
    end
end
y=y(elems{:});