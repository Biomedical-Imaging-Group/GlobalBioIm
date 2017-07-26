function sizeOK = checkSize(x, sz)
% Check if the size of matrix x matches what is specified in sz, returning
% true if so and false if not. Handles the special cases where the length
% of size(x) does not match the length of sz.
%

xSize = size(x);
if length(xSize) < length(sz)
	xSize(end+1:length(sz)) = 1;
elseif length(sz) < length(xSize)
	sz(end+1:length(xSize)) = 1;
end

sizeOK = isequal(xSize,sz);

end