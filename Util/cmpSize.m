function sizeOK = cmpSize(sz1, sz2)
% Check if the size sz1 matches sz2, returning
% true if so and false if not. Handles the special cases where the length
% of sz1 does not match the length of sz.
%

if length(sz1) < length(sz2)
	sz1(end+1:length(sz2)) = 1;
elseif length(sz2) < length(sz1)
	sz2(end+1:length(sz1)) = 1;
end

sizeOK = isequal(sz1,sz2);

end