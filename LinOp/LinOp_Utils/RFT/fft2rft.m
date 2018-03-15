% out=fft2rft(in) : deletes half of the already ffted data to convert it into an rft
function out=fft2rft(in, transformDirs)

if nargin < 2
    if ndims(in) > 1
        cutdirection=2;
    else
        cutdirection=1;
    end
else
    cutdirection = find(transformDirs,1,'first');
end

if nargin < 2
    if ndims(in) > 1
        cutdirection=2;
    else
        cutdirection=1;
    end
end
%safety unequal sizes
if isa(in,'dip_image') && (mod(size(in,cutdirection),2) && ~size(in,cutdirection)==1) || ((~ isa(in,'dip_image')) && (mod(size(in,1),2) && ~size(in,1)==1))
     error('The fft2rft function only accepts even size along dim 1/2 (matlab/dipImage), as the rift of the result would yield a different size');
end

subrange = repmat({':'}, 1, ndims(in));
if isa(in,'dip_image') 
    subrange(cutdirection)={0:size(in,cutdirection)/2};
else
    subrange(cutdirection)={1:size(in,cutdirection)/2+1};
end
out = in(subrange{:});
