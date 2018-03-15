% out=rft2fft(in) : Converts an rft half back to a full sized fft
function out=rft2fft(in, transformDirs)

if nargin < 2
    if ndims(in) > 1
        copydir=2;
    else
        copydir=1;
    end
else
    copydir = find(transformDirs,1,'first');
    if isempty(copydir)
        error('No direction to transform. This can be caused by transforming only singleton dimensions.');
    end
end

subrange = repmat({':'}, 1, ndims(in));  % default for the non-transformed directions
for d=1:ndims(in)
    if transformDirs(d)
        if isa(in,'dip_image')
            subrange(d)={[0,size(in,d)-1:-1:1]};
        else
            subrange(d)={[1,size(in,d):-1:2]};
        end
    end
end
if isa(in,'dip_image') 
    subrange(copydir)={size(in,copydir)-2:-1:1};
else
    subrange(copydir)={size(in,copydir)-1:-1:2};
end

out=cat(copydir,in,conj(in(subrange{:})));
