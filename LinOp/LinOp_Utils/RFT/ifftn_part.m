% out=ifftn_part(in, transformDirs) :partial fourier transform as given by the binary vector transformDirs
function out=ifftn_part(in, transformDirs)
if (ndims(in) ~= numel(transformDirs))
    error('fftn_part: Number of transform dimensions has to agree to number of dimensions of input.')
end
for d=1:ndims(in)
    if transformDirs(d)
        in=ifft(in,[],d);
    end
end
out=in;