% out=rift(in) : Simulates an inverse rft with a complex valued half complex array as input
function out=rift(in,transformDirs)
global FFTW_Threads
global FFTW_WisdomFilename
global FFTW_ForceWisdom   % if this exists and is set to 1, the wisdom mechanism is invoked in every call
if isempty(FFTW_Threads)
    FFTW_Threads=8;
end
if isempty(FFTW_WisdomFilename)
    mp=userpath();
    if mp(end)==';' || mp(end)==':'    % Windows and Linux
        mp=mp(1:end-1);
    end
    FFTW_WisdomFilename=[mp filesep 'FFTW_wisdom.txt'];
end

wasdip=isa(in,'dip_image');
if nargin <2
    transformDirs=ones(1,ndims(in)) .* (size(in) > 1);
else
    transformDirs=transformDirs .* (size(in) > 1);
end

if exist('fftw_rft','file')
    if isa(in,'dip_image')
        tmp=transformDirs(1);transformDirs(1)=transformDirs(2);transformDirs(2)=tmp;
    end
    TDir=-1; % inverse transform
    if ~isempty(FFTW_ForceWisdom) && FFTW_ForceWisdom
        TDir=TDir*2;
    end
    if isempty(FFTW_WisdomFilename)
        out=fftw_rft(single(in),TDir,transformDirs,FFTW_Threads);
    else
        out=fftw_rft(single(in),TDir,transformDirs,FFTW_Threads,FFTW_WisdomFilename);
    end
    Fac = transformDirs .* size(out);
    Fac(Fac==0)=[]; Fac=prod(Fac);
    if isa(in,'dip_image')
        out=dip_image(out/sqrt(Fac));
    else
        out = out/Fac;
    end
    return;
else

in=double(in);
in=rft2fft(in,transformDirs);
if wasdip
    if nargin < 2
        out=real(ifftn(in))*sqrt(prod(size(in)));
    else
        if wasdip && numel(transformDirs) > 1
            tmp=transformDirs(1);transformDirs(1)=transformDirs(2);transformDirs(2)=tmp;
        end
        out=real(ifftn_part(in,transformDirs));
        validSizes=size(in) .* transformDirs; validSizes(validSizes==0)=[];
        out=out*sqrt(prod(validSizes));
    end
    out=dip_image(out);
else
    if nargin < 2
        out=real(ifftn(in));
    else
        out=real(ifftn_part(in,transformDirs));
    end
end
end