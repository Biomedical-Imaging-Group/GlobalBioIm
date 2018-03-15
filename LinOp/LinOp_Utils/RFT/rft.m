% out=rft(in) : Simulates an rft with a dipimage (or matlab array) as an input object
function out=rft(in, transformDirs)
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

if nargin <2
    transformDirs=ones(1,ndims(in)) .* (size(in) > 1);
else
    transformDirs=transformDirs .* (size(in) > 1);
end

if exist('fftw_rft','file')
    if isa(in,'dip_image')
        Fac = transformDirs .* size(in);
        Fac(Fac==0)=[]; Fac=prod(Fac);
        tmp=transformDirs(1);transformDirs(1)=transformDirs(2);transformDirs(2)=tmp;
    end
    TDir=1;
    if ~isempty(FFTW_ForceWisdom) && FFTW_ForceWisdom
        TDir=TDir*2;
    end
    if isempty(FFTW_WisdomFilename)
        out=fftw_rft(single(in),TDir,transformDirs,FFTW_Threads);
    else
        out=fftw_rft(single(in),TDir,transformDirs,FFTW_Threads,FFTW_WisdomFilename);
    end
    if isa(in,'dip_image')
        out=dip_image(out/sqrt(Fac));
    end
    return;
else
    
firstTrans = find(transformDirs,1,'first');

%safety unequal sizes
if mod(size(in,firstTrans),2) ~= 0
     error('The rft function only accepts even size along dim 1/2 (matlab/dipImage), as the rift of the result would yield a different size');
end
% if ndims(in) > 3
%     error('rft only defined for 2d and 3d arrays.');
% end

% in=readim('orka'); %test object (display the ft in log stretch)
if nargin <2
    b=fftn(double(in)); % dip_image Fourier Transform: whole plane
    transformDirs2=transformDirs;
else
    transformDirs2=transformDirs;
    if isa(in,'dip_image') && numel(transformDirs) > 1
         tmp=transformDirs2(1);transformDirs2(1)=transformDirs2(2);transformDirs2(2)=tmp;
    end
    b=fftn_part(double(in),transformDirs2); % partial Fourier transformations
end

if isa(in,'dip_image') 
    if nargin < 1
        b=b/sqrt(prod(size(b)));
    else
        validSizes=size(b) .* transformDirs2;
        validSizes(validSizes==0)=[];
        b=b/sqrt(prod(validSizes));
    end
end

out=fft2rft(b,transformDirs2);
if isa(in,'dip_image')
    out=dip_image(out);
end

end
