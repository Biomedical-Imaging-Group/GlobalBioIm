classdef LinOpConv <  LinOp
    % LinOpConv: Convolution operator
    %
    % :param mtf: Fourier transform of Point Spread Function
    % :param isReal: if true (default) the result of the convolution should be real
    % :param index: dimensions along which the convolution is performed (the MTF/PSF must have a comformable size
    % :param 'MTF': keyword to provide MTF (default)
    % :param 'PSF': keyword to provide PSF instead of MTF
    % :param 'Centered': if the PSF is centered in the image
    % :param 'Pad': if the PSF must be padded to the size SZ with the value padvalue (default 0) 
    % :param 'useRFT': keyword to use the real-to-half-complex fourier transformation (works with a given 'PSF', not 'MTF')
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % **Example** H=LinOpConv(mtf,isReal,index)
    %
    % **Example** H=LinOpConv('MTF',mtf,isReal,index)
    %
    % **Example** H=LinOpConv('PSF', psf,isReal,index)
    %
    % **Example** H=LinOpConv('PSF', psf,isReal,index,'useRFT')
    %
    % **Example** H=LinOpConv('PSF', psf,isReal,index,'Centered')
    %
    % **Example** H=LinOpConv('PSF', psf,isReal,index,'Pad',sz,padvalue)
    %
    % **Example** H=LinOpConv('PSF', psf,isReal,index,'Centered','Pad',sz,padvalue)
    %
    % **Note** the order of keywords do not matters. For example the two following lines are equivalent
    %
    % H=LinOpConv('PSF', psf,isReal,index,'Pad',sz,padvalue,'useRFT')
    %
    % H=LinOpConv('PSF', psf,isReal,index,'useRFT','Pad',sz,padvalue)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %% GUI-Header
    % GUInotation-C-
    % GUIcall-LinOpConv('PSF',PSF,1,index,'Centered','Pad',InputSize,0)-
    % GUIparam-InputSize-vecInt-[]-Input size of the convolution operator (e.g. [512 512]). If empty, the size of the PSF will be used.
    % GUIparam-PSF-file-[]-Select PSF file
    % GUIparam-index-vecInt-[]-Used to specify the dimensions along which the convolution is performed. If empty (default), all dimensions are used. For instance, setting [1 2] in 3D will instanciate a LinOpConv that applies a 2D convolution to each plane of the input 3D volume.
    
    %%    Copyright (C) 2015
    %     F. Soulez ferreol.soulez@epfl.ch
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

    properties (SetAccess = protected,GetAccess = public)
        mtf;       % Fourier transform of the PSF
        index;     % Dimensions along which the convolution is performed
        isReal;    % true (default) if the result of the convolution should be real
        useRFT=0;  % true if the real-to-half-complex fourier transformation is used rather than the complex-to-complex FFT
                   % Default : false
        Notindex;  % Remaining dimensions
        ndms;      % number of dimensions
    end
    
    %% Constructor
    methods
        function this = LinOpConv(varargin)

            %====  Read arguments
            ispsf = false;
            centering = false;
            pad = false;
            padvalue=0;           
            if isnumeric(varargin{1})
                mtf = varargin{1};
                ndms = ndims(mtf);
                if nargin == 1
                    index =  1:ndms;
                    isReal=true;
                elseif nargin<3
                    isReal = varargin{2};
                    index =  1:ndms;
                elseif nargin==3
                    isReal = varargin{2};
                    index = varargin{3};
                end        
            else switch varargin{1}
                    case('PSF')
                        ispsf = true;
                        psf = varargin{2};
                        assert(isnumeric(psf),'The PSF should be numeric');
                        ndms = ndims(psf);
                    case('MTF')
                        mtf = varargin{2};
                        assert(isnumeric(mtf),'The PSF should be numeric');
                        ndms = ndims(mtf);
                    otherwise
                        error('Unknown keyword.');
                end
                if nargin == 2
                    index =  1:ndms;
                    isReal=true;
                end
                if nargin<4
                    isReal = varargin{3};
                    index =  1:ndms;
                end
                if nargin>3
                    isReal = varargin{3};
                    index = varargin{4};
                    if isempty(index)
                        index= 1:ndms;
                    end
                end
            end
            c=5;
            while c<=length(varargin)
                switch varargin{c}
                    case('Centered')
                        centering = true;
                    case('Pad')
                        if (nargin>c+1) && isnumeric(varargin{c+1})
                            pad = true;
                            sz = varargin{c+1};
                            c = c+1;
                        end
                        if  isscalar(varargin{c+1})
                            padvalue = varargin{c+1};
                            c = c+1;
                        end
                    case('useRFT')
                        this.useRFT= true;
                        assert(isReal==1,'RTF can only be used for real data. isReal needs to be one.');
                        assert(ispsf,'To use RFT a PSF should be given (not MTF)');
                    otherwise
                        error('Unknown keyword.');
                end
                c=c+1;
            end

            if ispsf
                if pad
                    if ~centering
                        for n=1:ndims(psf)
                            if any(index==n)
                                psf = fftshift(psf,n);
                            end
                        end
                        centering=true;
                    end
                    psfsize = size(psf);
                    assert(length(sz)==length(psfsize),['The input size [',num2str(sz),'] and the size of the PSF [',num2str(psfsize),'] must have the same number of dimensions']);
                    padsz= (sz - psfsize)/2.;                    
                    assert(all(padsz>=0),['The input size [',num2str(sz),'] must be larger than the size of the PSF [',num2str(psfsize),']']);
                    for n=1:ndims(psf)
                        if any(index==n)
                            padsize = zeros(size(psfsize));
                            padsize(n) = floor(gather(padsz(n)));
                            psf = padarray(psf,padsize,padvalue,'post');
                            padsize(n) = ceil(gather(padsz(n)));
                            psf = padarray(psf,padsize,padvalue,'pre');
                        end
                    end
                end
                
                if centering
                    for n=1:ndims(psf)
                        if any(index==n)
                            psf = ifftshift(psf,n);
                            
                        end
                    end
                end
                
                this.sizeout =size(psf);
            else
                this.sizeout =size(mtf);
            end

            %===== Set name, dimensions etc...
            this.name ='LinOpConv';
            this.isInvertible=false;
            this.isDifferentiable=true;            
            this.isReal= isReal;          
            this.sizein = this.sizeout;
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.ndms==2 && (this.sizein(2) ==1 || this.sizein(1) ==1)
                this.ndms = 1;
            end           
            if (~isempty(index))
                assert(isvector(index),'Parameter index should be a vector');
                assert(length(index)<= this.ndms,'The lenght of vector index must comply with the length of the input size (i.e. the number of dimensions of the input)');
                assert(max(index)<= this.ndms,'The larger value in vector index should not exceed the number of dimensions of the input (i.e. the length of the input size)');
                this.index = index;
                dim = 1:this.ndms;
                Iidx = true(this.ndms,1);
                Iidx(index) = 0;
                this.Notindex = dim(Iidx);
            else
                this.index = 1:this.ndms;
                this.Notindex = [];
            end

            %====== Set mtf and compute it if psf is given
            if ispsf      
                if this.useRFT
                    this.mtf = Srft(psf, this.Notindex);
                else
                    this.mtf = Sfft(psf, this.Notindex);
                end
            else
                this.mtf = mtf;
            end
            
            if sum(this.mtf(:)==0)==0
                this.isInvertible=true;
            else
                this.isInvertible=false;
            end
            
            % -- Norm of the operator
            this.norm=max(abs(this.mtf(:)));
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if (this.useRFT)
                y = iSrft( this.mtf .* Srft(x, this.Notindex), this.Notindex );
            else
                y = iSfft( this.mtf .* Sfft(x, this.Notindex), this.Notindex );
                if (this.isReal) && isreal(x)
                    y = real(y);
                end
            end
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if (this.useRFT)
                y = iSrft( conj(this.mtf) .* Srft(x, this.Notindex), this.Notindex );
            else
                y = iSfft( conj(this.mtf) .* Sfft(x, this.Notindex), this.Notindex );
                if (this.isReal)&&isreal(x)
                    y = real(y);
                end
            end
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if (this.useRFT)
                y = iSrft( (real(this.mtf).^2 + imag(this.mtf).^2) .* Srft(x, this.Notindex), this.Notindex );
            else
                y = iSfft( (real(this.mtf).^2 + imag(this.mtf).^2) .* Sfft(x, this.Notindex), this.Notindex );
                if (this.isReal)&&isreal(x)
                    y = real(y);
                end
            end
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=this.applyHtH_(x);
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if (this.useRFT)
                    y = iSrft( 1./this.mtf .* Srft(x, this.Notindex), this.Notindex );
                else
                    y = iSfft( 1./this.mtf .* Sfft(x, this.Notindex), this.Notindex );
                    if (this.isReal)&&isreal(x)
                        y = real(y);
                    end
                end
            else
                y = applyInverse_@LinOp(this,x);
            end
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if (this.useRFT)
                    y = iSrft( 1./conj(this.mtf) .* Srft(x, this.Notindex), this.Notindex );
                else
                    y = iSfft( 1./conj(this.mtf) .* Sfft(x, this.Notindex), this.Notindex );
                    if (this.isReal)&&isreal(x)
                        y = real(y);
                    end
                end
            else
                y = applyAdjointInverse_@LinOp(this,x);
            end
        end
        function M = plus_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            if isa(G,'LinOpDiag') && G.isScaledIdentity
                if (this.useRFT)
                    M=LinOpConv('PSF',iSrft(G.diag+this.mtf,this.Notindex),this.isReal,this.index,'useRFT');
                else
                    M=LinOpConv(G.diag+this.mtf,this.isReal,this.index);
                end
            elseif isa(G,'LinOpConv')
                if this.useRFT==G.useRFT
                    if  this.useRFT
                        M=LinOpConv('PSF',iSrft(this.mtf+G.mtf,this.Notindex),this.isReal,this.index,'useRFT');
                    else
                        M=LinOpConv(this.mtf+G.mtf,this.isReal,this.index);
                    end
                else
                    warning('Adding two LinOpConv with different flags useRFT: results with useRFT=false.')
                    if this.useRFT
                        newmtf=G.mtf + Sfft(iSrft(this.mtf,this.Notindex),this.Notindex);
                    else
                        newmtf=this.mtf + Sfft(iSrft(G.mtf,G.Notindex),G.Notindex);
                    end
                    M=LinOpConv(newmtf,this.isReal,this.index);
                end
            else
                M=plus_@LinOp(this,G);
            end
        end
% Not needed anymore : will goes at the level of Map and use sum with a wight (-1) which is ok because of the
% case "isa(H,'LinOpDiag') && H.isScaledIdentity" is makeComposition_ below. 
% Should be removed if everithing goes well without it 
%         function M = minus_(this,G)
%             % Reimplemented from parent class :class:`LinOp`.
%             if isa(G,'LinOpDiag')  && G.isScaledIdentity
%                 M=LinOpConv(this.mtf-G.diag,this.isReal,this.index);
%             elseif isa(G,'LinOpConv')
%                 if this.useRFT==G.useRFT
%                     M=LinOpConv(this.mtf-G.mtf,this.isReal,this.index);
%                 else
%                     warning('Subtracting two LinOpConv with different flags useRFT: results with useRFT=false.')
%                     if this.useRFT
%                         newmtf=Sfft(iSrft(this.mtf,this.Notindex),this.Notindex)-G.mtf;
%                     else
%                         newmtf=this.mtf - Sfft(iSrft(G.mtf,G.Notindex),G.Notindex);
%                     end
%                     M=LinOpConv(newmtf,this.isReal,this.index);
%                 end
%             else
%                 M=minus_@LinOp(this,G);
%             end
%         end
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if (this.useRFT)
                M=LinOpConv('PSF',iSrft(conj(this.mtf),this.Notindex),this.isReal,this.index,'useRFT');
            else
                M=LinOpConv(conj(this.mtf),this.isReal,this.index);
            end
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if (this.useRFT)
                M=LinOpConv('PSF',iSrft(complex(abs(this.mtf).^2),this.Notindex),this.isReal,this.index,'useRFT');
            else
                M=LinOpConv(abs(this.mtf).^2,this.isReal,this.index);
            end
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=this.makeHtH_();
        end
        function M = makeInversion_(this)
            % Reimplemented from parent class :class:`LinOp`
            
            if this.isInvertible
                if (this.useRFT)
                    M=LinOpConv('PSF',iSrft(1./this.mtf,this.Notindex),this.isReal,this.index,'useRFT');
                else
                    M=LinOpConv(1./this.mtf,this.isReal,this.index);
                end
            else
                M=makeInversion_@LinOp(this);
            end
        end
        function G = makeComposition_(this, H)
            % Reimplemented from parent class :class:`LinOp`
            if isa(H, 'LinOpConv')
                if this.useRFT==H.useRFT
                    if  this.useRFT
                        G=LinOpConv('PSF',iSrft(this.mtf.*H.mtf,this.Notindex),this.isReal,this.index,'useRFT');
                    else
                        G=LinOpConv(this.mtf.*H.mtf,this.isReal,this.index);
                    end
                else
                    warning('Composing two LinOpConv with different flags useRFT: results with useRFT=false.')
                    if this.useRFT
                        newmtf=H.mtf.*Sfft(iSrft(this.mtf,this.Notindex),this.Notindex);
                    else
                        newmtf=this.mtf.*Sfft(iSrft(H.mtf,H.Notindex),H.Notindex);
                    end
                    G=LinOpConv(newmtf,this.isReal,this.index);
                end
            elseif isa(H,'LinOpDiag')
                sz=ones(size(H.sizein));   % to consider singleton dimension that are at the end
                sz(1:length(size(H.diag)))=size(H.diag);
                if ( H.isScaledIdentity || (all(sz(this.index)==1)) )
                    if  this.useRFT
                        G = LinOpConv('PSF',iSrft(H*this.mtf,this.Notindex),this.isReal,this.index,'useRFT');
                    else
                        G = LinOpConv(H*this.mtf,this.isReal,this.index);
                    end
                else
                    G = makeComposition_@LinOp(this, H);
                end
            else
                G = makeComposition_@LinOp(this, H);
            end
        end
    end
    
end
