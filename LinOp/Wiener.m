classdef Wiener <  LinOp   
    %% Wiener : Wiener operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = Wiener(LinOp1, alpha1, LinOp2, alpha2,...)
    % Build the wiener linear operator such:
    % Obj*x = (alpha1.*Linop1.HtH +alpha2.*Linop2.HtH)^-1 *( alpha1.*Linop1' *x alpha2.*Linop2' *x)
    %
    % It is only valid if the Linop are diagonal in Fourier domain, namely:
    % Identity(), Convolution(), Fresnel()
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
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
   
    properties (SetAccess = protected,GetAccess = public)
        alpha;
        imtf2;
        Notindex;
        index;
        nlinop;
        ALinOp;
    end
    methods
        function this = Wiener(varargin)
            assert(isa(varargin{1}(1),'LinOp'),'First input should be a LinOp');
            this.name='Wiener';
            this.nlinop=0;
            this.alpha=[];
            this.ALinOp=[];
            this.Notindex=[];
            this.index=[];
            c = 1;
            while c<=length(varargin)
                
                switch(varargin{c}(1).name)
                    case('OneToMany')
                this.ALinOp =cat(2,this.ALinOp,num2cell(varargin{c}.ALinOp));
                 this.alpha =cat(1, this.alpha, varargin{c+1} .* varargin{c}.alpha);
                this.nlinop = this.nlinop + varargin{c}.numLinOp;
                        c = c+2;
                    case('Convolution')
                 this.Notindex =    varargin{c}(1).Notindex;
                 this.index =   varargin{c}(1).index;
                this.ALinOp =cat(2, this.ALinOp,num2cell(varargin{c}));
                 this.alpha =cat(1, this.alpha,varargin{c+1});
                 this.nlinop =  this.nlinop + 1;
                        c = c+2;
                        
                    case('Identity')
                this.ALinOp =cat(2, this.ALinOp,{num2cell(varargin{c})});
                 this.alpha =cat(1, this.alpha,varargin{c+1});
                this.nlinop =  this.nlinop + 1;
                        c = c+2;
                    case('Fresnel')
                this.ALinOp =cat(2, this.ALinOp,num2cell(varargin{c}));
                 this.alpha =cat(2, this.alpha,varargin{c+1});
                 this.nlinop =  this.nlinop + 1;
                        c = c+2;
                    case('DFT')
                this.ALinOp =cat(2, this.ALinOp,num2cell(varargin{c}));
                 this.alpha =cat(2, this.alpha,varargin{c+1});
                 this.nlinop =  this.nlinop + 1;
                        c = c+2;
                    otherwise
                        error('LinOp should only be OneToMany, Convolution, Identity and Fresnel');
                end
                
            end
            
            this.sizein =  this.ALinOp{1}{1}.sizein;
            this.sizeout =  this.ALinOp{1}{1}.sizeout;
            mtf2 =0; 
            for n=1:this.nlinop;
                    switch(this.ALinOp{n}{1}.name)
                    case('Convolution')
                 this.Notindex =    this.ALinOp{n}{1}.Notindex;
                 this.index =    this.ALinOp{n}{1}.index;
                        mtf2 = mtf2 + this.alpha(n) .* ( real(this.ALinOp{n}{1}.mtf).^2 + imag(this.ALinOp{n}{1}.mtf).^2);
                    case('Identity')
                        mtf2 = mtf2 + this.alpha(n) ;
                    case('Fresnel')
                        mtf2 = mtf2 + this.alpha(n) .* ( real(this.ALinOp{n}{1}.F).^2 + imag(this.ALinOp{n}{1}.F).^2);
                    case('DFT')
                        mtf2 = mtf2 + this.alpha(n) .* this.ALinOp{n}{1}.N ;
                    end
            end
            this.imtf2 = 1./mtf2;
        end
        function y = Apply(this,x)
        htx = 0;
            for n=1:this.nlinop;
                    switch(this.ALinOp{n}{1}.name)
                    case('Convolution')
                        htx = htx + this.alpha(n) .* conj(this.ALinOp{n}{1}.mtf).*Sfft(x{n},this.Notindex);
                    case('Identity')
                        htx = htx + this.alpha(n) .* Sfft(x{n},this.Notindex) ;
                    case('Fresnel')
                        htx = htx + this.alpha(n) .*conj(this.ALinOp{n}{1}.F .*fftn(x{n})); 
                    end
            end
            y =   iSfft(this.imtf2 .* htx,this.Notindex);
        end
    end
end