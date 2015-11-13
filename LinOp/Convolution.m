classdef Convolution <  LinOp
    %% Convolution : Convolution operator
    %  Matlab Linear Operator Library
    %
    % Example:
    % Obj = Convolution(psf, index)
    % Convolution operator with  PSF psf along the dimension
    % indexed in INDEX (all by default)
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp DFT Sfft iSFFT
    
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
    
    properties (SetAccess = protected, GetAccess = public)
        psf
        mtf
        index
        Notindex
        ndms
    end
    methods
        function this = Convolution(psf,index)
            if nargin == 1
                index = [];
            end
            this.name ='Convolution';
            this.isinvertible=false;
            this.issquare = true;
            
            assert(isnumeric(psf),'The psf should be a');
            this.psf = psf;
            if isreal(psf)
            this.iscomplex= false;
            end
            this.sizeout =size( this.psf);
            
            this.sizein = this.sizeout;
            
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.ndms==2 && (this.sizein(2) ==1 || this.sizein(1) ==1)
                this.ndms = 1;
            end
            
            if (~isempty(index))
                assert(isvector(index) && length(index)<= this.ndms && max(index)<= this.ndms,'The index should be a conformable  to sz');
                this.index = index;
                dim = 1:this.ndms;
                Iidx = true(this.ndms,1);
                Iidx(index) = 0;
                this.Notindex = dim(Iidx);
            else
                this.index = 1:this.ndms;
                this.Notindex = [];
            end
            
            this.mtf = Sfft(this.psf, this.Notindex);
            
                if all(this.mtf)
                    this.isinvertible=true;
                else
                    this.isinvertible=false;
                end
            
            
        end
        function y = Apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d]',this.sizein);
            y = iSfft( this.mtf .* Sfft(x, this.Notindex), this.Notindex );
            if (~this.iscomplex)&&isreal(x)
                y = real(y);
            end
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y = iSfft( conj(this.mtf) .* Sfft(x, this.Notindex), this.Notindex );
            if (~this.iscomplex)&&isreal(x)
                y = real(y);
            end
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = iSfft( (real(this.mtf).^2 + imag(this.mtf).^2) .* Sfft(x, this.Notindex), this.Notindex );
            if (~this.iscomplex)&&isreal(x)
                y = real(y);
            end
        end
        function y=Inverse(this,x) % Apply the inverse
            if this.isinvertible
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d]',this.sizein);
            y = iSfft( 1./this.mtf .* Sfft(x, this.Notindex), this.Notindex );       
            else
                error('Operator not invertible');
            end
        end
        function AdjointInverse(this,~) % Apply the inverse
            if this.isinvertible
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y = iSfft( 1./conj(this.mtf) .* Sfft(x, this.Notindex), this.Notindex );
            else
                error('Operator not invertible');
            end
        end
    end
end
