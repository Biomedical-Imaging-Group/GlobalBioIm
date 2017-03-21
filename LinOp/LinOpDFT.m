classdef LinOpDFT <  LinOp
    %% LinOpDFT : Discrete Fourier operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = LinOpDFT(sz,pad)
    % Compute the dicrete fourier transfom padded by PAD (optionnal) such
    % y = Obj.apply(x)  <=> y = fftn(x,pad)
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp fftn LinOpSDFT ifftn
    
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
        N
        pad  % size of of the fourier (default all)
        real % 1 if real complex Fourier transform (default 0).  FIXME: NOT IMPLEMENTED
        unitary = false
    end
    methods
        function this = LinOpDFT(sz, pad)
            
            this.name ='Complex DFT';
            this.iscomplex= true;
            this.isinvertible=true;
            this.issquare = true;
            
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein = sz;`
            if nargin>1
                if issize(pad)
                    this.pad= pad  ;
                else
                    error('PAD should be conformable to size() output');
                end
            end
            
        end
        function y = apply(this,x)
            this.sizein = size(x);
            this.sizeout = size(x);
            this.N=numel(x);
            if this.unitary
                y = 1./sqrt(this.N) * fftn(x,this.pad);
            else
                y =  fftn(x,this.pad);
                
            end
        end
        function y = adjoint(this,x)
            this.N=numel(x);
            if this.unitary
                y = sqrt(this.N) * ifftn(x,this.pad);
            else
                y = this.N * ifftn(x,this.pad);
            end
        end
        function y = inverse(this,x)
            y = ifftn(x,this.pad);
        end
        function y = HtH(this,x)%Beware of unitary fft of Matlab
            this.N=numel(x);
            
            if this.unitary
                y = x;
            else
                y =this.N *x;
            end
        end
        function y = adjointInverse(this,x)%Beware of unitary fft of Matlab
            this.N=numel(x);
            y = 1/this.N *fftn(x,this.pad);
        end
    end
end
