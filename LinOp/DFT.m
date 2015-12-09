classdef DFT <  LinOp
    %% DFT : Discrete Fourier operator
    %  Matlab Linear Operator Library 
    % 
    % Example
    % Obj = DFT(pad)
    % Compute the dicrete fourier transfom padded by PAD (optionnal) such
    % y = Obj.Apply(x)  <=> y = fftn(x,pad)
    % 
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp fftn SDFT ifftn
    
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
        function this = DFT(varargin)
            
            this.name ='DFT';
            this.iscomplex= true;   
            this.isinvertible=true;
            this.issquare = true;
            for n=1:length(varargin)
                if n<2
                    if issize(varargin{n})
                        this.pad= varargin{n}  ; 
                    else
                        error('PAD should be conformable to size() output');
                    end
                end
            end
        end
        function y = Apply(this,x)
            this.sizein = size(x);
            this.sizeout = size(x);
            this.N=numel(x);
            if this.unitary               
            y = 1./sqrt(this.N) * fftn(x,this.pad); 
            else
            y =  fftn(x,this.pad); 
            
            end
        end
        function y = Adjoint(this,x)
            this.N=numel(x);
            if this.unitary
            y = sqrt(this.N) * ifftn(x,this.pad); 
            else
            y = this.N * ifftn(x,this.pad); 
            end
        end
        function y = Inverse(this,x)
            y = ifftn(x,this.pad);
        end
        function y = HtH(this,x)%Beware of unitary fft of Matlab :-(
            this.N=numel(x);
            
            if this.unitary
                y = x;
            else
            y =this.N *x;
            end
        end
        function y = AdjointInverse(this,x)%Beware of unitary fft of Matlab :-(
            this.N=numel(x);
            y = 1/this.N *fftn(x,this.pad);
        end
    end
end
