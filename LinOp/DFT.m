classdef DFT <  LinOp
    %% DFT : Discrete Fourier operator
    %  Matlab Linear Operator Library 
    %
    % Obj = FourierTransform(pad,dim, real)
    % Discrete Fourier transform operator
    % 
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % see also LinOp
    
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
        pad  % size of of the fourier (default all)
        dim  % Fourier transform will be apply across the dimension DIM (default all). NOT IMPLEMENTED
        real % 1 if real complex Fourier transform (default 0). NOT IMPLEMENTED
    end
    methods
        function this = Fourier(varargin)
            
            this.name ='FourierTransform';
            this.iscomplex= true;   
            this.isinvertible=true;
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
            y = fftn(x,this.pad);
        end
        function y = Adjoint(this,x)
            y = ifftn(x,this.pad);
        end
        function y = Inverse(this,x)
            y = ifftn(x,this.pad);
        end
        function y = Gram(~,x)
            y =x;
        end
        function y = AdjointInverse(this,x)
            y = fftn(x,this.pad);
        end
    end
end
