classdef LinOpCpx <  LinOp
    %% LinOpCpx : Complex representation  operator
    %  Matlab Linear Operator Library
    %
    % Example:
    % Obj = LinOpCpx(sz)
    %
    % Build the operator transforming a complex vector as a 2D real vector with [Real, Im] the imaginary part of the input vector
    %
    % HANDLE WITH CARE ; it is not exactly a linear operator as the dot
    % product is not well defined
    % Please refer to the LinOp superclass for documentation
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
        nbDim
    end
    methods
        function this = LinOpCpx(sz)
            this.name ='LinOp Complrx';
            this.issquare = true;
            this.iscomplex = false;
            this.isinvertible = true;
            
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizeout = [sz 2];
            this.sizein = sz;
            this.nbDim = numel(sz);
            
            
        end
        function y = apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein); 
            y = cat(this.nbDim+1, real(x),imag(x));
        end
        function y = adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizein); 
            y = reshape(x,[],2);
            y = complex(y(:,1),y(:,2));
            y = reshape(y,this.sizein);
        end
        function y = HtH(~,x)
            y =x;
        end
        function y = HHt(~,x)
            y =x;
        end
        
        function y = inverse(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizein); 
            y = complex(x(:,:,1),x(:,:,2));
        end
    end
end
