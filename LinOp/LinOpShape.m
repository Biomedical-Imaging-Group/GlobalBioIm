classdef LinOpShape <  LinOp
    %% LinOpShape : reshaping operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = LinOpShape(sizein,sizeout)
    % Shape operator
    % Reshape an array of size sizein in a array of size sizeout
    % 
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
    
    end
    methods
        function this = LinOpShape(sizein, sizeout)
            this.name ='LinOp Shape';
            this.iscomplex= false;
            this.isinvertible=true;
            this.issquare = false;
			
			this.norm = 1;
            
            assert(issize(sizein),'The input size sizein should be a conformable  to a size ');
            this.sizein = sizein;
            assert(issize(sizeout),'The input size sizeout should be a conformable  to a size ');
            this.sizeout = sizeout;
            
            assert(prod(sizeout)==prod(sizein),'The inputs sizein and sizeout should be a conformable');
            
        end
        function y = apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = reshape(x, this.sizeout);
        end
        function y = adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y = reshape(x, this.sizein);
        end
        function y = HHt(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y=x;
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y=x;
        end        
        function y = inverse(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
             y = reshape(x, this.sizein);
        end
        function y = adjointInverse(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y = reshape(x, this.sizeout);
        end
    end
end
