classdef Inverse < LinOp
    %% Inverse : overload of Inverse function for LinOp
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = Inverse(LinOp)
    % Obj is the Inverse of the LinOp
    %
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
        TLinOp     % linop
    end
    
    methods 
        function this = inverse(TLinOp)
            this.name ='inverse';
            
            
            assert(isa(TLinOp,'LinOp'),'Input should be a  LinOp');
            assert(TLinOp.isinvertible,'Input should be a  invertible');
            this.TLinOp = TLinOp;
            this.iscomplex= this.TLinOp.iscomplex;
            this.isinvertible=this.TLinOp.isinvertible;
            this.sizein =  this.TLinOp.sizein;
            this.sizeout =  this.TLinOp.sizeout;
            
          end
        
        function y = apply(this,x) % apply the operator
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);            
            y =this.TLinOp.inverse(x);
        end
        function y = adjoint(this,x) % apply the adjoint
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
            y =this.TLinOp.adjointInverse(x);
        end
        function y = inverse(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
            y =this.TLinOp.apply(x);
        end
        function y = adjointInverse(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);   
            y =this.TLinOp.adjoint(x); 
        end
    end
end

