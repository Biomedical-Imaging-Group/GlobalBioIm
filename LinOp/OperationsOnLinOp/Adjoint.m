classdef Adjoint < LinOp
    %% Adjoint : overload of adjoint function for LinOp
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = Adjoint(LinOp)
    % Obj is the adjoint of the LinOp
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
        function this = Adjoint(TLinOp)
            this.name ='Adjoint';           
            assert(isa(TLinOp,'LinOp'),'Input should be a  LinOp');
            this.TLinOp = TLinOp;
            this.isComplex= this.TLinOp.isComplex;
            this.isInvertible=this.TLinOp.isInvertible;
            this.sizein =  this.TLinOp.sizeout;
            this.sizeout =  this.TLinOp.sizein;			
			this.norm = this.TLinOp.norm;        
          end
        
        function y = apply(this,x) % apply the operator         
            y =this.TLinOp.adjoint(x);
        end
        function y = adjoint(this,x) % apply the adjoint
            y =this.TLinOp.apply(x);
        end
        function y = HtH(this,x)
            y =this.TLinOp.HHt(x);
        end
        function y = HHt(this,x)
            y =this.TLinOp.HtH(x);
        end
        function y = inverse(this,x)
            y =this.TLinOp.adjointInverse(x);
        end
        function y = adjointInverse(this,x)
            y =this.TLinOp.adjointInverse(x); 
        end
    end
end

