classdef LinOpComposition < MapComposition
    %% LinOpComposition : Composition of LinOps
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = MulLinOp(LinOp1,LinOp2)
    % Multiplication of LinOps:
    % Obj = LinOp1 * LinOp2
    %
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
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
		isHTH = 0;
		isHHt = 0;
    end
    
    %% Constructor
    methods
        function this = LinOpComposition(H1,H2)
            this@MapComposition(H1,H2);
            this.name='LinOpComposition';          
            assert((isa(H1,'LinOp') || isscalar(H1)),'H1 have to be a LinOp object or a scalar');
            assert(isa(H2,'LinOp'),'H2 have to be a LinOp');
            
            %TODO computing the IsHtH and isHHt
        end
    end
        
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - applyInverse_(this,y)
    methods (Access = protected)
		function y = apply_(this,x) % apply the operator
			if this.isHTH
				y = this.LinOp2.HtH(x);
			elseif this.isHHt
				y = this.LinOp1.HHt(x);
			elseif this.isnum
				y = this.LinOp1.*this.LinOp2.apply(x);
			else
				y = this.LinOp1.apply( this.LinOp2.apply(x));
			end
		end
		function y = applyAdjoint_(this,x) % apply the adjoint
			if this.isHTH || this.isHHt
				y = this.apply(x); % because self-adjoint
			elseif this.isnum
				y = this.LinOp2.adjoint(this.LinOp1.*x);
			else
				y = this.LinOp2.adjoint(this.LinOp1.adjoint(x));
			end
		end
		function y = applyHtH_(this,x)
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isnum
				y = this.LinOp2.adjoint(this.LinOp1.^2.*( this.LinOp2.apply(x)));
			else
				y = this.LinOp2.adjoint(this.LinOp1.HtH( this.LinOp2.apply(x)));
			end
        end
		function y = applyHHt_(this,x)
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isnum
				y = this.LinOp1.*(this.LinOp2.HHt( this.LinOp1.*x));
			else
				y = this.LinOp1.apply(this.LinOp2.HHt( this.LinOp1.adjoint(x)));
			end
		end
		function y = applyInverse_(this,x) % apply the inverse
			if this.isinvertible
				y = this.LinOp2.inverse(this.LinOp1.inverse(x));
			else
				error('Operator not invertible');
			end
		end
        function y = applyAdjointInverse_(this,x) % apply the inverse
            if this.isinvertible             
                y = this.LinOp2.adjointInverse(this.LinOp1.adjointInverse(x));
            else
                error('Operator not invertible');
            end
        end
    end
end

