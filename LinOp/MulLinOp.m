classdef MulLinOp < LinOp
    %% MulLinop : Multiplication of linear operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = SumLinop(LinOp1,LinOp2)
    % Multiplication of LinOps:
    % Obj = LinOp1 * LinOp2
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
        LinOp1
        LinOp2
    end
    
    methods
        function this = MulLinOp(LinOp1, LinOp2)
            
            if isnumeric(LinOp1)
                if isscalar(LinOp1)
             LinOp1 = Scaling(LinOp1);
                else
             LinOp1 = Diagonal(LinOp1);                    
                end
            end
            assert(isa(LinOp1,'LinOp'),'MulLinOp: First input should be a LinOp');
            assert(isa(LinOp2,'LinOp'),'MulLinOp: Second input should be a LinOp');
            
            
            this.name ='MulLinOp';
            this.LinOp1 = LinOp1;
            this.LinOp2 = LinOp2;
            assert(isempty(LinOp1.sizein) || isequal(LinOp1.sizein, LinOp2.sizeout),'size of LinOp not conformable');
            this.sizein = LinOp2.sizein;
            this.sizeout = LinOp1.sizeout;
            if LinOp1.iscomplex || LinOp2.iscomplex
                this.iscomplex= true;
            else
                this.iscomplex= false;
            end
            
            if LinOp1.isinvertible && LinOp2.isinvertible
                this.isinvertible= true;
            else
                this.isinvertible= false;
            end
            if LinOp1.issquare && LinOp2.issquare
                this.issquare= true;
            else
                this.issquare= false;
            end
            
            
        end
        
        function y = Apply(this,x) % Apply the operator
            y = this.LinOp1.Apply( this.LinOp2.Apply(x));
        end
        function y = Adjoint(this,x) % Apply the adjoint
            y = this.LinOp2.Adjoint(this.LinOp1.Adjoint(x));
        end
        function y = HtH(this,x)
            y = this.LinOp2.Adjoint(this.LinOp1.HtH( this.LinOp2.Apply(x)));
        end
        function y = HHt(this,x)
            y = this.LinOp1.Apply(this.LinOp2.HHt( this.LinOp1.Adjoint(x)));
        end
        function y = Inverse(this,x) % Apply the inverse
            if this.isinvertible
                y = this.LinOp2.Inverse(this.LinOp1.Inverse(x));
            else
                error('Operator not invertible');
            end
        end
        function y = AdjointInverse(this,x) % Apply the inverse
            if this.isinvertible             
                y = this.LinOp1.AdjointInverse(this.LinOp1.AdjointInverse(x));
            else
                error('Operator not invertible');
            end
        end
    end
end

