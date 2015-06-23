classdef SumLinOp < LinOp
    %% SumLinop : Summation of linear operator
    %  Matlab Linear Operator Library
    %
    % Obj = SumLinop(LinOp1,LinOp2,alpha1, alpha2)
    % Element wise sum  of LinOps:
    % Obj = alpha1 LinOp1 + alpha2 LinOp2
    %
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
        LinOp1
        LinOp2
        alpha1
        alpha2
    end
    
    methods
        function this = SumLinOp(LinOp1, LinOp2, alpha1, alpha2)
            this.name ='SumLinOp';
            if nargin == 2
                alpha1 = 1;
                alpha2 = 1;
            end
            if nargin == 3
                alpha2 = 1;
            end
            
            
            assert(isa(LinOp1,'LinOp'),'First input should be a LinOp');
            assert(isa(LinOp2,'LinOp'),'Second input should be a LinOp');
            assert(isscalar(alpha1),'Third input should be a scalar');
            assert(isscalar(alpha2),'Fourth input should be a scalar');
            
            this.alpha1 = alpha1;
            this.LinOp1 = LinOp1;
            this.alpha2 = alpha2;
            this.LinOp2 = LinOp2;
            if LinOp1.iscomplex || LinOp2.iscomplex
                this.iscomplex= true;
            else
                this.iscomplex= false;
            end
            
            this.isinvertible=false;
            
            
        end
        
        function y = Apply(this,x) % Apply the operator
            y = this.alpha1 * this.LinOp1.Apply(x) + this.alpha2 * this.LinOp2.Apply(x);
        end
        function y = Adjoint(this,x) % Apply the adjoint
            y = this.alpha1 * this.LinOp1.Adjoint(x) + this.alpha2 * this.LinOp2.Adjoint(x);
        end
    end
end

