classdef CostL2 < Cost
    % CostL2: Weighted L2 norm cost function
    % $$C(\\mathrm{x}) := \\frac12\\|\\mathrm{x} - \\mathrm{y}\\|^2_W $$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param W: weighting :class:`LinOpDiag` object or scalar (default :class:`LinOpIdentity`)
    %
    % **Example** C=CostL2(sz,y,wght)
    %
    % See also :class:`Map`, :class:`Cost`,

    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch 
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
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        W=1;    % weight matrix      
    end

    %% Constructor
    methods
        function this = CostL2(sz,y,wght)
            if nargin<2, y=0; end
            this@Cost(sz,y);
            this.name='CostL2';
            if nargin==3
                assert((isnumeric(wght) && isscalar(wght))||isa(wght,'LinOpDiag'),'weight WGHT must be scalar or Diagonal LinOp');
                this.W=wght;
            end
            if isnumeric(this.W)
                this.lip=this.W;
            else
                this.lip=this.W.norm;
            end
            this.isConvex=true;
            this.isDifferentiable=true;
        end        
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,x,alpha)
	methods (Access = protected)
        function y=apply_(this,x)
        	% Reimplemented from parent class :class:`Cost`.       
            if (isscalar(this.y)&&(this.y==0))
                wr=this.W*(x);
                y=0.5*dot(x(:),wr(:));
            else
                wr=this.W*(x-this.y);
                y=0.5*dot(x(:)-this.y(:),wr(:));
            end
        end
        function g=applyGrad_(this,x)
        	% Reimplemented from parent class :class:`Cost`.
        	% $$ \\nabla C(\\mathrm{x}) = \\mathrm{W (x - y)} $$
        	% It is L-Lipschitz continuous with \\( L \\leq \\|\\mathrm{W}\\|\\).      	
            if(isscalar(this.y)&&(this.y==0))
                g=this.W*x;
            else
                g=this.W*(x-this.y);
            end
        end
        function y=applyProx_(this,x,alpha)
        	% Reimplemented from parent class :class:`Cost`        	
            if  isscalar(this.W)
                y=(x+alpha*this.W*this.y)./(1+alpha.*this.W);
            else
                y=(x+alpha*this.W*this.y)./(1+alpha.*this.W.diag);
            end
        end
        function M=makeComposition_(this,G)
            % Reimplemented from parent class :class:`Cost`  
            M = CostL2Composition(this,G);
        end
    end
end
