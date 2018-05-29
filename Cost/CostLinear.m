classdef CostLinear < Cost
    % CostLinear: Linear cost function
    % $$C(\\mathrm{x}) := \\mathrm{x}^T\\mathrm{y}$$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param W: weighting :class:`LinOpDiag` object or scalar (default 1)
    %
    % **Example** C=CostL2(sz,y)
    %
    % See also :class:`Map`, :class:`Cost`
    %%    Copyright (C) 2018
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
    
    %% Constructor
    methods
        function this = CostLinear(sz,y)
            this@Cost(sz,y);
            this.isConvex=false;
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
                y=dot(conj(x(:)),this.y(:));
        end
        function g=applyGrad_(this,x)
        	% Reimplemented from parent class :class:`Cost`. 	
                g=this.y;
          
        end
        function y=applyProx_(this,x,alpha)
        	% Reimplemented from parent class :class:`Cost` 
            
                y = x - alpha.*this.y;
            
        end
    end
end
