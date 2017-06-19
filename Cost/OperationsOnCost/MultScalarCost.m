classdef MultScalarCost < Cost
    %% MultScalarCost : Multiply a Cost by a scalar
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % F = MultScalarCost(Cost,s)
    % Multiply the Cost by the scalar s
    %
    % See also Cost
	%
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
        cost;      % Cost
        s;         % scalar factor
    end
    
    methods 
    	%% Constructor
        function this = MultScalarCost(cost,s)
            warning('MultScalarCost deprecated -> Use MulCost or operator * instead');
            this.name='Multiply Cost by Scalar';
			this.cost = cost;
			assert(isscalar(s),'s must be a scalar');
			this.s=s;
			this.isconvex=cost.isconvex; 
			if cost.lip~=-1
				this.lip=cost.lip*s;
			end
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.s*this.cost.eval(x);
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			g=this.s*this.cost.grad(x);
        end
        %% Proximity operator of the Functional
        function y=prox(this,x,alpha)
        	assert(isscalar(alpha),'alpha must be a scalar');
			y = this.cost.prox(x,this.s*alpha);
        end
    end
end
