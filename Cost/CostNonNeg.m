classdef CostNonNeg < Cost
    %% CostNonNeg : Non negativity indicator functionnal
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implement the indicator over positive vector function:
    % $$ \phi(Hx) = 0 \textrm{ if } Hx \ge 0 textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    % -- Example
    % F = CostNonNeg();
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Cost, LinOp
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
    
    methods 
    	%% Constructor
        function this = CostNonNeg(H)   
        	 if nargin==0 || isempty(H)
        	 	H=LinOpIdentity();
        	 end
        	 this.set_H(H);
        	 this.name='Func NonNegativity';
			 
			 if isa(this.H,'LinOpIdentity')
				 this.isconvex = true;
			 end
			 
    	end
    	%% Evaluation of the Functional
		function y=eval(this,x)
			y = 0; % We should put -> (norm(min(this.H*x,0.))>0)*realmax; But for some algorithm there is small negative residuals
			% which would cause unreadability of the evolution of the cost function along iterates
			y = inf * any(reshape(this.H*x, [], 1) < 0);
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
        	y=[];
        	if isa(this.H,'LinOpIdentity')
				y = max(x,0.);
        	end
        	if isempty(y)
        		error('Prox not implemented');
        	end
        end
    end
end
