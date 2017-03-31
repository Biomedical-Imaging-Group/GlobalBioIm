classdef SumCost < Cost
    %% SumCost : Sum of Costtionnals
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % F = SumCost(ACost,alpha)
    % Sum the all Cost contained in vector AFUNC weighted by ALPHA (default 1)
    % F  sum_n alpha(n) * ACost(n)
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
        costs      % cell of Cost
        numcosts   % number of Cost
        alpha      % scalar factor (array)
    end
    
    methods 
    	%% Constructor
        function this = SumCost(costs,alpha)
            if nargin == 1, alpha = ones(size(costs));end
            this.name='Sum of Costs';
            this.numcosts = numel(costs);
            assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numcosts))),'Second input should be a scalar or an array of scalar of the same size as the first input');
            allcosts = all( cellfun(@(x)(isa(x, 'Cost')), costs) );
			assert(iscell(costs) && allcosts, 'First input should be a cell array Cost');
			if  isscalar(alpha)
				this.alpha = repmat(alpha, 1, this.numcosts) ;
			else
				this.alpha = alpha;
			end
			this.costs = costs;
			this.isconvex=costs{1}.isconvex; 
            for n =2:this.numcosts
                this.isconvex = this.isconvex & costs{1}.isconvex; 
            end
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.alpha(1)*this.costs{1}.eval(x);
			for n=2:this.numcosts
				y=y+this.alpha(n)*this.costs{n}.eval(x);
			end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			g=this.alpha(1)*this.costs{1}.grad(x);
			for n=2:this.numcosts
				g=g+this.alpha(n)*this.costs{n}.grad(x);
			end
        end
    end
end
