classdef ComposeLinOpCost < Cost
    %% ComposeLinOpCost : Compose a Functional with a linear operator
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % G =  ComposeLinOpCost(cost,Hcomp)
    % where cost is a COST object and Hcomp a LINOP one
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
		cost;      % Functional
    end
    
    methods 
    	%% Constructor
        function this = ComposeLinOpCost(cost,Hcomp)
            this.name='ComposeLinOpCost';
            this.cost=cost;
            this.set_H(Hcomp);
			this.isconvex= cost.isconvex; 			
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.cost.eval(this.H.apply(x));
        end
    end
end
