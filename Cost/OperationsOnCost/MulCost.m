classdef MulCost < Cost
    %% MulCost : Multiplication of Costs
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % F = MulCost(Cost1,Cost2)
    % Multiplication of Costs
    % F = Cost1 * Cost2
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
        cost1;      
        cost2;
        isnum;
    end
    
     methods 
    	%% Constructor
        function this = MulCost(cost1,cost2)
            this.name='Multiply Costs';
			this.cost1 = cost1;
            this.cost2 = cost2;
            if isnumeric(cost1)
                this.isnum =1;
                this.isconvex=cost2.isconvex;
                if cost2.lip~=-1
                    this.lip=cost2.lip*cost1;
                end
            else               
                this.isconvex=0;  % It can be but we cannot conclude in a generic way ...
                this.lip=-1;     
            end		          
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
            if this.isnum
                y=this.cost1*this.cost2.eval(x);
            else
                y=this.cost1.eval(x)*this.cost2.eval(x);
            end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
            if this.isnum
                g=this.cost1*this.cost2.grad(x);
            else
                error('MulCost: Grad not implemented');
            end
        end
        %% Proximity operator of the Functional
        function y=prox(this,x,alpha)
        	assert(isscalar(alpha),'alpha must be a scalar');
            if this.isnum
                y = this.cost2.prox(x,this.cost1*alpha);
            else
                error('MulCost: Prox not implemented');
            end
        end
    end
end