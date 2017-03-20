classdef CostIndicator < Cost
    %% CostIndicator : Indicator function  
    %  Matlab Inverse Problems Library
    % 
    % Obj = CostIndicator(xmin, xmax,H,y):
    % Implement the indicator function of the retangle set defined by xmin and xmax for the 
    % $$ \phi(x) = 0 \textrm{ if } xmin <= (x-y) <= xmax textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    %
    %% Properties
    % * |xmin|         - minimum value (default 0)
    % * |xmax|         - minimum value (default +inf)
    %
    %%
    
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
        xmin=0
        xmax=+inf;
    end
    
    methods
        function this = CostIndicator(xmin, xmax,H,y)
          
            this.name='Cost Indicator';
            
             % -- Set entries
            if nargin<4
                y=[];
            end
            if nargin<3
                H=[];
            end
            set_H(this,H,y);
            
            if nargin<2
                xmax=[];
            end 
            if nargin==0
                xmin =0;
            end
            
            if isempty(xmax)
                xmax = +inf;
            end
            this.xmin = xmin;
            this.xmax = xmax;
            
        end
        function z = prox(this,x,~) % Apply the operator
            
        	if this.H.isinvertible
				z = this.H.Inverse(min(max(this.H.Apply(x)-this.y, this.xmin),this.xmax)+this.y);
            else
        		error('Prox not implemented');
        	end
            
          
                
        end
    end
end
