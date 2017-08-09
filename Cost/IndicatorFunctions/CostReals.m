classdef CostReals < CostRectangle
    %% CostReals : Reals Indicator function
    %  Matlab Inverse Problems Library
    %
    % Obj = CostReals(xmin, xmax,H,y):
    % Implement the indicator function of the Real set defined by xmin and xmax for the
    % $$ \phi(x) = 0 \textrm{ if } xmin <= (x-y) <= xmax textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    %
    %% Properties
    % * |xmin|         - minimum value (default -inf)
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
    end
    
    methods
        function this = CostReals(xmin, xmax,H,y)
            
            this.name='Cost Reals';
            
            % -- Set entries
            if nargin<4
                y=0;
            end
            if nargin<3
                H=[];
            end
            set_y(this,y);
            set_H(this,H);
            
            if nargin<2
                xmax=[];
            end
            if nargin==0
                xmin =-inf;
            end
            
            if isempty(xmax)
                xmax = +inf;
            end
            this.xmin = xmin;
            this.xmax = xmax;
            
        this.isComplex = false;
            
        end
    end
end
