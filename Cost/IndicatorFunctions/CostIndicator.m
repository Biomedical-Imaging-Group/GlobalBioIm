classdef  CostIndicator < Cost
    % CostIndicator: Indicator cost function
    % $$ C(\\mathrm{x}) = \\left\\lbrace \\begin{array}[l] \\text{0~if } \\mathrm{x -y} \\in \\mathrm{C},  
    % \\newline + \\infty \\text{ otherwise.} \\end{array} \\right. $$
    % where \\(\\mathrm{C} \\subset \\mathrm{X} \\) is a constraint set.
    %
    % All attributes of parent class :class:`Cost` are inherited 
    %
    %
    % **Note** :class:`CostIndicator` is an generic class for all
    % indicator cost functions
    %
    % See also :class:`Map`, :class:`Cost`

    %%    Copyright (C) 2017 
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
        function this = CostIndicator(sz,y)   
            if nargin<2, y=0; end
            this@Cost(sz,y);
            this.isDifferentiable=false;
        end
    end
            
    %% Core Methods containing implementations (Protected)
    % - applyGrad_(this,x)
    methods (Access = protected)
        function g = applyGrad_(this,x)
            error('Cannot evaluate gradient of an indicator function');
        end
    end
end
