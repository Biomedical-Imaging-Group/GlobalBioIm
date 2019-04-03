classdef CostComplexCircle < CostComplexRing
    % CostComplexCircle: Implements the indicator on the complex circle of radius c
    % $$ C(x) = \\left\\lbrace \\begin{array}[l]
    % \\text{0~if } \\vert \\mathrm{x-y}\\vert =c  \\newline
    % + \\infty \\text{ otherwise.} \\end{array} \\right. $$
    %
    % All attributes of parent class :class:`CostComplexRing` are inherited 
    %
    % 
    % :param radius: radius of the circle (default 1)
    % 
    % **Example** C=CostComplexCircle(radius,sz,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostIndicator`, :class:`CostComplexRing`
    
    %%    Copyright (C) 2015 
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
        function this = CostComplexCircle(sz,radius,y)
            if nargin<3, y=0; end
            if nargin >0
                assert(isnumeric(radius),'The radius argument should be numeric');
            else
                radius = 1;
            end
            this@CostComplexRing(sz,radius,radius,y);
            this.name='CostComplexCircle';
			this.isConvex= false;    
        end
    end
end
