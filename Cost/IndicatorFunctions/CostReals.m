classdef CostReals < CostRectangle
    % CostReals: Reals Indicator function
    % $$ C(x) = \\left\\lbrace \\begin{array}[l]
    % \\text{0~if } \\mathrm{xmin} \\leq \\mathrm{x-y} \\leq \\mathrm{xmax} \\newline
    % + \\infty \\text{ otherwise.} \\end{array} \\right. $$
    %
    % :param xmin: minimum value (default -inf)
    % :param xmax: maximum value (default +inf)
    %
    % All attributes of parent class :class:`CostRectangle` are inherited 
    %
    % **Example** C=CostReals(sz,xmin,xmax,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostIndicator`, :class:`CostRectangle`
    
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
        function this = CostReals(sz,xmin,xmax,y)
            if nargin<4, y=0; end
            assert(all(isreal(xmin))&&all(isreal(xmax)),'In CostReal bounds xmin and xmax must be real');
            this@CostRectangle(sz,xmin,xmax,y);
            this.name='CostReals';   
        end
    end
end
