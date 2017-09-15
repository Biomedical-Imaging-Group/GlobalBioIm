classdef CostNonNeg < CostReals
    % CostNonNeg : Non negativity indicator 
    % $$ C(x) = \\left\\lbrace \\begin{array}[l]
    % \\text{0~if } \\mathrm{x-y} \\geq 0 \\newline
    % + \\infty \\text{ otherwise.} \\end{array} \\right. $$
    %
    % All attributes of parent class :class:`CostRectangle` are inherited 
    %
    % **Example** C=CostNonNeg(sz,y) 
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostIndicator`, :class:`CostRectangle`
    % :class:`CostReals` 
	
    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch
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
        function this = CostNonNeg(sz,y) 
            if nargin<2, y=0; end
            this@CostReals(sz,0,[],y);
            this.name='CostNonNeg';
        end
    end
end
