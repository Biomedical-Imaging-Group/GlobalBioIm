classdef CostSummation <  MapSummation & Cost
    % CostSummation : Sum of :class:`Costs`
    % $$C(\\mathrm{x}) = \\sum_i \\alpha_i C_i(\\mathrm{x}) $$
    %
    % :param costs:  cell of :class:`Costs`
    % :param alpha:  array of coefficients
    %
    % **Example** F = SumCost(ACost,alpha)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`MapOpSummation`
	
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
    
    %% Constructor
    methods 
        function this = CostSummation(costs,alpha)
            this@MapSummation(costs,alpha); 
            this.name='CostSummation';
            allcosts = all( cellfun(@(x)(isa(x, 'Cost')), costs) );
			assert(iscell(costs) && allcosts, 'First input should be a cell array Cost');
			this.isConvex=costs{1}.isConvex; 
            for n =2:this.numMaps
                this.isConvex = this.isConvex & costs{1}.isConvex; 
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x) 
    % - applyHHt_(this,y) 
    % - makeAdjoint_(this)
    methods (Access = protected)    
        function g=applyGrad_(this,x)
            % Reimplemented from :class:`Cost`
			g=this.alpha(1)*this.costs{1}.applyGrad(x);
			for n=2:this.numMaps
				g=g+this.alpha(n)*this.costs{n}.applyGrad(x);
			end
        end
        % the function reimplementations below is needed because of
        % the multiple inheritance to specifies which method to use from
        % parent classes
        function M = makeComposition_(this,G)
            % Reimplemented from class :class:`MapSummation`
            M=makeComposition_@MapSummation(this,G);
        end
    end
end
