classdef LinOpSummation < MapSummation &  LinOp
    % LinOpSummation: Sum of linear operators
    % $$ \\mathrm{H}(\\mathrm{x}) = \\sum_i \\alpha_i \\mathrm{H}_i(\\mathrm{x}) $$
    %
    % :param LinOps:  cell of :class:`LinOp`
    % :param alpha:  array of coefficients
    %
    % **Example** L=LinOpSummation(LinOps,alpha)
    %
    % See also :class:`Map`, :class:`LinOp`, :class:`MapOpSummation`
    
    %%    Copyright (C) 2017
    %     F. Soulez ferreol.soulez@epfl.ch
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
        function this = LinOpSummation(LinOps,alpha)
            this@MapSummation(LinOps,alpha); 
            this.name ='LinOpSummation';			
			allLinOps = all( cellfun(@(x)(isa(x, 'LinOp')), LinOps) );
			assert(iscell(LinOps) && allLinOps, 'First input should be a cell array LinOp'); 
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x) 
    % - applyHHt_(this,y) 
    % - makeAdjoint_(this)
    methods (Access = protected)      		
        function x = applyAdjoint_(this,y) 
            % Reimplemented from :class:`LinOp` 
            x =  zeros(this.sizein);
            for n = 1:this.numMaps
                x = x + this.alpha(n) .* this.mapsCell{n}(1).applyAdjoint(y);
            end
        end	
		function y = applyHtH_(this,x) 
            % Reimplemented from :class:`LinOp` 
			y =  zeros(this.sizein);
			for n = 1:this.numMaps
				y = y + this.alpha(n) .* this.mapsCell{n}.applyHtH(x);
			end
        end		
		function x = applyHHt_(this,y) 
            % Reimplemented from :class:`LinOp` 
            x =  zeros(this.sizeout);
            for n = 1:this.numMaps
                x = x + this.alpha(n) .* this.mapsCell{n}.applyHHt(y);
            end
		end
        function M = makeAdjoint_(this)
            % Reimplemented from :class:`LinOp`
            adjointCell = cellfun(@(x)x',this.mapsCell,'UniformOutput',false);
            M=LinOpSummation(adjointCell,this.alpha);
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

