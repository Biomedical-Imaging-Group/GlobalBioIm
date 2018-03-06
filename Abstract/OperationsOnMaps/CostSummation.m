classdef CostSummation <  MapSummation & Cost
    % CostSummation : Sum of :class:`Cost`
    % $$C(\\mathrm{x}) = \\sum_i \\alpha_i C_i(\\mathrm{x}) $$
    %
    % :param costs:  cell of :class:`Cost`
    % :param alpha:  array of coefficients
    %
    % **Example** F = CostSummation(ACost,alpha)
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
            this.lip=costs{1}.lip;
            for n =2:length(costs)
                this.isConvex = this.isConvex & costs{n}.isConvex;
                if costs{n}.lip~=-1
                    this.lip = this.lip + costs{n}.lip;
                else
                    this.lip=-1;
                    break;
                end
            end
        end
        function M = makePartialSummation(this,Lsub)
            % Instanciation of :class:`CostPartialSummation`.
            % 
            % :param Lsub: number of :class:`Cost` used for computation

            M = CostPartialSummation(this.mapsCell,this.alpha,Lsub);
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function g=applyGrad_(this,x)
            % Reimplemented from :class:`Cost`
            g=this.alpha(1)*this.mapsCell{1}.applyGrad(x);
            for n=2:this.numMaps
                g=g+this.alpha(n)*this.mapsCell{n}.applyGrad(x);
            end
        end
        % the function reimplementations below is needed because of
        % the multiple inheritance to specifies which method to use from
        % parent classes
        function M = makeComposition_(this,G)
            M=makeComposition_@MapSummation(this,G);
        end
    end
end
