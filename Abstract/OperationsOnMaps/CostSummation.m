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
            % Call superclass constructor
            this@MapSummation(costs,alpha);
            % Set properties 
            this.name='CostSummation';            
            % Initialize
            this.initialize('CostSummation');
        end
        function M = makePartialSummation(this,Lsub)
            % Instanciation of :class:`CostPartialSummation`.
            % 
            % :param Lsub: number of :class:`Cost` used for computation

            M = CostPartialSummation(this.mapsCell,this.alpha,Lsub);
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`MapSummation` and :class:`Cost`
            
            % Call superclass methods
            updateProp@MapSummation(this,prop);
            updateProp@Cost(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'mapsCell') ||  strcmp(prop,'all')
                allcosts = all( cellfun(@(x)(isa(x, 'Cost')), this.mapsCell) );
                assert(iscell(this.mapsCell) && allcosts, 'Property mapsCell should be a cell array Cost');
                this.isConvex=this.mapsCell{1}.isConvex;
                this.isSeparable=this.mapsCell{1}.isSeparable;
                for n =2:length(this.mapsCell)
                    this.isConvex = this.isConvex & this.mapsCell{n}.isConvex;
                    this.isSeparable = this.isSeparable & this.mapsCell{n}.isSeparable;
                end
                this.lip=this.mapsCell{1}.lip;
                if this.lip~=-1
                    for n =2:length(this.mapsCell)
                        if this.mapsCell{n}.lip~=-1
                            this.lip = this.lip + this.mapsCell{n}.lip;
                        else
                            this.lip=-1;
                            break;
                        end
                    end
                end
            end
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
        function x=applyProx_(this,z,alpha)
            % Reimplemented from :class:`Cost` in the case of the sum
            % between a :class:`CostRectangle` \\(i_C \\) and a
            % :class:`Cost` \\(f \\) which is separable [1]
            % $$ \\mathrm{prox}_{\\alpha(i_C +f)}(z) = \\mathrm{prox}_{i_c} \\circ \\mathrm{prox}_{\\alpha f}(z) $$
            %
            % **Reference**
            %
            % [1] "A Douglas?Rachford splitting approach to nonsmooth convex variational signal recovery"
            % P. L. Combettes, and J.C. Pesquet, Journal of Selected Topics in Signal Processing, 1(4), 564-574, 2007
            if this.numMaps==2
                if (isa(this.mapsCell{1},'CostRectangle') && this.mapsCell{2}.isSeparable)
                    x=this.mapsCell{1}.applyProx_(this.mapsCell{2}.applyProx_(z,alpha),1);
                elseif (isa(this.mapsCell{2},'CostRectangle') && this.mapsCell{1}.isSeparable)
                    x=this.mapsCell{2}.applyProx_(this.mapsCell{1}.applyProx_(z,alpha),1);
                else
                    x=applyProx_@Cost(this,z,alpha);
                end             
            else
                x=applyProx_@Cost(this,z,alpha);
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
