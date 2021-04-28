classdef CostPartialSummation <  CostSummation
    % CostPartialSummation : Sum of :class:`Cost` with apply, applyGrad,...
    % computed from a subset of Cost
    % $$C(\\mathrm{x}) = \\sum_i \\alpha_i C_i(\\mathrm{x}) $$
    %
    % :param costs:  cell of :class:`Cost`
    % :param alpha:  array of coefficients
    % :param Lsub:  number of :class:`Cost` used for computation
    % :param partialGrad: parameter for subset selection (0: no partial
    %  gradient; 1: stochastic gradient descent; 2: equally spaced indices)
    %
    % **Example** F = CostPartialSummation(ACost,alpha,Lsub)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`MapOpSummation`
    
    %%    Copyright (C) 2017
    %     T. Pham thanh-an.pham@epfl.ch
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
        counter=0; % counter for subset update
        subset; % current subset of angles used to compute F grad
        Lsub;   % Number of costs used in partial gradient
        fixedorder;
    end
    properties (SetAccess = public,GetAccess = public)
        partialGrad=1; % activate partial gradient option, 0 : no; 1 : stochastic; 2 : equally spaced;
    end
    
    %% Constructor
    methods
        function this = CostPartialSummation(costs,alpha,Lsub)
            this@CostSummation(costs,alpha);
            this.name='CostPartialSummation';
            this.setLsub(Lsub);
            
        end
        function setLsub(this,Lsub)
            % Set Lsub parameter
            this.Lsub = Lsub;
            this.counter = 0;
            if this.partialGrad==3
                this.fixedorder = randperm(this.numMaps);
            end
            this.updateSubset();
        end
        function setPartialGrad(this,newpartialGrad)
            this.partialGrad = newpartialGrad;
            if this.partialGrad == 3
                this.fixedorder = randperm(this.numMaps);
            end
        end
    end
    methods (Access = protected)
        function updateSubset(this)
            switch this.partialGrad
                case 0
                    this.subset = 1:this.numMaps;
                case 1
                    this.subset = sort(randi(this.numMaps,this.Lsub,1));
                case 2
                    this.subset = sort(1 + mod(round(this.counter...
                        + (1:this.numMaps/this.Lsub:this.numMaps)),this.numMaps));
                    this.counter = this.counter + 1;
                case 3
                    this.subset = this.fixedorder(1 + this.Lsub*this.counter:min(this.Lsub*(1+this.counter),end));
                    if this.Lsub*(1+this.counter) >=this.numMaps
                        this.counter = 0;
                    else
                        this.counter = this.counter + 1;
                    end
                    
                otherwise
                    error('Non existing partialGrad');
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from :class:`Cost`
            if this.partialGrad > 0
                y = 0;
                for kk = 1:length(this.subset)
                    ind = this.subset(kk);
                    y = y + this.alpha(ind)*this.mapsCell{ind}*x;
                end
                y = y/length(this.subset);%to decide whether it is kept.
            else
                y = apply_@CostSummation(this,x);
            end
        end
        function g=applyGrad_(this,x)
            % Reimplemented from :class:`Cost`
            if this.partialGrad > 0
                g = zeros_(size(x));
                for kk = 1:length(this.subset)
                    ind = this.subset(kk);
                    g = g + this.alpha(ind)*this.mapsCell{ind}.applyGrad(x);
                end
                g = g/length(this.subset);%to decide whether it is kept. real AD HOC
                this.updateSubset();%because of hierarchical ADMM, I need it after
            else
                g = applyGrad_@CostSummation(this,x);
            end
        end
    end
end
