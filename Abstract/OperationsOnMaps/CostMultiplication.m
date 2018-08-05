classdef CostMultiplication <  Cost
    % CostMultiplication: Multiplication of Costs
    % $$C(\\mathrm{x}) = C_1(\\mathrm{x}) \\times C_1(\\mathrm{x}) $$
    %
    % :param C1: a :class:`Cost` object or a scalar
    % :param C2: a :class:`Cost` object
    %
    % **Example** F = MulCost(Cost1,Cost2)
    %
    % See also :class:`Map`, :class:`Cost`
    
    % TODO: Write a MapMultiplication (pointwise) from which CostMultiplication will
    % derive. Also overload the .* operation in Map to perform such a
    % MapMultiplication.
    %
    
	
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
    
    %% Properties
    % - Public 
    properties (SetObservable, AbortSet)
        cost1;
        cost2;
    end
    % - Private
    properties (SetAccess = protected,GetAccess = public)
        isnum;
    end
    
    %% Constructor 
    methods
        function this = CostMultiplication(C1,C2)
            % Call superclass constructor
            this@Cost(C2.sizein);
            % Basic properties
            this.name='CostMultiplications';
            this.cost1 = C1;
            this.cost2 = C2;
            % Initialize
            this.initialize('CostMultiplication');
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`Cost`
            
            % Call superclass method
            updateProp@Cost(this,prop);
            % Update current-class specific properties
            if (strcmp(prop,'cost1') ||  strcmp(prop,'all')) && isnumeric(this.cost1)  && isscalar(this.cost1)
                this.cost1=LinOpDiag(this.sizeout,this.cost1);
                this.isnum =1;
                if this.cost2.lip~=-1, this.lip=this.cost1*this.cost2.lip; end
            end
            if strcmp(prop,'cost2') ||  strcmp(prop,'all')
                if this.isnum
                    this.isConvex=this.cost2.isConvex;
                    this.isSeparable=this.cost2.isSeparable;
                    this.isDifferentiable=this.cost2.isDifferentiable;
                    if this.cost2.lip~=-1, this.lip=this.cost1*this.cost2.lip; end
                else
                    this.isConvex=0;  % It can be but we cannot conclude in a generic way ...
                    this.isDifferentiable=this.cost1.isDifferentiable && this.cost2.isDifferentiable;
                    % TODO: set the lip parameter properly ?
                end
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - makeComposition_(this,G)
     methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from :class:`Cost`
            if this.isnum
                y=this.cost1*this.cost2.apply(x);
            else
                y=this.cost1.apply(x)*this.cost2.apply(x);
            end
        end
        function g=applyGrad_(this,x)
            % Reimplemented from :class:`Cost`
            if this.isDifferentiable
                if this.isnum
                    g=this.cost1.diag*this.cost2.applyGrad(x);
                else
                    g=this.cost1.apply(x)*this.cost2.applyGrad(x) + this.cost1.applyGrad(x)*this.cost2.apply(x);
                end
            else
                g=applyGrad_@Cost(this,x);
            end
        end
        function y=applyProx_(this,x,alpha)
            % Reimplemented from :class:`Cost`
            if this.isnum
                y = this.cost2.applyProx(x,this.cost1*alpha);
            else
                y=applyProx_@Cost(this,x,alpha);
            end
        end
        function M=makeComposition_(this,G)
            % Reimplemented from :class:`Cost`
            if this.isnum
                M=CostMultiplication(this.cost1.diag,this.cost2*G);
            else
                M=makeComposition_@Cost(this,G);
            end
        end
    end
end