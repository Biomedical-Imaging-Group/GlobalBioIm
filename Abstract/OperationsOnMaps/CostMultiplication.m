classdef CostMultiplication < Cost
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
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        cost1;
        cost2;
        isnum;
    end
    
    %% Constructor
    methods
        function this = CostMultiplication(C1,C2)
            this@Cost(C2.sizein);
            this.name='CostMultiplication';
            this.cost1 = C1;
            this.cost2 = C2;
            if isnumeric(C1)
                this.isnum =1;
                this.isConvex=this.cost2.isConvex;
                this.isSeparable=this.cost2.isSeparable;
                this.isDifferentiable=this.cost2.isDifferentiable;
                if this.cost2.lip~=-1
                    this.lip=this.cost2.lip*this.cost1;
                end
            else
                this.isConvex=0;  % It can be but we cannot conclude in a generic way ...
                this.isDifferentiable=this.cost1.isDifferentiable && this.cost2.isDifferentiable;
                % TODO: set the lip parameter properly ?
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
                    g=this.cost1*this.cost2.applyGrad(x);
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
                M=CostMultiplication(this.cost1,this.cost2*G);
            else
                M=makeComposition_@Cost(this,G);
            end
        end
     end
    
    methods (Access = protected)
         %% Copy
         function this = copyElement(obj)
             this = copyElement@Cost(obj);
            this.cost1 = copy(obj.cost1);
            this.cost2 =  copy(obj.cost1);
      end
    end
end