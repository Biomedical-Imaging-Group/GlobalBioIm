classdef CostL32 < Cost
	% CostL32: L3/2 norm cost function
    % $$ C(\\mathrm{x}) := \\|\\mathrm{x}\\|_\\frac32^\\frac32 $$
	%
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **References**
    % 
    % [1] C. Chaux, P. L. Combettes, J.-C. Pesquet, and V. R. Wajs. A variational formulation for 
    % frame-based inverse problems. Inverse Problems, 23(4):1495-1518, June 2007.
    %
    % **Example** C=CostL32(sz)
    %
    % See also :class:`Map` :class:`Cost`
    
    %%    Copyright (C) 2025
    %     A. Floquet arthur.floquet@irit.fr
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
        function this = CostL32(sz,y)
            if nargin<2, y=0; end 
            this@Cost(sz,y);
            this.name='CostL32';
            this.isConvex=true;
            this.isDifferentiable=false;
            this.isSeparable=true;      
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            y = sum(abs(x).^(3/2), "all");
        end

        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{x} + \\frac98 \\alpha^2 \\mathrm{sign(x)} \\left( 1 - \\sqrt(1 + \\frac{16|x|}{9\\alpha^2} \\right) $$
            y = x + (9/8).*(alpha^2).*sign(x).*(1 - sqrt(1 + (16*abs(x)./(9*alpha^2))));
        end    
      

    end
end
