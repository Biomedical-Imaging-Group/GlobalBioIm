classdef CostL43 < Cost
	% CostL43: L4/3 norm cost function
    % $$ C(\\mathrm{x}) := \\|\\mathrm{x}\\|_\\frac43^\\frac43 $$
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
        function this = CostL43(sz,y)
            if nargin<2, y=0; end 
            this@Cost(sz,y);
            this.name='CostL43';
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
            y = sum(abs(x).^(4/3), "all");
        end

        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{x} + \\frac{4\\alpha}{3\\sqrt[3]{2}} \\left( \\sqrt[3]{\\xi-\\mathrm{x}} - \\sqrt[3]{\\xi+\\mathrm{x}} \\right)$$
            % where $$\\xi = \\sqrt{\\mathrm{x}^2 + 256 \\alpha^3 / 729}$$.
            xi = (x.^2 + 256*alpha^3/729).^(1/2);
            y  = x + (4*alpha)/(3*2^(1/3)) * ((xi-x).^(1/3) - (xi+x).^(1/3));
        end    
      

    end
end
