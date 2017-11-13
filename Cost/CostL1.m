classdef CostL1 < Cost
    % CostL1: L1 norm cost function
    % $$C(x) := \\|\\mathrm{x} - \\mathrm{y}\\|_1 $$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % **Example** C=CostL1(sz,y)
    %
    % See also :class:`Map` :class:`Cost`, :class:`LinOp`

    %%    Copyright (C) 2015 
    %     F. Soulez ferreol.soulez@epfl.ch
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
        function this = CostL1(sz,y)
            if nargin<2, y=0; end
            this@Cost(sz,y);
            this.name='CostL1';
            this.isConvex=true;
            this.isDifferentiable=false;      
        end
    end

    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
	methods (Access = protected)
        function y=applyProx_(this,x,alpha) 
            % Reimplemented from parent class :class:`Cost`.
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{sign(x-y)} \\mathrm{max}(\\vert x-y \\vert- \\alpha,0)+  \\mathrm{y} $$
            if(isscalar(this.y)&&(this.y==0))
                y =  sign(x) .* max( abs( x) - alpha,0);
            else
                tmp = x-this.y ;
                y =  sign(tmp) .* max( abs( tmp) - alpha,0)+this.y;
            end
        end
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            if(isscalar(this.y)&&(this.y==0))
                y=sum(abs(x(:)));
            else
                y=sum(abs(x(:)-this.y(:)));
            end
        end
    end
end
