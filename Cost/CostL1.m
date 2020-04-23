classdef CostL1 < Cost
    % CostL1: L1 norm cost function
    % $$C(x) := \\|\\mathrm{x} - \\mathrm{y}\\|_1 $$
    %
    % If nonneg is set to true, it adds positivity constraint on x:
    % $$C(x) := \\|\\mathrm{x} - \\mathrm{y}\\|_1 + i(\\mathrm{x} - \\mathrm{y}) $$
    % where i() is the indicator function of $\\mathbb{R}^{+}$
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Example** C=CostL1(sz,y,nonneg)
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
    properties
        nonneg = false;
    end
    methods
        function this = CostL1(sz,y,nonneg)
            if nargin<2
                y=0;
                nonneg=false;
            else if nargin<3
                    nonneg=false;
                else
                    nonneg = nonneg;
                    warning('nonneg parameter will be removed in future releases. Use the sum of a CostL1 and a CostNonNeg instead.');
                end
            end
            this@Cost(sz,y);
            this.name='CostL1';
            this.isConvex=true;
            this.isSeparable=true;
            this.isDifferentiable=false;
            this.nonneg = nonneg;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{sign(x-y)} \\mathrm{max}(\\vert x-y \\vert- \\alpha,0)+  \\mathrm{y} $$
            if this.nonneg
                if(isscalar(this.y)&&(this.y==0))
                    y =    max( x - alpha,0);
                else
                    tmp = x-this.y ;
                    y =    max( tmp - alpha,0)+this.y;
                end
            else
                if(isscalar(this.y)&&(this.y==0))
                    y =  sign(x) .* max( abs( x) - alpha,0);
                else
                    tmp = x-this.y ;
                    y =  sign(tmp) .* max( abs( tmp) - alpha,0)+this.y;
                end
            end
        end
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.

            if(isscalar(this.y)&&(this.y==0))
                if this.nonneg && any(x(:)<0)
                    y = + inf;
                else
                    y=sum(abs(x(:)));
                end
            else
                tmp = x(:)-this.y(:);
                if this.nonneg && any(tmp(:)<0)
                    y = + inf;
                else
                    y=sum(abs(tmp));
                end
            end
        end
        function y = applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            % Subgradient of CostL1
            if ~this.nonneg
                y = sign(x-this.y);
            else
                this.applyGrad_@Cost(this,x);
            end
        end
    end
end
