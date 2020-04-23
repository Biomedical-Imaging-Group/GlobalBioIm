classdef CostLInf < Cost
    % CostLInf: $L_\infty$ norm cost function
    % $$C(x) := \\|\\mathrm{x} - \\mathrm{y}\\|_\infty $$
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Example** C=CostLInf(sz,y)
    %
    % See also :class:`Map` :class:`Cost`, :class:`LinOp`
    
    %%    Copyright (C) 2020
    %     T. Debarre thomas.debarre@epfl.ch
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
    end
    methods
        function this = CostLInf(sz,y)
            if nargin<2
                y=0;
            end
            this@Cost(sz,y);
            this.name='CostLInf';
            this.isConvex=true;
            this.isSeparable=false;
            this.isDifferentiable=false;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y = applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{x} - \\alpha \\mathrm{proj_{\\Vert \\cdot \\Vert_1 \\leg 1}} ((\\mathrm{x} - \\mathrm{y})/\\alpha) $$
            % This code is adapted from DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher
            % https://www.epfl.ch/labs/lions/technology/decopt/
            % Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
            % EPFL Lausanne, 1015-Lausanne, Switzerland.
            tmp = (x(:) - this.y(:)) / alpha;
            sorted   = sort(abs(nonzeros(tmp)), 'descend');
            cum_sorted  = cumsum(sorted);
            nidx = find( cum_sorted - (1:numel(sorted))'.*[sorted(2:end); 0] >= 1 ...
                + 2*eps, 1);
            if ~isempty(nidx)
                dx = ( cum_sorted(nidx) - 1 ) /nidx;
                l1_proj = tmp.*( 1 - dx./ max(abs(tmp), dx) );
            else
                l1_proj = tmp;
            end
            y = reshape(x(:) - alpha*l1_proj, this.sizein);
        end
        
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            y = max(abs(x(:) - this.y(:)));
        end
        function y = applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            % Subgradient of CostLInf
            y = zeros(size(x));
            [~, idx] = max(abs(x(:) - this.y(:)));
            if isscalar(this.y)
                y(idx) = sign(x(idx) - this.y);
            else
                y(idx) = sign(x(idx) - this.y(idx));
            end
        end
    end
end
