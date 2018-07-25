classdef CostMixNorm21NonNeg < Cost
    % CostMixNorm21: Mixed norm 2-1 with non negativity constraints cost function
    % $$C(\\mathrm{x}) :=  \\left\\lbrace \\begin{array}{ll} \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{x}-y)_{k,l}^2} & \\text{ if } \\mathrm{x-y} \\geq 0 \\newline
    % + \\infty & \\text{ otherwise.} \\end{array} \\right. $$
    %
    % :param index: dimensions along which the l2-norm will be applied (inner sum over l)
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Example** C=CostMixNorm21(sz,index,y)
    %
    % See also :class:`Map` :class:`Cost`, :class:`LinOp`
    
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
        index;    % dimensions along which the l2-norm will be applied
        kerdims
        imdims
    end
    
    %% Constructor
    methods
        function this = CostMixNorm21NonNeg(sz,index,y)
            if nargin<3, y=0; end
            this@Cost(sz,y);
            this.name='CostMixNorm21NonNeg';
            assert(isnumeric(index)&&isvector(index),'The index should be a vector of integers');
            this.index=index;
            this.isConvex=true;
            this.isDifferentiable=false;
            
            
            ndms = length(this.sizein);
            T = true(ndms,1);
            T(this.index)=false;
            this.kerdims = this.sizein; this.kerdims(T)=1;
            this.imdims = this.sizein; this.imdims(~T)=1;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            if(isscalar(this.y)&&(this.y==0))
                if any(x(:)<0)
                    y = + inf;
                    return
                else
                    u=abs(x).^2;
                end
            else
                tmp = x(:)-this.y(:);
                if   any(tmp(:)<0)
                    y = + inf;
                    return
                else
                    u=abs(tmp).^2;
                end
            end
            % Computes the l2-norm along the dimensions given by index
            for n=1:length(this.index)
                u = sum(u,this.index(n));
            end
            u = sqrt(u);
            y=sum(u(:));
        end
        function z=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}_{k\\cdot} ) = \\left\\lbrace
            % \\begin{array}{ll}
            % \\max(\\mathrm{x}_{k\\cdot} - y_{k\\cdot} ,0)
            % \\left(1-\\frac{\\alpha}{\\Vert(\\mathrm{x}-y)_{k\\cdot}\\Vert_2}
            % \\right) +y_{k\\cdot} cd .. & \\; \\mathrm{if } \\;
            % \\Vert (\\mathrm{x-y})_{k\\cdot}\\Vert_2 > \\alpha,
            % \\newline
            % 0 & \\; \\mathrm{otherwise},
            % \\end{array}\\right. \\; \\forall \\, k $$
            % where the division is component-wise.
            
            % Computes the l2-norm along the dimensions given by index
            if(isscalar(this.y)&&(this.y==0))
                sx = max(x,0.).^2;
            else
                sx = max(x-this.y,0).^2;
            end
            for n=1:length(this.index)
                sx = sum(sx,this.index(n));
            end
            sx = sqrt(sx);
            
            % Computes the prox
            t = sx > alpha;
            b = zeros_(size(sx));
            
            b(t) = 1-alpha./sx(t);
            if(isscalar(this.y)&&(this.y==0))
                z = reshape(repmat(reshape(b ,this.imdims),this.kerdims),this.sizein).*x;
            else
                z = reshape(repmat(reshape(b ,this.imdims),this.kerdims),this.sizein).*x+this.y;
            end
            % result:
            % x(||x|| <= alpha) = 0
            % x(||x|| > alpha) = x(||x|| > alpha) - ...
            %      x(||x|| > alpha) / ||x|| * alpha
        end
        
    end
end
