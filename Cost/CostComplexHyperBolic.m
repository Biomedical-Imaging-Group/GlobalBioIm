classdef CostComplexHyperBolic < Cost
    % CostComplexHyperBolic: Hyperbolic cost function
    % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{x}-y)_{k,l}^2 + \\varepsilon^2}$$
    %
    % :param index: dimensions along which the l2-norm will be applied (inner sum over l)
    % :param epsilon: \\(\\in \\mathbb{R}_+\\) smoothing parameter (default
    %                 \\(10^{-3}\\))
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % **Example** C=CostComplexHyperBolic(sz,epsilon,index,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`
    
    %%    Copyright (C) 2017 
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
    
    %% Properties
    properties
        epsilon;
        index;
        sumOp;
    end
    
    %% Constructor
    methods
        function this = CostComplexHyperBolic(sz,epsilon,index,y)
            if nargin<4, y=0; end
            this@Cost(sz,y);
            this.name='CostComplexHyperBolic';
            this.isConvex=true;
            this.isDifferentiable=true;
           
            if nargin<3|| isempty(index) %no sum
                index=0;
            end
            if nargin<2|| isempty(epsilon)
                epsilon=1e-3;
            end
            this.epsilon = epsilon;
            this.index = index;
            
            if index~=0
                this.sumOp = LinOpSum(sz,index);
            else
                this.sumOp = LinOpDiag(this.H.sizeout);
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.            
            if(isscalar(this.y)&&(this.y==0))
                u=abs(x).^2;
            else
                u=abs(x-this.y).^2;
            end
            R = this.sumOp.apply(u);
            
            F = sqrt(R + this.epsilon.^2);
            y = sum(F(:)) - numel(F).*this.epsilon;
        end    
        function g=applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.             
            if(isscalar(this.y)&&(this.y==0))
                u=abs(x).^2;
            else
                u=abs(x-this.y).^2;
            end
            R = this.sumOp.apply(u);
            
            F = sqrt(R + this.epsilon.^2);
            g = x.*this.sumOp.applyAdjoint(1./F);  
        end
    end
end