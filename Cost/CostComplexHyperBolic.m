classdef CostComplexHyperBolic < Cost
    % Hyperbolic cost function
    % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{Hx}-y)_{k,l}^2 + \\varepsilon^2}$$
    %
    % :param index: dimensions along which the l2-norm will be applied (inner sum over l)
    % :param epsilon: \\(\\in \\mathbb{R}_+\\) smoothing parameter (default
    % \\(10^{-3}\\))
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % See also :class:`Cost` :class:`LinOp`
    
    %     Copyright (C) 2017 F. Soulez ferreol.soulez@epfl.ch
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
    
    
    properties
        epsilon;
        index;
        sumOp;
    end
    
    methods
        function this = CostComplexHyperBolic(epsilon,index,H,y)
            this.name='CostComplexHyperBolic';
            
            this.isconvex=true;
            % -- Set entries
            if nargin<4
                y=0;
            end
            if nargin<3
                H=[];
            end           
            set_y(this,y);
            set_H(this,H);
            
            if nargin<2|| isempty(epsilon) %no sum
                index=0;
            end
            if nargin<1|| isempty(epsilon)
                epsilon=1e-3;
            end
            this.epsilon = epsilon;
            this.index = index;
            
            if index~=0
                this.sumOp = LinOpSum(this.H.sizeout,index);
            else
                this.sumOp = LinOpIdentity(this.H.sizeout);
            end
        end
        
        function cost  = eval(this,x) 
            % Reimplemented from parent class :class:`Cost`.
            
            if(isscalar(this.y)&&(this.y==0))
                z = this.H.apply(x);
            else
                z = this.H.apply(x) -this.y;
            end
            R = this.sumOp.apply(abs(z).^2);
            
            F = sqrt(R + this.epsilon.^2);
            cost = sum(F(:)) - numel(F).*this.epsilon;
        end
        
        function gradient = grad(this,x)
            % Reimplemented from parent class :class:`Cost`. 
            
            if(isscalar(this.y)&&(this.y==0))
                z = this.H.apply(x);
            else
                z = this.H.apply(x) -this.y;
            end
            R = this.sumOp.apply(abs(z).^2);
            
            F = sqrt(R + this.epsilon.^2);
            gradient = this.H.adjoint(z.*this.sumOp.adjoint(1./F)+ this.y);           
        end
        
        function [cost , gradient] = eval_grad(this,x) 
            % Reimplemented from parent class :class:`Cost`.
            
            if(isscalar(this.y)&&(this.y==0))
                z = this.H.apply(x);
            else
                z = this.H.apply(x) -this.y;
            end
            R = this.sumOp.apply(abs(z).^2);
            
            F = sqrt(R + this.epsilon.^2);
            cost = sum(F(:)) - numel(F).*this.epsilon;
            
            gradient = this.H.adjoint(z.*this.sumOp.adjoint(1./F)+ this.y);          
        end
    end
end