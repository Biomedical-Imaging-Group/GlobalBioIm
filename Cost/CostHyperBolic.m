classdef CostHyperBolic < Cost
    % CostHyperBolic: Hyperbolic cost function
    % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{x}-y)_{k,l}^2 + \\varepsilon^2}$$
    %
    % :param index: dimensions along which the l2-norm will be applied (inner sum over l)
    % :param epsilon: \\(\\in \\mathbb{R}_+\\) smoothing parameter (default
    %                 \\(10^{-3}\\))
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % **Example** C=CostHyperBolic(sz,epsilon,index,y)
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
    % - Public
    properties (SetObservable, AbortSet)
        epsilon;
        index;
    end
    % - Private
    properties (SetAccess = protected,GetAccess = protected)
        sumOp;
    end
    
    %% Constructor and handlePropEvents method
    methods
        function this = CostHyperBolic(sz,epsilon,index,y)
            % Default values
            if nargin<4, y=0; end
            if nargin<3|| isempty(index), index=0; end
            if nargin<2|| isempty(epsilon), epsilon=1e-3; end
            % Call superclass constructor
            this@Cost(sz,y);
            % Listeners to PostSet events
            addlistener(this,'index','PostSet',@this.handlePropEvents);
            addlistener(this,'epsilon','PostSet',@this.handlePropEvents);
            % Set properties
            this.name='CostHyperBolic';
            this.isConvex=true;
            this.isDifferentiable=true;
            this.epsilon = epsilon;
            this.index = index;          
            this.memoizeOpts.computeF=true;
            this.memoCache.computeF= struct('in',[],'out', []);
        end
        function handlePropEvents(this,src,~)
            % Reimplemented from superclass :class:`Cost`
            switch src.Name
                case 'index'
                    if this.index~=0
                        this.sumOp = LinOpSum(this.sizein,this.index);
                    else
                        this.sumOp = LinOpDiag(this.sizein);
                    end
            end
            % Call superclass method (important to ensure the right execution order)
            handlePropEvents@Cost(this,src);
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.            
            F = this.memoize('computeF',@this.computeF_,x);
            y = sum(F(:)) - numel(F).*this.epsilon;
        end    
        function g=applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.             
            F = this.memoize('computeF',@this.computeF_,x);
            g = x.*this.sumOp.applyAdjoint(1./F);  
        end
         function F = computeF_(this,x)
            % Reimplemented from parent class :class:`Cost`.             
            if(isscalar(this.y)&&(this.y==0))
                u=abs(x).^2;
            else
                u=abs(x-this.y).^2;
            end
            
            F = sqrt(this.sumOp*u + this.epsilon.^2);
        end
    end
end