classdef CostTV < CostComposition
    % CostTV : Composition of a :class:`CostTV` with a
    % :class:`Map`
    % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{H_2 x}-y)_{k,l}^2}= \\sum_{k=1}^K \\Vert (\\mathrm{H_2 x-y})_{k\\cdot} \\Vert_2$$
    %
    % :param H1: :class:`CostMixNorm21` object
    % :param H2: :class:`LinOpGrad` object
    %
    % All attributes of parent :class:`CostComposition` are inherited.
    %
    % **Example** C=CostTV(sz)
    %
    % **Example** C=CostTV(H1,H2); with H1=CostMixNorm21(sz,index,y); and H2=LinOpGrad(sz);
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostL2`, :class:`CostComposition`, :class:`LinOp`
    
    %%    Copyright (C) 2017
    %     T-a Pham  thanh-an.pham@epfl.ch &
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
        bounds   % Bounds for set constraint
        maxiter
        xtol
    end
    % - Protected
    properties (SetAccess = protected,GetAccess = protected)
        optim;
    end
    % Emmanuel : We should put optim public and remove the other properties
    
    %% Constructor 
    methods
        function this = CostTV(varargin)
            % Default values
            if nargin==1
                assert(issize(varargin{1}),'argument must be conformable to a size');
                H2=LinOpGrad(varargin{1});
                H1=CostMixNorm21(H2.sizeout,numel(H2.sizeout));
            elseif nargin==2
                assert(isa(varargin{1},'CostMixNorm21'),'First argument must be a CostMixNorm21 object');
                assert(isa(varargin{2},'LinOpGrad'),'Second argument must be a LinOpGrad object');
                H1=varargin{1};H2=varargin{2};
            end    
            % Call superclass constructor
            this@CostComposition(H1,H2);
            % Set properties
            this.name='CostTV';
            this.bounds = [-inf,inf];% Bounds for set constraint
            this.maxiter = 20;
            this.xtol=1e-5;
            % Initialize
            this.initialize('CostTV');
        end
        function setProxAlgo(this,bounds,maxiter,xtol,Outop)
            % Set the parameters of :class:`OptiFGP` used to compute the proximity
            % operator. 
            
            this.bounds = bounds;
            this.maxiter = maxiter;
            this.xtol=xtol;
            this.optim.OutOp=Outop;
            warning('setProxAlgo will be removed in future releases. Instead set directly this.myParam=new_value');
        end       
    end
        %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`CostComposition`
            
            % Call superclass method
            updateProp@CostComposition(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'bounds') ||  strcmp(prop,'all')
                LS = CostL2(this.sizein);
                this.optim = OptiFGP(LS,this,this.bounds);
                this.optim.ItUpOut = 0;
                this.optim.verbose = 0;
                this.optim.OutOp=OutputOpti(0);
            end
            if strcmp(prop,'xtol') ||  strcmp(prop,'all')
                this.optim.CvOp = TestCvgStepRelative(this.xtol);
            end
            if strcmp(prop,'maxiter') ||  strcmp(prop,'all')
                this.optim.maxiter = this.maxiter;
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyProx_(this,z,alpha)
    methods (Access = protected)
        function y = applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`CostComposition`.
            % Computed using the iterative :class:`OptiFGP` 
            
            % If y==0
            if this.y==0                
                this.optim.F0.y=x;
                if isa(this.optim.TV,'CostTV')
                    this.optim.TV = alpha*this.optim.TV;
                else
                    this.optim.TV = alpha*this.optim.TV.cost2;
                end
                this.optim.run(x);
                y = this.optim.xopt;
            else
                error('Proximity operator not implemented for the composition CostMixNorm21 and %s with non-null this.y',this.H2.name);
            end
        end
        
        
        function M = sum_(this,H)
            % Reimplemented from parent class :class:`CostComposition`.
            
            if isa(H,'CostRectangle')
                this.bounds = [max(H.xmin,this.bounds(1)), min(H.xmax,this.bounds(end))]; %Combine existing bounds with new one
                if diff(this.bounds) < 0
                    error('There is no intersection between the combined set constraints. If you want to replace the already existing bounds, change this.bounds property');
                end
                M = this;
            else
                M=sum_@CostComposition(H);
            end
        end
    end
end
