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
    
    properties  (SetAccess = protected,GetAccess = protected)
        warn=0;
    end
    
    properties
        optim;
        bounds = [-inf,inf];% Bounds for set constraint
        maxiter = 20;
    end
    
    %% Constructor
    methods
        function this = CostTV(varargin)
            if nargin==1
                assert(issize(varargin{1}),'argument must be conformable to a size');
                H2=LinOpGrad(varargin{1});
                H1=CostMixNorm21(H2.sizeout,numel(H2.sizeout));
            elseif nargin==2
                assert(isa(varargin{1},'CostMixNorm21'),'First argument must be a CostMixNorm21 object');
                assert(isa(varargin{2},'LinOpGrad'),'Second argument must be a LinOpGrad object');
                H1=varargin{1};H2=varargin{2};
            end    
            this@CostComposition(H1,H2);
            this.name='CostTV';
            LS = CostL2(this.sizein);
            %lambda is always 1 here. It can be set differently by multiplying this cost with a scalar
            this.optim = OptiFGP(LS,this,this.bounds);
            this.optim.maxiter = this.maxiter;
            this.optim.ItUpOut = 0;
            this.optim.verbose = 0;
        end
        
        function setProxAlgo(this,bounds,maxiter,xtol,Outop)
            % Set the parameters of :class:`OptiFGP` used to compute the proximity
            % operator. 
            
            if nargin<=1 || isempty(bounds),bounds = this.bounds;end
            if nargin<=2 || isempty(maxiter),maxiter = this.maxiter;end
            if nargin<=3 || isempty(xtol),xtol = [];end
            if nargin<=4 || isempty(Outop),Outop = OutputOpti(0);end
            
            this.bounds = bounds;
            this.maxiter = maxiter;
      
            LS = CostL2(this.sizein);
            %lambda is always 1 here. It can be set differently by multiplying this cost with a scalar
            this.optim = OptiFGP(LS,this,this.bounds);
            this.optim.OutOp=Outop;
            this.optim.maxiter = this.maxiter;
            this.optim.ItUpOut = 0;
            this.optim.CvOp = TestCvgStepRelative(xtol);
            this.optim.verbose = 0;
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
                if ~this.warn % To raise the warning only once
                   warnStruct = warning('off','backtrace'); warning('The prox in CostTV is iterative (OptiFGP): This may lead to slow computations.');this.warn=1;warning(warnStruct);
                end
                %for k = 1:length(this.optim.Fn)
                %    this.optim.Fn{k} = alpha*this.optim.Fn{k};%what is that
                %end
                this.optim.F0.set_y(x);
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
