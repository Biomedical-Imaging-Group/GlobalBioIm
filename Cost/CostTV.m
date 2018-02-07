classdef CostTV < CostComposition
    % CostTV : Composition of a :class:`CostTV` with a
    % :class:`Map`
    % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{Hx}-y)_{k,l}^2}= \\sum_{k=1}^K \\Vert (\\mathrm{Hx-y})_{k\\cdot} \\Vert_2$$
    %
    % :param H1: :class:`CostMixNorm21` object
    % :param H2:  :class:`Map` object
    %
    % All attributes of parent :class:`CostComposition` are inherited.
    %
    % **Example** C=CostMixNorm21(sz,index,y)*LinOpGrad(sz);
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostL2`, :class:`CostComposition`, :class:`LinOp`
    
    %%    Copyright (C) 2017
    %     T-A. Pham thanh-an.pham@epfl.ch
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
    
    properties (Access=private)
        OpSumP;  % Operator used to compute the prox when H2 is the combination of a LinOpDownsample with a LinOpConv
        LLt;     % averaged convolution kernel used to compute the prox when H2 is the combination of a LinOpDownsample with a LinOpConv
    end
    properties
        optimizer = 'FGP';  % flag used to choose the computation of the proximity operator (e.g. 'FGP' for total variation if H2 is a LinOpGrad). Note a set constraint (e.g. nonnegativity) can be imposed.
        optim;
        bounds = [-inf,inf];% Bounds for set constraint
        maxiter = 20;
        %Parameters related to ADMM
    end
    
    %% Constructor
    methods
        function this = CostTV(H1,H2)
            assert(isa(H1,'CostMixNorm21'),'First argument must be a CostMixNorm21 object');
            this@CostComposition(H1,H2);
            %if this.H2.norm>=0
            %    this.lip=this.H1.lip*this.H2.norm^2;
            %end
            this.name=sprintf('CostMixNorm21Composition ( %s )', H2.name);
            LS = CostL2(this.sizein);
            switch this.optimizer
                case 'FGP'
                    this.optim = OptiFGP(LS,this,this.bounds,OutputOpti(0));
                    this.optim.Nesterov = true;
                    this.optim.maxiter = this.maxiter;
                    this.optim.ItUpOut = 0;
                    
                case 'ADMM'
                    NN = this.bounds;
                    this.optim = OptiADMM(LS,{NN, this.L12},...
                        this.Hn,this.rhon, [], []);
                    this.optim.maxiter = this.maxiter;
                    this.optim.xtol = this.xtol;
            end
        end
        
        function setProxAlgo(this,optiType,bounds,maxiter,par,xtol,Outop)
            this.optimizer = optiType;
            if nargin<=2 || isempty(bounds),bounds = this.bounds;end
            if nargin<=3 || isempty(maxiter),maxiter = this.maxiter;end
            if nargin<=4 || isempty(par), par = [];end
            if nargin<=5 || isempty(xtol),xtol = [];end
            if nargin<=6 || isempty(Outop),Outop = OutputOpti(0);end
            
            this.bounds = bounds;
            this.maxiter = maxiter;
            
            switch this.optimizer
                case 'FGP'
                    LS = CostL2(this.sizein);
                    this.optim = OptiFGP(LS,this,this.bounds,Outop);
                    this.optim.Nesterov = true;
                    if ~isempty(par)
                        this.optim.L = par;
                    end
                    this.optim.maxiter = this.maxiter;
                    this.optim.ItUpOut = 0;
                    if ~isempty(xtol)
                        this.optim.xtol = xtol;
                    end
                case 'ADMM'
                    LS = CostL2(this.sizein);
                    if isempty(par)
                        par = 4e-1*[0.5,1];
                    end
                    Hn = {LinOpIdentity(this.sizein), this.H2};%circular boundary
                    NN = CostRectangle(this.sizein,this.bounds(1),this.bounds(2));
                    
                    this.optim = OptiADMM(LS,{NN, this.H1},Hn,par, [], Outop);
                    this.optim.maxiter = this.maxiter;
                    this.optim.ItUpOut = 0;
                    if ~isempty(xtol)
                        this.optim.xtol = xtol;
                    end
            end
        end
        
%         function setZ(this,new_z)
%             switch this.optimizer
%                 case 'FGP'
%                     this.optim.F0.set_y(new_z);
%                 case 'ADMM'
%                     this.optim.F0.set_y(new_z);
%             end
%         end
        
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this, x)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function y = apply_(this, x)
            % Reimplemented from parent class :class:`CostComposition`.
            % $$C(\\mathrm{x}) := \\sum_{k=1}^K \\sqrt{\\sum_{l=1}^L (\\mathrm{Hx}-y)_{k,l}^2}
            y = apply_@CostComposition(this, x);
            % Should the constraint set be added ? (i.e. cost indicator)
        end
        function y = applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`CostComposition`.
            % Implemented
            %  - if the operator \mathrm{H} is a LinOpGrad (i.e., this class is Total Variation)
            %    optimizer 'FGP', 'ADMM'
            
            % If the composed operator is a LinOpGrad & y==0
            if this.y==0
                switch this.optimizer
                    case 'FGP'
                        this.optim.setLambda(alpha);
                    case 'ADMM'
                        for k = 1:length(this.optim.Fn)
                            this.optim.Fn{k} = alpha*this.optim.Fn{k};
                        end
                end
                this.optim.F0.set_y(x);
                this.optim.run(x);
                y = this.optim.xopt;
            else
                error('Proximity operator not implemented for the composition CostMixNorm21 and %s with non-null this.y',this.H2.name);
            end
        end
        
        
        function M = sum_(this,H)
            if isa(H,'CostRectangle')
                this.bounds = [max(H.xmin,this.bounds(1)), min(H.xmax,this.bounds(end))]; %Combine existing bounds with new one
                if diff(this.bounds) < 0
                    error('There is no intersection between the combined set constraints. If you want to replace the already existing bounds, change this.bounds property');
                end
                M = this;
            end
        end
    end
end
