classdef OptiVMLMB<Opti
    % Variable Metric Limited Memory Bounded (VMLMB) [1] algorithm that
    % minimizes a cost \\(C(\\mathrm{x})\\) which is differentiable with bound
    % constraints and/or preconditioning.
    %
    % :param C: minimized cost
    % :param xmin: min bound (optional)
    % :param xmax: max bound (optional)
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    % **Note**
    % This Optimizer has many other variables that are set by
    % default to reasonable values. See the function m_vmlmb_first.m in the
    % MatlabOptimPack folder for more details.
    %
    % **Reference**
    %
    % [1] Eric Thiebaut, "Optimization issues in blind deconvolution algorithms",
    % SPIE Conf. Astronomical Data Analysis II, 4847, 174-183 (2002). See
    % `here <https://github.com/emmt/OptimPackLegacy>`_.
    %
    % **Example** VMLMB=OptiVMLMB(C,xmin,xmax,OutOp)
    %
    % See also :class:`Opti`, :class:`OptiConjGrad` :class:`OutputOpti`, :class:`Cost`
    
    %%    Copyright (C) 2017
    %     Ferreol Soulez ferreol.soulez@univ-lyon1.fr
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
    
    properties (Constant)
        OP_TASK_START  = 0; % first entry, start search
        OP_TASK_FG     = 1; % computation of F and G requested
        OP_TASK_NEWX   = 2; % new improved solution available for inspection
        OP_TASK_CONV   = 3; % search has converged
        OP_TASK_WARN   = 4; % search aborted with warning
        OP_TASK_ERROR  = 5; % search aborted with error
    end
    % Full public properties
    properties
        m=3;
        gtol=0;
        fatol=0.0;
        frtol=1e-8;
        sftol=0.001;
        sgtol=0.9;
        sxtol=0.1;
        epsilon=0.01;
        costheta=0.4;
        
        doFullGradient = 1; % boolean. If false, F gradient is computed from a subset of "angles" (requires F = CostSummation)
        stochastic_gradient = 0; % boolean, stochastic gradient descent rule (requires F = CostSummation and doFullGradient = false)
        
        Lsub; % Number of "angles" used if doFullGradient==0
        subset; % current subset of angles used to compute F grad
        nrepeat; % number of repetition of one subset
    end
    properties (SetAccess = protected,GetAccess = public)
        nparam;
        task;
        xmin=[];
        xmax=[];
        neval;
        bounds =0;
        csave;
        isave;
        dsave;
        
        set; %indices of mapsCell (CostSummation) to consider as an "angle" (e.g. F is composed of 100 CostL2 (i.e., angles) + 1 CostHyperbolic)
        nonset; %indices of mapsCell (CostSummation) not to consider as an "angle"
        L; % Total number of available "angles"
        counter; % counter for subset update
        count4rep; % counter for subset repeat
    end
    properties (SetAccess = protected,GetAccess = protected)
        nbeval;
        x;
    end
    methods
        function this = OptiVMLMB(C,xmin,xmax,OutOp)
            this.name='OptiVMLMB';
            if(nargin>1)
                if(~isempty(xmin))
                    this.bounds=1;
                    this.xmin = xmin;
                end
                if(~isempty(xmax))
                    this.bounds=bitor(this.bounds,2);
                    this.xmax = xmax;
                end
                if nargin==4 && ~isempty(OutOp),this.OutOp=OutOp;end
            end
            this.cost=C;
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            this.nparam =numel(x0);
            if isscalar(this.xmin)
                this.xmin=ones(size(x0))*this.xmin;
            end
            if isscalar(this.xmax)
                this.xmax=ones(size(x0))*this.xmax;
            end
            [this.csave, this.isave, this.dsave] = m_vmlmb_first(this.nparam, this.m, this.fatol, this.frtol,...
                this.sftol, this.sgtol, this.sxtol, this.epsilon, this.costheta);
            this.task =  this.isave(3);
            this.task = this.OP_TASK_FG;
            this.nbeval=0;
            this.x = x0;
            this.x(1)= x0(1); %Warning : side effect on x0 if x=x0 (due to the fact that x is passed-by-reference in the mexfiles)
        end
        
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].
            
            flag=0;
            if (this.task == this.OP_TASK_FG)
                % apply bound constraints
                % op_bounds_apply(n, x, xmin, xmax);
                if(bitand(this.bounds,1))
                    test = (this.x<this.xmin);
                    if any(test(:)), this.x(test) = this.xmin(test); end
                end
                if (bitand(this.bounds,2))
                    test = (this.x>this.xmax);
                    if any(test(:)), this.x(test) = this.xmax(test); end
                end
                
                
                %cost = this.cost.apply(x);
                %grad = this.cost.applyGrad(x);
                
                cost = this.computeCost(this.x);
                grad = this.computeGrad(this.x);
                
                %[cost,grad] = this.cost.eval_grad(x);     % evaluate the function and its gradient at X;
                normg= sum(grad(:).^2);
                
                this.nbeval=this.nbeval+1;
                if (normg< this.gtol)
                    fprintf('Convergence: normg < gtol \n %d\t%d\t%7.2e\t%6.2g\t\t%d\t%d \n',this.niter,this.nbeval,cost,normg,this.task,this.isave(4));
                    this.time=toc(tstart);
                    this.ending_verb();
                    flag=1;
                end
            elseif (this.task == this.OP_TASK_NEWX)
                this.niter = this.niter +1;
                this.xopt = this.x;
                if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end% New successful step: the approximation X, function F, and
                % gradient G, are available for inspection.
            else
                % Convergence, or error, or warning
                fprintf('Convergence, or error, or warning : %d  , %s\n',this.task,this.csave);
                this.time=toc(tstart);
                this.ending_verb();
                flag=1;
            end
            if ( (this.nbeval==1) || (this.task == this.OP_TASK_NEWX))
                % Computes set of active parameters :
                % op_bounds_active(n, active, x, g, xmin, xmax);
                switch(this.bounds)
                    case 0
                        active = [];
                    case 1
                        active = int32( (this.x>this.xmin) + (grad<0) );
                    case 2
                        active = int32( (grad>0) + (this.x<this.xmax) );
                    case 3
                        active = int32( ( (this.x>this.xmin) + (grad<0) ).*( (this.x<this.xmax) + (grad>0) ) );
                end
            end
            % Computes next step:
            [this.task, this.csave]= m_vmlmb_next(this.x,cost,grad,active,this.isave,this.dsave);
            
            %                 if (this.niter==this.maxiter)
            %                     this.ending_verb();
            %                 end
        end
        
        
        function grad = computeGrad(this,x)
            
            if this.Lsub < this.L
                %                 if isempty(this.subset)
                %                     this.updateSubset(sort(1 + mod(round(this.counter...
                %                         + (1:this.L/this.Lsub:this.L)),this.L)));
                %                 end
                grad = zeros(size(x));
                for kk = 1:this.Lsub
                    ind = this.set(this.subset(kk));
                    grad = grad + this.cost.alpha(ind)*this.cost.mapsCell{ind}.applyGrad(x);
                end
                grad = real(grad)/this.Lsub;%ad hoc
                
                for kk = 1:length(this.nonset)
                    grad = grad + this.cost.alpha(this.nonset(kk))*this.cost.mapsCell{this.nonset(kk)}.applyGrad(x);
                end
                
                if mod(this.count4rep,this.nrepeat) == this.nrepeat - 1
                    if this.stochastic_gradient
                        this.updateSubset(randi(this.L,this.Lsub,1));
                    else
                        this.updateSubset(sort(1 + mod(round(this.counter...
                            + (1:this.L/this.Lsub:this.L)),this.L)));
                        this.counter = this.counter + 1;
                    end
                    this.count4rep = 0;
                else
                    this.count4rep = this.count4rep + 1;
                end
            else
                grad = this.cost.applyGrad(x);
            end
        end
        
        function cost = computeCost(this,x)
            
            if this.Lsub < this.L
                if isempty(this.subset)
                    this.updateSubset(sort(1 + mod(round(this.counter...
                        + (1:this.L/this.Lsub:this.L)),this.L)));
                end
                cost = 0;
                for kk = 1:this.Lsub
                    ind = this.set(this.subset(kk));
                    cost = cost + this.cost.alpha(ind)*this.cost.mapsCell{ind}.apply(x);
                end
                cost = real(cost)/this.Lsub;%ad hoc
                
                for kk = 1:length(this.nonset)
                    cost = cost + this.cost.alpha(this.nonset(kk))*this.cost.mapsCell{this.nonset(kk)}.apply(x);
                end
            else
                cost = this.cost.apply(x);
            end
        end
        
        function updateSet(this,new_set)
            this.set = new_set;
            this.L = length(new_set);
            this.nonset = find(~ismember(1:this.cost.numMaps,this.set));
        end
        
        function updateSubset(this,subset)
            this.subset = subset;
            this.Lsub = length(subset);
        end
        
        function reset(this)
            this.counter = 0;
            this.count4rep = 0;
            if this.stochastic_gradient
                this.updateSubset(randi(this.L,this.Lsub,1));
            else
                this.updateSubset(unique(round(1:this.L/this.Lsub:this.L)));
            end
        end
    end
end
