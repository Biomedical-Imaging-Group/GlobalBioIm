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
        function Init(this)    
            [this.csave, this.isave, this.dsave] = m_vmlmb_first(this.nparam, this.m, this.fatol, this.frtol,...
                this.sftol, this.sgtol, this.sxtol, this.epsilon, this.costheta);
            this.task =  this.isave(3);
        end
        function  run(this,x0)
            % Reimplementation from :class:`Opti`. For details see [1].

            this.nparam =numel(x0);
            if isscalar(this.xmin)
                this.xmin=ones(size(x0))*this.xmin;
            end
            if isscalar(this.xmax)
                this.xmax=ones(size(x0))*this.xmax;
            end
            this.Init();
            this.xopt = x0;
            x = x0;
            x(1)= x0(1); % Warning : side effect on x0 if x=x0 (due to the fact that x is passed-by-reference in the mexfiles)
            this.task = this.OP_TASK_FG;
            tstart=tic;
            this.OutOp.init();
            nbeval=0;
            this.niter =1;
            starting_verb(this);
            while(this.niter< this.maxiter)
                if (this.task == this.OP_TASK_FG)
                    % apply bound constraints
                    % op_bounds_apply(n, x, xmin, xmax);
                    if(bitand(this.bounds,1))
                        test = (x<this.xmin);
                        if any(test(:)), x(test) = this.xmin(test); end
                    end
                    if (bitand(this.bounds,2))
                        test = (x>this.xmax);
                        if any(test(:)), x(test) = this.xmax(test); end
                    end
                    [cost,grad] = this.cost.eval_grad(x);     % evaluate the function and its gradient at X;
                    normg= sum(grad(:).^2);
                    
                    nbeval=nbeval+1;
                    if (normg< this.gtol)
                        fprintf('Convergence: normg < gtol \n %d\t%d\t%7.2e\t%6.2g\t\t%d\t%d \n',this.niter,nbeval,cost,normg,this.task,this.isave(4));
                        this.time=toc(tstart);
                        this.ending_verb();
                        break;
                    end
                elseif (this.task == this.OP_TASK_NEWX)
                    this.niter = this.niter +1;
                    this.xopt = x;
                    if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end% New successful step: the approximation X, function F, and
                    % gradient G, are available for inspection.
                else
                    % Convergence, or error, or warning
                    fprintf('Convergence, or error, or warning : %d  , %s\n',this.task,this.csave);                    
                    this.time=toc(tstart);
                    this.ending_verb();
                    break;
                end
                if ( (nbeval==1) || (this.task == this.OP_TASK_NEWX))
                    % Computes set of active parameters :
                    % op_bounds_active(n, active, x, g, xmin, xmax);
                    switch(this.bounds)
                        case 0
                            active = [];
                        case 1
                            active = int32( (x>this.xmin) + (grad<0) );
                        case 2
                            active = int32( (grad>0) + (x<this.xmax) );
                        case 3
                            active = int32( ( (x>this.xmin) + (grad<0) ).*( (x<this.xmax) + (grad>0) ) );
                    end
                end
                % Computes next step:
                [this.task, this.csave]= m_vmlmb_next(x,cost,grad,active,this.isave,this.dsave);
                
                if (this.niter==this.maxiter)
                    this.ending_verb();
                end
            end
            
        end
    end
end
