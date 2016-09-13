classdef VMLMB<handle
    properties (Constant)
        OP_TASK_START  = 0; % first entry, start search
        OP_TASK_FG     = 1; % computation of F and G requested
        OP_TASK_NEWX   = 2; % new improved solution available for inspection
        OP_TASK_CONV   = 3; % search has converged
        OP_TASK_WARN   = 4; % search aborted with warning
        OP_TASK_ERROR  = 5; % search aborted with error
    end
    properties
        m=3;
        fatol=0.0;
        frtol=1e-8;
        sftol=0.001;
        sgtol=0.9;
        sxtol=0.1;
        epsilon=0.01;
        costheta=0.4;
        nbitermax=10;
        nbevalmax;
        verb=0;
    end
    properties (SetAccess = protected,GetAccess = public)
        nparam;
        task;
        xmin=[];
        xmax=[];
        neval;
        niter;
        bounds =0;
        csave;
        isave;
        dsave;
    end
    methods
        function this = VMLMB(nparam,xmin,xmax)
            this.nparam =nparam;
            if(nargin>1)
                if(~isempty(xmin))
                    this.bounds=1;
                    this.xmin = xmin;
                end
                if(~isempty(xmax))
                    this.bounds=bitor(this.bounds,2);
                    this.xmax = xmax;
                end
            end
            [this.csave, this.isave, this.dsave] = m_vmlmb_first(this.nparam, this.m, this.fatol, this.frtol,...
                this.sftol, this.sgtol, this.sxtol, this.epsilon, this.costheta);
            this.task =  this.isave(3);
        end
        function ReInit(this)
            [this.csave, this.isave, this.dsave] = m_vmlmb_first(this.nparam, this.m, this.fatol, this.frtol,...
                this.sftol, this.sgtol, this.sxtol, this.epsilon, this.costheta);
            this.task =  this.isave(3);
        end
        function bestx =Optimize(this,F,x0)
            bestx = x0;
            x = x0;
            x(1)= x0(1); % Warning : side effect on x0 if x=x0 (due to the fact that x is passed-by-reference in the mexfiles)
            this.task = this.OP_TASK_FG;
            nbeval=0;
            iter =1;
            if(this.verb)
                fprintf('it\t nbeval\t cost\t\t  normg\t\t task\tstage\n');
            end
            while(iter< this.nbitermax)
                if (this.task == this.OP_TASK_FG)
                    % apply bound constraints
                    % op_bounds_apply(n, x, xmin, xmax);
                    if(bitand(this.bounds,1))
                        test = (x<this.xmin);
                        if any(test), x(test) = this.xmin(test); end
                    end
                    if (bitand(this.bounds,2))
                        test = (x>this.xmax);
                        if any(test), x(test) = this.xmax(test); end
                    end
                    cost = F.GetCost(x);     % evaluate the function at X;
                    grad = F.GetGradient();   % evaluate the gradient of F at X;
                    normg= sum(grad(:).^2);
                    nbeval=nbeval+1;
                elseif (this.task == this.OP_TASK_NEWX)
                    iter = iter +1;
                    bestx = x;
                    F.UpdateLagrangians();
                    if (mod(iter,this.verb)==0)
                        fprintf('%d\t%d\t%7.2e\t%6.2g\t\t%d\t%d \n',iter,nbeval,cost,normg,this.task,this.isave(4));
                    end
                    % New successful step: the approximation X, function F, and
                    % gradient G, are available for inspection.
                else
                    % Convergence, or error, or warning
                    fprintf('Convergence, or error, or warning : %d  , %s\n',this.task,this.csave);
                    break;
                end
                if ( (nbeval==1) || (this.task == this.OP_TASK_NEWX))
                    % Computes set of active parameters :
                    % op_bounds_active(n, active, x, g, xmin, xmax);
                    switch(this.bounds)
                        case 0
                            active = [];
                        case 1
                            active = int32( (grad>0) + (x<this.xmax) );
                        case 2
                            active = int32( (x>this.xmin) + (grad<0) );
                        case 3
                            active = int32( ( (x>this.xmin) + (grad<0) ).*( (x<this.xmax) + (grad>0) ) );
                    end
                end
                % Computes next step:
                [this.task, this.csave]= m_vmlmb_next(x,cost,grad,active,this.isave,this.dsave);
                
            end
            
        end
    end
end