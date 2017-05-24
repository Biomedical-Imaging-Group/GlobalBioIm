classdef OptiVMLMB<Opti
        %% OptiVMLMB : VMLMB Optimizer
    %  Matlab inverse Problems Library
    %
    % Variable Metric Limited Memory Bounded (VMLMB) algorithm by Éric Thiébaut.
    % It is a limited memory BFGS (variable metric) method possibly with bound 
    % constraints and/or preconditioning.
    % See https://github.com/emmt/OptimPackLegacy
    %
    %
    % -- References
    % Éric Thiébaut, "Optimization issues in blind deconvolution algorithms", SPIE Conf. Astronomical Data Analysis II, 4847, 174-183 (2002).
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OutputOpti
    %
    %     Copyright (C) 2017 Ferreol Soulez ferreol.soulez@univ-lyon1.fr
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
        function this = OptiVMLMB(nparam,xmin,xmax)
            
    		this.name='OptiVMLMB'; 
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
        function bestx =run(this,F,x0)
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
                        if any(test(:)), x(test) = this.xmin(test); end
                    end
                    if (bitand(this.bounds,2))
                        test = (x>this.xmax);
                        if any(test(:)), x(test) = this.xmax(test); end
                    end
                    this.cost = F.eval(x);     % evaluate the function at X;
                    grad = F.grad(x);   % evaluate the gradient of F at X;
                    normg= sum(grad(:).^2);

                    nbeval=nbeval+1;
                    if (normg< this.gtol)
                        fprintf('Convergence: normg < gtol \n %d\t%d\t%7.2e\t%6.2g\t\t%d\t%d \n',iter,nbeval,this.cost,normg,this.task,this.isave(4));
                        break;
                    end
                elseif (this.task == this.OP_TASK_NEWX)
                    iter = iter +1;
                    bestx = x;
                    if (mod(iter,this.verb)==0)
                        fprintf('%d\t%d\t%7.2e\t%6.2g\t\t%d\t%d \n',iter,nbeval,this.cost,normg,this.task,this.isave(4));
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
                            active = int32( (x>this.xmin) + (grad<0) );
                        case 2
                            active = int32( (grad>0) + (x<this.xmax) );
                        case 3
                            active = int32( ( (x>this.xmin) + (grad<0) ).*( (x<this.xmax) + (grad>0) ) );
                    end
                end
                % Computes next step:
                [this.task, this.csave]= m_vmlmb_next(x,this.cost,grad,active,this.isave,this.dsave);
                
            end
            
        end
    end
end
