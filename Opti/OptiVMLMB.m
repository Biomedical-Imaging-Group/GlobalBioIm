classdef OptiVMLMB<Opti
    % Variable Metric Limited Memory Bounded (VMLMB) from OptimPackLegacy [1].
    % This algorithm
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
    % SPIE Conf. Astronomical Data Analysis II, 4847, 174-183 (2002).
    % See OptimPackLegacy `repository <https://github.com/emmt/OptimPackLegacy>`_.
    %
    % **Example** VMLMB=OptiVMLMB(C,xmin,xmax)
    %
    % See also :class:`Opti`, :class:`OptiConjGrad` :class:`OutputOpti`, :class:`Cost`
    
    %% GUI-Header
    % GUInotation-
    % GUIcall-OptiVMLMB(C,xmin,xmax)-
    % GUIparam-xmin-vecInt-[]-min bound
    % GUIparam-xmax-vecInt-[]-max bound
    
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
        OPL_TASK_START  = 0; % first entry, start search
        OPL_TASK_FG     = 1; % computation of F and G requested
        OPL_TASK_FREEVARS  = 2; % caller has to determine the free variables
        OPL_TASK_NEWX      = 3; % new variables available for inspection
        OPL_TASK_CONV      = 4; % search has converged
        OPL_TASK_WARN      = 5; % search aborted with warning
        OPL_TASK_ERROR     = 6; % search aborted with error
    end
    % Full public properties
    properties
        m=3;                %  M is the number of correction pairs to remember in order to compute the limited memory variable metric (BFGS) approximation of the inverse of the Hessian.  For large problems, M = 3 to 5 gives good results.  For small problems, M should be less or equal N.  The larger is M (and N) the more computer memory will be needed to store the workspace WS.
        fatol=0.0;          % absolute error desired in the function (e.g. FATOL=0.0). Convergence occurs if the estimate of the absolute error between F(X) and F(XSOL), where XSOL is a local minimizer, is less or equal FATOL. FATOL must have a non-negative floating point value.
        frtol=0.;           % relative error desired in the function (e.g.  FRTOL=1e-9). Convergence occurs if the estimate of the relative error between F(X) and F(XSOL), where XSOL is a local minimizer, is less or equal FRTOL. FRTOL must have a non-negative floating point value.
        gtol=0;             % Convergence occurs if the norm of gradient is lower than GTOL
        % Tolerance for the line search function
        sftol=0.001;
        sgtol=0.9;
        sxtol=0.1;
        epsilon=0.01;       % a small value, in the range [0,1), equals to the cosine of the maximum angle between the search direction and the anti-gradient. The BFGS recursion is restarted, whenever the search direction is not sufficiently "descending".
        delta=0.1;          %   DELTA is a small nonegative value used to compute a small initial step.
        xmin=[];
        xmax=[];
        active;
    end
    properties (SetAccess = protected,GetAccess = public)
        nparam;
        task;
        neval;
        bounds =0;
        ws;
    end
    properties (SetAccess = protected,GetAccess = protected)
        nbeval;
        %active;
        grad;
        cc;
    end
    methods
        function this = OptiVMLMB(C,xmin,xmax)
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
            end
            this.cost=C;
            if (exist('m_opl_vmlmb_create')~=3)||(exist('m_opl_vmlmb_restore')~=3)||(exist('m_opl_vmlmb_iterate')~=3)||(exist('m_opl_vmlmb_get_reason')~=3)
                installOptimPack();
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            this.nparam =numel(x0);
            if isscalar(this.xmin)
                this.xmin=ones_(size(x0))*this.xmin;
            end
            if isscalar(this.xmax)
                this.xmax=ones_(size(x0))*this.xmax;
            end
            this.ws = m_opl_vmlmb_create(this.nparam, this.m, this.fatol, this.frtol,...
                this.sftol, this.sgtol, this.sxtol, this.epsilon, this.delta);
            
            this.task = this.OPL_TASK_FG;
            this.nbeval=0;
            this.xopt = x0;
            this.xopt(1)= x0(1); %Warning : side effect on x0 if x=x0 (due to the fact that x is passed-by-reference in the mexfiles)
            
            % apply bound constraints
            % op_bounds_apply(n, x, xmin, xmax);
            if(bitand(this.bounds,1))
                test = (this.xopt<this.xmin);
                if any(test(:)), this.xopt(test) = this.xmin(test); end
            end
            if (bitand(this.bounds,2))
                test = (this.xopt>this.xmax);
                if any(test(:)), this.xopt(test) = this.xmax(test); end
            end
            
            
            this.cc = gather(this.cost.apply(this.xopt));
            this.grad = gather(real(this.cost.applyGrad(this.xopt)));
            
            this.nbeval=this.nbeval+1;
            
        end
        
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].
            
            % Computes next step:
            this.xopt(1) = this.xopt(1); %Warning : side effect on x0 if x=x0 (due to the fact that x is passed-by-reference in the mexfiles)
            this.task = m_opl_vmlmb_iterate(this.ws,gather(this.xopt),this.cc,this.grad,this.active);
            
            flag=this.OPTI_REDO_IT;
            if (this.task == this.OPL_TASK_FG)
                % apply bound constraints
                % op_bounds_apply(n, x, xmin, xmax);
                if(bitand(this.bounds,1))
                    test = (this.xopt<this.xmin);
                    if any(test(:)), this.xopt(test) = this.xmin(test); end
                end
                if (bitand(this.bounds,2))
                    test = (this.xopt>this.xmax);
                    if any(test(:)), this.xopt(test) = this.xmax(test); end
                end
                
                
                this.cc = gather(this.cost.apply(this.xopt));
                this.grad = gather(real(this.cost.applyGrad(this.xopt)));
                
                
                this.nbeval=this.nbeval+1;
                
                if this.gtol>0
                    normg= sum(this.grad(:).^2);
                    if (normg< this.gtol)
                        this.endMessage = ['Convergence: normg < gtol ',this.niter,this.nbeval,this.cc,normg,this.task];
                        flag=this.OPTI_STOP;
                    end
                    
                end
            elseif (this.task == this.OPL_TASK_NEWX)
                flag=this.OPTI_NEXT_IT;
                
            elseif (this.task == this.OPL_TASK_FREEVARS)
                % Computes set of active parameters :
                % op_bounds_active(n, active, x, g, xmin, xmax);
                switch(this.bounds)
                    case 0
                        this.active = [];
                    case 1
                        this.active = int32( (this.xopt>this.xmin) + (this.grad<0) );
                    case 2
                        this.active = int32( (this.grad>0) + (this.xopt<this.xmax) );
                    case 3
                        this.active = int32( ( (this.xopt>this.xmin) + (this.grad<0) ).*( (this.xopt<this.xmax) + (this.grad>0) ) );
                end
                this.active = gather(this.active);
            else
                % Convergence, or error, or warning
                this.endingMessage = ['Convergence, or error, or warning : ',this.task,m_opl_vmlmb_get_reason(this.ws)];
                
                flag=this.OPTI_STOP;
                this.task = m_opl_vmlmb_restore(this.ws,gather(this.xopt),this.cc,this.grad);
                
                return;
            end
        end
    end
    
end
