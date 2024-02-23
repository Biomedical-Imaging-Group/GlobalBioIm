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
    
    %% Properties common to mex and matlab implementations
    properties
        xmin=[];
        xmax=[];
        m=3;                % M is the number of correction pairs to remember in order to compute the limited memory variable metric (BFGS) approximation of the inverse of the Hessian.  For large problems, M = 3 to 5 gives good results.  For small problems, M should be less or equal N.  The larger is M (and N) the more computer memory will be needed to store the workspace WS.
        epsilon=0.01;       % a small value, in the range [0,1), equals to the cosine of the maximum angle between the search direction and the anti-gradient. The BFGS recursion is restarted, whenever the search direction is not sufficiently "descending".
        gtol=0.;            % Convergence occurs if the norm of gradient is lower than GTOL
        impl = 'mat'
    end
    properties (SetAccess = protected,GetAccess = protected)
        nbeval = 0;       % number of calls to `cost`    
        grad = [];        % gradient
        cc ;              % cost value
        fatol=0.0;        % absolute error desired in the function (e.g. FATOL=0.0). Convergence occurs if the estimate of the absolute error between F(X) and F(XSOL), where XSOL is a local minimizer, is less or equal FATOL. FATOL must have a non-negative floating point value.
        frtol=0.;         % relative error desired in the function (e.g.  FRTOL=1e-9). Convergence occurs if the estimate of the relative error between F(X) and F(XSOL), where XSOL is a local minimizer, is less or equal FRTOL. FRTOL must have a non-negative floating point value.
        task;
        bounds=false;         % flag indicating the constrained case
    end
    %% Properties specific to each implementation
    % -- mex implementation
    properties (Constant)
        OPL_TASK_START  = 0; % first entry, start search
        OPL_TASK_FG     = 1; % computation of F and G requested
        OPL_TASK_FREEVARS  = 2; % caller has to determine the free variables
        OPL_TASK_NEWX      = 3; % new variables available for inspection
        OPL_TASK_CONV      = 4; % search has converged
        OPL_TASK_WARN      = 5; % search aborted with warning
        OPL_TASK_ERROR     = 6; % search aborted with error
    end
    properties
        % Tolerance for the line search function
        sftol=0.001;
        sgtol=0.9;
        sxtol=0.1;
        delta=0.1;          %   DELTA is a small nonegative value used to compute a small initial step.
        active;
    end
    properties (SetAccess = protected,GetAccess = public)
        nparam;
        ws;
    end
    % -- matlab implementation
    properties (Constant)
        % Status indicating the termination of the algorithm
        OPL_STATUS_TOO_MANY_EVALUATIONS =  1;
        OPL_STATUS_TOO_MANY_ITERATIONS  =  2;
        OPL_STATUS_FTEST_SATISFIED =  3;
        OPL_STATUS_XTEST_SATISFIED =  4;
        OPL_STATUS_GTEST_SATISFIED =  5;
        OPL_STATUS_NOT_POSITIVE_DEFINITE = -1 ;

        % Stage flags
        OPL_STAGE_INIT = 0 ;                % initial task
        OPL_STAGE_LNSRCH_IN_PROGRESS = 1 ;  % line-search in progress
        OPL_STAGE_LNSRCH_CONV = 2 ;         % line-search has converged
    end
    properties   
        maxeval = Inf;      % Maximum number of evaluations. By default, it is unlimited.
        ftol = 1e-8;
        xtol = 1e-6;        % specify tolerances for deciding the
                            % convergence of the algorithm
                            
        lnsrch = [];        % specify line-search settings different than the
                            % default (see `optm_new_line_search`).
        verb = 0;	    % flag for "homemade" verbose mode
        f2nd = nan();
        fmin = nan();       
        dxrel = nan();
        dxabs = 1.0;        % These settings may be used to determine the step
                            % length along the steepest descent.

        blmvm = false;      % (false by default) specifies whether to use BLMVM trick to
                            % account for the bound constraints in the L-BFGS model of the Hessian.
    end
    properties (SetAccess = protected,GetAccess = protected)
        bounded=false;         % flag indicating the constrained case
        lbfgs ;          % lbfgs context   
        x0;
        f0 = Inf;       % function value at start of line-search
        g0 = [];         % gradient at start of line-search
        d = [];          % search direction
        s = [];          % effective step
        pg = [];         % projected gradient
        pg0 = [];        % projected gradient at start of line search
        gnorm = 0.0;     % Euclidean norm of the (pojected) gradient
        gtest;		 
        		 % Tolerances        
        gatol;		 % -
        grtol;		 % -
        xatol;		 % -
        xrtol;           % -
        		 % Tolerances   
        alpha = 0.0;     % step length
        amin = -Inf;     % first step length threshold
        amax = +Inf;     % last step length threshold
        iters = 0;       % number of iterations
        projs = 0;       % number of projections onto the feasible set
        rejects = 0;     % number of search direction rejections
        status = 0;      % non-zero when algorithm is about to terminate
        best_f = +Inf;   % function value at `best_x`
        best_g = [];     % gradient at `best_x`
        best_x = [];     % best solution found so far
        best_gnorm = -1; % norm of projected gradient at `best_x` (< 0 if unknown)
        best_alpha =  0; % step length at `best_x` (< 0 if unknown)
        best_evals = -1; % number of calls to `fg` at `best_x`
        last_evals = -1; % number of calls to `fg` at last iterate
        last_print = -1; % iteration number for last print
        freevars = [];   % subset of free variables (not yet known)
        t0;
    end
   
    %% Methods
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
            
        end
                
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.

            initialize@Opti(this,x0);
            if strcmp(this.impl,'mex')
                initialize_mex(this,x0);
            elseif strcmp(this.impl,'mat')
                initialize_mat(this,x0);
            end
        end
        function initialize_mat(this,x0)
            if ((exist('optm_new_line_search.m')~=2)...
                    ||(exist('optm_clamp.m')~=2)...
                    ||(exist('optm_new_line_search.m')~=2)...
                    ||(exist('optm_iterate_line_search.m')~=2)...
                    ||(exist('optm_norm2.m')~=2)...
                    ||(exist('optm_unblocked_variables.m')~=2)...
                    ||(exist('optm_update_lbfgs.m')~=2)...
                    ||(exist('optm_apply_lbfgs.m')~=2)...
                    ||(exist('optm_inner.m')~=2)...
                    ||(exist('optm_status.m')~=2))
                installVMLMB();
            end
            if ~this.bounds
                this.blmvm = false; % no needs to use BLMVM trick in the unconstrained case
            end
            if isempty(this.lnsrch)
                this.lnsrch = optm_new_line_search();
            end
            this.lbfgs = optm_new_lbfgs(this.m);
            if this.verb > 0
                time = @() 86400E3*now(); % yields number of milliseconds
                this.t0 = time();
            end

            % Tolerances.  Most of these are forced to be nonnegative to simplify tests.
            if isscalar(this.ftol)
                this.fatol = -Inf;
                this.frtol = max(0.0, this.ftol);
            else
                this.fatol = max(0.0, this.ftol(1));
                this.frtol = max(0.0, this.ftol(2));
            end
            if isscalar(this.gtol)
                this.gatol = 0.0;
                this.grtol = max(0.0, this.gtol);
            else
                this.gatol = max(0.0, this.gtol(1));
                this.grtol = max(0.0, this.gtol(2));
            end
            if isscalar(this.xtol)
                this.xatol = 0.0;
                this.xrtol = max(0.0, this.xtol);
            else
                this.xatol = max(0.0, this.xtol(1));
                this.xrtol = max(0.0, this.xtol(2));
            end

            this.x0 = x0 ;
            this.task = this.OPL_STAGE_INIT ;
        end
        function initialize_mex(this,x0)

            if (exist('m_opl_vmlmb_create')~=3)||(exist('m_opl_vmlmb_restore')~=3)||(exist('m_opl_vmlmb_iterate')~=3)||(exist('m_opl_vmlmb_get_reason')~=3)
                installOptimPack();
            end

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
            if strcmp(this.impl,'mex')
                flag = doIteration_mex(this);
            elseif strcmp(this.impl,'mat')
                flag = doIteration_mat(this);
            end
        end
        function flag=doIteration_mat(this)
            % Reimplementation from :class:`Opti`. For details see [1].

            flag=this.OPTI_REDO_IT;
            %-- Make the variables feasible.
            if this.bounds
                % In principle, we can avoid projecting the variables whenever
                % `alpha ≤ amin` (because the feasible set is convex) but rounding
                % errors could make this wrong.  It is safer to always project the
                % variables.  This cost O(n) operations which are probably
                % negligible compared to, say, computing the objective function
                % and its gradient.
                this.xopt = optm_clamp(this.xopt, this.xmin, this.xmax);
                this.projs = this.projs + 1;
            end
            %-- Compute objective function and its gradient.
            this.cc = gather(this.cost.apply(this.xopt));
            this.grad = gather(real(this.cost.applyGrad(this.xopt)));

            this.nbeval = this.nbeval + 1;
            if this.cc < this.best_f || this.nbeval == 1
                % Save best solution so far.
                this.best_f = this.cc;
                this.best_g = this.grad;
                this.best_x = this.xopt;
                this.best_gnorm = -1; % must be recomputed
                this.best_alpha = this.alpha;
                this.best_evals = this.nbeval;
            end
            if this.task ~= this.OPL_STAGE_INIT
                % Line-search in progress, check for line-search convergence.
                this.lnsrch = optm_iterate_line_search(this.lnsrch, this.cc);
                this.task = this.lnsrch.stage;
                if this.task ==  this.OPL_STAGE_LNSRCH_CONV
                    % Line-search has converged, `x` is the next iterate.
                    this.iters = this.iters + 1;
                    this.last_evals = this.nbeval;
                    flag=this.OPTI_NEXT_IT;
                elseif this.task == this.OPL_STAGE_LNSRCH_IN_PROGRESS
                    % Line-search has not converged, peek next trial step.
                    this.alpha = this.lnsrch.step;
                else
                    % Convergence, or error, or warning
                    this.endingMessage = ['something is wrong in line search !'];
                    flag=this.OPTI_STOP;
                    return;
                end
            end
            if this.task ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                %%Initial or next iterate after convergence of line-search.
                if this.bounds
                    % Determine the subset of free variables and compute the norm
                    % of the projected gradient (needed to check for convergence).
                    this.freevars = optm_unblocked_variables(this.xopt, this.xmin, this.xmax, this.grad);
                    this.pg = this.freevars .* this.grad;
                    this.gnorm = optm_norm2(this.pg);
                    if ~this.blmvm
                        % Projected gradient no longer needed, free some memory.
                        this.pg = [];
                    end
                else
                    % Just compute the norm of the gradient.
                    this.gnorm = optm_norm2(this.grad);
                end
                if this.nbeval == this.best_evals
                    % Now we know the norm of the (projected) gradient at the best
                    % solution so far.
                    this.best_gnorm = this.gnorm;
                end
                % Check for algorithm convergence or termination.
                if this.nbeval == 1
                    % Compute value for testing the convergence in the gradient.
                    this.gtest = max(this.gatol, this.grtol*this.gnorm);
                end
                if this.gnorm <= this.gtest
                    this.endingMessage = ['Gradient tolerance reached'];
                    % Convergence in gradient.
                    this.status = this.OPL_STATUS_GTEST_SATISFIED ; % optm_status('GTEST_SATISFIED');
                    flag=this.OPTI_STOP;
                    return;
                end
                if this.task == this.OPL_STAGE_LNSRCH_CONV
                    % Check convergence in relative function reduction.
                    if this.cc <= this.fatol || abs(this.cc - this.f0) <= this.frtol*max(abs(this.cc), abs(this.f0))
                        this.endingMessage = ['Relative function reduction tolerance reached'];
                        this.status = this.OPL_STATUS_FTEST_SATISFIED ; % optm_status('FTEST_SATISFIED');
                        flag=this.OPTI_STOP;
                        return;
                    end
                    % Compute the effective change of variables.
                    this.s = this.xopt - this.x0;
                    snorm = optm_norm2(this.s);
                    % Check convergence in variables.
                    if snorm <= this.xatol || (this.xrtol > 0 && snorm <= this.xrtol*optm_norm2(this.xopt))
                        this.endingMessage = ['Relative norm variable reduction tolerance reached'];
                        this.status = this.OPL_STATUS_XTEST_SATISFIED; % optm_status('XTEST_SATISFIED');
                        flag=this.OPTI_STOP;
                        return;
                    end
                end
                if this.iters >= this.maxiter
                    this.status = this.OPL_STATUS_TOO_MANY_ITERATIONS; % optm_status('TOO_MANY_ITERATIONS');
                    flag=this.OPTI_STOP;
                    return;
                end
            end
            if this.nbeval >= this.maxeval
                this.endingMessage = ['Max number of evaluations reached'];
                this.status = this.OPL_STATUS_TOO_MANY_EVALUATIONS; % optm_status('TOO_MANY_EVALUATIONS');
                flag=this.OPTI_STOP;
                return;
            end
            if this.task ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                % Possibly print iteration information.
                if this.verb > 0 && mod(this.iters, this.verb) == 0
                    this.print_iteration();
                    this.last_print = this.iters;
                end
                if this.task ~= this.OPL_STAGE_INIT
                    % At least one step has been performed, L-BFGS approximation
                    % can be updated.
                    if this.blmvm
                        this.lbfgs = optm_update_lbfgs(this.lbfgs, this.s, this.pg - this.pg0);
                    else
                        this.lbfgs = optm_update_lbfgs(this.lbfgs, this.s, this.grad - this.g0);
                    end
                end
                % Determine a new search direction `d`.  Parameter `dir` is set to:
                %   0 if `d` is not a search direction,
                %   1 if `d` is unscaled steepest descent,
                %   2 if `d` is scaled sufficient descent.
                dir = 0;
                % Use L-BFGS approximation to compute a search direction and check
                % that it is an acceptable descent direction.
                if this.blmvm
                    [this.d, scaled] = optm_apply_lbfgs(this.lbfgs, -this.pg);
                    this.d = this.d .* this.freevars;
                else
                    [this.d, scaled] = optm_apply_lbfgs(this.lbfgs, -this.grad, this.freevars);
                end
                dg = optm_inner(this.d, this.grad);
                if ~scaled
                    % No exploitable curvature information, `d` is the unscaled
                    % steepest feasible direction, that is the opposite of the
                    % projected gradient.
                    dir = 1;
                else
                    % Some exploitable curvature information were available.
                    dir = 2;
                    if dg >= 0
                        % L-BFGS approximation does not yield a descent direction.
                        dir = 0; % discard search direction
                        if ~this.bounds
                            this.status = this.OPL_STATUS_NOT_POSITIVE_DEFINITE ; % optm_status('NOT_POSITIVE_DEFINITE');
                            this.endingMessage = ['L-BFGS approximation is not positive definite !'];
                            flag=this.OPTI_STOP;
                            return;
                        end
                    elseif this.epsilon > 0
                        % A more restrictive criterion has been specified for
                        % accepting a descent direction.
                        if dg > -this.epsilon*optm_norm2(this.d)*this.gnorm
                            dir = 0; % discard search direction
                        end
                    end
                end
                if dir == 0
                    % No exploitable information about the Hessian is available or
                    % the direction computed using the L-BFGS approximation failed
                    % to be a sufficient descent direction.  Take the steepest
                    % feasible descent direction.
                    if this.bounds
                        this.d = -this.grad .* this.freevars;
                    else
                        this.d = -this.grad;
                    end
                    dg = -this.gnorm^2;
                    dir = 1; % scaling needed
                end
                if dir ~= 2 && this.iters > 0
                    this.rejects = this.rejects + 1;
                end
                % Determine the length `alpha` of the initial step along `d`.
                if dir == 2
                    % The search direction is already scaled.
                    this.alpha = 1.0;
                else
                    % Find a suitable step size along the steepest feasible
                    % descent direction `d`.  Note that `gnorm`, the Euclidean
                    % norm of the (projected) gradient, is also that of `d` in
                    % that case.
                    this.alpha = optm_steepest_descent_step(this.xopt, this.gnorm, this.cc, this.f2nd, this.fmin, ...
                        this.dxrel, this.dxabs);

                end
                if this.bounds
                    % Safeguard the step to avoid searching in a region where
                    % all bounds are overreached.
                    [this.amin, this.amax] = optm_line_search_limits(this.xopt, this.xmin, this.xmax, ...
                        this.d, this.alpha);
                    this.alpha = min(this.alpha, this.amax);
                end
                % Initialize line-search.
                this.lnsrch = optm_start_line_search(this.lnsrch, this.cc, dg, this.alpha);
                this.task = this.lnsrch.stage;
                if this.task ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                    % Convergence, or error, or warning
                    this.endingMessage = ['something is wrong in line search !'];
                    flag=this.OPTI_STOP;
                    return;
                end
                % Save iterate at start of line-search.
                this.f0 = this.cc;
                this.g0 = this.grad;
                this.x0 = this.xopt;
                if this.blmvm
                    this.pg0 = this.pg;
                end
            end
            % Compute next iterate.
            if this.alpha == 1
                this.xopt = this.x0 + this.d;
            else
                this.xopt = this.x0 + this.alpha*this.d;
            end
        end
        function flag=doIteration_mex(this)
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
                        this.endingMessage = ['Convergence: normg < gtol ',this.niter,this.nbeval,this.cc,normg,this.task];
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
