classdef OptiVMLMB_v2<Opti
    % Variable Metric Limited Memory Bounded (VMLMB) from VMLMB [1].
    % This algorithm
    % minimizes a cost \\(C(\\mathrm{x})\\) which is differentiable with bound
    % constraints and/or preconditioning.
    %
    % :param C: minimized cost
    % :param lower: min bound (optional)
    % :param upper: max bound (optional)
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    % **Note**
    %
    % This Optimizer is a new version of the initial VMLMB Optimizer using 
    % the OptimPackLegacy package (see OptiVMLMB). This one use a 
    % "purely-Matlab" implementation available in 
    % the VMLMB `repository <https://github.com/emmt/VMLMB>`_.
    % This implementation is based on the whole VMLMB algorithm that is found 
    % in the function optm_vmlmb.m in the VMLMB/matlab/src
    %
    % This Optimizer has many other variables that are set by
    % default to reasonable values. See the function optm_vmlmb.m in the
    % VMLMB/matlab/src folder for more details.
    %
    % **Reference**
    %
    % [1] Eric Thiebaut, "Optimization issues in blind deconvolution algorithms",
    % SPIE Conf. Astronomical Data Analysis II, 4847, 174-183 (2002).
    % See VMLMB `repository <https://github.com/emmt/VMLMB>`_.
    %
    % **Example** VMLMB=OptiVMLMB_v2(C,lower,upper)
    %
    % See also :class:`Opti`, :class:`OptiVMLMB`, :class:`OptiConjGrad`, 
    %	       :class:`OutputOpti`, :class:`Cost`

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
        % Status indicating the termination of the algorithm
        OPL_STATUS_TOO_MANY_EVALUATIONS =  1;
        OPL_STATUS_TOO_MANY_ITERATIONS  =  2;
        OPL_STATUS_FTEST_SATISFIED =  3;
        OPL_STATUS_XTEST_SATISFIED =  4;
        OPL_STATUS_GTEST_SATISFIED =  5;
        OPL_STATUS_NOT_POSITIVE_DEFINITE = -1 ;

        % Stage flags
        OPL_STAGE_INIT = 0 ;                % initial stage
        OPL_STAGE_LNSRCH_IN_PROGRESS = 1 ;  % line-search in progress
        OPL_STAGE_LNSRCH_CONV = 2 ;         % line-search has converged
    end
    % Full public properties
    properties
        lower = [];
        upper = []; % lower and/or an upper bounds
                    % for the variables.  If unspecified or set to an empty array, a given bound
                    % is considered as unlimited.  Bounds must be conformable with the
                    % variables.
        mem = 5;            % specifies the memory used by the algorithm, that is the
                            % number of previous steps memorized to approximate the Hessian of the
                            % objective function.  With `mem=0`, the algorithm behaves as a steepest
                            % descent method.  The default is `mem=5`.
        maxeval = Inf;      % Maximum number of evaluations. By default, it is unlimited.
        ftol = 1e-8;
        gtol = 1e-5;
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
        epsilon = 0.0;      % specifies a threshold for a sufficient descent
                            % condition.
        blmvm = false;      % (false by default) specifies whether to use BLMVM trick to
                            % account for the bound constraints in the L-BFGS model of the Hessian.
    end
    properties (SetAccess = protected,GetAccess = protected)
        %% Other initialization.
        stage;           % current saage of the algorithm
        bounded=false;         % flag indicating the constrained case
        f ;              % cost value
        lbfgs ;          % lbfgs context   
        g = [];          % gradient
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
        fatol;		 % -	
        frtol;		 % -
        gatol;		 % -
        grtol;		 % -
        xatol;		 % -
        xrtol;           % -
        		 % Tolerances   
        alpha = 0.0;     % step length
        amin = -Inf;     % first step length threshold
        amax = +Inf;     % last step length threshold
        evals = 0;       % number of calls to `cost`
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
    methods
        function this = OptiVMLMB_v2(C,lower,upper,varargin)

            this.name = 'OptiVMLMB_v2';
            this.cost = C;

            % Install VMLMB sources if not available
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

            if(nargin>1)
                if(~isempty(lower))
                    this.bounded=1;
                    this.lower = lower;
                end
                if(~isempty(upper))
                    this.bounded=bitor(this.bounded,2);
                    this.upper = upper;
                end
            end
            if ~this.bounded
                this.blmvm = false; % no needs to use BLMVM trick in the unconstrained case
            end

            if (nargin>3 && ~isempty(varargin{1}))
                this.mem = varargin{1} ;
            end
            
            if (nargin>4 && ~isempty(varargin{2}))
                this.maxiter = varargin{2} ;
            end

            if (nargin>5 && ~isempty(varargin{3}))
                this.maxeval = varargin{3} ;
            end

            %%
            %% TODO: MANAGE THE OTHER ALGORITHM PARAMETERS ? => LET TO DEFAULT FOR THE MOMENT
            %%

            if isempty(this.lnsrch)
                this.lnsrch = optm_new_line_search();
            end
            this.lbfgs = optm_new_lbfgs(this.mem);
            if this.verb > 0
                time = @() 86400E3*now(); % yields number of milliseconds
                this.t0 = time();
            end

            %% Tolerances.  Most of these are forced to be nonnegative to simplify
            %% tests.
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

        end

        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.      
            initialize@Opti(this,x0);
            this.x0 = x0 ;
            this.stage = this.OPL_STAGE_INIT ;
        end
        
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].

            flag=this.OPTI_REDO_IT;
            %% Make the variables feasible.
            if this.bounded
                %% In principle, we can avoid projecting the variables whenever
                %% `alpha ≤ amin` (because the feasible set is convex) but rounding
                %% errors could make this wrong.  It is safer to always project the
                %% variables.  This cost O(n) operations which are probably
                %% negligible compared to, say, computing the objective function
                %% and its gradient.
                this.xopt = optm_clamp(this.xopt, this.lower, this.upper);
                this.projs = this.projs + 1;
            end
            %% Compute objective function and its gradient.
            this.f = gather(this.cost.apply(this.xopt));
            this.g = gather(real(this.cost.applyGrad(this.xopt)));

            this.evals = this.evals + 1;
            if this.f < this.best_f || this.evals == 1
                %% Save best solution so far.
                this.best_f = this.f;
                this.best_g = this.g;
                this.best_x = this.xopt;
                this.best_gnorm = -1; % must be recomputed
                this.best_alpha = this.alpha;
                this.best_evals = this.evals;
            end
            if this.stage ~= this.OPL_STAGE_INIT
                %% Line-search in progress, check for line-search convergence.
                this.lnsrch = optm_iterate_line_search(this.lnsrch, this.f);
                this.stage = this.lnsrch.stage;
                if this.stage ==  this.OPL_STAGE_LNSRCH_CONV
                    %% Line-search has converged, `x` is the next iterate.
                    this.iters = this.iters + 1;
                    this.last_evals = this.evals;
                    flag=this.OPTI_NEXT_IT;
                elseif this.stage == this.OPL_STAGE_LNSRCH_IN_PROGRESS
                    %% Line-search has not converged, peek next trial step.
                    this.alpha = this.lnsrch.step;
                else
                    % Convergence, or error, or warning
                    this.endingMessage = ['something is wrong in line search !'];
                    flag=this.OPTI_STOP;
                    return;
                end
            end
            if this.stage ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                %% Initial or next iterate after convergence of line-search.
                if this.bounded
                    %% Determine the subset of free variables and compute the norm
                    %% of the projected gradient (needed to check for convergence).
                    this.freevars = optm_unblocked_variables(this.xopt, this.lower, this.upper, this.g);
                    this.pg = this.freevars .* this.g;
                    this.gnorm = optm_norm2(this.pg);
                    if ~this.blmvm
                        %% Projected gradient no longer needed, free some memory.
                        this.pg = [];
                    end
                else
                    %% Just compute the norm of the gradient.
                    this.gnorm = optm_norm2(this.g);  
                end
                if this.evals == this.best_evals
                    %% Now we know the norm of the (projected) gradient at the best
                    %% solution so far.
                    this.best_gnorm = this.gnorm;
                end
                %% Check for algorithm convergence or termination.
                if this.evals == 1
                    %% Compute value for testing the convergence in the gradient.
                    this.gtest = max(this.gatol, this.grtol*this.gnorm);
                end
                if this.gnorm <= this.gtest

                    %% Convergence in gradient.
                    this.status = this.OPL_STATUS_GTEST_SATISFIED ; % optm_status('GTEST_SATISFIED');
                    flag=this.OPTI_STOP;
                    return;
                end
                if this.stage == this.OPL_STAGE_LNSRCH_CONV
                    %% Check convergence in relative function reduction.
                    if this.f <= this.fatol || abs(this.f - this.f0) <= this.frtol*max(abs(this.f), abs(this.f0))
                        this.status = this.OPL_STATUS_FTEST_SATISFIED ; % optm_status('FTEST_SATISFIED');
                        flag=this.OPTI_STOP;
                        return;
                    end
                    %% Compute the effective change of variables.
                    this.s = this.xopt - this.x0;
                    snorm = optm_norm2(this.s);
                    %% Check convergence in variables.
                    if snorm <= this.xatol || (this.xrtol > 0 && snorm <= this.xrtol*optm_norm2(this.xopt))
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
            if this.evals >= this.maxeval
                this.status = this.OPL_STATUS_TOO_MANY_EVALUATIONS; % optm_status('TOO_MANY_EVALUATIONS');
                flag=this.OPTI_STOP;
                return;
            end
            if this.stage ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                %% Possibly print iteration information.
                if this.verb > 0 && mod(this.iters, this.verb) == 0
                    this.print_iteration();
                    this.last_print = this.iters;
                end
                if this.stage ~= this.OPL_STAGE_INIT
                    %% At least one step has been performed, L-BFGS approximation
                    %% can be updated.
                    if this.blmvm
                        this.lbfgs = optm_update_lbfgs(this.lbfgs, this.s, this.pg - this.pg0);
                    else
                        this.lbfgs = optm_update_lbfgs(this.lbfgs, this.s, this.g - this.g0);
                    end
                end
                %% Determine a new search direction `d`.  Parameter `dir` is set to:
                %%   0 if `d` is not a search direction,
                %%   1 if `d` is unscaled steepest descent,
                %%   2 if `d` is scaled sufficient descent.
                dir = 0;
                %% Use L-BFGS approximation to compute a search direction and check
                %% that it is an acceptable descent direction.
                if this.blmvm
                    [this.d, scaled] = optm_apply_lbfgs(this.lbfgs, -this.pg);
                    this.d = this.d .* this.freevars;
                else
                    [this.d, scaled] = optm_apply_lbfgs(this.lbfgs, -this.g, this.freevars);
                end
                dg = optm_inner(this.d, this.g);
                if ~scaled
                    %% No exploitable curvature information, `d` is the unscaled
                    %% steepest feasible direction, that is the opposite of the
                    %% projected gradient.
                    dir = 1;
                else
                    %% Some exploitable curvature information were available.
                    dir = 2;
                    if dg >= 0
                        %% L-BFGS approximation does not yield a descent direction.
                        dir = 0; % discard search direction
                        if ~this.bounded
                            this.status = this.OPL_STATUS_NOT_POSITIVE_DEFINITE ; % optm_status('NOT_POSITIVE_DEFINITE');
                            this.endingMessage = ['L-BFGS approximation is not positive definite !'];
                            flag=this.OPTI_STOP;
                            return;
                        end
                    elseif this.epsilon > 0
                        %% A more restrictive criterion has been specified for
                        %% accepting a descent direction.
                        if dg > -this.epsilon*optm_norm2(this.d)*this.gnorm
                            dir = 0; % discard search direction
                        end
                    end
                end
                if dir == 0
                    %% No exploitable information about the Hessian is available or
                    %% the direction computed using the L-BFGS approximation failed
                    %% to be a sufficient descent direction.  Take the steepest
                    %% feasible descent direction.
                    if this.bounded
                        this.d = -this.g .* this.freevars;
                    else
                        this.d = -this.g;
                    end
                    dg = -this.gnorm^2;
                    dir = 1; % scaling needed
                end
                if dir ~= 2 && this.iters > 0
                    this.rejects = this.rejects + 1;
                end
                %% Determine the length `alpha` of the initial step along `d`.
                if dir == 2
                    %% The search direction is already scaled.
                    this.alpha = 1.0;
                else
                    %% Find a suitable step size along the steepest feasible
                    %% descent direction `d`.  Note that `gnorm`, the Euclidean
                    %% norm of the (projected) gradient, is also that of `d` in
                    %% that case.
                    this.alpha = optm_steepest_descent_step(this.xopt, this.gnorm, this.f, this.f2nd, this.fmin, ...
                        this.dxrel, this.dxabs);

                end
                if this.bounded
                    %% Safeguard the step to avoid searching in a region where
                    %% all bounds are overreached.
                    [this.amin, this.amax] = optm_line_search_limits(this.xopt, this.lower, this.upper, ...
                        this.d, this.alpha);
                    this.alpha = min(this.alpha, this.amax);
                end
                %% Initialize line-search.
                this.lnsrch = optm_start_line_search(this.lnsrch, this.f, dg, this.alpha);
                this.stage = this.lnsrch.stage;
                if this.stage ~= this.OPL_STAGE_LNSRCH_IN_PROGRESS
                    % Convergence, or error, or warning
                    this.endingMessage = ['something is wrong in line search !'];
                    flag=this.OPTI_STOP;
                    return;
                end
                %% Save iterate at start of line-search.
                this.f0 = this.f;
                this.g0 = this.g;
                this.x0 = this.xopt;
                if this.blmvm
                    this.pg0 = this.pg;
                end
            end
            %% Compute next iterate.
            if this.alpha == 1
                this.xopt = this.x0 + this.d;
            else
                this.xopt = this.x0 + this.alpha*this.d;
            end
        end

	%% Method for the "homemade" verbose mode found in optm_vmlmb.m
        function print_iteration(this)
            time = @() 86400E3*now(); % yields number of milliseconds
            if this.iters < 1
                fprintf('%s%s\n%s%s\n', ...
                    '# Iter.   Time (ms)    Eval. Reject.', ...
                    '       Obj. Func.           Grad.       Step', ...
                    '# ----------------------------------', ...
                    '-----------------------------------------------');
            end
            fprintf('%7d %11.3f %7d %7d %23.15e %11.3e %11.3e\n', ...
                this.iters, time()-this.t0, this.evals, this.rejects, this.f, this.gnorm, this.alpha);
        end

    end
    
end
