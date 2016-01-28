function [x, J , list_x] = optim_cont(fun, grad_fun, x_init, x_min, x_max, option);

% Function optim_cont : bound constraint optimisation with gradient
%    [x [, J , list_x]] = optim_cont(fun, grad_fun, x_init [, x_min, x_max, option]);
%
%    Use the VMLMB - limited memory variable metric method (BFGS) - of the
%    OptimPack package to minimise function fun with/without bound constraints 
%    Require the OptimPack package and the mexfiles of the programs m_vmlmb_first.c
%    m_vmlmb_next.c compiled with :         
%           mex m_vmlmb_first.c -loptimpack
%           mex m_vmlmb_next.c -loptimpack
%
%   INPUT :
%       - fun : function to be optimized
%       - grad_fun : gradient of function fun
%       - x_init : initial value for the parameters
%       - x_min : lower bound for parameters x 
%               (same size than x_init or empty vector)
%       - x_max : upper bound for parameters x 
%               (same size than x_init or empty vector)
%       - option : 
%         [ m, fatol, frtol, sftol, sgtol, sxtol, epsilon, costheta, nbitermax, verb]
%               vector of options for the algorithm 
%               (see optimpack.h for more details on parameters m, fatol, frtol, sftol, 
%               sgtol, sxtol, epsilon, costheta)
%               * default values : [ 3; 0.0; 1e-8; 0.001; 0.9; 0.1; 0.01; 0.4; 10000; 1];
%               For large problems, m = 3 to 5 gives good results.  
%               For  small problems, m should  be less or equal to n=length(x_init).  
%               The larger is  m (and n)  the more computer  memory will be  needed to
%               store the workspaces (see DSAVE).
%               * nbitermax : maximum number of iterations
%               * verb : 0, 1 or 2, verbose mode
%   OUTPUT :
%       - x : solution
%       - J : list of explored function values (not mendatory)
%       - list_x : list of explored parameters values (not mendatory)

    % Default values for the option vector
    if nargin<6,
        option = [3; 0.0; 1e-8; 0.001; 0.9; 0.1; 0.01; 0.4; 10000; 1];
    end
    if (length(option)<8),
        help optim_cont
        error('vectors option should be of length 8 to 10')
    elseif (length(option)<9),
        nbitermax = 1000;
        verb=1;
    elseif (length(option)<10),
        nbitermax=option(9);
        verb=1;
    else
        nbitermax=option(9);
        verb=option(10);
    end
    n = length(x_init);
    m = int32(option(1));

    fatol=option(2); frtol=option(3); sftol=option(4); sgtol=option(5); 
    sxtol=option(6); epsilon=option(7); costheta=option(8); 

    n = length(x_init);    
    x(1:n) = x_init; % Warning : side effects if x=x_init; is used ==> x_init value changed with x !
    x = x(:);
    % Control of the x_min and x_max parameters
    if ( isempty(x_max) && isempty(x_min) )
        active = [];
    else
        if ( (~isempty(x_min)) && (~isempty(x_max)) )
            if ( ~(length(x_min)==length(x)) || ~(length(x_max)==length(x)) )
                help optim_cont
                error('Vectors x_min and x_max should be empty or with the same length than x_init')
            elseif ( sum(x_min>x_max) )
                help optim_cont
                error('Incompatible bounds constraints : x_min > x_max for some parameters')
            end
        end
    end

    % Initialization of the ouptut parameters
    if nargout>1,
        J = zeros(1,nbitermax);
        if nargout==3,
            list_x = zeros(length(x),nbitermax);
        end
    end
    
    % Constants used to controle the algorithms tasks
    OP_TASK_START  = 0; % first entry, start search 
    OP_TASK_FG     = 1; % computation of F and G requested 
    OP_TASK_NEWX   = 2; % new improved solution available for inspection 
    OP_TASK_CONV   = 3; % search has converged 
    OP_TASK_WARN   = 4; % search aborted with warning 
    OP_TASK_ERROR  = 5; % search aborted with error 

    % Initialization with m_vmlmb_first
    [csave, isave, dsave] = m_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, sxtol, epsilon, costheta);
    task = isave(3);
    f = fun(x)+1e8;
    task = OP_TASK_FG;
    nbeval=0;
    if verb==2
        % fprintf('Minimisation with vmlmb\n');    
        fprintf('it   nbeval\t f\t\t df\t\t normg\t task stage\n');
    end
    % iterations
    for iter=1:nbitermax
        if (task == OP_TASK_FG)
            % apply bound constraints
            % op_bounds_apply(n, x, xmin, xmax);
            if ~isempty(x_min) 
                test = (x<x_min);
                if sum(test), x(test) = x_min(test); end
            end
            if ~isempty(x_max)
                test = (x>x_max);
                if sum(test), x(test) = x_max(test); end
            end
            fold = f;
            f = fun(x);         % evaluate the function at X; store in F
            g = grad_fun(x);    % evaluate the gradient of F at X; store in G	
            normg= sum(g.^2);
            nbeval=nbeval+1;
        elseif (task == OP_TASK_NEWX) 
            % New successful step: the approximation X, function F, and
            % gradient G, are available for inspection.
        else
            % Convergence, or error, or warning
            % fprintf('Convergence, or error, or warning :\n');
            break;
        end
        if ( (nbeval==1) || (task == OP_TASK_NEWX)) 
            % Computes set of active parameters :
            % op_bounds_active(n, active, x, g, xmin, xmax);
            if (isempty(x_min) && isempty(x_max))
                active = [];
            elseif isempty(x_min),
                active = int32( (g>0) + (x<x_max) );
            elseif isempty(x_max),
                active = int32( (x>x_min) + (g<0) );
            else
                active = int32( ( (x>x_min) + (g<0) ).*( (x<x_max) + (g>0) ) );
            end
        end
        % Computes next step:
        [task, csave]= m_vmlmb_next( x, f, g, active,isave, dsave);
        if nargout>1,
            J(iter) = f;
            if nargout==3,
                list_x(:,iter) = x;
            end
        end
        if verb==2
            fprintf('%d\t%d\t%6.2e\t%6.2e\t%6.2e    %d    %d\n',iter,nbeval,f,fold-f,normg,task,isave(4));
        end
    end
    if (iter==nbitermax)
        fprintf('vmlmb : maximum number of iteration %d reached\n' ,iter);
    else
        iter = iter-1;
        if nargout>1,
            J = J(1:iter);
            if nargout==3,
                list_x = list_x(:,1:iter);
            end
        end
        if verb,
              fprintf('vmlmb for %d iterations : %s\n' ,iter, csave);
        end
    end