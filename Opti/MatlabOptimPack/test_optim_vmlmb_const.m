 % Matlab function similar to the C program test_optim_vmlmb_const.c
 % to test the Matlab interface to the OptimPack optimisation Package
 % for a constrained minimization of the rosenbrock function
 %
 % The Optimpack functions can be called through the 
 % m_vmlmb_first and m_vmlmb_next functions which are interface between Matlab 
 % and the OptimPack functions op_vmlmb_first and op_vmlmb_next 
 % or through the optim_vmlmb.m optimisation function
 % 
 % The results are identical to the test_optim_vmlmb_const.c ones.
  clc
  clear all

  % Fonction to minimize and its gradient
  fun = @(x) ros(x);
  grad_fun = @(x) grad_ros(x);

  % Initial value
  x_init = [ -1 ; 1];
  nbitermax=1000;  

  % lower and upper bounds
  xmin = [-1 ; 0] ; 
  xmax = [ 1 ; .75];

  
%% Test with calls to the m_vmlmb_first.c and m_vmlmb_next.c mex_files 
  
  n=length(x_init); m=3; fatol=0.0; frtol=1e-12; sftol=0.001; sgtol=0.9; sxtol=0.1; epsilon=0.01; costheta=0.4;
  [csave, isave, dsave] = m_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, sxtol, epsilon, costheta);
  task = isave(3);

  OP_TASK_START  = 0; % first entry, start search 
  OP_TASK_FG     = 1; % computation of F and G requested 
  OP_TASK_NEWX   = 2; % new improved solution available for inspection 
  OP_TASK_CONV   = 3; % search has converged 
  OP_TASK_WARN   = 4; % search aborted with warning 
  OP_TASK_ERROR  = 5; % search aborted with error 

  % Initialisaton
  x(1:n) = x_init; % Warning : side effect on x_init if x=x_init (due to the fact that x is passed-by-reference in the mexfiles)
  x = x(:);
  fprintf('Minimization with VMLMB\n')
  fprintf('Initialization x= [%f ; %f]\n',x(1),x(2));
  fprintf('Lower bounds xmin = [%f ; %f]\n',xmin(1),xmin(2));
  fprintf('Upper bounds xmax = [%f ; %f]\n',xmax(1),xmax(2));
  f = fun(x)+1e8;
  task = OP_TASK_FG;
  nbeval=0;
  list_x = zeros(length(x),nbitermax);
  list_J = zeros(1,nbitermax);
  fprintf('it\t nbeval\t x \t\t\t f\t df\t normg\t task\tstage\n');
  % Iterations
  for iter=1:nbitermax
      if (task == OP_TASK_FG)
        % apply bound constraints
        % op_bounds_apply(n, x, xmin, xmax);
        if ~isempty(xmin) 
            test = (x<xmin);
            if sum(test), x(test) = xmin(test); end
        end
        if ~isempty(xmax)
            test = (x>xmax);
            if sum(test), x(test) = xmax(test); end
        end
        % x = x.*(x>=xmin).*(x<=xmax) + xmin.*(x<xmin) + xmax.*(x>xmax);
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
        %fprintf('Convergence, or error, or warning :\n');
        break;
      end
      if ( (nbeval==1) || (task == OP_TASK_NEWX)) 
        % Computes set of active parameters :
        % op_bounds_active(n, active, x, g, xmin, xmax);
        if (isempty(xmin) && isempty(xmax))
            active = [];
        elseif isempty(xmin),
            active = int32( (g>0) + (x<xmax) );
        elseif isempty(xmax),
            active = int32( (x>xmin) + (g<0) );
        else
            active = int32( ( (x>xmin) + (g<0) ).*( (x<xmax) + (g>0) ) );
        end
      end
      % Computes next step:
      [task, csave]= m_vmlmb_next(x,f,g,active,isave,dsave);
      list_x(:,iter) = x;
      list_J(iter) = f;
      fprintf('%d\t%d\t[%7.2e ; %7.2e]\t%6.2g\t%6.2g\t%6.2g\t    %d\t%d\n',iter,nbeval,x(1),x(2),f,fold-f,normg,task,isave(4));
  end
  iter = iter-1;
  list_x = list_x(:,1:iter);
  list_J = list_J(1:iter);
  
  fprintf('vmlmb for %d iterations : %s\n' ,iter, csave);
  fprintf('function = %.8f, grad = [%.8f ; %.8f]\n' , f, g(1), g(2));
  fprintf('Solution : x= [%.8f ; %.8f]\n', x(1),x(2));
  
  % Display the results
  figure(1)
  xx = linspace(min(list_x(1,:)),max(list_x(1,:)),100);
  yy = linspace(min(list_x(2,:)),max(list_x(2,:)),100);
  [XX,YY] = meshgrid(xx,yy);
  ZZ = (1.0-XX).^2 + 100.0*(YY-XX.^2).^2;
  imagesc(xx,yy,ZZ); hold on; axis xy
  plot(list_x(1,:),list_x(2,:),'*-r',x(1),x(2),'ow')
  if (~isempty(xmin) && ~isempty(xmax))
    plot([xmin(1) xmax(1) xmax(1) xmin(1) xmin(1)],[xmin(2) xmin(2) xmax(2) xmax(2) xmin(2)],'g'); 
  elseif (~isempty(xmin) && isempty(xmax))
    plot([xmin(1) max(list_x(1,:)) max(list_x(1,:)) xmin(1) xmin(1)],[xmin(2) xmin(2) max(list_x(2,:)) max(list_x(2,:)) xmin(2)],'g'); 
  elseif (isempty(xmin) && ~isempty(xmax))
    plot([min(list_x(1,:)) xmax(1) xmax(1) min(list_x(1,:)) min(list_x(1,:))],[min(list_x(2,:)) min(list_x(2,:)) xmax(2) xmax(2) min(list_x(2,:))],'g'); 
  end
  hold off
  
  
  %% Test with call to the optim_vmlmb.m function
 
  option = [3; 0.0; 1e-8; 0.001; 0.9; 0.1; 0.01; 0.4; 10000; 2];

  fprintf('Initialization x= [%f ; %f]\n',x(1),x(2));
  fprintf('Lower bounds xmin = [%f ; %f]\n',xmin(1),xmin(2));
  fprintf('Upper bounds xmax = [%f ; %f]\n',xmax(1),xmax(2));
  [x, J , list_x] = optim_vmlmb(fun,grad_fun,x_init,xmin,xmax,option);
  f = fun(x);
  g = grad_fun(x);
  fprintf('function = %.8f, grad = [%.8f ; %.8f]\n' , f, g(1), g(2));
  fprintf('Solution : x= [%.8f ; %.8f]\n', x(1),x(2));
  
  % Display the results
  figure(2)
  xx = linspace(min(list_x(1,:)),max(list_x(1,:)),100);
  yy = linspace(min(list_x(2,:)),max(list_x(2,:)),100);
  [XX,YY] = meshgrid(xx,yy);
  ZZ = (1.0-XX).^2 + 100.0*(YY-XX.^2).^2;
  imagesc(xx,yy,ZZ); hold on; axis xy
  plot(list_x(1,:),list_x(2,:),'*-r',x(1),x(2),'ow'); 
  if (~isempty(xmin) && ~isempty(xmax))
    plot([xmin(1) xmax(1) xmax(1) xmin(1) xmin(1)],[xmin(2) xmin(2) xmax(2) xmax(2) xmin(2)],'g'); 
  elseif (~isempty(xmin) && isempty(xmax))
    plot([xmin(1) max(list_x(1,:)) max(list_x(1,:)) xmin(1) xmin(1)],[xmin(2) xmin(2) max(list_x(2,:)) max(list_x(2,:)) xmin(2)],'g'); 
  elseif (isempty(xmin) && ~isempty(xmax))
    plot([min(list_x(1,:)) xmax(1) xmax(1) min(list_x(1,:)) min(list_x(1,:))],[min(list_x(2,:)) min(list_x(2,:)) xmax(2) xmax(2) min(list_x(2,:))],'g'); 
  end
  hold off
   
  