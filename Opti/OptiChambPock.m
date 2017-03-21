classdef OptiChambPock < Opti
    %% OptiChambPock : Chambolle-Pock optimization algorithm
    %  Matlab inverse Problems Library
    %
    % -- Description
    % Implements the Chambolle-Pock algorithm [1] to minimize a function of the form:
	%      F(Hx) + G(x)
	% where F and G are two functionals (Func) where F and G have an implementation for
	% the proximity operator (.prox). H is a linear operator (LinOp).
	%
	% Note: In fact, F needs only the prox of its fenchel transform (which is implemented 
	%       as soon as F as an implementation of the prox, see FUNC superclass).
	%
    % -- Example
    % OptiCP = OptiChambPock(F,H,G,OutOp)
    % where F and G are FUNC objects, H LINOP object and OutOp a OutputOpti object 
    %
    % -- Properties
    % * |tau|    parameter of the algorithm (see the note below) default 1
    % * |sig|    parameter of the algorithm (see the note below) default computed according to
    %            the inequality below (if H.norm is implemented)
    % * |gam|    if non-empty, then the accelerated version of the algorithm is used (see [1])
    %            Hence G of F^* is uniformly with Grad(G^*) or Grad(F) is 1/gam-Lipschitz. 
	%            If G is uniformly convex then set the parameter var to 1
	%            If F^* is uniformly convex then set the parameter var to 2
    % * |var|    select the "bar" variable of the algorithm (see [1]):
	%              - if 1 (default) then the primal variable xbar = x_n + theta(x_n - x_{n-1}) is used 
	%              - if 2 then the dual variable ybar = y_n + theta(y_n - y_{n-1}) is used
    %   
    % Notes: 1- to ensure convergence (see [1]), parameters sig and tau have to verify 
	%             sig*tau*||H||^2 < 1
	%           where ||H|| denotes the norm of the linear operator H.
	%        2- when the accelerated version is used (i.e. parameter gam is non-empty), 
	%           sig and tau will be updated at each iteration and the initial ones (given by user) 
	%           have to verify
	%             sig*tau*||K||^2 <= 1
    %
    % -- References
    % [1] Chambolle, Antonin, and Thomas Pock. "A first-order primal-dual algorithm for convex problems with 
	%     applications to imaging." Journal of Mathematical Imaging and Vision 40.1, pp 120-145 (2011).	
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OutputOpti
    %
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		F;  % Func F
		G;  % Func G
		H;  % LinOp H
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		theta=1; % parameter of the algorithm (fixed to 1 but can could be set as a user defined 
		         % parameter for future upgrades of the Library)
    end
    % Full public properties
    properties
		sig=[];   % parameter of the algorithm
		tau=1;    % parameter of the algorithm
		gam=[];   % parameter of the algorithm for accelerated version
		var=1;    % select the bar variable in the algorith (primal 1/ dual 2)
    end
    
    methods
    	%% Constructor
    	function this=OptiChambPock(F,H,G,OutOp)
    		this.name='Opti Chambolle-Pock';
    		this.cost=F.o(H)+G;
    		this.F=F;
    		this.G=G;
    		this.H=H;
    		if nargin==4 && ~isempty(OutOp)
    			this.OutOp=OutOp;
    		end
    		if this.H.norm>=0
    			this.sig=1/(this.tau*this.H.norm^2)-eps;
    		end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
        	% Check parameters
        	assert(~isempty(this.sig),'parameter sig is not setted');
        	assert(~isempty(this.tau),'parameter tau is not setted');
        	if ~isempty(x0) % To restart from current state if wanted
				this.xopt=x0;
			end
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			y=this.H.Apply(this.xopt);
			if this.var==1
				xbar=this.xopt;
				Kxbar=y;
			else
				ybar=y;
				KTybar=this.H.adjoint(ybar);
				KTy=KTybar;
			end
			Kxopt=y;
			sig=this.sig; tau=this.tau; theta=this.theta; gam=this.gam;
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;			
				if this.var==1 % === using xbar
					Kxold=Kxopt;
					% - Algorithm iteration
					y=this.F.prox_fench(y+sig*Kxbar,sig);
					this.xopt=this.G.prox(this.xopt-tau*this.H.adjoint(y),tau);
					Kxopt=this.H.Apply(this.xopt);
					if ~isempty(gam) % acceleration => uodate theta, tau and sig according to [1]
						theta=1/sqrt(1+2*gam*tau);
						tau=theta*tau;			 
						sig=sig/theta;           
					end
					xbar=this.xopt+theta*(this.xopt - xold);    
					Kxbar=Kxopt+theta*(Kxopt - Kxold);
				else % === using ybar
					yold=y;
					KTyold=KTy;
					% -- Algorithm iteration
					this.xopt=this.G.prox(this.xopt-tau*KTybar,tau);
					Kxopt=this.H.Apply(this.xopt);
					y=this.F.prox_fench(y+sig*Kxopt,sig);
					KTy=this.H.adjoint(y);
					if ~isempty(gam) % acceleration => uodate theta, tau and sig according to [1]
						theta=1/sqrt(1+2*gam*sig); 
						tau=tau/theta;		 
						sig=sig*theta;        
					end
					ybar=(1+theta)*y - theta*yold;    
					KTybar=(1+theta)*KTy - theta*KTyold;
				end				
				% - Convergence test
				if this.test_convergence(xold), break; end
				% - Call OutputOpti object
				if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end
			end 
			this.time=toc(tstart);
			this.ending_verb();
        end
	end
end
