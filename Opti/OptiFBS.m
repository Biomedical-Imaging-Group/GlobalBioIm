classdef OptiFBS < Opti
    %% OptiFBS : Forward-Backward Splitting optimization algorithm
    %  Matlab Inverse Problems Library
    %  Implements the Forward-Backward Splitting algorithm [1] to minimize a function of the form:
    %         F(x) + w*G(x)
    % where  G has  an implementation for the proximity operator (.prox) and F is differentiable 
    % (i.e. has gradient (.grad))
    %
    % -- Example
    % OptiGD=OptiFBS(F,G,w,verbup)
    % where F and G are FUNC object, w a positive real and verbup a VerbUpdate object 
    % 
    % -- Properties
    % * |name|      - name of the optimization algorithm (inherited from parent Opti class)
    % * |cost|      - functional to minimize (inherited from parent Opti class,should have
    %                 an implementation of the gradient)
    % * |gam|       - descent step (public to be setted by the user if necessary)
    % * |fista|     - boolean true if the accelerated version FISTA [3] is used (default false)
	%
	% Note: when the functional are convex and F has a Lipschitz continuous gradient, convergence is
	%       ensured by taking gam in (0,2/L] where L is the Lipschitz constant of grad(F) (see [1]).
	%       When FISTA is used [3], gam should be in (0,1/L]. For nonconvex functions [2] take gam in (0,1/L].    
    %       If F.lip is known (i.e. different from -1), parameter gam is automatically setted to 1/L
    %
    % -- Methods
    % * |run(x0)|   - run the algorithm from the initial point x0. If x0=[], restarts from the current state
    %
    % -- References 
	% [1] P.L. Combettes and V.R. Wajs, "Signal recovery by proximal forward-backward splitting", SIAM Journal on
	%     Multiscale Modeling & Simulation, vol 4, no. 4, pp 1168-1200, (2005).
	% [2] Hedy Attouch, Jerome Bolte and Benar Fux Svaiter "Convergence of descent methods for semi-algebraic and 
	%     tame problems: proximal algorithms, forward-backward splitting, and regularized gaussiedel methods." 
	%     Mathematical Programming, 137 (2013).
	% [3] Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems",
	%     SIAM Journal on Imaging Science, vol 2, no. 1, pp 182-202 (2009)
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti
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
		w;  % weight w
    end
    % Full public properties
    properties
    	gam;           % descent step
    	fista=false;   % FISTA option [3]
    end
    
    methods
    	%% Constructor
    	function this=OptiFBS(F,G,w,verbup)
    		this.name='Opti FBS';
    		this.cost=F+FuncMultScalar(G,w);
    		this.F=F;
    		this.G=G;
    		this.w=w;
    		if F.lip~=-1
    			this.gam=1/F.lip;
    		end
    		if nargin==4
    			this.verbup=verbup;
    		end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
        	if isempty(this.gam), error('Parameter gam is not setted'); end
			tstart=tic;
			this.verbup.init();
			this.xopt=x0;
			this.niter=1;
			this.starting_verb();
			if this.fista,tk=1; y=this.xopt; end
			while (this.niter<this.maxiter)
				xold=this.xopt;
				% - Algorithm iteration
				if this.fista  % if fista
					this.xopt=this.G.prox(y - this.gam*this.F.grad(y),this.gam*this.w);
					told=tk;
					tk=0.5*(1+sqrt(1+4*tk^2));
					y=this.xopt + (told-1)/tk*(this.xopt-xold);
				else 
					this.xopt=this.G.prox(this.xopt - this.gam*this.F.grad(this.xopt),this.gam*this.w);
				end
				% - Convergence test
				if this.test_convergence(xold), break; end
				% - Call VerbUpdate object
				if (mod(this.niter,this.verb)==0),this.verbup.exec(this);end
				this.niter=this.niter+1;
			end 
			this.time=toc(tstart);
			this.ending_verb();
        end
	end
end
