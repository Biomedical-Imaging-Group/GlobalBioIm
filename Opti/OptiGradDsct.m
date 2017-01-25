classdef OptiGradDsct < handle
    %% OptiGradDsct : Gradient Descent optimization algorithm
    %  Matlab Inverse Problems Library
    %  The Opti meta class implements generic methods for all optimization algorithms
    %
    % -- Properties
    % * |name|      - name of the optimization algorithm (inherited from parent Opti class)
    % * |cost|      - functional to minimize (inherited from parent Opti class,should have
    %                 an implementation of the gradient)
    % * |gam|       - descent step 
    %
    % Note: If the Functional F is gradient Lipschitz, gam has to be lower than 2/L where
    %       L is the Lipschitz constant of the gradient. The optimal choice is 1/L (see [1]).
    %       If F.lip is known (i.e. different from -1), parameter gam is automatically setted to 1/L
    %
    % -- Methods
    % * |run(x0)|   - run the algorithm from the initial point x0. If x0=[], restarts from the current state
    %
    % -- References 
    % [1] Nesterov, Yurii. "Introductory lectures on convex programming." Lecture Notes (1998): 119-120.
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
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'Opti Gradient Descent'   % name of the optimization algorithm
    end
    properties
    	gam;      % descent step
    end
    
    methods
    	%% Constructor
    	function this=OptiGradDsct(F,verbup)
    		this.cost=F;
    		if F.lip~=-1
    			gam=1/F.lip;
    		end
    		if nargin==2
    			this.verbup=verbup;
    		end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			tstart=tic;
			starting_verb();
			while (this.niter<this.maxiter)
				this.verbup.init();
				xold=this.xopt;
				% - Algorithm iteration
				this.xopt=this.xopt-this.gam*this.cost.grad(this.xopt);
				% - Convergence test
				if test_convergence(xold), break; end
				% - Call VerbUpdate object
				if (mod(this.niter,this.verb)==0),this.verbup.exec(this);end
				this.niter=this.niter+1;
			end 
			ending_verb();
			this.time=toc(tstart);
        end
	end
end
