classdef OptiGradDsct < Opti
    % Gradient Descent optimization algorithm to minimize a differentiable :class:`Cost` \\(C(\\mathrm{x})\\)
    %
    % :param C: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`grad`).
    % :param gam: descent step
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Note** If the cost \\(C\\) is gradient Lipschitz, convergence is ensured by taking 
    % \\(\\gamma \\in (0,2/L] \\) where \\(L\\) is the Lipschitz constant of \\(\\nabla C\\) (see [1]).
	% The optimal choice is \\(\\gamma = 1/L \\) (see [1]). If \\(L\\) is known (i.e. F.lip different from -1), 
    % parameter \\(\\gamma\\) is automatically set to \\(1/L\\).
    %
    % **Reference**
    %
    % [1] Nesterov, Yurii. "Introductory lectures on convex programming." Lecture Notes (1998): 119-120.
    %
    % See also :class:`Opti` :class:`OutputOpti` :class:`Cost`
    
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

	% Full public properties
    properties
    	gam=[];      % descent step
    end
    
    methods
    	%% Constructor
    	function this=OptiGradDsct(F,OutOp)
    		this.name='Opti Gradient Descent';
    		this.cost=F;
    		if F.lip~=-1
    			this.gam=1/F.lip;
    		end
    		if nargin==2 && ~isempty(OutOp)
    			this.OutOp=OutOp;
    		end
    	end 
    	%% Run the algorithm
        function xopt = run(this,x0) 
            % Reimplementation from :class:`Opti`.
            
        	if isempty(this.gam), error('Parameter gam is not setted'); end
			if ~isempty(x0),this.xopt=x0;end;  % To restart from current state if wanted
			assert(~isempty(this.xopt),'Missing starting point x0');
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			while (this.niter<=this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;
				% - Algorithm iteration
				this.xopt=this.xopt-this.gam*this.cost.grad(this.xopt);
				% - Convergence test
				if this.test_convergence(xold), break; end
				% - Call OutputOpti object
				if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end
			end 
			this.time=toc(tstart);
			this.ending_verb();
			xopt = this.xopt;
        end
	end
end
