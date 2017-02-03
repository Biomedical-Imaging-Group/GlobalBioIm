classdef OptiRichLucy < Opti
    %% OptiRichLucy : Richardson-Lucy algorithm
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implements the Richardson-Lucy algorithm [1] to minimize the functional (FuncKullLeib):
	% $$\sum_n -d_n log((H*x)_n + bet) + (Hx)_n$$
    % where H is a LinOp object (default LinOpIdentity), d are the data
    % and bet is a small scalar (default 1e-3) to smooth the function at zero.
    %
    % -- Example
    % OptiRL=OptiRichLucy(F,OutOp)
    % where F is a FuncKullLeib object and OutOp a OutputOpti object.
    %
    % -- References
	% [1] Lucy, Leon B. "An iterative technique for the rectification of observed distributions" The astronomical journal (1974)
	% [2] Richardson, William Hadley. "Bayesian-based iterative method of image restoration." JOSA (1972): 55-59.
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, FuncKullLeib, OutputOpti
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
   
    methods
    	%% Constructor
    	function this=OptiRichLucy(F,OutOp)
    		this.name='Opti Richardson-Lucy';
    		assert(isa(F,'FuncKullLeib'), 'The minimized functional should be the FuncKullLeib');
    		this.cost=F;
    		if nargin==2 && ~isempty(OutOp)
    			this.OutOp=OutOp;
    		end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			if ~isempty(x0),this.xopt=x0;end;  % To restart from current state if wanted
			assert(~isempty(this.xopt),'Missing starting point x0');
        	data=this.cost.data;
        	He1=this.cost.H.Adjoint(ones(this.cost.sizein));
        	bet=this.cost.bet;
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;
				% - Algorithm iteration
				this.xopt=this.xopt./He1.*this.cost.H.Adjoint(data./(this.cost.H.Apply(this.xopt)+bet));			
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
