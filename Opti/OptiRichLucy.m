classdef OptiRichLucy < Opti
    %% OptiRichLucy : Richardson-Lucy algorithm
    %  Matlab inverse Problems Library
    %
    % -- Description
    % Implements the Richardson-Lucy algorithm [1] to minimize the functional (CostKullLeib):
	% $$\sum_n -d_n log((H*x)_n + bet) + (Hx)_n$$
    % where H is a LinOp object (default LinOpIdentity), d are the data
    % and bet is a small scalar (default 1e-3) to smooth the function at zero.
    %
    % -- Example
    % OptiRL=OptiRichLucy(F,TV,lamb,OutOp)
    % where F is a CostKullLeib object, TV a boolean to activate the TV regularized version or not
    % (default false), lamb the regularization parameter used when TV and OutOp a OutputOpti object.
    %
    % -- Properties
    % *|lamb|   regularization parameter when TV is activated
    % *|epsl|   smoothing parameter for TV
    %
    % -- References
	% [1] Lucy, Leon B. "An iterative technique for the rectification of observed distributions" The astronomical journal (1974)
	% [2] Richardson, William Hadley. "Bayesian-based iterative method of image restoration." JOSA (1972): 55-59.
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, CostKullLeib, OutputOpti
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
    properties (SetAccess = public,GetAccess = public)
		epsl=1e-10; % smoothing parameter for TV
    end
    % Protected set public read    
    properties (SetAccess = public,GetAccess = public)
		TV=false;   % boolean true if RL-TV version activated (default false)
		lamb=1e-2;  % regularization parameter
	end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
		G;   % gradient operator (used if TV activated)
		Fkl; % Kullback-Leibler divergence
    end
   
    methods
    	%% Constructor
    	function this=OptiRichLucy(F,TV,lamb,OutOp)
    		this.name='Opti Richardson-Lucy';
    		assert(isa(F,'CostKullLeib'), 'The minimized functional should be the FuncKullLeib');
    		this.cost=F;
    		if nargin==4 && ~isempty(OutOp)
    			this.OutOp=OutOp;
    		end
    		if nargin>=2 && ~isempty(TV), this.TV=TV; end
    		this.Fkl=F;
    		if nargin>=3 && ~isempty(lamb), this.lamb=lamb; end
    		if this.TV
    			this.G=LinOpGrad(this.Fkl.H.sizein);
    			if length(this.Fkl.H.sizein)==2  % 2D
    				this.cost=this.cost + MultScalarCost(CostMixNorm12([3],this.G),this.lamb);
    			elseif length(this.Fkl.H.sizein)==3
    				this.cost=this.cost + MultScalarCost(CostMixNorm12([4],this.G),this.lamb);
    			end
    		end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			if ~isempty(x0),this.xopt=x0;end;  % To restart from current state if wanted
			assert(~isempty(this.xopt),'Missing starting point x0');
        	data=this.Fkl.y;
        	He1=this.Fkl.H.adjoint(ones(this.Fkl.H.sizein));
        	bet=this.Fkl.bet;
            if bet==0, error('Smoothing parameter beta has to be different from 0 (see constructor of CostKullLeib)'); end;
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;
				% - Algorithm iteration
				if ~this.TV
					this.xopt=this.xopt./He1.*this.Fkl.H.adjoint(data./(this.Fkl.H.apply(this.xopt)+bet));	
				else
					tmp=this.G.apply(this.xopt);
					if length(size(tmp))==2     % 1D
						nor=sqrt(tmp.^2+this.epsl);
					elseif length(size(tmp))==3 % 2D
						nor=repmat(sqrt(sum(tmp.^2,3)+this.epsl),[1,1,size(tmp,3)]);
					elseif length(size(tmp))==4 % 3D
						nor=repmat(sqrt(sum(tmp.^2,4)+this.epsl),[1,1,1,size(tmp,4)]);
					end
					gradReg=this.G.adjoint(tmp./nor);
					this.xopt=this.xopt./(He1 + this.lamb*gradReg).*this.Fkl.H.adjoint(data./(this.Fkl.H.apply(this.xopt)+bet));
					if sum(this.xopt(:)<0)~=0
  						warning('Violation of the positivity of the solution (the regularization parameter should be decreased).');
  					end
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
