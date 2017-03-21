classdef OptiPrimalDualCondat < Opti
    %% OptiPrimalDualCondat : Primal-Dual algorithm proposed by Condat [1]
    %  Matlab inverse Problems Library
    %
    % -- Description
    % Implements the primal-dual algorithm proposed by Condat [1] to minimizes a function
    % of the form
	%	     F0(x) + G(x) sum_n F_n(H_n(x))
	% where F0 has an implementation of the gradient, G and F_n have an implementation of the
	% proximity operator and H_n are linear operators (LinOp).
	%
    % -- Example
    % Op=OptiPrimalDualCondat(F0,G,Fn,Hn,OutOp)
    % where F0 and G are FUNC objects, Fn a cell of FUNC and Hn a cell of LINOP (same length as Fn)
    % Finally OutOp a OutputOpti object.
    %
    % -- Properties
    % * |tau|    parameter of the algorithm (see the note below) 
    % * |sig|    parameter of the algorithm (see the note below) 
    % * |rho|    parameter of the algorithm (see the note below) 
    %
	% Note: 1- when F=0, parameters sig and tau have to verify
	%             sig*tau*||sum_n Hn*Hn|| <= 1
	%          and rho must be in ]0,2[, to ensure convergence (see [1, Theorem 5.3]).
	%       2- otherwise, when F~=0, parameters sig and tau have to verify
	%             1/tau - sig*||sum_n Hn*Hn|| >= bet/2
	%          where bet is the Lipschitz constant of the gradient of F and rho must be in ]0,delta[ with
	%             delta = 2 - bet/2*(1/tau - sig*||sum_n Hn*Hn||)^{-1} \in [1,2[
	%          to ensure convergence (see [1, Theorem 5.1]).
    %
    % -- References
	%    [1] Laurent Condat, "A Primal-Dual Splitting Method for Convex Optimization Involving Lipchitzian, Proximable and Linear 
	%        Composite Terms", Journal of Optimization Theory and Applications, vol 158, no 2, pp 460-479 (2013).
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OutputOpti, LinOp
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
		F0;    % function F0
		G;     % function G
		Fn;    % functions F_n (cell)
		Hn;    % associated LinOp (cell))
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		y;    % cell containing the dual variables
    end
    % Full public properties
    properties
		tau;    % parameter of the algorithm
		sig;    % parameter of the algorithm
		rho;    % parameter of the algorithm
    end
    
    methods
    	%% Constructor
    	function this=OptiPrimalDualCondat(F0,G,Fn,Hn,OutOp)
    		this.name='Opti Primal-Dual Condat';
    		if nargin==5 && ~isempty(OutOp), this.OutOp=OutOp; end
    		assert(length(Fn)==length(Hn),'Fn, Hn and rho_n must have the same length');
    		this.Fn=Fn;
    		this.Hn=Hn;
    		this.F0=F0;
    		this.G=G;
    		if ~isempty(F0), this.cost=F0;end
    		if ~isempty(G)
    			if isempty(this.cost), this.cost=G;
    			else, this.cost=this.cost + G; end
    		end
    		if ~isempty(Fn)
    			if isempty(this.cost), this.cost=Fn{1}.o(Hn{1});
    			else, this.cost=this.cost + Fn{1}.o(Hn{1}); end
    		end
    		for n=2:length(Fn)
    			this.cost=this.cost+Fn{n}.o(Hn{n});
			end
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			if ~isempty(x0)   % To restart from current state if wanted
				this.xopt=x0;
				% initialization of the dual variables y
				for n=1:length(this.Hn)
    				this.y{n}=this.Hn{n}.Apply(x0); 	 
				end
			end; 
			% Check parameters
        	assert(~isempty(this.sig),'parameter sig is not setted');
        	assert(~isempty(this.tau),'parameter tau is not setted');
        	assert(~isempty(this.rho),'parameter rho is not setted');
			assert(~isempty(this.xopt),'Missing starting point x0');
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;
				% - Algorithm iteration
				% Update xtilde
				if ~isempty(this.F0)
					temp=this.xopt-this.tau*this.F0.grad(this.xopt);
				else
					temp=this.xopt;
				end
				for n=1:length(this.Hn) 
					temp=temp-this.tau*this.Hn{n}.adjoint(this.y{n});
				end
				if ~isempty(this.G)
					xtilde=this.G.prox(temp,this.tau);
				else
					xtilde=temp;
				end
				% Update xopt
				this.xopt=this.rho*xtilde+(1-this.rho)*this.xopt;
				% Update ytilde and y
				for n=1:length(this.Fn) 
					ytilde{n}=this.Fn{n}.prox_fench(this.y{n}+this.sig*this.Hn{n}.Apply(2*xtilde-xold),this.sig);
					this.y{n}=this.rho*ytilde{n}+(1-this.rho)*this.y{n};
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
