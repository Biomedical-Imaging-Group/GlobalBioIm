classdef OptiDouglasRachford < Opti
	%% OptiDouglasRachford : Douglas Rachford splitting optimization algorithm
	%
	% Implements the Douglas Rachford splitting algorithm that solves:
	% x^+ = argmin_x( F1(x) + F2(L x))
	% with the proxes of F1 and F2 known
	% L a linear operator  such  L.L^T = nu I
	% y0 the initialization
	% maxiter the number of iteration
	% gamma \in [0,+\inf[
	% lambda\in ]0,2[ the relaxation parmeter.
	%
	
	properties
		nu = 1;
		useL = 0;
		F1;
		F2;
		lambda
		gamma
	end
	
	methods
		function this = OptiDouglasRachford(F1, F2, L, gamma, lambda, OutOp)
			% F1, F2
			this.F1 = F1;
			this.F2 = F2;
			
			% L
			if exist('L', 'var') && ~isempty(L)
				this.useL = 1;
				r = randn(L.sizeout);
				nu = r ./ L.HHt(r);
				assert(std(nu(:)) <1e-6, 'LLt != nu I');
				this.nu = real(mean(nu(:)));
				if this.nu==1
					this.useL = 2;
				end
			end
			
			% gamma
			if numel(gamma)==1
				gamma = [gamma, gamma];
			end
			this.gamma = gamma;
			
			% lambda
			this.lambda = lambda;
			
			this.OutOp = OutOp;
			
		end
		
		function xopt = run(this, x0)
			y = x0;
			this.xopt = x0;
			
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			
			while(this.niter < this.maxiter)
				xold=this.xopt;		
				if this.useL
					Ly = this.L*y;
					if this.useL==2
						this.xopt = L.Adjoint( this.F2.prox(Ly, this.gamma(2)));
					else
						this.xopt = y + (1./this.nu).* L.Adjoint( this.F2.prox(Ly, this.nu.*this.gamma(2)) - Ly);
					end
				else
					this.xopt = this.F2.prox(y, this.gamma(2));
				end
				y = y + this.lambda .* ( this.F1.prox(2.*this.xopt- y,this.gamma(1)) - this.xopt);
				
				if this.test_convergence(xold), break; end
				
				if (mod(this.niter,this.ItUpOut)==0)
					this.OutOp.update(this)
				end
				
				this.niter = this.niter + 1;
			end
			
			xopt = this.xopt;
			this.time=toc(tstart);
			this.ending_verb();
		end
		
	end
	
end
