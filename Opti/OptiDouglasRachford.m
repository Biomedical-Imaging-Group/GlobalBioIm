classdef OptiDouglasRachford < Opti
	% Douglas Rachford splitting optimization algorithm which minimizes
	% :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F_1(\\mathrm{x}) + F_2(\\mathrm{L x}) $$
    %
    % :param F_1: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param F_2: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param L: a :class:`LinOp` such that \\(\\mathrm{LL^T} = \\nu \\mathrm{I} \\)
	% :param gamma: \\(\\in [0,+\\inf[\\)
	% :param lambda: \\(\\in ]0,2[\\) the relaxation parmeter
    %
    % All attributes of parent class :class:`Opti` are inherited. 
	%
    % **Example** DR=OptiDouglasRachford(F1, F2, L, gamma, lambda, OutOp)
    %
	% See also :class:`Opti`, :class:`OutputOpti`, :class:`Cost`
    
    %%    Copyright (C) (2017)
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
            % Reimplementation from :class:`Opti`.
            
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
						this.xopt = L.adjoint( this.F2.prox(Ly, this.gamma(2)));
					else
						this.xopt = y + (1./this.nu).* L.adjoint( this.F2.prox(Ly, this.nu.*this.gamma(2)) - Ly);
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
