classdef OptiConjGrad < Opti
    % Conjugate gradient optimization algorithm which solves the linear
    % system \\(\\mathrm{Ax=b}\\) by minimizing
    % $$ C(\\mathrm{x})= \\frac12 \\mathrm{x^TAx - b^Tx} $$
    %
    % :param A: symmetric definite positive :class:`LinOp`
    % :param b: right-hand term
    % :param W: :class:`LinOpDiag`, if set then the algorithm solves \\(\\mathrm{A^TWAx=b}\\)
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Note**: In this algorithm the parameter cost is not fixed to the above functional
    % but to the squared resildual \\( 0.5\\Vert \\mathrm{Ax - b}\\Vert^2 \\)
    %
    % See also :class:`Opti`, :class:`OutputOpti` :class:`Cost`

	%     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
		A;  % Linear operator
		b;  % right hand side term
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		r; % residual
    end
    
    methods
    	%% Constructor
    	function this=OptiConjGrad(A,b,W,OutOp)
    		this.name='Opti Conjugate Gradient';
    		if nargin<=2, W=[]; end
    		if nargin==4 && ~isempty(OutOp),this.OutOp=OutOp;end  
    		if isempty(W)
    			this.A=A;
    		else
    			assert(isa(W,'LinOpDiag'),'W must be a LinOpDiag object');
    			this.A=A'*W*A; 			
    		end
    		this.cost=CostL2(this.A,b);
    		assert(isequal(this.A.sizeout,size(b)),'A sizeout and size of b must be equal');
    		this.b=b;
        end 
        %% Set data b
        function set_b(this,b)
            % Set the right-hand side \\(\\mathrm{b}\\)
        	
            assert(isequal(this.A.sizeout,size(b)),'A sizeout and size of b must be equal');
            this.b=b;
        end
    	%% Run the algorithm
        function run(this,x0) 
            % Reimplementation from :class:`Opti`.
            
			if ~isempty(x0) % To restart from current state if wanted
				this.xopt=x0;
				this.r= this.b - this.A.apply(this.xopt);
			end;  
			assert(~isempty(this.xopt),'Missing starting point x0');
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();			
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				% - Algorithm iteration
   				rho = dot(this.r(:),this.r(:));
    			if this.niter == 2 && ~isempty(x0)
        			p = this.r;
                else
                    beta = rho/rho_prec;
                    p = this.r + beta*p;
    			end
        		q = this.A*p;
    			alpha = rho/dot(p(:), q(:));
    			this.xopt = this.xopt + alpha*p;
    			this.r = this.r - alpha*q;
    			rho_prec = rho;
				% - Call OutputOpti object
				if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end
			end 
			this.time=toc(tstart);
			this.ending_verb();
        end
	end
end
