classdef OptiConjGrad < Opti
    %% OptiConjGrad : Conjugate gradient optimization algorithm
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Solves the linear system A*x=b by minimizing 
    %   $$ 0.5 x'*A*w - b'*x $$
    % using the Conjugate Gradient algorithm
    %
    % Note: In this algorithm the parameter cost is not fixed to the above functional
    %       but to the squered resildual:
    %         $$ 0.5||A*x - b||^2
    %
    % -- Example
    % CG = OptiConjGrad(A,b,W,OutOp)
    % where A is a LINOP object, b the data, W a LINOP so that the algorithm solves A'*W*Ax=b
    % and OutOp a OutputOpti object
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OutputOpti
    %
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
    			assert(isa(W,'LinOp'),'W must be a LinOp object');
    			this.A=A'*W*A; 			
    		end
    		this.cost=FuncLeastSquares(b,this.A);
    		assert(isequal(this.A.sizeout,size(b)),'A sizeout and size of b must be equal');
    		this.b=b;
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			if ~isempty(x0) % To restart from current state if wanted
				this.xopt=x0;
				this.r= this.b - this.A.Apply(this.xopt);
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
        			this.xtol = eps*rho;
        			p = this.r;
    			elseif rho <= this.xtol
        			break
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
