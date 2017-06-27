classdef OptiADMM < Opti
    %% OptiADMM : Alternating Direction Method of Multipliers algorithm
    %  Matlab inverse Problems Library
    %
    % -- Description
    % Implements the ADMM algorithm [1] to minimize:
    %    $$ F_0(x) + \sum_{n=1}^N F_n(y_n) $$
    % subject to:
    %    $$ H_n*x=y_n \forall n \in {1,...,N}$$
    % where the H_n are linear operators and F_n are functional with an implementation of the
    % proximity operator for n = 1,...,N (and not necessarily for n=0)
    %
    % In fact the algorithm aims to minimize the Lagrangian formulation of the above problem:
    % $$ L(x,y_1...y_n,w_1...w_) = F_0(x) + \sum_{n=1}^N 0.5*\rho_n*||H_n*x - y_n + w_n/rho_n||^2 + F_n(y_n)$$
    % where the \rho_n >0 n=1...N are the multipliers.
    %
    % -- Example
    %   ADMM= OptiADMM(F0,H0,Fn,Hn,rho_n,solver,OutOp)
    % where F0 is a FUNC object, H0 a LINOP object, Fn a cell of N COST, Hn a cell of N LINOP,
    % rho_n a vector of N nonnegative scalars and solver a function handle such that:
    %   solver(z_n,rho_n,x0)
    % where z_n is a cell of N elements, rho_n is as above.
    % The solver minimizes the following function starting from x0:
    %    $$ F_0(x) + \sum_{n=1}^N 0.5*\rho_n||H_n*x -z_n||^2 $$
    % Finally OutOp is an OutputOpti object.
    %
    % Note: If F0 not empty and F0 not CostL2, then solver is mandatory.
    %       Otherwise not and by default the ADMM algorithm will use the Conjugate Gradient algorithm
    %       (see OptiConjGrad) to make the minimization task of solver. However, if one has a
    %       faster method than applying a conjugate gradient to perform this step, it is
    %       recommended to provide a solver. If F0 is nonempty, then solver is MANDATORY.
    %
    % -- Properties
    % * |maxiterCG|   number of iteration for Conjugate Gradient (when used)
    % * |rho_n|       vector containing the multipliers
    % * |OutOpCG|     OutputOpti object for Conjugate Gradient (when used)
    % * |ItUpOutCG|   ItUpOut parameter for Conjugate Gracdient (when used, default 0)
    %
    % -- References
    % [1] Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction
    %     method of multipliers." Foundations and Trends in Machine Learning, 2011.
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OptiConjGrad, LinOp, Cost, OutputOpti
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
		F0=[];               % func F0
		Fn;                  % Cell containing the Cost Fn
		Hn;                  % Cell containing the LinOp Hn
		solver=[];           % solver for the last step of the algorithm
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
		yn;    % Internal parameters
		zn;
		wn;
		Hnx;
        b0=[];
		A;     % LinOp for conjugate gradient (if used)
    end
    % Full public properties
    properties
		rho_n;                 % vector containing the multipliers
        CG;                    % conjugate gradient algot (Opti Object, when used)
    end
    
    methods
    	%% Constructor
        function this=OptiADMM(F0,Fn,Hn,rho_n,solver,OutOp)
            this.name='Opti ADMM';
            if ~isempty(F0), this.F0=F0; end
            if nargin<=4, solver=[]; end
            if nargin==6 && ~isempty(OutOp),this.OutOp=OutOp;end
            assert(length(Fn)==length(Hn),'Fn, Hn and rho_n must have the same length');
            assert(length(Hn)==length(rho_n),'Fn, Hn and rho_n must have the same length');
            this.Fn=Fn;
            this.Hn=Hn;
            this.rho_n=rho_n;
            if ~isempty(F0)
                assert(~isempty(solver) || isa(F0,'CostL2'),'when F0 is nonempty and is not CostL2 a solver must be given (see help)');
                this.cost=F0 +Fn{1}.o(Hn{1});
            else
                this.cost=Fn{1}.o(Hn{1});
            end
            this.solver=solver;
            for n=2:length(Fn)
                this.cost=this.cost+Fn{n}.o(Hn{n});
            end
            if isempty(this.solver)
	    this.A=SumLinOp({ (this.Hn{1})' *this.Hn{1}},[this.rho_n(1)]);
                for n=2:length(this.Hn)
                    this.A=SumLinOp({ this.A, (this.Hn{n})' * this.Hn{n}},[1,this.rho_n(n)]);
                end
                if ~isempty(this.F0) && isa(this.F0,'CostL2')
                    this.A=SumLinOp({this.A,this.F0.H'*this.F0.H},[1,1]);
                    this.b0=this.F0.H'*this.F0.y;
                end
                this.CG=OptiConjGrad(this.A,zeros(this.A.sizeout),[],OutputOpti());
                this.CG.maxiter=20;
                this.CG.ItUpOut=0;
            end
        end
    	%% Run the algorithm
        function run(this,x0)
            if ~isempty(x0), % To restart from current state if wanted
                this.xopt=x0;
                for n=1:length(this.Hn)
                    this.yn{n}=this.Hn{n}.apply(this.xopt);
                    this.Hnx{n}=this.yn{n};
                    this.wn{n}=zeros(size(this.yn{n}));
                end
            end;
            assert(~isempty(this.xopt),'Missing starting point x0');
            tstart=tic;
            this.OutOp.init();
            this.niter=1;
            this.starting_verb();
            while (this.niter<this.maxiter)
                this.niter=this.niter+1;
                xold=this.xopt;
                % - Algorithm iteration
                for n=1:length(this.Fn)
                    this.yn{n}=this.Fn{n}.prox(this.Hnx{n} + this.wn{n}/this.rho_n(n),1/this.rho_n(n));
                    this.zn{n}=this.yn{n}-this.wn{n}/this.rho_n(n);
                end
                if isempty(this.solver)
                    b=this.rho_n(1)*this.Hn{1}.adjoint(this.zn{1});
                    for n=2:length(this.Hn)
                        b=b+this.rho_n(n)*this.Hn{n}.adjoint(this.zn{n});
                    end
                    if ~isempty(this.b0)
                        b=b+this.b0;
                    end
                    this.CG.set_b(b);
                    this.CG.run(this.xopt);
                    this.xopt=this.CG.xopt;
                else
                    this.xopt=this.solver(this.zn,this.rho_n, this.xopt);
                end
                for n=1:length(this.wn)
                    this.Hnx{n}=this.Hn{n}.apply(this.xopt);
                    this.wn{n}=this.wn{n} + this.rho_n(n)*(this.Hnx{n}-this.yn{n});
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
