classdef OptiADMM < Opti
    % Alternating Direction Method of Multipliers [1] algorithm which minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F_0(\\mathrm{x}) + \\sum_{n=1}^N F_n(\\mathrm{H_n x}) $$
    %
    % :param F_0: :class:`Cost` object
    % :param F_n: cell of N :class:`Cost` with an implementation of the :meth:`prox` for each one
    % :param H_n: cell of N :class:`LinOp`
    % :param rho_n: array of N positive scalars
    % :param solver: a handle function taking parameters solver(z_n,rho_n,x0) (see the note below)
    % :param maxiterCG: max number of iterations for conjugate-gradient (CG) (when used)
    % :param OutOpCG: :class:`OutputOpti` object for CG (when used)
    % :param ItUpOutCG: :attr:`ItUpOut` parameter for CG (when used, default 0)
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Principle** 
    % The algorithm aims at minimizing the Lagrangian formulation of the above problem:
    % $$ \\mathcal{L}(\\mathrm{x,y_1...y_n,w_1...w_n}) = F_0(\\mathrm{x}) + \\sum_{n=1}^N \\frac12\\rho_n\\Vert \\mathrm{H_nx - y_n + w_n/\\rho_n} \\Vert^2 + F_n(\\mathrm{y_n})$$
    % using an alternating minimization scheme [1].
    %
    % **Note** The minimization of \\(\\mathcal{L}\\) over \\(\\mathrm{x}\\), 
    % $$ F_0(\\mathrm{x}) + \\sum_{n=1}^N \\frac12\\rho_n\\Vert \\mathrm{H_nx -z_n}\\Vert^2, \\quad \\mathrm{z_n= y_n - w_n/\\rho_n} $$
    % is performed  either using the conjugate-gradient :class:`OptiConjGrad` algoriothm or the given solver
    %
    %  - If \\(F_0\\) is empty or is a :class:`CostL2`, then :class:`OptiConjGrad` is used by
    %    default if no more efficient solver is provided. 
    %
    %  - Otherwithe the solver is required.
    %
    % **Reference**
    %
    % [1] Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction
    % method of multipliers." Foundations and Trends in Machine Learning, 2011.
    %
    % See also :class:`Opti`, :class:`OptiConjGrad` :class:`OutputOpti`, :class:`Cost`
    
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
            % Reimplementation from :class:`Opti`.
            
            if ~isempty(x0) % To restart from current state if wanted
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
                        b=b+this.F0.H'*this.F0.y;
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
