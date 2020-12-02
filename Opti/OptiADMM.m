classdef OptiADMM < Opti
    % Alternating Direction Method of Multipliers [1] algorithm which minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F_0(\\mathrm{x}) + \\sum_{n=1}^N F_n(\\mathrm{H_n x}) $$
    %
    % :param F_0: :class:`Cost` object
    % :param F_n: cell of N :class:`Cost` with an implementation of the :meth:`prox` for each one
    % :param H_n: cell of N :class:`LinOp`
    % :param rho_n: array of N positive scalars
    % :param maxiterCG: maximal number of inner conjugate gradient (CG) iterations (when required, default 20)
    % :param solver: a handle function taking parameters solver(z_n,rho_n,x0) (see the note below)
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
    % is performed  either using the conjugate-gradient :class:`OptiConjGrad` algorithm, a direct inversion or the given solver
    %
    %  - If \\(F_0\\) is empty or is a :class:`CostL2`, then if the :class:`LinOp` \\(\\sum_{n=0}^N \\mathrm{H_n}^*\\mathrm{H_n}\\)
    %    is not invertible the :class:`OptiConjGrad` is used by default if no more efficient solver is provided.
    %    Here \\(\\mathrm{H_0}\\) is the :class:`LinOp` associated to \\(F_0\\).
    %
    %  - Otherwise the solver is required.
    %
    % **Reference**
    %
    % [1] Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction
    % method of multipliers." Foundations and Trends in Machine Learning, 2011.
    %
    % **Example** ADMM=OptiADMM(F0,Fn,Hn,rho_n,solver)
    %
    % See also :class:`Opti`, :class:`OptiConjGrad` :class:`OutputOpti`, :class:`Cost`

    %%    Copyright (C) 2017
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    properties (SetAccess = public,GetAccess = public)
        F0=[];               % func F0
        Fn;                  % Cell containing the Cost Fn
        Hn;                  % Cell containing the LinOp Hn
        solver=[];           % solver for the last step of the algorithm
        A;     % LinOp for conjugate gradient (if used)
    %end
    % Full protected properties
    %properties (SetAccess = protected,GetAccess = {?OutputOpti ,
    %?TestCvg})  % This makes problems for doc compilation. For the moment
    %I let all the attributes public
        yn;    % Internal parameters
        wn;
        Hnx;
        b0=0;
        
        
        % To manage the new maxiterCG property of ADMM (from v1.1.2)
        % without "breaking" script written with previous versions of
        % GlobalBioIm
        prevIterCG=20;
    end
    % Full public properties
    properties
        rho_n;                 % vector containing the multipliers
        CG;                    % conjugate gradient algo (Opti Object, when used)
        maxiterCG=20;
    end
    
    methods
        %% Constructor
        function this=OptiADMM(F0,Fn,Hn,rho_n,solver)
            this.name='Opti ADMM';
            if ~isempty(F0), this.F0=F0; end
            if nargin<=3, rho_n=ones(size(Hn)); end
            if nargin<=4, solver=[]; end
            if isscalar(rho_n)
                rho_n=rho_n*ones(size(Hn));
            end
            assert(length(Fn)==length(Hn),'Fn, Hn and rho_n must have the same length');
            assert(length(Hn)==length(rho_n),'Fn, Hn and rho_n must have the same length');
            this.Fn=Fn;
            this.Hn=Hn;
            this.rho_n=rho_n;
            if ~isempty(F0)
                assert(~isempty(solver) || isa(F0,'CostL2') || isa(F0,'CostL2Composition') || ...
                    (isa(F0,'CostMultiplication') && F0.isnum && (isa(F0.cost2,'CostL2') || isa(F0.cost2,'CostL2Composition'))) || ...
                    ( isa(F0,'CostSummation')  && all(cellfun(@(x) (isa(x,'CostL2Composition') || isa(x,'CostL2') || ...
                    (isa(x,'CostMultiplication') && x.isnum && (isa(x.cost2,'CostL2') || isa(x.cost2,'CostL2Composition'))) ), F0.mapsCell))), ...
                    'when F0 is nonempty and is not CostL2 or CostL2Composition (or a sum of them), a solver must be given (see help)');
                this.cost=F0 +Fn{1}*Hn{1};
            else
                this.cost=Fn{1}*Hn{1};
            end
            this.solver=solver;
            for n=2:length(Fn)
                this.cost=this.cost+Fn{n}*Hn{n};
            end
            if isempty(this.solver)
                this.A=this.rho_n(1)* ((this.Hn{1})' *this.Hn{1});
                for n=2:length(this.Hn)
                    this.A=this.A + this.rho_n(n)*((this.Hn{n})' * this.Hn{n});
                end
                if ~isempty(this.F0)
                    if isa(this.F0,'CostL2Composition')
                        this.A=this.A + this.F0.H2'*(this.F0.H1.W*this.F0.H2);
                        this.b0=this.F0.H2'*this.F0.H1.y;
                    elseif isa(this.F0,'CostL2')
                        this.A=this.A + LinOpDiag(this.A.sizein);
                        this.b0=this.F0.y;
                    elseif isa(this.F0,'CostMultiplication') && F0.isnum
                        if isa(F0.cost2,'CostL2')
                            this.A=this.A + this.F0.cost1 * LinOpDiag(this.A.sizein);
                            this.b0=this.F0.cost1 * this.F0.cost2.y;
                        elseif isa(F0.cost2,'CostL2Composition')
                            this.A=this.A + this.F0.cost1 * this.F0.cost2.H2'*(this.F0.cost2.H1.W*this.F0.cost2.H2);
                            this.b0=this.F0.cost1 * this.F0.cost2.H2'*this.F0.cost2.H1.y;
                        else
                            error('If F0 is not a CostL2 / CostL2Composition / CostSummation of them, a solver is required in ADMM');
                        end
                    elseif isa(this.F0,'CostSummation')
                        for n=1:length(this.F0.mapsCell)
                            if isa(this.F0.mapsCell{n},'CostL2Composition')
                                this.A=this.A   + this.F0.alpha(n)*this.F0.mapsCell{n}.H2'*this.F0.mapsCell{n}.H2;
                                this.b0=this.b0 + this.F0.alpha(n)*this.F0.mapsCell{n}.H2'*this.F0.mapsCell{n}.H1.y;
                            elseif isa(this.F0.mapsCell{n},'CostL2')
                                this.A=this.A + this.F0.alpha(n)*LinOpDiag(this.A.sizein);
                                this.b0=this.b0 + this.F0.alpha(n)*this.F0.mapsCell{n}.y;
                            elseif isa(this.F0.mapsCell{n}.cost2,'CostL2Composition') % If CostMultiplication with a scalar
                                this.A=this.A   + this.F0.alpha(n)*this.F0.mapsCell{n}.cost1*this.F0.mapsCell{n}.cost2.H2'*this.F0.mapsCell{n}.cost2.H2;
                                this.b0=this.b0 + this.F0.alpha(n)*this.F0.mapsCell{n}.cost1*this.F0.mapsCell{n}.cost2.H2'*this.F0.mapsCell{n}.cost2.H1.y;
                            elseif isa(this.F0.mapsCell{n}.cost2,'CostL2')
                                this.A=this.A + this.F0.alpha(n)*this.F0.mapsCell{n}.cost1*LinOpDiag(this.A.sizein);
                                this.b0=this.b0 + this.F0.alpha(n)*this.F0.mapsCell{n}.cost1*this.F0.mapsCell{n}.cost2.y;
                            end
                        end                      
                    else
                        error('If F0 is not a CostL2 / CostL2Composition / CostSummation of them, a solver is required in ADMM');
                    end
                end
                if ~this.A.isInvertible  % If A is non invertible -> intanciate a CG
                    this.CG=OptiConjGrad(this.A,zeros_(this.A.sizeout));
                    this.CG.verbose=0;
                    this.CG.maxiter=this.maxiterCG;
                    this.CG.ItUpOut=0;
                    warning('ADMM will use a Conjugate Gradient to compute the linear step: This may lead to slow computations.');
                end
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            % To manage the new maxiterCG property of ADMM (from v1.1.2)
            % without "breaking" script written with previous versions of
            % GlobalBioIm
            if ~isempty(this.CG)
                if this.CG.maxiter==this.prevIterCG && this.CG.maxiter~=this.maxiterCG
                    this.CG.maxiter=this.maxiterCG;
                    this.prevIterCG=this.maxiterCG;
                else
                    this.maxiterCG=this.CG.maxiter;
                    this.prevIterCG=this.CG.maxiter;
                end
            end
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                for n=1:length(this.Hn)
                    this.yn{n}=this.Hn{n}.apply(x0);
                    this.Hnx{n}=this.yn{n};
                    this.wn{n}=zeros_(size(this.yn{n}));
                end
            end
            % This is done here in case one change the y in F0 between two
            % runs
            if isempty(this.solver)
                if ~isempty(this.F0) && isa(this.F0,'CostL2Composition')
                    this.b0=this.F0.H2'*this.F0.H1.y;
                elseif ~isempty(this.F0) && isa(this.F0,'CostL2')
                    this.b0=this.F0.y;
                end
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].
            for n=1:length(this.Fn)
                this.yn{n}=this.Fn{n}.applyProx(this.Hnx{n} - this.wn{n},1/this.rho_n(n));
            end
            if isempty(this.solver)
                b=this.rho_n(1)*this.Hn{1}.applyAdjoint(this.yn{1}+this.wn{1});
                for n=2:length(this.Hn)
                    b=b+this.rho_n(n)*this.Hn{n}.applyAdjoint(this.yn{n}+this.wn{n});
                end
                if ~isempty(this.b0)
                    b=b+this.b0;
                end
                if this.A.isInvertible
                    this.xopt=this.A.applyInverse(b);
                else
                    this.CG.set_b(b);
                    this.CG.run(this.xopt);
                    this.xopt=this.CG.xopt;
                end
            else
                for n=1:length(this.Fn)
                    zn{n} = this.yn{n}+this.wn{n};
                end
                this.xopt=this.solver(zn,this.rho_n, this.xopt);
                clear zn;
            end
            for n=1:length(this.wn)
                this.Hnx{n}=this.Hn{n}.apply(this.xopt);
                this.wn{n}=this.wn{n} - (this.Hnx{n}-this.yn{n});
            end
            flag=0;
        end
        
        
    end
end
