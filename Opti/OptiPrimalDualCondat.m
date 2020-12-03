classdef OptiPrimalDualCondat < Opti
    % Primal-Dual algorithm proposed by L. Condat in [1] which minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x})= F_0(\\mathrm{x}) + G(\\mathrm{x}) + \\sum_n F_n(\\mathrm{H_nx}) $$
    %
    % :param F_0: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`grad`).
    % :param G: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param F_n: cell of N :class:`Cost` with an implementation of the :meth:`prox` for each one
    % :param H_n: cell of N :class:`LinOp`
    % :param tau: parameter of the algorithm (see the note below)
    % :param sig: parameter of the algorithm (see the note below)
    % :param rho: parameter of the algorithm (see the note below)
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    % **Note**:
    %
    %   - When \\(F_0=0\\), parameters sig and tau have to verify
    %     $$ \\sigma \\times \\tau \\Vert \\sum_n \\mathrm{H_n^*H_n} \\Vert \\leq 1 $$
    %     and \\(\\rho \\in ]0,2[\\), to ensure convergence (see [1, Theorem 5.3]).
    %
    %   - Otherwise, when \\(F_0\\neq 0\\), parameters sig and tau have to verify
    %     $$ \\frac{1}{\\tau} - \\sigma \\times \\Vert \\sum_n \\mathrm{H_n^*H_n} \\Vert \\geq \\frac{\\beta}{2} $$
    %     where \\(\\beta\\) is the Lipschitz constant of \\(\\nabla F\\) and we need \\(\\rho \\in ]0,\\delta[ \\) with
    %     $$ \\delta = 2 - \\frac{\\beta}{2}\\times\\left(\\frac{1}{\\tau}
    %     - \\sigma \\times \\Vert \\sum_n \\mathrm{H_n^*H_n}  \\Vert\\right)^{-1} \\in [1,2[ $$
    %     to ensure convergence (see [1, Theorem 5.1]).
    %
    % **Reference**
    %
    % [1] Laurent Condat, "A Primal-Dual Splitting Method for Convex Optimization Involving Lipchitzian, Proximable and Linear
    % Composite Terms", Journal of Optimization Theory and Applications, vol 158, no 2, pp 460-479 (2013).
    %
    % **Example** A=OptiPrimalDualCondat(F0,G,Fn,Hn)
    %
    % See also :class:`Opti`, :class:`OutputOpti`, :class:`Cost`
    
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
    properties (SetAccess = protected,GetAccess = public)
        F0;    % cost F0
        G;     % cost G
        Fn;    % costs F_n (cell)
        Hn;    % associated LinOp (cell))
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        y;    % cell containing the dual variables
    end
    % Full public properties
    properties
        tau;       % parameter of the algorithm
        sig;       % parameter of the algorithm
        rho=1.95;  % parameter of the algorithm
    end
    
    methods
        %% Constructor
        function this=OptiPrimalDualCondat(F0,G,Fn,Hn)
            this.name='Opti Primal-Dual Condat';
            assert(length(Fn)==length(Hn),'Fn, Hn and rho_n must have the same length');
            this.Fn=Fn;
            this.Hn=Hn;
            this.F0=F0;
            this.G=G;
            this.needxold = true;
            if ~isempty(F0), this.cost=F0;end
            if ~isempty(G)
                if isempty(this.cost), this.cost=G;
                else, this.cost=this.cost + G; end
            end
            if ~isempty(Fn)
                if isempty(this.cost), this.cost=Fn{1}*Hn{1};
                else, this.cost=this.cost + Fn{1}*Hn{1}; end
            end
            for n=2:length(Fn)
                this.cost=this.cost+Fn{n}*Hn{n};
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                % initialization of the dual variables y
                for n=1:length(this.Hn)
                    this.y{n}=zeros_(this.Hn{n}.sizeout);
                end
            end
            % Check parameters
            assert(~isempty(this.sig),'parameter sig is not set');
            assert(~isempty(this.tau),'parameter tau is not set');
            assert(~isempty(this.rho),'parameter rho is not set');
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].
            % Update xtilde
            
            if ~isempty(this.F0)
                temp=this.xopt-this.tau*this.F0.applyGrad(this.xopt);
            else
                temp=this.xopt;
            end
            for n=1:length(this.Hn)
                temp=temp-this.tau*this.Hn{n}.applyAdjoint(this.y{n});
            end
            if ~isempty(this.G)
                xtilde=this.G.applyProx(temp,this.tau);
            else
                xtilde=temp;
            end
            % Update xopt
            this.xopt=this.rho*xtilde+(1-this.rho)*this.xopt;
            % Update ytilde and y
            for n=1:length(this.Fn)
                ytilde=this.Fn{n}.applyProxFench(this.y{n}+this.sig*this.Hn{n}.apply(2*xtilde-this.xold),this.sig);
                this.y{n}=this.rho*ytilde +(1-this.rho)*this.y{n};
            end  

            flag=this.OPTI_NEXT_IT;
        end
    end
end
