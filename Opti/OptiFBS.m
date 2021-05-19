classdef OptiFBS < Opti
    % Forward-Backward Splitting optimization algorithm [1] which minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F(\\mathrm{x}) + G(\\mathrm{x}) $$
    %
    % :param F: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`applyGrad`).
    % :param G: a :class:`Cost` with an implementation of the :meth:`applyProx`.
    % :param gam: descent step
    % :param fista: boolean true if the accelerated version FISTA [3] is used (default false)
    % :param momRestart: boolean true if the moment restart strategy is used [4](default false)
    % :param updateGam: Rule for updating gamma (none : default, reduced : the parameter gam is decreased according to \\(\\gamma / \\sqrt{k} \\), backtracking : backtracking rule following [3])
    % :param eta: parameter greater than 1 that is used with backtracking (see [3])
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    % **Note**: When the functional are convex and \\(F\\) has a Lipschitz continuous gradient, convergence is
    % ensured by taking \\(\\gamma \\in (0,2/L] \\) where \\(L\\) is the Lipschitz constant of \\(\\nabla F\\) (see [1]).
    % When FISTA is used [3], \\(\\gamma \\) should be in \\((0,1/L]\\). For nonconvex functions [2] take \\(\\gamma \\in (0,1/L]\\).
    % If \\(L\\) is known (i.e. F.lip different from -1), parameter \\(\\gamma\\) is automatically set to \\(1/L\\).
    %
    % **References**:
    %
    % [1] P.L. Combettes and V.R. Wajs, "Signal recovery by proximal forward-backward splitting", SIAM Journal on
    % Multiscale Modeling & Simulation, vol 4, no. 4, pp 1168-1200, (2005).
    %
    % [2] Hedy Attouch, Jerome Bolte and Benar Fux Svaiter "Convergence of descent methods for semi-algebraic and
    % tame problems: proximal algorithms, forward-backward splitting, and regularized gaussiedel methods."
    % Mathematical Programming, 137 (2013).
    %
    % [3] Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear inverse Problems",
    % SIAM Journal on Imaging Science, vol 2, no. 1, pp 182-202 (2009)
    %
    % [4] Brendan O'donoghue and Emmanuel Candès. 2015. Adaptive Restart for Accelerated Gradient Schemes. 
    % Found. Comput. Math. 15, 3 (June 2015), 715-732.
    % 
    % **Example** FBS=OptiFBS(F,G)
    %
    % See also :class:`Opti` :class:`OutputOpti` :class:`Cost`
       
    %%     Copyright (C) 2017
    %     E. Soubies emmanuel.soubies@irit.fr
    %     T-A. Pham thanh-an.pham@epfl.ch
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
    
    % Full protected properties
    properties (SetAccess = protected,GetAccess = public)
        y;         % Internal parameters
        tk;
        F;                   % Cost F
        G;                   % Cost G
    end
    % Full public properties
    properties
        fista=false;         % FISTA option [3]
        momRestart = false;  % boolean true if the moment restart strategy is used [4]
        gam=[];              % descent step
        
        updateGam = 'none'; % reduce the step size (TODO : add the possibilty to chose the update rule)
        alpha = 1;   % see Kamilov paper
        eta = 1.1;  % parameter > 1 used in the backtracking rule
    end
    
    methods
        
        function this=OptiFBS(F,G)
            this.name='Opti FBS';
            this.cost=F+G;
            this.F=F;
            this.G=G;
            this.needxold = true;
            if F.lip~=-1
                this.gam=1/F.lip;
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                assert(~isempty(this.gam),'parameter gam is not setted');
                if this.fista
                    this.tk=1;
                    this.y=x0;
                end
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1-3].
            
            if strcmp(this.updateGam,'backtracking')
                if isa(this.F,'CostPartialSummation')
                    orig_subset = this.F.subset;
                end
                if this.fista
                    tmp=this.y;
                    this.xopt=this.G.applyProx(this.y - this.gam.*this.F.applyGrad(this.y),this.gam);
                else
                    tmp=this.xopt;
                    this.xopt=this.G.applyProx(this.xopt - this.gam.*this.F.applyGrad(this.xopt),this.gam);
                end
                g = this.F*tmp;
                grad = this.F.applyGrad(tmp);
                if isa(this.F,'CostPartialSummation')
                    this.F.subset = orig_subset;
                end
                t = this.gam;
                
                satisfied = this.F*this.xopt ...
                    <= g + dot(grad(:),(this.xopt(:) - tmp(:))) + norm(this.xopt(:) - tmp(:))^2/(2*t);
                
                while ~gather(satisfied)
                    t = t/this.eta;
                    this.xopt = this.G.applyProx(tmp - t*grad,t);%p_L(yk)
                    satisfied = this.F*this.xopt ...
                        <= g + dot(grad(:),(this.xopt(:) - tmp(:))) + norm(this.xopt(:) - tmp(:))^2/(2*t);
                end
                this.gam=t;
            else
                if this.fista
                    this.xopt=this.G.applyProx(this.y - this.gam.*this.F.applyGrad(this.y),this.gam);
                else
                    this.xopt=this.G.applyProx(this.xopt - this.gam.*this.F.applyGrad(this.xopt),this.gam);
                end
                if strcmp(this.updateGam,'reduced')
                    this.gam = this.gam/sqrt(this.niter+1)*max(sqrt(this.niter),1);
                end
            end
            if this.fista  % if fista
                if this.momRestart && dot(this.y(:) - this.xopt(:),this.xopt(:) - this.xold(:)) > 0
                    fprintf('Restarting\n');
                    this.y = this.xold;
                    this.tk = 1;
                else
                    told=this.tk;
                    this.tk=0.5*(1+sqrt(1+4*this.tk^2));
                    this.y=this.xopt + this.alpha*(told-1)/this.tk*(this.xopt - this.xold);
                end
            end
            
            flag=this.OPTI_NEXT_IT;            
        end
    end
end
