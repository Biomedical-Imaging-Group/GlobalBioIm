classdef OptiFBS < Opti
    % Forward-Backward Splitting optimization algorithm [1] which minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F(\\mathrm{x}) + G(\\mathrm{x}) $$
    %
    % :param F: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`applyGrad`).
    % :param G: a :class:`Cost` with an implementation of the :meth:`applyProx`.
    % :param gam: descent step
    % :param fista: boolean true if the accelerated version FISTA [3] is used (default false)
    % :param partialGrad: when F is a :class:`CostSummation`, uses a subset of this sum of Costs to compute the gradient. 0: desactivated (default), 1: circular shift of uniform subset, 2: stochastic subset selection
    % :param set: when partialGrad>0, subset of F.mapsCells on which the partial gradient is applied (remaining cost will always be used in the gradient computation).
    % :param Lsub: number of costs used to compute the partial gradient (Lsub < length(set))
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
    % **Example** FBS=OptiFBS(F,G,OutOp)
    %
    % See also :class:`Opti` :class:`OutputOpti` :class:`Cost`
    
    
    %%     Copyright (C) 2017
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    properties (SetAccess = protected,GetAccess = protected)
        y;         % Internal parameters
        tk;
        counter=0; % counter for subset update
        set;       %indices of mapsCell (CostSummation) to consider as an "angle" (e.g. F is composed of 100 CostL2 (i.e., angles) + 1 CostHyperbolic)
        nonset;    %indices of mapsCell (CostSummation) not to consider as an "angle"
        L;         % Total number of available "angles"
        subset; % current subset of angles used to compute F grad
    end
    % Full public properties
    properties
        F;  % Cost F
        G;  % Cost G
        
        partialGrad=0; % activate partial gradient option
        gam=[];        % descent step
        fista=false;   % FISTA option [3]
        
        reducedStep = false; % reduce the step size (TODO : add the possibilty to chose the update rule)
        mingam;
        alpha = 1; % see Kamilov paper
        
        Lsub;   % Number of costs used in partial gradient
    end
    
    methods
        
        function this=OptiFBS(F,G,OutOp)
            this.name='Opti FBS';
            this.cost=F+G;
            this.F=F;
            this.G=G;
            if F.lip~=-1
                this.gam=1/F.lip;
            end
            if nargin==3 && ~isempty(OutOp)
                this.OutOp=OutOp;
            end
            if isa(F,'CostSummation')
                this.L = F.numMaps;
                this.Lsub=this.L;
                this.updateSet(1:this.L);
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                assert(~isempty(this.gam),'parameter gam is not setted');
                this.counter=0;
                if this.fista
                    this.tk=1;
                    this.y=x0;
                end
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1-3].
            
            if this.fista  % if fista
                this.xopt=this.G.applyProx(this.y - this.gam*this.computeGrad(this.y),this.gam);
                told=this.tk;
                this.tk=0.5*(1+sqrt(1+4*this.tk^2));
                this.y=this.xopt + this.alpha*(told-1)/this.tk*(this.xopt - this.xold);
            else
                this.xopt=this.G.applyProx(this.xopt - this.gam*this.computeGrad(this.xopt),this.gam);
            end
            flag=0;
        end   
        function updateSet(this,new_set)
            % Changes the attribute set
            
            this.set = new_set;
            this.L = length(new_set);
            this.nonset = find(~ismember(1:this.F.numMaps,this.set));
        end   
        
        function grad = computeGrad(this,x)
            if this.partialGrad >0
                switch this.partialGrad
                    case 1
                        this.updateSubset(randi(this.L,this.Lsub,1));
                    case 2
                        this.updateSubset(sort(1 + mod(round(this.counter...
                            + (1:this.L/this.Lsub:this.L)),this.L)));
                        this.counter = this.counter + 1;
                end
                grad = zeros(size(x));
                for kk = 1:this.Lsub
                    ind = this.set(this.subset(kk));
                    grad = grad + this.F.alpha(ind)*this.F.mapsCell{ind}.applyGrad(x);
                end
                grad = real(grad)/this.Lsub;%ad hoc               
                for kk = 1:length(this.nonset)
                    grad = grad + this.F.alpha(this.nonset(kk))*this.F.mapsCell{this.nonset(kk)}.applyGrad(x);
                end
            else
                grad = this.F.applyGrad(x);
            end
        end
        
        function updateSubset(this,subset)
            this.subset = subset;
            this.Lsub = length(subset);
        end        
    end
end
